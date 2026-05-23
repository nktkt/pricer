// pricer_server.cpp — a small REST pricing service.
//
// A self-contained HTTP/1.1 server built only on POSIX sockets (no third-party
// dependencies), exposing the pricer library over HTTP:
//   GET /health
//   GET /price?type=call&S=100&K=100&r=0.05&sigma=0.2&T=1   -> price + Greeks
//   GET /impliedvol?type=call&price=10.45&S=100&K=100&r=0.05&T=1
//   GET /mc?type=call&S=100&K=100&r=0.05&sigma=0.2&T=1&paths=2000000
//   GET /submit?...&paths=...   -> {"job_id":N}   (async Monte Carlo job)
//   GET /job?id=N               -> {"status":"running|done", ...}
//   GET /metrics                -> Prometheus text exposition (observability)
// Responses are JSON (except /metrics, which is Prometheus text). Every request
// is also logged to stderr as a structured JSON line. POSIX sockets =>
// Linux/macOS; not built on Windows.
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

#include "pricer/black_scholes.hpp"
#include "pricer/implied_vol.hpp"
#include "pricer/parallel.hpp"

using namespace pricer;
using Params = std::map<std::string, std::string>;

// --- tiny request helpers ---
static Params parse_query(const std::string& q) {
    Params p;
    size_t i = 0;
    while (i < q.size()) {
        size_t amp = q.find('&', i);
        if (amp == std::string::npos) amp = q.size();
        const std::string kv = q.substr(i, amp - i);
        const size_t eq = kv.find('=');
        if (eq != std::string::npos) p[kv.substr(0, eq)] = kv.substr(eq + 1);
        i = amp + 1;
    }
    return p;
}
static double num(const Params& p, const char* k, double def = 0.0) {
    auto it = p.find(k);
    if (it == p.end()) return def;
    try { return std::stod(it->second); } catch (...) { return def; }
}
static OptionType opt_type(const Params& p) {
    auto it = p.find("type");
    return (it != p.end() && (it->second == "put" || it->second == "Put")) ? OptionType::Put
                                                                           : OptionType::Call;
}

// --- async Monte Carlo job registry ---
struct Job {
    std::atomic<bool> done{false};
    double result{0.0};
    long paths{0};
};
static std::map<int, std::shared_ptr<Job>> g_jobs;
static std::mutex g_jobs_mtx;
static int g_next_id = 1;

// --- observability: counters, per-endpoint latency, structured request log ---
struct Metrics {
    std::atomic<long> requests_total{0};
    std::atomic<long> responses_2xx{0};
    std::atomic<long> responses_4xx{0};
    std::atomic<long> jobs_submitted{0};
    std::atomic<long> jobs_completed{0};
    std::mutex mtx;                                  // guards the per-path maps
    std::map<std::string, long> path_count;
    std::map<std::string, double> path_latency_us;  // cumulative microseconds
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
};
static Metrics g_metrics;

// Record one finished request: bump counters and accumulate per-endpoint latency.
static void record_request(const std::string& path, int code, double us) {
    g_metrics.requests_total.fetch_add(1, std::memory_order_relaxed);
    (code < 400 ? g_metrics.responses_2xx : g_metrics.responses_4xx)
        .fetch_add(1, std::memory_order_relaxed);
    std::lock_guard<std::mutex> lk(g_metrics.mtx);
    g_metrics.path_count[path] += 1;
    g_metrics.path_latency_us[path] += us;
}

// Structured (JSON) access log line on stderr — one per request, machine-parseable.
static void log_request(const std::string& path, int code, double us, std::size_t bytes) {
    std::fprintf(stderr, "{\"ts\":%lld,\"path\":\"%s\",\"status\":%d,\"dur_us\":%.1f,\"bytes\":%zu}\n",
                 static_cast<long long>(std::time(nullptr)), path.c_str(), code, us, bytes);
}

// --- route handlers (return a JSON body) ---
static std::string h_price(const Params& p) {
    const Greeks g = black_scholes_greeks(opt_type(p), num(p, "S", 100), num(p, "K", 100),
                                          num(p, "r", 0.05), num(p, "sigma", 0.2), num(p, "T", 1),
                                          num(p, "q", 0.0));
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "{\"price\":%.8f,\"delta\":%.8f,\"gamma\":%.8f,\"vega\":%.8f,"
                  "\"theta\":%.8f,\"rho\":%.8f}",
                  g.price, g.delta, g.gamma, g.vega, g.theta, g.rho);
    return buf;
}
static std::string h_impliedvol(const Params& p) {
    const double iv = implied_vol(opt_type(p), num(p, "price"), num(p, "S", 100), num(p, "K", 100),
                                  num(p, "r", 0.05), num(p, "T", 1), 1e-8, 100, num(p, "q", 0.0));
    char buf[128];
    std::snprintf(buf, sizeof buf, "{\"implied_vol\":%.8f}", iv);
    return buf;
}
static double mc_price(const Params& p, long paths) {
    const double K = num(p, "K", 100);
    auto payoff = (opt_type(p) == OptionType::Call)
                      ? std::function<double(double)>([K](double ST) { return ST > K ? ST - K : 0.0; })
                      : std::function<double(double)>([K](double ST) { return ST < K ? K - ST : 0.0; });
    return mc::price_terminal_parallel(payoff, num(p, "S", 100), num(p, "r", 0.05),
                                       num(p, "sigma", 0.2), num(p, "T", 1), paths);
}
static std::string h_mc(const Params& p) {
    char buf[128];
    std::snprintf(buf, sizeof buf, "{\"mc_price\":%.8f}", mc_price(p, (long)num(p, "paths", 1000000)));
    return buf;
}
static std::string h_submit(const Params& p) {
    const long paths = (long)num(p, "paths", 10000000);
    auto job = std::make_shared<Job>();
    job->paths = paths;
    int id;
    {
        std::lock_guard<std::mutex> lk(g_jobs_mtx);
        id = g_next_id++;
        g_jobs[id] = job;
    }
    g_metrics.jobs_submitted.fetch_add(1, std::memory_order_relaxed);
    Params copy = p;
    std::thread([job, copy, paths]() {
        const double r = mc_price(copy, paths);
        job->result = r;
        job->done.store(true, std::memory_order_release);
        g_metrics.jobs_completed.fetch_add(1, std::memory_order_relaxed);
    }).detach();
    char buf[64];
    std::snprintf(buf, sizeof buf, "{\"job_id\":%d}", id);
    return buf;
}

// Prometheus text exposition format (text/plain). Counters and gauges that a
// scraper (or `curl`) can read to monitor the service.
static std::string h_metrics(const Params&) {
    const double uptime =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - g_metrics.start).count();
    std::ostringstream o;
    o << "# HELP pricer_requests_total Total HTTP requests received.\n"
      << "# TYPE pricer_requests_total counter\n"
      << "pricer_requests_total " << g_metrics.requests_total.load() << "\n"
      << "# HELP pricer_responses_total HTTP responses by status class.\n"
      << "# TYPE pricer_responses_total counter\n"
      << "pricer_responses_total{class=\"2xx\"} " << g_metrics.responses_2xx.load() << "\n"
      << "pricer_responses_total{class=\"4xx\"} " << g_metrics.responses_4xx.load() << "\n"
      << "# HELP pricer_endpoint_requests_total Requests per endpoint.\n"
      << "# TYPE pricer_endpoint_requests_total counter\n";
    {
        std::lock_guard<std::mutex> lk(g_metrics.mtx);
        for (const auto& kv : g_metrics.path_count)
            o << "pricer_endpoint_requests_total{path=\"" << kv.first << "\"} " << kv.second << "\n";
        o << "# HELP pricer_endpoint_duration_seconds_sum Cumulative handler latency per endpoint.\n"
          << "# TYPE pricer_endpoint_duration_seconds_sum counter\n";
        for (const auto& kv : g_metrics.path_latency_us)
            o << "pricer_endpoint_duration_seconds_sum{path=\"" << kv.first << "\"} "
              << (kv.second / 1e6) << "\n";
    }
    o << "# HELP pricer_jobs_submitted_total Async Monte Carlo jobs submitted.\n"
      << "# TYPE pricer_jobs_submitted_total counter\n"
      << "pricer_jobs_submitted_total " << g_metrics.jobs_submitted.load() << "\n"
      << "# HELP pricer_jobs_completed_total Async Monte Carlo jobs completed.\n"
      << "# TYPE pricer_jobs_completed_total counter\n"
      << "pricer_jobs_completed_total " << g_metrics.jobs_completed.load() << "\n"
      << "# HELP pricer_uptime_seconds Seconds since the server started.\n"
      << "# TYPE pricer_uptime_seconds gauge\n"
      << "pricer_uptime_seconds " << uptime << "\n";
    return o.str();
}
static std::string h_job(const Params& p) {
    const int id = (int)num(p, "id", -1);
    std::shared_ptr<Job> job;
    {
        std::lock_guard<std::mutex> lk(g_jobs_mtx);
        auto it = g_jobs.find(id);
        if (it != g_jobs.end()) job = it->second;
    }
    if (!job) return "{\"error\":\"unknown job id\"}";
    char buf[128];
    if (job->done.load(std::memory_order_acquire))
        std::snprintf(buf, sizeof buf, "{\"status\":\"done\",\"mc_price\":%.8f,\"paths\":%ld}",
                      job->result, job->paths);
    else
        std::snprintf(buf, sizeof buf, "{\"status\":\"running\",\"paths\":%ld}", job->paths);
    return buf;
}

// --- minimal HTTP plumbing ---
static void send_response(int fd, int code, const char* status, const std::string& body,
                          const char* content_type) {
    char header[256];
    const int n = std::snprintf(header, sizeof header,
                                "HTTP/1.1 %d %s\r\nContent-Type: %s\r\n"
                                "Content-Length: %zu\r\nConnection: close\r\n\r\n",
                                code, status, content_type, body.size());
    ::send(fd, header, n, 0);
    ::send(fd, body.data(), body.size(), 0);
}

static void handle_connection(int fd, const std::map<std::string, std::string (*)(const Params&)>& routes) {
    const auto t0 = std::chrono::steady_clock::now();
    std::string req;
    char buf[4096];
    // Read until the end of the request headers.
    while (req.find("\r\n\r\n") == std::string::npos) {
        const ssize_t n = ::recv(fd, buf, sizeof buf, 0);
        if (n <= 0) break;
        req.append(buf, n);
        if (req.size() > 65536) break;
    }
    // Parse the request line: "GET /path?query HTTP/1.1".
    const size_t sp1 = req.find(' ');
    const size_t sp2 = (sp1 == std::string::npos) ? std::string::npos : req.find(' ', sp1 + 1);
    if (sp1 == std::string::npos || sp2 == std::string::npos) { ::close(fd); return; }
    const std::string target = req.substr(sp1 + 1, sp2 - sp1 - 1);
    const size_t qm = target.find('?');
    const std::string path = target.substr(0, qm);
    const Params params = (qm == std::string::npos) ? Params{} : parse_query(target.substr(qm + 1));

    int code = 200;
    const char* status = "OK";
    std::string body;
    auto it = routes.find(path);
    if (it == routes.end()) {
        code = 404;
        status = "Not Found";
        body = "{\"error\":\"unknown endpoint\"}";
    } else {
        body = it->second(params);
    }
    // /metrics is Prometheus text; everything else is JSON.
    const char* ctype = (path == "/metrics") ? "text/plain; version=0.0.4" : "application/json";
    send_response(fd, code, status, body, ctype);
    ::close(fd);

    const double us =
        std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
    record_request(path, code, us);
    log_request(path, code, us, body.size());
}

int main(int argc, char** argv) {
    const int port = (argc > 1) ? std::atoi(argv[1]) : 8080;

    const std::map<std::string, std::string (*)(const Params&)> routes = {
        {"/health", [](const Params&) { return std::string("{\"status\":\"ok\"}"); }},
        {"/price", h_price},
        {"/impliedvol", h_impliedvol},
        {"/mc", h_mc},
        {"/submit", h_submit},
        {"/job", h_job},
        {"/metrics", h_metrics},
    };

    const int srv = ::socket(AF_INET, SOCK_STREAM, 0);
    if (srv < 0) { std::perror("socket"); return 1; }
    int yes = 1;
    ::setsockopt(srv, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof yes);
    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(static_cast<uint16_t>(port));
    if (::bind(srv, reinterpret_cast<sockaddr*>(&addr), sizeof addr) < 0) {
        std::perror("bind");
        return 1;
    }
    ::listen(srv, 16);
    std::printf("pricer server listening on http://127.0.0.1:%d\n", port);
    std::fflush(stdout);

    for (;;) {
        const int c = ::accept(srv, nullptr, nullptr);
        if (c < 0) continue;
        std::thread(handle_connection, c, std::cref(routes)).detach();
    }
}
