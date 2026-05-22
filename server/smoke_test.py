#!/usr/bin/env python3
"""Self-contained smoke test for the pricer REST server.

Spawns the server binary, waits for it to come up, exercises the endpoints, and
exits non-zero on any failure. Used by CI and runnable locally:

    python server/smoke_test.py build-srv/server/pricer_server [port]
"""
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request


def get(base, path):
    with urllib.request.urlopen(base + path, timeout=5) as r:
        return json.loads(r.read().decode())


def main():
    if len(sys.argv) < 2:
        print("usage: smoke_test.py <server_binary> [port]")
        return 1
    binary = sys.argv[1]
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 8139
    base = f"http://127.0.0.1:{port}"

    proc = subprocess.Popen([binary, str(port)])
    try:
        # Wait for /health (up to ~10s).
        for _ in range(100):
            try:
                if get(base, "/health").get("status") == "ok":
                    break
            except Exception:
                time.sleep(0.1)
        else:
            print("FAIL: server did not become healthy")
            return 1

        ok = True

        price = get(base, "/price?type=call&S=100&K=100&r=0.05&sigma=0.2&T=1")
        if abs(price["price"] - 10.450584) > 1e-3:
            print("FAIL price:", price); ok = False
        else:
            print("ok   price:", price["price"])

        iv = get(base, "/impliedvol?type=call&price=10.450584&S=100&K=100&r=0.05&T=1")
        if abs(iv["implied_vol"] - 0.20) > 1e-3:
            print("FAIL iv:", iv); ok = False
        else:
            print("ok   implied_vol:", iv["implied_vol"])

        mc = get(base, "/mc?type=call&S=100&K=100&r=0.05&sigma=0.2&T=1&paths=2000000")
        if abs(mc["mc_price"] - 10.450584) > 0.1:
            print("FAIL mc:", mc); ok = False
        else:
            print("ok   mc_price:", mc["mc_price"])

        # Async job API.
        jid = get(base, "/submit?type=call&S=100&K=100&r=0.05&sigma=0.2&T=1&paths=5000000")["job_id"]
        result = None
        for _ in range(200):
            st = get(base, "/job?id=%d" % jid)
            if st["status"] == "done":
                result = st["mc_price"]; break
            time.sleep(0.05)
        if result is None or abs(result - 10.450584) > 0.1:
            print("FAIL job:", st); ok = False
        else:
            print("ok   job mc_price:", result)

        # Unknown endpoint must be a 404.
        try:
            get(base, "/nope")
            print("FAIL: /nope did not 404"); ok = False
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print("ok   404 on unknown endpoint")
            else:
                print("FAIL 404:", e.code); ok = False

        print("ALL PASSED" if ok else "SMOKE FAILED")
        return 0 if ok else 1
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            proc.kill()


if __name__ == "__main__":
    sys.exit(main())
