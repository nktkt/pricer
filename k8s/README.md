# Kubernetes deployment

Manifests to run the pricer REST service on Kubernetes:

- `deployment.yaml` — replicas of the stateless server (non-root, read-only root
  filesystem, dropped capabilities, `/health` liveness/readiness probes, and
  `prometheus.io/*` scrape annotations for the `/metrics` endpoint).
- `service.yaml` — a `ClusterIP` service exposing port 80 → container 8080.
- `hpa.yaml` — CPU-based autoscaling (2–10 replicas; Monte Carlo pricing is
  CPU-bound). Requires `metrics-server` in the cluster.

## Deploy

Build the image and make it available to your cluster, then apply the manifests:

```sh
docker build -t pricer-server:latest .
# kind:  kind load docker-image pricer-server:latest
# or push to a registry and edit `image:` in deployment.yaml to <registry>/pricer-server:<tag>

kubectl apply -f k8s/
kubectl rollout status deploy/pricer-server
```

## Use

```sh
kubectl port-forward svc/pricer-server 8080:80
curl 'http://127.0.0.1:8080/price?type=call&S=100&K=100&r=0.05&sigma=0.2&T=1'
curl http://127.0.0.1:8080/metrics
```

The HPA needs `metrics-server` (so `kubectl top pods` works). The
`prometheus.io/*` pod annotations let a Prometheus scraper discover `/metrics`.

> These manifests are validated offline (YAML well-formedness + structure). Apply
> them to a real cluster (e.g. kind / minikube) for full API-schema and runtime
> validation.
