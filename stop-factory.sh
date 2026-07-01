#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════════
#  HCLS AI Factory — Stop All Services
# ═══════════════════════════════════════════════════════════════════════════════

echo ""
echo "Stopping HCLS AI Factory services..."
echo ""

# Stop Python/Streamlit services by port
for port in 8080 5000 5001 8501 8505 8510; do
    pid=$(ss -tlnp "sport = :$port" 2>/dev/null | grep -oP 'pid=\K[0-9]+' | head -1)
    if [ -n "$pid" ]; then
        name=$(ps -p "$pid" -o args= 2>/dev/null | head -c 60)
        kill "$pid" 2>/dev/null && echo "  Stopped :$port (PID $pid) $name"
    fi
done

echo ""
echo "Note: Docker services (Milvus, Grafana, Prometheus, Ollama) left running."
echo "To stop those too: docker stop milvus-standalone grafana prometheus node-exporter dcgm-exporter ollama-compose attu"
echo ""
