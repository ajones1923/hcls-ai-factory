/**
 * RAG Chat Pipeline Portal JavaScript
 */

// Global state
let eventSource = null;
let gpuMemoryChart = null;
let gpuMemoryData = [];
const MAX_DATA_POINTS = 60;

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    initializeEventSource();
    updateStatus();
    loadConfig();
    loadCurrentModel();
    initCharts();

    // Set up periodic updates
    setInterval(updateStatus, 5000);
    setInterval(updateLlmMetrics, 3000);

    // Set up event listeners
    setupEventListeners();
});

function setupEventListeners() {
    // Workflow step clicks
    document.querySelectorAll('#workflow-steps .list-group-item').forEach(item => {
        item.addEventListener('click', (e) => {
            e.preventDefault();
            const step = item.dataset.step;
            runStep(step);
        });
    });

    // Config form
    document.getElementById('config-form').addEventListener('submit', (e) => {
        e.preventDefault();
        saveConfig();
    });

    // Stop button
    document.getElementById('stop-btn').addEventListener('click', stopStep);
}

function initializeEventSource() {
    if (eventSource) {
        eventSource.close();
    }

    eventSource = new EventSource('/api/stream');

    eventSource.onmessage = (event) => {
        const data = JSON.parse(event.data);

        if (data.type === 'status') {
            updatePipelineStatus(data.data);
        } else if (data.type === 'log') {
            appendLog(data.data);
        }
    };

    eventSource.onerror = () => {
        console.log('SSE connection error, reconnecting...');
        setTimeout(initializeEventSource, 3000);
    };
}

function updatePipelineStatus(status) {
    const statusBadge = document.getElementById('pipeline-status-badge');
    const stopBtn = document.getElementById('stop-btn');
    const currentStepCard = document.getElementById('current-step-card');
    const summaryStatus = document.getElementById('summary-status');

    // Update status badge
    let badgeClass = 'bg-secondary';
    let icon = 'bi-circle-fill text-secondary';

    if (status.status === 'running') {
        badgeClass = 'bg-warning text-dark';
        icon = 'bi-circle-fill text-warning';
        stopBtn.style.display = 'block';
        currentStepCard.style.display = 'block';
    } else if (status.status === 'success') {
        badgeClass = 'bg-success';
        icon = 'bi-check-circle-fill text-success';
        stopBtn.style.display = 'none';
        setTimeout(() => { currentStepCard.style.display = 'none'; }, 3000);
    } else if (status.status === 'error') {
        badgeClass = 'bg-danger';
        icon = 'bi-x-circle-fill text-danger';
        stopBtn.style.display = 'none';
    } else {
        stopBtn.style.display = 'none';
        currentStepCard.style.display = 'none';
    }

    statusBadge.className = `badge ${badgeClass}`;
    statusBadge.innerHTML = `<i class="bi ${icon}"></i> ${capitalize(status.status)}`;
    summaryStatus.textContent = capitalize(status.status);

    // Update current step info
    if (status.current_step) {
        document.getElementById('current-step-title').innerHTML =
            `<span class="spinner-border spinner-border-sm me-2"></span> ${status.current_step}`;
        document.getElementById('step-start-time').textContent =
            status.start_time ? new Date(status.start_time).toLocaleTimeString() : '-';

        const stepStatusText = document.getElementById('step-status-text');
        stepStatusText.textContent = capitalize(status.status);
        stepStatusText.className = `badge ${status.status === 'running' ? 'bg-primary' : status.status === 'success' ? 'bg-success' : 'bg-danger'}`;
    }

    // Update workflow step icons
    updateWorkflowStepIcons(status);
}

function updateWorkflowStepIcons(status) {
    document.querySelectorAll('#workflow-steps .list-group-item').forEach(item => {
        const step = item.dataset.step;
        const icon = item.querySelector('.step-icon');

        item.classList.remove('running', 'completed', 'error');

        if (status.current_step && status.current_step.toLowerCase().includes(step.replace('-', ' '))) {
            if (status.status === 'running') {
                item.classList.add('running');
                icon.className = 'bi bi-arrow-repeat step-icon text-warning';
            } else if (status.status === 'success') {
                item.classList.add('completed');
                icon.className = 'bi bi-check-circle-fill step-icon text-success';
            } else if (status.status === 'error') {
                item.classList.add('error');
                icon.className = 'bi bi-x-circle-fill step-icon text-danger';
            }
        }
    });
}

function appendLog(logEntry) {
    const console = document.getElementById('console-output');
    const timestamp = new Date(logEntry.timestamp).toLocaleTimeString();
    const message = logEntry.message;

    // Color code based on content
    let color = '#d4d4d4';
    if (message.includes('ERROR') || message.includes('error')) {
        color = '#f44336';
    } else if (message.includes('SUCCESS') || message.includes('complete') || message.includes('Complete')) {
        color = '#4caf50';
    } else if (message.includes('WARNING') || message.includes('warning')) {
        color = '#ff9800';
    } else if (message.includes('INFO')) {
        color = '#2196f3';
    }

    console.innerHTML += `<span style="color: #888;">[${timestamp}]</span> <span style="color: ${color};">${escapeHtml(message)}</span>\n`;
    console.scrollTop = console.scrollHeight;
}

function runStep(step) {
    // Clear console
    document.getElementById('console-output').textContent = 'Starting step...\n';

    fetch(`/api/run/${step}`)
        .then(response => response.json())
        .then(data => {
            if (data.error) {
                alert('Error: ' + data.error);
            } else {
                console.log('Step started:', data.step);
            }
        })
        .catch(error => {
            alert('Failed to start step: ' + error);
        });
}

function stopStep() {
    fetch('/api/stop')
        .then(response => response.json())
        .then(data => {
            console.log('Step stopped');
        })
        .catch(error => {
            alert('Failed to stop: ' + error);
        });
}

function updateStatus() {
    fetch('/api/status')
        .then(response => response.json())
        .then(data => {
            updateServices(data.services);
            updateSystemStatus(data.system);
            updateCollectionStatus(data.collection);
            updateDataStatus(data.data);
            updateConfigModel(data);
        })
        .catch(error => {
            console.error('Failed to get status:', error);
        });
}

function updateServices(services) {
    ['milvus', 'ollama'].forEach(service => {
        const el = document.getElementById(`service-${service}`);
        if (el) {
            if (services[service]) {
                el.className = 'badge bg-success';
                el.textContent = 'Online';
            } else {
                el.className = 'badge bg-danger';
                el.textContent = 'Offline';
            }
        }
    });
}

function updateSystemStatus(system) {
    // GPU
    if (system.gpu_utilization && system.gpu_utilization.available && system.gpu_utilization.devices.length > 0) {
        const gpu = system.gpu_utilization.devices[0];
        document.getElementById('gpu-percent').textContent = `${gpu.utilization}%`;
        document.getElementById('gpu-progress').style.width = `${gpu.utilization}%`;
        document.getElementById('gpu-details').textContent =
            `${gpu.name} - ${gpu.memory_used_gb}/${gpu.memory_total_gb} GB, ${gpu.temperature}C`;

        // Update chart
        updateGpuChart(gpu.memory_percent);
    }

    // Memory
    if (system.memory_utilization) {
        const mem = system.memory_utilization;
        document.getElementById('memory-percent').textContent = `${mem.percent}%`;
        document.getElementById('memory-progress').style.width = `${mem.percent}%`;
        document.getElementById('memory-details').textContent =
            `${mem.used_gb}/${mem.total_gb} GB used`;
    }

    // Disk
    if (system.disk_space) {
        document.getElementById('disk-status').textContent =
            `${system.disk_space.available} free (${system.disk_space.use_percent} used)`;
    }
}

function updateCollectionStatus(collection) {
    if (collection.exists) {
        document.getElementById('collection-name').textContent = collection.name;
        document.getElementById('collection-name').className = 'badge bg-success';
        document.getElementById('collection-count').textContent = collection.num_entities.toLocaleString();
        document.getElementById('summary-variants').textContent = collection.num_entities.toLocaleString();
    } else {
        document.getElementById('collection-name').textContent = 'Not created';
        document.getElementById('collection-name').className = 'badge bg-secondary';
        document.getElementById('collection-count').textContent = '0';
        document.getElementById('summary-variants').textContent = '0';
    }
}

function updateDataStatus(data) {
    if (data.vcf_exists) {
        document.getElementById('vcf-status').textContent = data.vcf_size;
        document.getElementById('vcf-status').className = 'badge bg-success';
    } else {
        document.getElementById('vcf-status').textContent = 'Not found';
        document.getElementById('vcf-status').className = 'badge bg-danger';
    }
}

function updateConfigModel(data) {
    // Update model display in header
    const model = document.getElementById('summary-model');
    // Try to get from config or show default
    model.textContent = 'Llama 70B';
}

function loadVcfPreview() {
    const tbody = document.getElementById('vcf-table-body');
    tbody.innerHTML = '<tr><td colspan="6" class="text-center"><span class="spinner-border spinner-border-sm"></span> Loading...</td></tr>';

    fetch('/api/vcf-preview?limit=100')
        .then(response => response.json())
        .then(data => {
            if (data.error) {
                tbody.innerHTML = `<tr><td colspan="6" class="text-center text-danger">${data.error}</td></tr>`;
                return;
            }

            document.getElementById('vcf-path-display').textContent = `VCF Path: ${data.path} (${data.count} variants shown)`;

            if (data.variants.length === 0) {
                tbody.innerHTML = '<tr><td colspan="6" class="text-center text-muted">No variants found</td></tr>';
                return;
            }

            tbody.innerHTML = data.variants.map(v => `
                <tr>
                    <td><span class="badge bg-secondary">${v.chrom}</span></td>
                    <td>${v.pos.toLocaleString()}</td>
                    <td><code>${v.ref}</code></td>
                    <td><code>${v.alt}</code></td>
                    <td>${v.qual}</td>
                    <td><span class="badge ${v.filter === 'PASS' ? 'bg-success' : 'bg-warning'}">${v.filter}</span></td>
                </tr>
            `).join('');
        })
        .catch(error => {
            tbody.innerHTML = `<tr><td colspan="6" class="text-center text-danger">Error: ${error}</td></tr>`;
        });
}

function updateLlmMetrics() {
    fetch('/api/llm-metrics')
        .then(response => response.json())
        .then(data => {
            // Update TTFT
            const ttftEl = document.getElementById('metric-ttft');
            if (ttftEl) {
                ttftEl.textContent = data.ttft !== null ? data.ttft.toFixed(0) : '--';
            }

            // Update tokens/sec
            const tpsEl = document.getElementById('metric-tps');
            if (tpsEl) {
                tpsEl.textContent = data.tokens_per_sec !== null ? data.tokens_per_sec.toFixed(1) : '--';
            }

            // Update total queries
            const queriesEl = document.getElementById('metric-queries');
            if (queriesEl) {
                queriesEl.textContent = data.total_queries !== undefined ? data.total_queries : '0';
            }

            // Update avg latency
            const latencyEl = document.getElementById('metric-latency');
            if (latencyEl) {
                latencyEl.textContent = data.avg_latency !== null ? `${data.avg_latency} ms` : '--';
            }

            // Update cache hits
            const cacheHitsEl = document.getElementById('metric-cache-hits');
            if (cacheHitsEl) {
                cacheHitsEl.textContent = data.cache_hits !== undefined ? data.cache_hits : '0';
            }

            // Update cache misses
            const cacheMissesEl = document.getElementById('metric-cache-misses');
            if (cacheMissesEl) {
                cacheMissesEl.textContent = data.cache_misses !== undefined ? data.cache_misses : '0';
            }

            // Calculate and update cache hit rate
            const hits = data.cache_hits || 0;
            const misses = data.cache_misses || 0;
            const total = hits + misses;
            const hitRate = total > 0 ? (hits / total * 100) : 0;

            const hitRateEl = document.getElementById('cache-hit-rate');
            if (hitRateEl) {
                hitRateEl.textContent = `${hitRate.toFixed(1)}%`;
            }

            const hitProgressEl = document.getElementById('cache-hit-progress');
            if (hitProgressEl) {
                hitProgressEl.style.width = `${hitRate}%`;
            }

            // Update KV Cache metrics
            const kvUsageEl = document.getElementById('metric-kv-usage');
            if (kvUsageEl) {
                if (data.kv_cache_usage !== undefined && data.kv_cache_usage !== null) {
                    kvUsageEl.textContent = `${data.kv_cache_usage.toFixed(1)}%`;
                } else {
                    kvUsageEl.textContent = '--';
                }
            }

            const prefixHitsEl = document.getElementById('metric-prefix-hits');
            if (prefixHitsEl) {
                prefixHitsEl.textContent = data.prefix_cache_hits !== undefined ? data.prefix_cache_hits : '--';
            }

            // Update Ollama status
            const ollamaStatusEl = document.getElementById('metric-ollama-status');
            if (ollamaStatusEl) {
                if (data.ollama && data.ollama.available) {
                    ollamaStatusEl.textContent = `${data.ollama.models_count} models`;
                    ollamaStatusEl.className = 'badge bg-success';
                } else {
                    ollamaStatusEl.textContent = 'Offline';
                    ollamaStatusEl.className = 'badge bg-danger';
                }
            }
        })
        .catch(error => {
            console.error('Failed to get LLM metrics:', error);
        });
}

function loadConfig() {
    fetch('/api/config')
        .then(response => response.json())
        .then(config => {
            const form = document.getElementById('config-form');
            for (const [key, value] of Object.entries(config)) {
                const input = form.querySelector(`[name="${key}"]`);
                if (input) {
                    if (input.type === 'select-one') {
                        input.value = value;
                    } else {
                        input.value = value;
                    }
                }
            }
        })
        .catch(error => {
            console.error('Failed to load config:', error);
        });
}

function saveConfig() {
    const form = document.getElementById('config-form');
    const formData = new FormData(form);
    const config = {};

    formData.forEach((value, key) => {
        config[key] = value;
    });

    fetch('/api/config', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify(config),
    })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                alert('Configuration saved!');
            } else {
                alert('Failed to save configuration');
            }
        })
        .catch(error => {
            alert('Error saving configuration: ' + error);
        });
}

function initCharts() {
    const ctx = document.getElementById('gpuMemoryChart');
    if (!ctx) return;

    gpuMemoryChart = new Chart(ctx.getContext('2d'), {
        type: 'line',
        data: {
            labels: [],
            datasets: [{
                label: 'GPU Memory %',
                data: [],
                borderColor: '#0d6efd',
                backgroundColor: 'rgba(13, 110, 253, 0.1)',
                fill: true,
                tension: 0.3,
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 100,
                    title: {
                        display: true,
                        text: 'Memory %'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Time'
                    }
                }
            },
            plugins: {
                legend: {
                    display: false
                }
            }
        }
    });
}

function updateGpuChart(memPercent) {
    if (!gpuMemoryChart) return;

    const now = new Date().toLocaleTimeString();

    gpuMemoryData.push({
        time: now,
        value: memPercent
    });

    if (gpuMemoryData.length > MAX_DATA_POINTS) {
        gpuMemoryData.shift();
    }

    gpuMemoryChart.data.labels = gpuMemoryData.map(d => d.time);
    gpuMemoryChart.data.datasets[0].data = gpuMemoryData.map(d => d.value);
    gpuMemoryChart.update('none');
}

// Utility functions
function capitalize(str) {
    if (!str) return '';
    return str.charAt(0).toUpperCase() + str.slice(1);
}

function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

function loadTargets() {
    const tbody = document.getElementById('targets-table-body');
    tbody.innerHTML = '<tr><td colspan="7" class="text-center"><span class="spinner-border spinner-border-sm"></span> Loading...</td></tr>';

    fetch('/api/targets')
        .then(response => response.json())
        .then(data => {
            if (!data.success) {
                tbody.innerHTML = `<tr><td colspan="7" class="text-center text-danger">${data.error}</td></tr>`;
                return;
            }

            // Update summary
            document.getElementById('targets-total').textContent = data.summary.total || 0;
            document.getElementById('targets-high-priority').textContent = data.summary.high_priority || 0;
            document.getElementById('targets-validated').textContent = data.summary.by_status?.validated || 0;
            document.getElementById('targets-selected').textContent = data.summary.by_status?.selected || 0;

            if (!data.targets || data.targets.length === 0) {
                tbody.innerHTML = '<tr><td colspan="7" class="text-center text-muted">No targets saved yet. Add targets via the Streamlit Chat UI.</td></tr>';
                return;
            }

            // Sort by priority
            const sorted = data.targets.sort((a, b) => b.priority - a.priority);

            tbody.innerHTML = sorted.map(t => {
                const priorityBadge = t.priority >= 4 ? 'bg-danger' : t.priority >= 2 ? 'bg-warning text-dark' : 'bg-secondary';
                const statusBadge = t.status === 'selected' ? 'bg-success' : t.status === 'validated' ? 'bg-info' : t.status === 'rejected' ? 'bg-danger' : 'bg-secondary';

                return `
                    <tr>
                        <td><span class="badge ${priorityBadge}">${t.priority}/5</span></td>
                        <td><strong>${t.gene}</strong></td>
                        <td>${t.protein || '--'}</td>
                        <td>${t.confidence}</td>
                        <td><span class="badge ${statusBadge}">${t.status}</span></td>
                        <td>${t.variant_count || 0}</td>
                        <td>${t.pdb_ids && t.pdb_ids.length > 0 ? t.pdb_ids.join(', ') : '--'}</td>
                    </tr>
                `;
            }).join('');
        })
        .catch(error => {
            tbody.innerHTML = `<tr><td colspan="7" class="text-center text-danger">Error: ${error}</td></tr>`;
        });
}

function exportTargets() {
    fetch('/api/targets/export')
        .then(response => response.json())
        .then(data => {
            if (!data.success) {
                alert('Export failed: ' + data.error);
                return;
            }

            // Create download
            const blob = new Blob([JSON.stringify(data.data, null, 2)], { type: 'application/json' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'targets_for_phase5.json';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);

            alert(`Exported ${data.data.count} targets for Phase 5 (Cryo-EM)`);
        })
        .catch(error => {
            alert('Export failed: ' + error);
        });
}

function getModelDisplayName(model) {
    const names = {
        'llama3.1:8b-fast-ctx2k': 'Llama 8B Fast',
        'llama3.1:8b-instruct-q4_0': 'Llama 8B Q4_0',
        'llama3.1:8b': 'Llama 8B',
        'llama3.1:70b': 'Llama 70B',
        'claude-sonnet-4-20250514': 'Claude Sonnet',
        'claude-opus-4-20250514': 'Claude Opus',
    };
    return names[model] || model;
}

function loadCurrentModel() {
    fetch('/api/model')
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                const selector = document.getElementById('model-selector');
                if (selector) {
                    selector.value = data.model;
                    updateModelInfo(data.model);
                }
                // Update summary display in header
                const summaryModel = document.getElementById('summary-model');
                if (summaryModel) {
                    summaryModel.textContent = getModelDisplayName(data.model);
                }
            }
        })
        .catch(error => {
            console.error('Failed to load current model:', error);
        });
}

function switchModel(model) {
    const statusDiv = document.getElementById('model-switch-status');
    const selector = document.getElementById('model-selector');

    // Show loading state
    if (statusDiv) statusDiv.style.display = 'block';
    if (selector) selector.disabled = true;

    fetch('/api/model', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ model: model }),
    })
        .then(response => response.json())
        .then(data => {
            if (data.success) {
                updateModelInfo(model);
                // Update summary display in header
                const summaryModel = document.getElementById('summary-model');
                if (summaryModel) {
                    summaryModel.textContent = getModelDisplayName(model);
                }
                // Show success message briefly
                if (statusDiv) {
                    statusDiv.innerHTML = '<i class="bi bi-check-circle text-success"></i> <small class="text-success">Model switched!</small>';
                    setTimeout(() => {
                        statusDiv.style.display = 'none';
                    }, 2000);
                }
            } else {
                alert('Failed to switch model: ' + data.error);
                // Revert selector
                loadCurrentModel();
            }
        })
        .catch(error => {
            alert('Error switching model: ' + error);
            loadCurrentModel();
        })
        .finally(() => {
            if (selector) selector.disabled = false;
        });
}

function updateModelInfo(model) {
    const infoText = document.getElementById('model-status-text');
    if (infoText) {
        const modelInfo = {
            'llama3.1:8b-fast-ctx2k': 'Fast: Q4_0, 2K context, optimized',
            'llama3.1:8b-instruct-q4_0': 'Q4_0: smaller quantization',
            'llama3.1:8b': '8B: standard Q4_K_M',
            'llama3.1:70b': '70B: ~30-60s, best local quality',
            'claude-sonnet-4-20250514': 'Sonnet: ~2-5s, fast cloud',
            'claude-opus-4-20250514': 'Opus: ~5-10s, best quality',
        };
        infoText.textContent = modelInfo[model] || 'Unknown model';
    }
}

function getProviderForModel(model) {
    if (model.startsWith('claude-')) {
        return 'anthropic';
    } else if (model.startsWith('gpt-')) {
        return 'openai';
    } else {
        return 'ollama';
    }
}

// Make runStep global for onclick handlers
window.runStep = runStep;
window.stopStep = stopStep;
window.loadVcfPreview = loadVcfPreview;
window.loadTargets = loadTargets;
window.exportTargets = exportTargets;
window.switchModel = switchModel;
