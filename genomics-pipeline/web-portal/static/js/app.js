// Genomics Pipeline Portal - Main JavaScript

let eventSource = null;
let currentConfig = {};
let throughputChart = null;
let throughputData = [];
let runtimeInterval = null;
let pipelineStartTime = null;
let completedStages = new Set();
let currentSubStep = 0;
let lastBasesPerMin = 0;

// Initialize on page load
document.addEventListener('DOMContentLoaded', function() {
    initializeApp();
    setupEventListeners();
    startStatusPolling();
    loadConfiguration();
    initializeThroughputChart();
    startThroughputSimulation();
});

function initializeApp() {
    console.log('Genomics Pipeline Portal initialized');
    // Reset bases metric on page load
    document.getElementById('summary-bases-per-min').textContent = '0.00/min';
    updateStatus();
    // Load any saved completed state
    loadSavedState();
    // Load file list
    loadFileList();
}

async function loadSavedState() {
    try {
        const response = await fetch('/api/metrics');
        const metrics = await response.json();

        // If there's a completed or error state with runtime, display it
        if (metrics.status === 'success' || metrics.status === 'error') {
            const summaryStatus = document.getElementById('summary-status');
            const summaryRuntime = document.getElementById('summary-runtime');

            if (metrics.status === 'success') {
                summaryStatus.textContent = 'Complete';
                summaryStatus.style.color = '#48bb78';  // Green
            } else if (metrics.status === 'error') {
                summaryStatus.textContent = 'Error';
                summaryStatus.style.color = '#dc3545';  // Red
            }

            // Display the saved runtime
            if (metrics.runtime_seconds > 0) {
                const hours = Math.floor(metrics.runtime_seconds / 3600);
                const minutes = Math.floor((metrics.runtime_seconds % 3600) / 60);
                const seconds = metrics.runtime_seconds % 60;
                const formatted = [
                    hours.toString().padStart(2, '0'),
                    minutes.toString().padStart(2, '0'),
                    seconds.toString().padStart(2, '0')
                ].join(':');
                summaryRuntime.textContent = formatted;
            }

            // Update stage display
            const summaryStage = document.getElementById('summary-stage');
            if (metrics.current_step) {
                summaryStage.textContent = metrics.current_step + (metrics.status === 'success' ? ' ✓' : ' ✗');
            }
        }
    } catch (error) {
        console.error('Error loading saved state:', error);
    }
}

function setupEventListeners() {
    // Workflow step clicks
    document.querySelectorAll('#workflow-steps .list-group-item').forEach(item => {
        item.addEventListener('click', function(e) {
            e.preventDefault();
            const step = this.dataset.step;
            runStep(step);
        });
    });

    // Stop button (in current step card)
    document.getElementById('stop-btn').addEventListener('click', function() {
        stopProcess();
    });

    // Stop All button (in top navbar)
    document.getElementById('stop-all-btn').addEventListener('click', function() {
        stopAllProcesses();
    });

    // Configuration form
    document.getElementById('config-form').addEventListener('submit', function(e) {
        e.preventDefault();
        saveConfiguration();
    });

    // Log file selector
    document.getElementById('log-select').addEventListener('change', function() {
        const logType = this.value;
        if (logType) {
            loadLogFile(logType);
        }
    });
}

function startStatusPolling() {
    // Poll status every 3 seconds
    updateStatus();
    setInterval(updateStatus, 3000);

    // Start SSE connection for real-time logs
    connectEventSource();
}

function connectEventSource() {
    if (eventSource) {
        eventSource.close();
    }

    eventSource = new EventSource('/api/stream');

    eventSource.onmessage = function(event) {
        const data = JSON.parse(event.data);

        if (data.type === 'status') {
            updatePipelineStatus(data.data);
        } else if (data.type === 'log') {
            console.log('Received log:', data.data.message);  // Debug
            appendLog(data.data);
        }
    };

    eventSource.onerror = function(error) {
        console.error('SSE Error:', error);
        // Reconnect after 5 seconds
        setTimeout(connectEventSource, 5000);
    };
}

async function updateStatus() {
    try {
        const response = await fetch('/api/status');
        const data = await response.json();

        // Update system status
        updateSystemStatus(data.system);

        // Update data files status
        updateDataStatus(data.data, data.file_sizes);

    } catch (error) {
        console.error('Error fetching status:', error);
    }
}

function updateSystemStatus(system) {
    // Docker status
    const dockerBadge = document.getElementById('docker-status');
    if (system.docker_installed) {
        dockerBadge.className = 'badge bg-success float-end';
        dockerBadge.textContent = 'Available';
    } else {
        dockerBadge.className = 'badge bg-danger float-end';
        dockerBadge.textContent = 'Not Found';
    }

    // GPU status
    const gpuBadge = document.getElementById('gpu-status');
    if (system.gpu_available) {
        gpuBadge.className = 'badge bg-success float-end';
        gpuBadge.textContent = 'Available';
    } else {
        gpuBadge.className = 'badge bg-danger float-end';
        gpuBadge.textContent = 'Not Found';
    }

    // Disk space
    const diskBadge = document.getElementById('disk-status');
    const diskSpace = system.disk_space;
    diskBadge.textContent = `${diskSpace.available} free`;

    const usePercent = parseInt(diskSpace.use_percent);
    if (usePercent > 90) {
        diskBadge.className = 'badge bg-danger float-end';
    } else if (usePercent > 75) {
        diskBadge.className = 'badge bg-warning float-end';
    } else {
        diskBadge.className = 'badge bg-success float-end';
    }

    // CPU Utilization
    if (system.cpu_utilization) {
        const cpu = system.cpu_utilization;
        updateUtilizationBar('cpu', cpu.percent, cpu.count, cpu.frequency);
    }

    // Memory Utilization
    if (system.memory_utilization) {
        const mem = system.memory_utilization;
        updateMemoryBar(mem);
    }

    // GPU Utilization
    if (system.gpu_utilization) {
        updateGPUBars(system.gpu_utilization);
    }
}

function updateUtilizationBar(type, percent, count, frequency) {
    const percentElem = document.getElementById(`${type}-percent`);
    const progressBar = document.getElementById(`${type}-progress`);
    const detailsElem = document.getElementById(`${type}-details`);

    percentElem.textContent = `${percent}%`;
    progressBar.style.width = `${percent}%`;

    // Color based on utilization
    progressBar.className = 'progress-bar';
    if (percent > 90) {
        progressBar.classList.add('bg-danger');
    } else if (percent > 75) {
        progressBar.classList.add('bg-warning');
    } else if (percent > 50) {
        progressBar.classList.add('bg-info');
    } else {
        progressBar.classList.add('bg-success');
    }

    // Add details
    if (count && frequency) {
        detailsElem.textContent = `${count} cores @ ${frequency} MHz`;
    }
}

function updateMemoryBar(mem) {
    const percentElem = document.getElementById('memory-percent');
    const progressBar = document.getElementById('memory-progress');
    const detailsElem = document.getElementById('memory-details');

    percentElem.textContent = `${mem.percent}%`;
    progressBar.style.width = `${mem.percent}%`;

    // Color based on utilization
    progressBar.className = 'progress-bar';
    if (mem.percent > 90) {
        progressBar.classList.add('bg-danger');
    } else if (mem.percent > 75) {
        progressBar.classList.add('bg-warning');
    } else if (mem.percent > 50) {
        progressBar.classList.add('bg-info');
    } else {
        progressBar.classList.add('bg-success');
    }

    // Add details
    detailsElem.textContent = `${mem.used_gb} / ${mem.total_gb} GB used`;
}

function updateGPUBars(gpuData) {
    const container = document.getElementById('gpu-utilization-container');

    if (!gpuData.available || gpuData.devices.length === 0) {
        container.innerHTML = '';
        return;
    }

    let html = '';
    gpuData.devices.forEach((gpu, index) => {
        html += `
            <div class="mb-3">
                <div class="d-flex justify-content-between align-items-center mb-1">
                    <small class="text-muted"><i class="bi bi-gpu-card"></i> ${gpu.name}</small>
                    <small class="fw-bold">${gpu.utilization}%</small>
                </div>
                <div class="progress" style="height: 8px;">
                    <div class="progress-bar ${getUtilizationColor(gpu.utilization)}"
                         role="progressbar"
                         style="width: ${gpu.utilization}%"></div>
                </div>
                <small class="text-muted" style="font-size: 0.7rem;">
                    ${gpu.temperature > 0 ? `${gpu.temperature}°C` : ''}
                    ${gpu.temperature > 0 && gpu.power_watts > 0 ? ' | ' : ''}${gpu.power_watts > 0 ? `${gpu.power_watts}W` : ''}
                </small>
            </div>
        `;
    });

    container.innerHTML = html;
}

function getUtilizationColor(percent) {
    if (percent > 90) return 'bg-danger';
    if (percent > 75) return 'bg-warning';
    if (percent > 50) return 'bg-info';
    return 'bg-success';
}

function updateDataStatus(data, fileSizes) {
    updateFileBadge('file-fastq-r1', data.fastq_r1, fileSizes.fastq_r1);
    updateFileBadge('file-fastq-r2', data.fastq_r2, fileSizes.fastq_r2);
    updateFileBadge('file-reference', data.reference, fileSizes.reference);
    updateFileBadge('file-chr20', data.chr20_vcf, fileSizes.chr20_vcf);
    updateFileBadge('file-genome', data.genome_vcf, fileSizes.genome_vcf);
}

function updateFileBadge(elementId, exists, size) {
    const badge = document.getElementById(elementId);
    if (exists) {
        badge.className = 'badge bg-success float-end';
        badge.textContent = size;
    } else {
        badge.className = 'badge bg-secondary float-end';
        badge.textContent = 'Not Found';
    }
}

function updatePipelineStatus(status) {
    const statusBadge = document.getElementById('pipeline-status-badge');
    const currentStepCard = document.getElementById('current-step-card');
    const currentStepTitle = document.getElementById('current-step-title');
    const stepStartTime = document.getElementById('step-start-time');
    const stepStatusText = document.getElementById('step-status-text');

    // Map step names to data-step attributes
    const stepNameToAttr = {
        'Prerequisites Check': 'check',
        'NGC Login': 'login',
        'Data Download': 'download',
        'Reference Setup': 'reference',
        'Chr20 Test': 'test',
        'Full Genome': 'full'
    };

    // Helper function to highlight a step
    function highlightStep(stepName, color) {
        const stepAttr = stepNameToAttr[stepName];
        if (stepAttr) {
            const stepItem = document.querySelector(`#workflow-steps .list-group-item[data-step="${stepAttr}"]`);
            if (stepItem) {
                stepItem.style.backgroundColor = color;
            }
        }
    }

    // Only reset highlights when idle (no active step)
    if (status.status === 'idle' || !status.current_step) {
        document.querySelectorAll('#workflow-steps .list-group-item').forEach(item => {
            item.style.backgroundColor = '';
        });
    }

    // Get stop button references
    const stopBtn = document.getElementById('stop-btn');
    const stopAllBtn = document.getElementById('stop-all-btn');

    // Update main status badge
    if (status.status === 'running') {
        statusBadge.innerHTML = '<span class="badge bg-warning">Running</span>';
        currentStepCard.style.display = 'block';
        currentStepTitle.innerHTML = `<span class="spinner-border spinner-border-sm me-2"></span>${status.current_step}`;
        stepStartTime.textContent = new Date(status.start_time).toLocaleTimeString();
        stepStatusText.className = 'badge bg-warning';
        stepStatusText.textContent = 'Running';
        // Show stop buttons while running
        if (stopBtn) stopBtn.style.display = 'inline-block';
        if (stopAllBtn) stopAllBtn.style.display = 'inline-block';

        // Highlight the running step with light green
        highlightStep(status.current_step, '#d4edda');
    } else if (status.status === 'success') {
        // Keep step highlighted light green on success
        highlightStep(status.current_step, '#d4edda');
        statusBadge.innerHTML = '<span class="badge bg-success">Success</span>';
        if (status.current_step) {
            currentStepCard.style.display = 'block';
            currentStepTitle.innerHTML = `<i class="bi bi-check-circle me-2"></i>${status.current_step}`;
            stepStatusText.className = 'badge bg-success';
            stepStatusText.textContent = 'Completed';
            // Hide stop buttons when completed
            if (stopBtn) stopBtn.style.display = 'none';
            if (stopAllBtn) stopAllBtn.style.display = 'none';
        }
    } else if (status.status === 'error') {
        statusBadge.innerHTML = '<span class="badge bg-danger">Error</span>';
        if (status.current_step) {
            currentStepCard.style.display = 'block';
            currentStepTitle.innerHTML = `<i class="bi bi-exclamation-triangle me-2"></i>${status.current_step}`;
            stepStatusText.className = 'badge bg-danger';
            stepStatusText.textContent = 'Failed';
            // Hide stop buttons when failed (nothing to stop)
            if (stopBtn) stopBtn.style.display = 'none';
            if (stopAllBtn) stopAllBtn.style.display = 'none';
        }
    } else {
        statusBadge.innerHTML = '<span class="badge bg-secondary">Idle</span>';
        currentStepCard.style.display = 'none';
        // Hide Stop All button when idle (show individual stop buttons for next run)
        if (stopBtn) stopBtn.style.display = 'inline-block';
        if (stopAllBtn) stopAllBtn.style.display = 'none';
    }

    // Update executive dashboard
    updateExecutiveSummary(status);
}

function appendLog(logEntry) {
    const consoleOutput = document.getElementById('console-output');
    const timestamp = new Date(logEntry.timestamp).toLocaleTimeString();
    const logLine = `[${timestamp}] ${logEntry.message}\n`;

    consoleOutput.textContent += logLine;

    // Auto-scroll to bottom
    consoleOutput.scrollTop = consoleOutput.scrollHeight;

    // Parse log message for sub-step updates
    parseSubStepFromLog(logEntry.message);
}

// Sub-steps tracking functions
function showSubStepsTracker(show) {
    const tracker = document.getElementById('sub-steps-tracker');
    if (tracker) {
        tracker.style.display = show ? 'block' : 'none';
    }
}

function resetSubSteps() {
    currentSubStep = 0;
    document.querySelectorAll('.sub-step').forEach(step => {
        step.classList.remove('running', 'completed', 'error');
        step.querySelector('.sub-step-icon').innerHTML = '<i class="bi bi-hourglass"></i>';
        step.querySelector('.sub-step-status').className = 'sub-step-status badge bg-secondary';
        step.querySelector('.sub-step-status').textContent = 'Pending';
    });
}

function updateSubStep(stepNum, status) {
    const step = document.querySelector(`.sub-step[data-substep="${stepNum}"]`);
    if (!step) return;

    // Remove all status classes
    step.classList.remove('running', 'completed', 'error');

    const icon = step.querySelector('.sub-step-icon');
    const statusBadge = step.querySelector('.sub-step-status');

    if (status === 'running') {
        step.classList.add('running');
        icon.innerHTML = '<i class="bi bi-arrow-repeat"></i>';
        statusBadge.className = 'sub-step-status badge bg-warning';
        statusBadge.textContent = 'Running';
        currentSubStep = stepNum;
    } else if (status === 'completed') {
        step.classList.add('completed');
        icon.innerHTML = '<i class="bi bi-check-circle-fill"></i>';
        statusBadge.className = 'sub-step-status badge bg-success';
        statusBadge.textContent = 'Complete';
    } else if (status === 'error') {
        step.classList.add('error');
        icon.innerHTML = '<i class="bi bi-x-circle-fill"></i>';
        statusBadge.className = 'sub-step-status badge bg-danger';
        statusBadge.textContent = 'Error';
    }
}

function parseSubStepFromLog(message) {
    // Check for step start patterns (new format with box headers)
    if (message.includes('STEP 1/4') || message.includes('Create Interval BED')) {
        // Mark step 1 as running
        updateSubStep(1, 'running');
    } else if (message.includes('STEP 2/4') || message.includes('fq2bam (FASTQ')) {
        // Mark step 1 as complete, step 2 as running
        updateSubStep(1, 'completed');
        updateSubStep(2, 'running');
    } else if (message.includes('STEP 3/4') || message.includes('BAM Indexing & QC')) {
        // Mark step 2 as complete, step 3 as running
        updateSubStep(2, 'completed');
        updateSubStep(3, 'running');
    } else if (message.includes('STEP 4/4') || message.includes('DeepVariant (Variant Calling)')) {
        // Mark step 3 as complete, step 4 as running
        updateSubStep(3, 'completed');
        updateSubStep(4, 'running');
    }
    // Check for step completion patterns
    else if (message.includes('CHR20 TEST PIPELINE COMPLETE') || message.includes('FULL GENOME PIPELINE COMPLETE')) {
        // Mark all steps as complete
        updateSubStep(1, 'completed');
        updateSubStep(2, 'completed');
        updateSubStep(3, 'completed');
        updateSubStep(4, 'completed');
    } else if (message.includes('STEP 1/4 COMPLETE') || message.includes('Interval BED created')) {
        updateSubStep(1, 'completed');
    } else if (message.includes('STEP 2/4 COMPLETE') || message.includes('fq2bam completed successfully')) {
        updateSubStep(2, 'completed');
    } else if (message.includes('STEP 3/4 COMPLETE') || message.includes('BAM indexed and QC complete')) {
        updateSubStep(3, 'completed');
    } else if (message.includes('STEP 4/4 COMPLETE') || message.includes('DeepVariant completed successfully')) {
        updateSubStep(4, 'completed');
    }
    // Check for GPU activity indicators
    else if (message.includes('Starting GPU alignment') || message.includes('Parabricks BWA')) {
        updateSubStep(2, 'running');
    } else if (message.includes('Starting DeepVariant')) {
        updateSubStep(4, 'running');
    }

    // Parse Parabricks metrics (bases/GPU/minute and total bases processed)
    parseParabricksMetrics(message);
}

function parseParabricksMetrics(message) {
    // Match pattern: "202798327139 bases/GPU/minute: 1371853674.0"
    const metricsMatch = message.match(/(\d+)\s+bases\/GPU\/minute:\s+([\d.]+)/);
    if (metricsMatch) {
        const basesPerMin = parseFloat(metricsMatch[2]);

        // Convert to billions
        const billionBasesPerMin = (basesPerMin / 1e9).toFixed(2);

        // Update display
        document.getElementById('summary-bases-per-min').textContent = billionBasesPerMin + '/min';

        // Store for other uses
        lastBasesPerMin = basesPerMin;
    }

    // Parse variant counts from DeepVariant output
    // Common patterns: "Found X variants", "X total variants", "variants: X"
    const variantPatterns = [
        /Found\s+(\d+)\s+variants/i,
        /(\d+)\s+total\s+variants/i,
        /variants:\s*(\d+)/i,
        /(\d+)\s+variant\s+records/i,
        /VCF contains\s+(\d+)/i
    ];

    for (const pattern of variantPatterns) {
        const variantMatch = message.match(pattern);
        if (variantMatch) {
            const variantCount = parseInt(variantMatch[1]);
            const variantElem = document.getElementById('metric-variants');
            if (variantElem && variantCount > 0) {
                // Format with K or M suffix for large numbers
                let formatted;
                if (variantCount >= 1000000) {
                    formatted = (variantCount / 1000000).toFixed(1) + 'M';
                } else if (variantCount >= 1000) {
                    formatted = (variantCount / 1000).toFixed(1) + 'K';
                } else {
                    formatted = variantCount.toString();
                }
                variantElem.textContent = formatted;
                variantElem.classList.add('active');
            }
            break;
        }
    }
}

async function runStep(step) {
    // Confirm for long-running steps
    if (step === 'download' || step === 'full') {
        const stepNames = {
            'download': 'Data Download (~200GB, several hours)',
            'full': 'Full Genome Pipeline (~30-110 minutes)'
        };

        if (!confirm(`Are you sure you want to start: ${stepNames[step]}?`)) {
            return;
        }
    }

    // Clear console and show initial message
    const consoleOutput = document.getElementById('console-output');
    const stepNames = {
        'check': 'Prerequisites Check',
        'login': 'NGC Login',
        'download': 'Data Download',
        'reference': 'Reference Setup',
        'test': 'Chr20 Test',
        'full': 'Full Genome Pipeline'
    };

    consoleOutput.textContent = `Starting ${stepNames[step] || step}...\n`;
    consoleOutput.textContent += `Launching script...\n`;
    consoleOutput.textContent += `Waiting for output stream...\n`;

    // Show/hide sub-steps tracker based on step type
    if (step === 'test' || step === 'full') {
        resetSubSteps();
        showSubStepsTracker(true);
    } else {
        showSubStepsTracker(false);
    }

    // Reset and start runtime timer immediately
    stopRuntimeTimer();  // Stop any existing timer first
    pipelineStartTime = new Date();
    document.getElementById('summary-runtime').textContent = '00:00:00';
    startRuntimeTimer();

    // Reset bases metrics
    lastBasesPerMin = 0;
    document.getElementById('summary-bases-per-min').textContent = '0.00/min';

    // Update status display
    document.getElementById('summary-status').textContent = 'Starting...';
    document.getElementById('summary-status').style.color = '#ffc107';

    try {
        const response = await fetch(`/api/run/${step}`);
        const data = await response.json();

        if (data.error) {
            consoleOutput.textContent += `\nERROR: ${data.error}\n`;
            stopRuntimeTimer();
            document.getElementById('summary-status').textContent = 'Error';
            document.getElementById('summary-status').style.color = '#dc3545';
            alert(`Error: ${data.error}`);
        } else {
            consoleOutput.textContent += `Script started successfully. Streaming output...\n\n`;
            console.log(`Started ${step}:`, data);
        }
    } catch (error) {
        console.error('Error running step:', error);
        consoleOutput.textContent += `\nERROR: Failed to start step\n`;
        stopRuntimeTimer();
        document.getElementById('summary-status').textContent = 'Error';
        document.getElementById('summary-status').style.color = '#dc3545';
        alert('Failed to start step. Check console for details.');
    }
}

async function stopProcess() {
    if (!confirm('Are you sure you want to stop the current process?')) {
        return;
    }

    try {
        const response = await fetch('/api/stop');
        const data = await response.json();

        if (data.success) {
            document.getElementById('console-output').textContent += '\n[SYSTEM] Process stopped by user\n';
            // Stop the runtime timer
            stopRuntimeTimer();
            document.getElementById('summary-status').textContent = 'Stopped';
            document.getElementById('summary-status').style.color = '#dc3545';
        }
    } catch (error) {
        console.error('Error stopping process:', error);
    }
}

async function stopAllProcesses() {
    if (!confirm('Are you sure you want to stop ALL pipeline processes?\n\nThis will terminate all running genome pipeline scripts.')) {
        return;
    }

    try {
        const response = await fetch('/api/stop-all');
        const data = await response.json();

        if (data.success) {
            document.getElementById('console-output').textContent += '\n[SYSTEM] All pipeline processes stopped by user\n';
            if (data.killed_count > 0) {
                document.getElementById('console-output').textContent += `[SYSTEM] Terminated ${data.killed_count} process(es)\n`;
            }
            // Stop the runtime timer
            stopRuntimeTimer();
            document.getElementById('summary-status').textContent = 'Stopped';
            document.getElementById('summary-status').style.color = '#dc3545';
        } else {
            document.getElementById('console-output').textContent += `\n[SYSTEM] ${data.message || 'No processes to stop'}\n`;
        }
    } catch (error) {
        console.error('Error stopping all processes:', error);
        document.getElementById('console-output').textContent += '\n[SYSTEM] Error stopping processes\n';
    }
}

async function loadConfiguration() {
    try {
        const response = await fetch('/api/config');
        const config = await response.json();
        currentConfig = config;

        // Populate form
        const form = document.getElementById('config-form');
        Object.keys(config).forEach(key => {
            const input = form.querySelector(`[name="${key}"]`);
            if (input) {
                input.value = config[key];
            }
        });
    } catch (error) {
        console.error('Error loading configuration:', error);
    }
}

async function saveConfiguration() {
    const form = document.getElementById('config-form');
    const formData = new FormData(form);
    const config = {};

    for (let [key, value] of formData.entries()) {
        config[key] = value;
    }

    try {
        const response = await fetch('/api/config', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(config)
        });

        const data = await response.json();

        if (data.success) {
            alert('Configuration saved successfully!');
            currentConfig = config;
        } else {
            alert('Error saving configuration');
        }
    } catch (error) {
        console.error('Error saving configuration:', error);
        alert('Failed to save configuration');
    }
}

async function loadLogFile(logType) {
    const logContent = document.getElementById('log-content');
    logContent.textContent = 'Loading...';

    try {
        const response = await fetch(`/api/logs/${logType}`);
        const data = await response.json();

        if (data.error) {
            logContent.textContent = `Error: ${data.error}`;
        } else if (data.content) {
            logContent.textContent = data.content;
            // Scroll to bottom
            logContent.scrollTop = logContent.scrollHeight;
        }
    } catch (error) {
        console.error('Error loading log:', error);
        logContent.textContent = 'Error loading log file';
    }
}

// Cleanup on page unload
window.addEventListener('beforeunload', function() {
    if (eventSource) {
        eventSource.close();
    }
    if (runtimeInterval) {
        clearInterval(runtimeInterval);
    }
});


// ============================================
// PIPELINE STATE RESET
// ============================================

async function resetPipelineState() {
    if (!confirm('Reset pipeline status to Idle? This will clear the current status display.')) {
        return;
    }

    try {
        const response = await fetch('/api/reset', { method: 'POST' });
        const data = await response.json();

        if (data.success) {
            // Reset UI
            document.getElementById('summary-status').textContent = 'Idle';
            document.getElementById('summary-status').style.color = '#6c757d';
            document.getElementById('summary-runtime').textContent = '00:00:00';
            document.getElementById('summary-stage').textContent = 'Ready to begin';
            stopRuntimeTimer();
            console.log('Pipeline state reset');
        } else {
            alert('Error resetting state: ' + data.error);
        }
    } catch (error) {
        console.error('Error resetting state:', error);
        alert('Error resetting pipeline state');
    }
}


// ============================================
// FILE MANAGER FUNCTIONS
// ============================================

async function loadFileList() {
    const directory = document.getElementById('file-directory').value;
    const fileList = document.getElementById('file-list');

    fileList.innerHTML = '<div class="text-center py-2"><div class="spinner-border spinner-border-sm"></div></div>';

    try {
        const response = await fetch(`/api/files/${directory}`);
        const data = await response.json();

        if (data.error) {
            fileList.innerHTML = `<div class="text-danger p-2"><small>${data.error}</small></div>`;
            return;
        }

        if (data.files.length === 0) {
            fileList.innerHTML = '<div class="text-muted text-center py-3"><i class="bi bi-folder2"></i><br><small>No files found</small></div>';
            return;
        }

        let html = '';
        for (const file of data.files) {
            const icon = getFileIcon(file.type);
            html += `
                <div class="file-item d-flex align-items-center justify-content-between p-1 border-bottom" style="font-size: 0.8rem;">
                    <div class="d-flex align-items-center" style="overflow: hidden; flex: 1;">
                        <i class="bi ${icon} me-1 text-muted"></i>
                        <span class="text-truncate" title="${file.name}">${file.name}</span>
                    </div>
                    <div class="d-flex align-items-center gap-1">
                        <span class="badge bg-light text-dark" style="font-size: 0.65rem;">${file.size_human}</span>
                        <button class="btn btn-sm btn-outline-primary p-0 px-1" onclick="downloadFile('${directory}', '${file.name}')" title="Download">
                            <i class="bi bi-download" style="font-size: 0.7rem;"></i>
                        </button>
                        <button class="btn btn-sm btn-outline-danger p-0 px-1" onclick="deleteFile('${directory}', '${file.name}')" title="Delete">
                            <i class="bi bi-trash" style="font-size: 0.7rem;"></i>
                        </button>
                    </div>
                </div>
            `;
        }
        fileList.innerHTML = html;

    } catch (error) {
        console.error('Error loading files:', error);
        fileList.innerHTML = '<div class="text-danger p-2"><small>Error loading files</small></div>';
    }
}

function refreshFileList() {
    loadFileList();
}

function getFileIcon(type) {
    const iconMap = {
        'fastq': 'bi-file-earmark-medical',
        'compressed': 'bi-file-earmark-zip',
        'bam': 'bi-file-earmark-binary',
        'index': 'bi-file-earmark-code',
        'vcf': 'bi-file-earmark-text',
        'fasta': 'bi-file-earmark-code',
        'bed': 'bi-file-earmark-spreadsheet',
        'text': 'bi-file-earmark-text',
        'json': 'bi-file-earmark-code',
        'file': 'bi-file-earmark'
    };
    return iconMap[type] || 'bi-file-earmark';
}

function downloadFile(directory, filename) {
    // Create a temporary link and click it
    const link = document.createElement('a');
    link.href = `/api/files/download/${directory}/${encodeURIComponent(filename)}`;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

async function deleteFile(directory, filename) {
    if (!confirm(`Are you sure you want to delete "${filename}"?`)) {
        return;
    }

    try {
        const response = await fetch(`/api/files/delete/${directory}/${encodeURIComponent(filename)}`, {
            method: 'DELETE'
        });
        const data = await response.json();

        if (data.success) {
            loadFileList();
            updateStatus();  // Refresh the Data Files status
        } else {
            alert(`Error: ${data.error}`);
        }
    } catch (error) {
        console.error('Error deleting file:', error);
        alert('Error deleting file');
    }
}

async function uploadFile() {
    const fileInput = document.getElementById('file-upload-input');
    const directory = document.getElementById('file-directory').value;
    const progressDiv = document.getElementById('upload-progress');
    const progressBar = document.getElementById('upload-progress-bar');
    const statusText = document.getElementById('upload-status');

    if (!fileInput.files || fileInput.files.length === 0) {
        alert('Please select a file to upload');
        return;
    }

    const file = fileInput.files[0];
    const formData = new FormData();
    formData.append('file', file);

    // Show progress
    progressDiv.style.display = 'block';
    progressBar.style.width = '0%';
    statusText.textContent = `Uploading ${file.name}...`;

    try {
        // Use XMLHttpRequest for progress tracking
        const xhr = new XMLHttpRequest();

        xhr.upload.addEventListener('progress', (e) => {
            if (e.lengthComputable) {
                const percent = Math.round((e.loaded / e.total) * 100);
                progressBar.style.width = `${percent}%`;
                statusText.textContent = `Uploading... ${percent}%`;
            }
        });

        xhr.addEventListener('load', () => {
            if (xhr.status === 200) {
                const data = JSON.parse(xhr.responseText);
                if (data.success) {
                    statusText.textContent = `Uploaded: ${data.filename} (${data.size_human})`;
                    progressBar.classList.remove('progress-bar-animated');
                    progressBar.classList.add('bg-success');
                    fileInput.value = '';
                    loadFileList();
                    updateStatus();
                    // Hide progress after 3 seconds
                    setTimeout(() => {
                        progressDiv.style.display = 'none';
                        progressBar.classList.add('progress-bar-animated');
                        progressBar.classList.remove('bg-success');
                    }, 3000);
                } else {
                    statusText.textContent = `Error: ${data.error}`;
                    progressBar.classList.add('bg-danger');
                }
            } else {
                statusText.textContent = 'Upload failed';
                progressBar.classList.add('bg-danger');
            }
        });

        xhr.addEventListener('error', () => {
            statusText.textContent = 'Upload failed';
            progressBar.classList.add('bg-danger');
        });

        xhr.open('POST', `/api/files/upload/${directory}`);
        xhr.send(formData);

    } catch (error) {
        console.error('Error uploading file:', error);
        statusText.textContent = 'Upload error';
        progressBar.classList.add('bg-danger');
    }
}

// ============================================
// EXECUTIVE DASHBOARD FUNCTIONS
// ============================================

function initializeThroughputChart() {
    const ctx = document.getElementById('throughputChart');
    if (!ctx) return;

    // Initialize with empty data
    const labels = [];
    const now = new Date();
    for (let i = 30; i >= 0; i--) {
        const time = new Date(now - i * 2000);
        labels.push(time.toLocaleTimeString('en-US', { hour12: false, hour: '2-digit', minute: '2-digit', second: '2-digit' }));
    }

    throughputChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: labels,
            datasets: [{
                label: 'Throughput (GB/hr)',
                data: new Array(31).fill(0),
                borderColor: 'rgba(72, 187, 120, 1)',
                backgroundColor: 'rgba(72, 187, 120, 0.1)',
                borderWidth: 2,
                fill: true,
                tension: 0.4,
                pointRadius: 0,
                pointHoverRadius: 4,
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: {
                intersect: false,
                mode: 'index'
            },
            plugins: {
                legend: {
                    display: false
                },
                tooltip: {
                    backgroundColor: 'rgba(26, 26, 46, 0.9)',
                    titleColor: '#fff',
                    bodyColor: '#48bb78',
                    borderColor: '#48bb78',
                    borderWidth: 1,
                    padding: 12,
                    displayColors: false,
                    callbacks: {
                        label: function(context) {
                            return `${context.parsed.y.toFixed(1)} GB/hr`;
                        }
                    }
                }
            },
            scales: {
                x: {
                    display: true,
                    grid: {
                        color: 'rgba(0, 0, 0, 0.05)',
                    },
                    ticks: {
                        maxTicksLimit: 6,
                        color: '#6c757d',
                        font: { size: 10 }
                    }
                },
                y: {
                    display: true,
                    beginAtZero: true,
                    grid: {
                        color: 'rgba(0, 0, 0, 0.05)',
                    },
                    ticks: {
                        color: '#6c757d',
                        font: { size: 10 },
                        callback: function(value) {
                            return value + ' GB/hr';
                        }
                    }
                }
            },
            animation: {
                duration: 300
            }
        }
    });
}

function updateThroughputChart(value) {
    if (!throughputChart) return;

    const now = new Date();
    const timeLabel = now.toLocaleTimeString('en-US', { hour12: false, hour: '2-digit', minute: '2-digit', second: '2-digit' });

    // Add new data point
    throughputChart.data.labels.push(timeLabel);
    throughputChart.data.datasets[0].data.push(value);

    // Remove old data (keep last 31 points)
    if (throughputChart.data.labels.length > 31) {
        throughputChart.data.labels.shift();
        throughputChart.data.datasets[0].data.shift();
    }

    throughputChart.update('none');

    // Update current throughput display
    document.getElementById('current-throughput').textContent = `${value.toFixed(1)} GB/hr`;
    document.getElementById('summary-throughput').textContent = `${value.toFixed(1)} GB/hr`;
}

function startThroughputSimulation() {
    // Update throughput every 2 seconds using real metrics
    setInterval(async () => {
        try {
            const response = await fetch('/api/metrics');
            const metrics = await response.json();

            let throughput = metrics.throughput_gb_hr || 0;

            // Add small variance for visual interest when running
            if (throughput > 0) {
                throughput += (Math.random() - 0.5) * 5;
                throughput = Math.max(0, throughput);
            }

            updateThroughputChart(throughput);

            // Update reads/sec estimate (approximately 1M reads = 0.3GB)
            const readsPerSec = throughput > 0 ? Math.round((throughput / 3600) * 3333333) : 0;
            const readsDisplay = readsPerSec > 0 ? `${(readsPerSec / 1000000).toFixed(1)}M reads/sec` : 'Reads/sec: --';
            document.getElementById('summary-reads').textContent = readsDisplay;

            // Update GPU metrics panel
            updateGPUMetricsPanel(metrics.gpu_metrics);

            // Update genomics metrics bar
            updateGenomicsMetricsBar(throughput, readsPerSec);

        } catch (error) {
            console.error('Error fetching metrics:', error);
            // Fallback to showing 0 on error
            updateThroughputChart(0);
            updateGPUMetricsPanel(null);
            updateGenomicsMetricsBar(0, 0);
        }
    }, 2000);
}

function updateGenomicsMetricsBar(throughput, readsPerSec) {
    const setGenomicsMetric = (id, value, isActive) => {
        const elem = document.getElementById(id);
        if (elem) {
            elem.textContent = value;
            if (isActive && value !== '--') {
                elem.classList.add('active');
            } else {
                elem.classList.remove('active');
            }
        }
    };

    // Bases per minute (use lastBasesPerMin from Parabricks parsing, or estimate from throughput)
    // 1 GB = ~3.3 billion bases (for compressed FASTQ)
    let basesPerMin = lastBasesPerMin > 0 ? (lastBasesPerMin / 1e9).toFixed(2) : '--';
    if (basesPerMin === '--' && throughput > 0) {
        // Estimate: throughput GB/hr * 3.3B bases/GB / 60 min = bases/min
        basesPerMin = ((throughput * 3.3) / 60).toFixed(2);
    }
    setGenomicsMetric('metric-bases-per-min', basesPerMin, throughput > 0);

    // Reads per second (in millions)
    const readsPerSecM = readsPerSec > 0 ? (readsPerSec / 1000000).toFixed(1) : '--';
    setGenomicsMetric('metric-reads-per-sec', readsPerSecM, readsPerSec > 0);

    // Estimated coverage (30x genome = ~100GB FASTQ, so GB processed / 3.3 = coverage)
    // This accumulates over time based on data processed
    const dataProcessedElem = document.getElementById('summary-data');
    if (dataProcessedElem) {
        const dataGB = parseFloat(dataProcessedElem.textContent) || 0;
        const coverage = dataGB > 0 ? (dataGB / 3.3).toFixed(1) + 'x' : '--';
        setGenomicsMetric('metric-coverage', coverage, dataGB > 0);
    }

    // Variants - this would be updated from log parsing when DeepVariant runs
    // For now, show -- unless we have variant data
    const variantsElem = document.getElementById('metric-variants');
    if (variantsElem && variantsElem.textContent === '--') {
        // Keep as -- until we parse variant counts from logs
    }
}

function updateGPUMetricsPanel(gpuMetrics) {
    const setMetric = (id, value, isActive) => {
        const elem = document.getElementById(id);
        if (elem) {
            elem.textContent = value;
            if (isActive && value !== '--') {
                elem.classList.add('active');
            } else {
                elem.classList.remove('active');
            }
        }
    };

    if (!gpuMetrics || (gpuMetrics.iops === 0 && gpuMetrics.power_draw === 0)) {
        // No active GPU work
        setMetric('metric-iops', '--', false);
        setMetric('metric-bandwidth', '--', false);
        setMetric('metric-sm', '--', false);
        setMetric('metric-tensor', '--', false);
        setMetric('metric-pcie', '--', false);
        setMetric('metric-power', '--', false);
        return;
    }

    // Update metrics with values
    setMetric('metric-iops', gpuMetrics.iops > 0 ? gpuMetrics.iops.toFixed(0) : '--', gpuMetrics.iops > 0);
    setMetric('metric-bandwidth', gpuMetrics.memory_bandwidth > 0 ? gpuMetrics.memory_bandwidth.toFixed(0) : '--', gpuMetrics.memory_bandwidth > 0);
    setMetric('metric-sm', gpuMetrics.sm_efficiency > 0 ? gpuMetrics.sm_efficiency.toFixed(0) : '--', gpuMetrics.sm_efficiency > 0);
    setMetric('metric-tensor', gpuMetrics.tensor_usage > 0 ? gpuMetrics.tensor_usage.toFixed(0) : '--', gpuMetrics.tensor_usage > 0);
    setMetric('metric-pcie', gpuMetrics.pcie_throughput > 0 ? gpuMetrics.pcie_throughput.toFixed(1) : '--', gpuMetrics.pcie_throughput > 0);
    setMetric('metric-power', gpuMetrics.power_draw > 0 ? gpuMetrics.power_draw.toFixed(0) : '--', gpuMetrics.power_draw > 0);
}

function updatePipelineFlow(currentStep, status) {
    const stepNameToStage = {
        'Prerequisites Check': 'check',
        'NGC Login': 'login',
        'Data Download': 'download',
        'Reference Setup': 'reference',
        'Chr20 Test': 'test',
        'Full Genome': 'full'
    };

    const stageOrder = ['check', 'login', 'download', 'reference', 'test', 'full'];

    // Update all stages
    document.querySelectorAll('.pipeline-stage').forEach(stage => {
        const stageId = stage.dataset.stage;
        const currentStageId = stepNameToStage[currentStep];

        // Remove all state classes
        stage.classList.remove('pending', 'running', 'completed', 'error');

        if (completedStages.has(stageId)) {
            stage.classList.add('completed');
        } else if (stageId === currentStageId) {
            if (status === 'running') {
                stage.classList.add('running');
            } else if (status === 'success') {
                stage.classList.add('completed');
                completedStages.add(stageId);
            } else if (status === 'error') {
                stage.classList.add('error');
            } else {
                stage.classList.add('pending');
            }
        } else {
            stage.classList.add('pending');
        }
    });
}

function updateExecutiveSummary(status) {
    const summaryStatus = document.getElementById('summary-status');
    const summaryStage = document.getElementById('summary-stage');
    const summaryEta = document.getElementById('summary-eta');

    if (status.status === 'running') {
        summaryStatus.textContent = 'Running';
        summaryStatus.style.color = '#48bb78';
        summaryStage.textContent = status.current_step || 'Processing...';

        // Start runtime timer if not already running
        if (!runtimeInterval && status.start_time) {
            pipelineStartTime = new Date(status.start_time);
            startRuntimeTimer();
        }

        // Update ETA based on current step
        const etaMap = {
            'Prerequisites Check': 'ETA: < 1 min',
            'NGC Login': 'ETA: < 1 min',
            'Data Download': 'ETA: 2-4 hours',
            'Reference Setup': 'ETA: 5-15 min',
            'Chr20 Test': 'ETA: 5-20 min',
            'Full Genome': 'ETA: 30-110 min'
        };
        summaryEta.textContent = etaMap[status.current_step] || 'Calculating...';

    } else if (status.status === 'success') {
        summaryStatus.textContent = 'Complete';
        summaryStatus.style.color = '#48bb78';  // Green
        summaryStage.textContent = status.current_step + ' ✓';
        summaryEta.textContent = 'Step completed successfully';
        // Stop timer when step completes but preserve the runtime display
        stopRuntimeTimer();
        document.getElementById('summary-bases-per-min').textContent = '0.00/min';
        // Don't reset the runtime - it should show the final value

    } else if (status.status === 'error') {
        summaryStatus.textContent = 'Error';
        summaryStatus.style.color = '#dc3545';
        summaryStage.textContent = status.current_step || 'Failed';
        summaryEta.textContent = 'Check logs for details';
        stopRuntimeTimer();
        document.getElementById('summary-bases-per-min').textContent = '0.00/min';

    } else {
        summaryStatus.textContent = 'Idle';
        summaryStatus.style.color = '#6c757d';
        summaryStage.textContent = 'Ready to begin';
        summaryEta.textContent = 'Click a workflow step to start';
        stopRuntimeTimer();
        document.getElementById('summary-runtime').textContent = '00:00:00';
        document.getElementById('summary-bases-per-min').textContent = '0.00/min';
    }

    // Update pipeline flow visualization
    updatePipelineFlow(status.current_step, status.status);
}

function startRuntimeTimer() {
    if (runtimeInterval) return;

    runtimeInterval = setInterval(() => {
        if (!pipelineStartTime) return;

        const now = new Date();
        const elapsed = Math.floor((now - pipelineStartTime) / 1000);

        const hours = Math.floor(elapsed / 3600);
        const minutes = Math.floor((elapsed % 3600) / 60);
        const seconds = elapsed % 60;

        const formatted = [
            hours.toString().padStart(2, '0'),
            minutes.toString().padStart(2, '0'),
            seconds.toString().padStart(2, '0')
        ].join(':');

        document.getElementById('summary-runtime').textContent = formatted;

        // Update data processed estimate based on runtime and throughput
        const throughputElem = document.getElementById('summary-throughput');
        const throughputValue = parseFloat(throughputElem.textContent) || 0;
        const dataProcessed = (throughputValue * elapsed / 3600).toFixed(1);
        document.getElementById('summary-data').textContent = `${dataProcessed} GB`;

    }, 1000);
}

function stopRuntimeTimer() {
    if (runtimeInterval) {
        clearInterval(runtimeInterval);
        runtimeInterval = null;
    }
}
