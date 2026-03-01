/**
 * AI Factory Landing Page - Interactive Features
 * With real-time service health checking
 */

document.addEventListener('DOMContentLoaded', () => {
    initParticles();
    initSmoothScroll();
    initNavbarScroll();
    initCounterAnimation();
    initScrollReveal();
    initServiceHealthChecks();
    initToastNotifications();
    initLaunchInterceptors();
});

// ============================================================================
// SERVICE HEALTH CHECKING
// ============================================================================

// Service mapping: CSS class -> service ID (matching server.py)
const SERVICE_MAP = {
    // Pipeline cards
    'genomics': 'genomics',
    'rag-portal': 'rag-portal',
    'rag-chat': 'rag-chat',
    'drug-main': 'drug-main',
    'drug-portal': 'drug-portal',
    // Intelligence Agent cards
    'cart-agent': 'cart-agent',
    'imaging-agent': 'imaging-agent',
    'onco-agent': 'onco-agent',
    // Monitor cards
    'grafana': 'grafana',
    'prometheus': 'prometheus',
    'node': 'node-exporter',
    'dcgm': 'dcgm'
};

// Cache for service status
let serviceStatusCache = {};
let lastStatusCheck = 0;
const STATUS_CACHE_TTL = 10000; // 10 seconds

/**
 * Initialize service health checking
 */
function initServiceHealthChecks() {
    // Initial check
    checkAllServicesStatus();

    // Periodic refresh every 30 seconds
    setInterval(checkAllServicesStatus, 30000);

    // Also check when page becomes visible again
    document.addEventListener('visibilitychange', () => {
        if (!document.hidden) {
            checkAllServicesStatus();
        }
    });
}

/**
 * Check all services status via API
 */
async function checkAllServicesStatus() {
    try {
        const response = await fetch('/api/check-services');
        if (!response.ok) throw new Error('API request failed');

        const data = await response.json();
        serviceStatusCache = data.services;
        lastStatusCheck = Date.now();

        // Update all status indicators
        updateAllStatusIndicators(data.services);

        // Update navbar status
        updateNavbarStatus(data.summary);

        console.log(`[Health Check] ${data.summary.online}/${data.summary.total} services online`);

    } catch (error) {
        console.error('[Health Check] Failed to check services:', error);
    }
}

/**
 * Check a single service status
 */
async function checkSingleServiceStatus(serviceId) {
    // Return cached result if recent
    if (serviceStatusCache[serviceId] && (Date.now() - lastStatusCheck < STATUS_CACHE_TTL)) {
        return serviceStatusCache[serviceId];
    }

    try {
        const response = await fetch(`/api/check-service/${serviceId}`);
        if (!response.ok) throw new Error('API request failed');

        const result = await response.json();
        serviceStatusCache[serviceId] = result;
        return result;

    } catch (error) {
        console.error(`[Health Check] Failed to check ${serviceId}:`, error);
        return { online: false, message: 'Check failed' };
    }
}

/**
 * Update all status indicators on the page
 */
function updateAllStatusIndicators(services) {
    // Update pipeline cards
    Object.entries(SERVICE_MAP).forEach(([cssClass, serviceId]) => {
        const service = services[serviceId];
        if (!service) return;

        // Find the card
        const pipelineCard = document.querySelector(`.pipeline-card.${cssClass}`);
        const monitorCard = document.querySelector(`.monitor-card.${cssClass}`);
        const card = pipelineCard || monitorCard;

        if (card) {
            updateCardStatus(card, service.online, service.name);
        }
    });
}

/**
 * Update a single card's status indicator
 */
function updateCardStatus(card, isOnline, serviceName) {
    const statusDiv = card.querySelector('.card-status');
    if (statusDiv) {
        if (isOnline) {
            statusDiv.className = 'card-status online';
            statusDiv.innerHTML = '<span class="status-indicator"></span><span>Online</span>';
        } else {
            statusDiv.className = 'card-status offline';
            statusDiv.innerHTML = '<span class="status-indicator"></span><span>Offline</span>';
        }
    }

    // Also update monitor card styling
    if (card.classList.contains('monitor-card')) {
        if (isOnline) {
            card.classList.remove('service-offline');
            card.classList.add('service-online');
        } else {
            card.classList.remove('service-online');
            card.classList.add('service-offline');
        }
    }
}

/**
 * Update navbar status indicator
 */
function updateNavbarStatus(summary) {
    const statusDot = document.querySelector('.nav-status .status-dot');
    const statusText = document.querySelector('.nav-status .status-text');

    if (statusDot && statusText) {
        if (summary.all_healthy) {
            statusDot.className = 'status-dot pulse';
            statusDot.style.background = '#22c55e';
            statusText.textContent = 'All Systems Operational';
            statusText.style.color = '#22c55e';
        } else if (summary.online > 0) {
            statusDot.className = 'status-dot pulse';
            statusDot.style.background = '#f59e0b';
            statusText.textContent = `${summary.online}/${summary.total} Services Online`;
            statusText.style.color = '#f59e0b';
        } else {
            statusDot.className = 'status-dot';
            statusDot.style.background = '#ef4444';
            statusText.textContent = 'Services Offline';
            statusText.style.color = '#ef4444';
        }
    }
}

// ============================================================================
// LAUNCH INTERCEPTORS
// ============================================================================

/**
 * Initialize click interceptors for service cards
 */
function initLaunchInterceptors() {
    // Intercept pipeline card clicks
    document.querySelectorAll('.pipeline-card').forEach(card => {
        card.addEventListener('click', handleServiceLaunch);
    });

    // Intercept monitor card clicks
    document.querySelectorAll('.monitor-card').forEach(card => {
        card.addEventListener('click', handleServiceLaunch);
    });
}

/**
 * Handle service launch with health check
 */
async function handleServiceLaunch(event) {
    const card = event.currentTarget;
    const href = card.getAttribute('href');

    // Find the service ID from card class
    let serviceId = null;
    for (const [cssClass, svcId] of Object.entries(SERVICE_MAP)) {
        if (card.classList.contains(cssClass)) {
            serviceId = svcId;
            break;
        }
    }

    // If we can't identify the service, allow normal navigation
    if (!serviceId) {
        return true;
    }

    // Prevent default navigation
    event.preventDefault();

    // Show checking indicator
    const originalContent = showCheckingState(card);

    // Check service status
    const status = await checkSingleServiceStatus(serviceId);

    // Restore card state
    restoreCardState(card, originalContent);

    if (status.online) {
        // Service is online - proceed with navigation
        window.open(href, '_blank');
    } else {
        // Service is offline - show notification
        showServiceOfflineNotification(status.name || serviceId, status.port, status.message);
    }
}

/**
 * Show checking state on card
 */
function showCheckingState(card) {
    const footer = card.querySelector('.card-footer, .monitor-port');
    const originalContent = footer ? footer.innerHTML : '';

    if (footer) {
        if (card.classList.contains('pipeline-card')) {
            footer.innerHTML = '<span class="checking-indicator">Checking...</span>';
        }
    }

    card.style.opacity = '0.7';
    card.style.pointerEvents = 'none';

    return originalContent;
}

/**
 * Restore card state after check
 */
function restoreCardState(card, originalContent) {
    const footer = card.querySelector('.card-footer, .monitor-port');
    if (footer && originalContent) {
        footer.innerHTML = originalContent;
    }

    card.style.opacity = '';
    card.style.pointerEvents = '';
}

// ============================================================================
// TOAST NOTIFICATIONS
// ============================================================================

let toastContainer = null;

/**
 * Initialize toast notification system
 */
function initToastNotifications() {
    // Create toast container
    toastContainer = document.createElement('div');
    toastContainer.className = 'toast-container';
    document.body.appendChild(toastContainer);

    // Add styles
    const style = document.createElement('style');
    style.textContent = `
        .toast-container {
            position: fixed;
            top: 100px;
            right: 20px;
            z-index: 10000;
            display: flex;
            flex-direction: column;
            gap: 10px;
            max-width: 400px;
        }

        .toast {
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            border: 1px solid rgba(239, 68, 68, 0.5);
            border-radius: 12px;
            padding: 16px 20px;
            box-shadow: 0 10px 40px rgba(0, 0, 0, 0.5);
            animation: slideIn 0.3s ease-out;
            display: flex;
            gap: 12px;
            align-items: flex-start;
        }

        .toast.toast-success {
            border-color: rgba(34, 197, 94, 0.5);
        }

        .toast.toast-warning {
            border-color: rgba(245, 158, 11, 0.5);
        }

        .toast.toast-error {
            border-color: rgba(239, 68, 68, 0.5);
        }

        .toast-icon {
            font-size: 24px;
            line-height: 1;
        }

        .toast-content {
            flex: 1;
        }

        .toast-title {
            font-weight: 600;
            color: #fff;
            margin-bottom: 4px;
            font-size: 14px;
        }

        .toast-message {
            color: #9ca3af;
            font-size: 13px;
            line-height: 1.4;
        }

        .toast-close {
            background: none;
            border: none;
            color: #6b7280;
            cursor: pointer;
            padding: 0;
            font-size: 18px;
            line-height: 1;
            transition: color 0.2s;
        }

        .toast-close:hover {
            color: #fff;
        }

        .toast-action {
            margin-top: 10px;
        }

        .toast-action button {
            background: rgba(118, 185, 0, 0.2);
            border: 1px solid rgba(118, 185, 0, 0.5);
            color: #76B900;
            padding: 6px 12px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 12px;
            font-weight: 500;
            transition: all 0.2s;
        }

        .toast-action button:hover {
            background: rgba(118, 185, 0, 0.3);
        }

        .checking-indicator {
            color: #76B900;
            font-size: 12px;
            animation: pulse 1s infinite;
        }

        @keyframes slideIn {
            from {
                transform: translateX(100%);
                opacity: 0;
            }
            to {
                transform: translateX(0);
                opacity: 1;
            }
        }

        @keyframes slideOut {
            from {
                transform: translateX(0);
                opacity: 1;
            }
            to {
                transform: translateX(100%);
                opacity: 0;
            }
        }

        @keyframes pulse {
            0%, 100% { opacity: 1; }
            50% { opacity: 0.5; }
        }

        /* Offline card styling */
        .card-status.offline {
            background: rgba(239, 68, 68, 0.1) !important;
        }

        .card-status.offline .status-indicator {
            background: #ef4444 !important;
            box-shadow: 0 0 10px rgba(239, 68, 68, 0.5) !important;
        }

        .card-status.offline span:last-child {
            color: #ef4444 !important;
        }

        .monitor-card.service-offline {
            opacity: 0.6;
            border-color: rgba(239, 68, 68, 0.3) !important;
        }
    `;
    document.head.appendChild(style);
}

/**
 * Show a toast notification
 */
function showToast(title, message, type = 'error', duration = 5000, actions = []) {
    const toast = document.createElement('div');
    toast.className = `toast toast-${type}`;

    const icons = {
        error: '⚠️',
        success: '✓',
        warning: '⚡',
        info: 'ℹ️'
    };

    let actionsHtml = '';
    if (actions.length > 0) {
        actionsHtml = '<div class="toast-action">' +
            actions.map(a => `<button onclick="${a.onclick}">${a.label}</button>`).join(' ') +
            '</div>';
    }

    toast.innerHTML = `
        <span class="toast-icon">${icons[type] || icons.info}</span>
        <div class="toast-content">
            <div class="toast-title">${title}</div>
            <div class="toast-message">${message}</div>
            ${actionsHtml}
        </div>
        <button class="toast-close" onclick="this.parentElement.remove()">×</button>
    `;

    toastContainer.appendChild(toast);

    // Auto-remove after duration
    if (duration > 0) {
        setTimeout(() => {
            toast.style.animation = 'slideOut 0.3s ease-out forwards';
            setTimeout(() => toast.remove(), 300);
        }, duration);
    }

    return toast;
}

/**
 * Show service offline notification
 */
function showServiceOfflineNotification(serviceName, port, reason) {
    showToast(
        `${serviceName} is Offline`,
        `The service on port ${port} is not responding. ${reason ? `(${reason})` : ''}<br><br>` +
        `To start services, run:<br><code style="background: #000; padding: 2px 6px; border-radius: 4px;">./start-services.sh</code>`,
        'error',
        8000
    );
}

// ============================================================================
// PARTICLES & ANIMATIONS
// ============================================================================

/**
 * Create floating particles in background
 */
function initParticles() {
    const container = document.getElementById('particles');
    if (!container) return;

    const particleCount = 30;

    for (let i = 0; i < particleCount; i++) {
        const particle = document.createElement('div');
        particle.className = 'particle';
        particle.style.left = `${Math.random() * 100}%`;
        particle.style.animationDuration = `${15 + Math.random() * 20}s`;
        particle.style.animationDelay = `${Math.random() * 15}s`;
        particle.style.width = `${2 + Math.random() * 4}px`;
        particle.style.height = particle.style.width;
        container.appendChild(particle);
    }
}

/**
 * Smooth scrolling for anchor links
 */
function initSmoothScroll() {
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function(e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                const navHeight = document.querySelector('.navbar').offsetHeight;
                const targetPosition = target.getBoundingClientRect().top + window.pageYOffset - navHeight - 20;

                window.scrollTo({
                    top: targetPosition,
                    behavior: 'smooth'
                });
            }
        });
    });
}

/**
 * Navbar background on scroll
 */
function initNavbarScroll() {
    const navbar = document.querySelector('.navbar');
    if (!navbar) return;

    let lastScroll = 0;

    window.addEventListener('scroll', () => {
        const currentScroll = window.pageYOffset;

        if (currentScroll > 100) {
            navbar.style.background = 'rgba(10, 10, 15, 0.95)';
            navbar.style.boxShadow = '0 4px 20px rgba(0, 0, 0, 0.3)';
        } else {
            navbar.style.background = 'rgba(10, 10, 15, 0.8)';
            navbar.style.boxShadow = 'none';
        }

        lastScroll = currentScroll;
    });
}

/**
 * Animate counter numbers when they come into view
 */
function initCounterAnimation() {
    const counters = document.querySelectorAll('.stat-value [data-value], .stat-value[data-value]');

    const observerOptions = {
        threshold: 0.5,
        rootMargin: '0px'
    };

    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                animateCounter(entry.target);
                observer.unobserve(entry.target);
            }
        });
    }, observerOptions);

    counters.forEach(counter => observer.observe(counter));
}

function animateCounter(element) {
    const target = parseFloat(element.dataset.value);
    const duration = 2000;
    const startTime = performance.now();
    const isDecimal = target % 1 !== 0;

    function update(currentTime) {
        const elapsed = currentTime - startTime;
        const progress = Math.min(elapsed / duration, 1);

        // Easing function
        const easeOut = 1 - Math.pow(1 - progress, 3);
        const current = target * easeOut;

        if (isDecimal) {
            element.textContent = current.toFixed(1);
        } else {
            element.textContent = Math.floor(current).toLocaleString();
        }

        if (progress < 1) {
            requestAnimationFrame(update);
        } else {
            if (isDecimal) {
                element.textContent = target.toFixed(1);
            } else {
                element.textContent = target.toLocaleString();
            }
        }
    }

    requestAnimationFrame(update);
}

/**
 * Reveal elements on scroll
 */
function initScrollReveal() {
    const revealElements = document.querySelectorAll('.pipeline-card, .stat-card, .monitor-card, .flow-step');

    const observerOptions = {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    };

    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.classList.add('active');
            }
        });
    }, observerOptions);

    revealElements.forEach(el => {
        el.classList.add('reveal');
        observer.observe(el);
    });
}

// ============================================================================
// KEYBOARD NAVIGATION
// ============================================================================

document.addEventListener('keydown', (e) => {
    // Press 1-5 to quickly navigate to pipeline cards
    if (e.key >= '1' && e.key <= '5' && !e.ctrlKey && !e.metaKey && !e.altKey) {
        const cards = document.querySelectorAll('.pipeline-card');
        const index = parseInt(e.key) - 1;
        if (cards[index] && document.activeElement.tagName !== 'INPUT') {
            cards[index].click();
        }
    }

    // Press 'g' to go to Grafana
    if (e.key === 'g' && !e.ctrlKey && !e.metaKey && !e.altKey) {
        const grafana = document.querySelector('.monitor-card.grafana');
        if (grafana && document.activeElement.tagName !== 'INPUT') {
            grafana.click();
        }
    }

    // Press 'r' to refresh service status
    if (e.key === 'r' && !e.ctrlKey && !e.metaKey && !e.altKey) {
        if (document.activeElement.tagName !== 'INPUT') {
            checkAllServicesStatus();
            showToast('Refreshing', 'Checking service status...', 'info', 2000);
        }
    }
});

// ============================================================================
// BUTTON EFFECTS
// ============================================================================

document.querySelectorAll('.btn').forEach(button => {
    button.addEventListener('click', function(e) {
        const ripple = document.createElement('span');
        const rect = this.getBoundingClientRect();

        ripple.style.cssText = `
            position: absolute;
            background: rgba(255, 255, 255, 0.3);
            border-radius: 50%;
            transform: scale(0);
            animation: ripple 0.6s ease-out;
            pointer-events: none;
            left: ${e.clientX - rect.left}px;
            top: ${e.clientY - rect.top}px;
            width: 100px;
            height: 100px;
            margin-left: -50px;
            margin-top: -50px;
        `;

        this.style.position = 'relative';
        this.style.overflow = 'hidden';
        this.appendChild(ripple);

        setTimeout(() => ripple.remove(), 600);
    });
});

// Add ripple animation to stylesheet
const rippleStyle = document.createElement('style');
rippleStyle.textContent = `
    @keyframes ripple {
        to {
            transform: scale(4);
            opacity: 0;
        }
    }
`;
document.head.appendChild(rippleStyle);

// ============================================================================
// CONSOLE WELCOME
// ============================================================================

console.log(`
%c AI Factory %c Precision Medicine to Drug Discovery
%c ═══════════════════════════════════════════════════════

 🧬 Genomics Pipeline    → Port 5000
 💬 RAG/Chat Pipeline    → Port 5001, 8501
 💊 Drug Discovery       → Port 8505, 8510

 🤖 Intelligence Agents
    CAR-T Agent          → Port 8521
    Imaging Agent        → Port 8525
    Oncology Agent       → Port 8526

 📊 Monitoring
    Grafana              → Port 3000
    Prometheus           → Port 9099

 ⌨️  Keyboard Shortcuts:
    1-5     Launch pipeline interfaces
    6-8     Launch intelligence agents
    g       Open Grafana
    r       Refresh service status

%c Powered by NVIDIA DGX Spark
`,
    'background: #76B900; color: #000; font-weight: bold; padding: 4px 8px; border-radius: 4px 0 0 4px;',
    'background: #1a1a24; color: #76B900; padding: 4px 8px; border-radius: 0 4px 4px 0;',
    'color: #666;',
    'color: #76B900; font-weight: bold;'
);
