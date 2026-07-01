# Web Portal Quick Start Guide

## 🎯 Launch the Portal (30 seconds)

```bash
cd genomics-pipeline/web-portal
./start-portal.sh
```

**Then open your browser to:** http://localhost:5000

That's it! The portal is running.

---

## 🖥️ First Time Setup (5 minutes)

### Step 1: Start the Portal
```bash
cd genomics-pipeline/web-portal
./start-portal.sh
```

Wait for: `Running on http://0.0.0.0:5000`

### Step 2: Open in Browser
Navigate to: **http://localhost:5000**

### Step 3: Run Prerequisites Check
1. Click **"Prerequisites Check"** in the left sidebar
2. Watch the Console Output tab for results
3. Verify all checks pass (green badges)

---

## 📱 Portal Interface Overview

```
┌─────────────────────────────────────────────────────────────┐
│  🧬 Genomics Pipeline Portal                    [Idle] 🟢    │
├──────────────┬──────────────────────────────────────────────┤
│              │                                              │
│  Workflow    │           Console Output                     │
│  Steps       │  ┌────────────────────────────────────────┐ │
│  ─────────   │  │ [2024-01-05 10:30:15] Starting...     │ │
│              │  │ [2024-01-05 10:30:16] Checking Docker │ │
│  1. Check ✓  │  │ [2024-01-05 10:30:17] ✓ Docker OK     │ │
│  2. Login    │  │ [2024-01-05 10:30:18] ✓ GPU OK        │ │
│  3. Download │  │                                        │ │
│  4. Reference│  └────────────────────────────────────────┘ │
│  5. Chr20    │                                              │
│  6. Full     │    [Configuration] [Logs] [Help]            │
│              │                                              │
│  System      │                                              │
│  Docker: ✓   │                                              │
│  GPU: ✓      │                                              │
│  Disk: 500GB │                                              │
└──────────────┴──────────────────────────────────────────────┘
```

---

## ⚡ Complete Workflow (Visual Guide)

### Run All Steps in Order:

#### 1️⃣ Prerequisites Check
**What it does:** Verifies Docker, GPU, and disk space
**Duration:** ~30 seconds
**Click:** "Prerequisites Check" button

**Expected output:**
```
✅ Docker version: Docker version 24.0.0
✅ Docker daemon is running
✅ NVIDIA Container Runtime is working
✅ Sufficient disk space available
```

---

#### 2️⃣ NGC Login
**What it does:** Authenticates with NVIDIA NGC
**Duration:** ~1 minute
**Click:** "NGC Login" button

**You'll need:**
- Username: `$oauthtoken`
- Password: Your NGC API key (from https://ngc.nvidia.com/setup/api-key)

**Expected output:**
```
Login Succeeded
✅ NGC authentication successful!
```

---

#### 3️⃣ Download Data
**What it does:** Downloads GIAB HG002 FASTQ files (~200GB)
**Duration:** 2-6 hours
**Click:** "Download Data" button

⚠️ **Warning:** This downloads 200GB of data!

**Progress indicators:**
```
Downloading R1 files...
[aria2] 15GB/100GB downloaded (15%)
```

**You can:**
- ☕ Take a break
- 📊 Monitor in "Console Output" tab
- 🔄 Close browser (portal keeps running)
- 🛑 Stop anytime (resume later)

---

#### 4️⃣ Setup Reference
**What it does:** Downloads and indexes GRCh38 genome
**Duration:** 5-15 minutes
**Click:** "Setup Reference" button

**Expected output:**
```
Downloading Parabricks sample bundle...
Extracting reference genome...
Building BWA index and FAI...
Creating sequence dictionary...
✅ Reference genome setup complete!
```

---

#### 5️⃣ Chr20 Test (Optional but Recommended)
**What it does:** Quick test on chromosome 20 only
**Duration:** 5-20 minutes
**Click:** "Chr20 Test" button

**Why run this:**
- ✓ Validates pipeline before full run
- ✓ Tests GPU acceleration
- ✓ Catches configuration issues early

**Expected output:**
```
Running fq2bam (chr20 only)...
Running DeepVariant (chr20 only)...
✅ Chr20 test complete!

Files created:
- HG002.chr20.bam
- HG002.chr20.vcf.gz
```

---

#### 6️⃣ Full Genome
**What it does:** Complete genome analysis
**Duration:** 30-110 minutes
**Click:** "Full Genome" button

**Progress stages:**
```
[1/4] Running fq2bam (whole genome)...      [20-75 min]
[2/4] Indexing and QC genome BAM...         [1-6 min]
[3/4] Running DeepVariant (whole genome)... [10-35 min]
[4/4] Indexing VCF...                       [1-2 min]
```

**Final output:**
```
✅ Full genome pipeline complete!

Output VCF: /data/output/HG002.genome.vcf.gz
Output VCF index: /data/output/HG002.genome.vcf.gz.tbi

This VCF is ready for Stage 2 (Evidence RAG/Chat)!
```

---

## 🎛️ Using the Tabs

### 📟 Console Output Tab
- **Real-time** command output
- **Auto-scrolls** to latest
- **Shows** all stdout/stderr

### ⚙️ Configuration Tab
Edit pipeline settings:
- **NUM_GPUS**: 1, 2, 4, etc.
- **LOW_MEMORY**: Enable if GPU OOM errors
- **PATIENT_ID**: Sample name
- **PB_IMG**: Container version

Click **"Save Configuration"** to apply.

### 📄 Pipeline Logs Tab
View historical logs:
1. Select log from dropdown
2. Browse complete output
3. Search with Ctrl+F

Available logs:
- Chr20 fq2bam
- Chr20 DeepVariant
- Full Genome fq2bam
- Full Genome DeepVariant

### ❓ Help Tab
- Quick start guide
- System requirements
- Expected timings
- Troubleshooting tips

---

## 🔍 Monitoring While Running

### Real-Time Indicators

**1. Status Badge (Top Right)**
- 🟢 Idle: No process running
- 🟡 Running: Step in progress
- ✅ Success: Completed successfully
- 🔴 Error: Failed (check logs)

**2. Current Step Card**
Shows when running:
```
⚠️ Running: Download Data
Started: 10:30:15 AM
Status: Running
[Stop] button
```

**3. Workflow Steps (Left Sidebar)**
- ⏳ Not started
- 🟡 Running (yellow highlight)
- ✅ Complete (green checkmark)
- ❌ Failed (red X)

**4. System Status Panel**
Live updates:
- Docker: Available / Not Found
- GPU: Available / Not Found
- Disk: 500GB free / 75% used

**5. Data Files Panel**
Track outputs:
- FASTQ R1: 45.2 GB
- FASTQ R2: 45.1 GB
- Reference: 3.1 GB
- Chr20 VCF: 125 MB
- Genome VCF: 2.3 GB

---

## 🛑 Stop a Running Process

1. Click **"Stop"** button (red, top of console)
2. Confirm the dialog
3. Process terminates gracefully
4. You can resume later

---

## 💡 Pro Tips

### Tip 1: Run in Background
```bash
# Start portal in background
cd web-portal
nohup ./start-portal.sh > portal.log 2>&1 &

# Access from any computer on your network
http://your-server-ip:5000
```

### Tip 2: Monitor GPU
Open a second terminal:
```bash
watch -n 1 nvidia-smi
```

### Tip 3: Check Disk Usage
```bash
df -h genomics-pipeline
```

### Tip 4: Resume After Disconnect
- Portal keeps running if you close browser
- Just open http://localhost:5000 again
- Console shows current status

### Tip 5: Parallel Monitoring
- Portal for high-level status
- Terminal for detailed `nvidia-smi`
- Logs tab for historical review

---

## 🚨 Common Issues & Fixes

### Issue: Port 5000 already in use
```bash
# Find what's using port 5000
lsof -i :5000

# Kill it or change portal port in server.py
```

### Issue: "Docker not found"
```bash
# Install Docker
curl -fsSL https://get.docker.com | sh

# Add user to docker group
sudo usermod -aG docker $USER
# Log out and back in
```

### Issue: "GPU not available"
```bash
# Check GPU driver
nvidia-smi

# Install nvidia-container-toolkit
# (see main README troubleshooting)
```

### Issue: Steps won't start
```bash
# Check script permissions
chmod +x ../scripts/*.sh

# Check paths in server.py
```

---

## 🎬 Complete Workflow Timeline

```
Time    Step                     Action
─────────────────────────────────────────────────────
0:00    Start Portal            ./start-portal.sh
0:01    Open Browser            http://localhost:5000
0:02    Prerequisites Check     Click → Wait 30s
0:03    NGC Login              Click → Enter credentials
0:05    Download Data          Click → Wait 2-6 hours
6:00    Setup Reference        Click → Wait 10 minutes
6:10    Chr20 Test             Click → Wait 15 minutes
6:25    Full Genome            Click → Wait 60 minutes
7:25    ✅ Complete!            HG002.genome.vcf.gz ready
```

---

## 🎯 You're Ready!

The web portal makes the entire genomics pipeline:
- ✅ **Visual** - See every step
- ✅ **Simple** - Click to run
- ✅ **Monitored** - Real-time status
- ✅ **Documented** - Built-in help

**Just open the portal and click through the steps!**

---

## 📚 Next Steps

After completing the pipeline:
1. **Find your VCF:** `data/output/HG002.genome.vcf.gz`
2. **Use for Stage 2:** Evidence RAG analysis
3. **Explore logs:** Review detailed execution logs
4. **Optimize:** Adjust GPU/memory settings if needed

---

**Questions?** Check the Help tab in the portal or see README.md
