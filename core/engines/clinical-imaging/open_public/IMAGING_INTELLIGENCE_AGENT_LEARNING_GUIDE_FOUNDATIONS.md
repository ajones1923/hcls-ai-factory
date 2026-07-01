# Imaging Intelligence Agent — Learning Guide (Foundations)

> **Purpose:** Explain every concept behind the HCLS Imaging Intelligence Agent at a high school reading level. Read this guide front to back to fill knowledge gaps so you can explain and demo the agent to any audience.
>
> **License:** Apache 2.0 | **Author:** Adam Jones | **Date:** February 2026

---

## Table of Contents

1. [What Is Medical Imaging?](#chapter-1--what-is-medical-imaging)
2. [How Computers See Images (AI/ML Basics)](#chapter-2--how-computers-see-images)
3. [The DICOM Standard](#chapter-3--the-dicom-standard)
4. [What Is a GPU and Why Does It Matter?](#chapter-4--what-is-a-gpu-and-why-does-it-matter)
5. [Containers and Docker](#chapter-5--containers-and-docker)
6. [Databases](#chapter-6--databases)
7. [The Four Clinical Workflows](#chapter-7--the-four-clinical-workflows)
8. [How the Agent Thinks (LangGraph and AI Reasoning)](#chapter-8--how-the-agent-thinks)
9. [Embeddings and Similarity Search](#chapter-9--embeddings-and-similarity-search)
10. [Talking to Hospital Systems (FHIR and Clinical Integration)](#chapter-10--talking-to-hospital-systems)
11. [RAG and Large Language Models](#chapter-11--rag-and-large-language-models)
12. [Monitoring and Observability](#chapter-12--monitoring-and-observability)
13. [The HCLS AI Factory (The Bigger Picture)](#chapter-13--the-hcls-ai-factory)
14. [Trust, Safety, and Regulation](#chapter-14--trust-safety-and-regulation)
15. [Test Yourself — Comprehensive Review](#test-yourself--comprehensive-review)
16. [Glossary](#glossary)

---

## Chapter 1 — What Is Medical Imaging?

### X-Rays: Shadows of Your Body

An X-ray works like shining a flashlight at your hand in a dark room. The light passes through your skin but gets blocked by your bones, so you see a shadow. In a real X-ray machine, the "flashlight" is a beam of X-ray radiation, and the "shadow" lands on a detector plate behind you. Dense things like bones appear white. Soft things like lungs appear dark. Air appears black.

The most common X-ray is a chest X-ray (often shortened to CXR). Doctors order it to check the lungs, heart, and ribcage. A single chest X-ray produces one flat image — like a photograph.

### CT Scans: Slicing a Loaf of Bread

CT stands for Computed Tomography. A CT scanner is a large doughnut-shaped machine. You lie on a table that slides through the hole. Inside the doughnut, an X-ray tube spins around you, taking hundreds of X-ray images from every angle. A computer then stacks those images together to build a 3D picture.

Think of it this way: if your body were a loaf of bread, each CT slice is like cutting one thin piece. You can look at any single slice, or you can stack all the slices to see the entire loaf in three dimensions. Radiologists scroll through these slices to find problems.

CT scans use X-ray radiation (like a regular X-ray) but take many more images. They are excellent for seeing bones, blood, and organs in detail. A single CT scan of the head can produce 200-300 individual slices.

### MRI: Magnets and Radio Waves

MRI stands for Magnetic Resonance Imaging. Unlike X-rays and CT, MRI uses no radiation at all. Instead, it uses a very powerful magnet and radio waves.

Here is the simplified version: your body is mostly water, and water contains hydrogen atoms. When you lie inside the MRI machine, the magnet lines up all the hydrogen atoms in your body. Then the machine sends radio wave pulses that knock those atoms out of alignment. When the atoms snap back into place, they release tiny signals. The machine detects those signals and uses them to build an image.

MRI is especially good at showing soft tissue — the brain, spinal cord, muscles, and ligaments. It produces detailed images without any radiation exposure. The tradeoff is that MRI scans take longer (often 30-60 minutes) and the machines are more expensive.

One important MRI sequence is called FLAIR (Fluid-Attenuated Inversion Recovery). FLAIR suppresses the signal from normal fluid in the brain so that abnormal spots — like damaged areas in multiple sclerosis — stand out brightly. The agent uses FLAIR images specifically for MS lesion tracking.

### What Doctors Order and Why

| Exam | Best For | Radiation? | Speed |
|---|---|---|---|
| X-Ray (CXR) | Lungs, bones, heart size | Yes (very low dose) | Seconds |
| CT Scan | Brain bleeds, lung nodules, trauma | Yes (moderate dose) | 1-5 minutes |
| MRI | Brain tissue, spinal cord, MS lesions | No | 30-60 minutes |

A doctor picks the imaging type based on what they are looking for and how quickly they need the answer. A suspected brain bleed needs a CT because it is fast. A suspected MS flare-up needs an MRI because it shows soft tissue detail.

### Studies and Series

When you get a scan, the hospital stores it as a **study**. A study is the complete collection of images from one scan session. For example, "CT Head performed on January 15, 2026" is one study.

Inside a study, there can be multiple **series**. Each series is a set of images acquired with the same settings. A CT head study might have a series taken without contrast dye and another series taken after injecting contrast dye. An MRI study might have a T1-weighted series, a T2-weighted series, and a FLAIR series — each one shows different tissue characteristics.

Think of it this way: a study is like an album of photos from a vacation. Each series is a group of photos taken at the same location with the same camera settings.

### What Radiologists Do

A radiologist is a doctor who specializes in reading medical images. Their daily workflow looks like this:

1. A scan arrives from the scanner.
2. It appears on the radiologist's **worklist** — their to-do list.
3. The radiologist opens the images, scrolls through them, and looks for abnormalities.
4. They write a **report** describing what they found.
5. The report goes back to the ordering doctor, who uses it to make treatment decisions.

A busy hospital radiologist may read hundreds of studies per day, each containing hundreds or thousands of images. That is millions of images per year. Fatigue and volume are real problems. AI can help by triaging urgent cases to the top of the worklist and flagging areas of concern — so the radiologist knows where to look first.

### Why It Matters

The Imaging Intelligence Agent exists because radiologists are overwhelmed. It does not replace them. It acts like an assistant that pre-reads every scan, flags emergencies, measures changes over time, and writes draft findings — so the radiologist can work faster and more accurately.

### Test Yourself — Chapter 1

1. What is the difference between a CT scan and an MRI?
2. What does FLAIR stand for, and why is it useful for MS?
3. What is the difference between a study and a series?
4. Why can't AI replace radiologists entirely?

---

## Chapter 2 — How Computers See Images

### Pixels and Voxels: Images Are Numbers

When you look at a photo on your phone, you see colors and shapes. But to a computer, that photo is just a grid of numbers. Each tiny square in the grid is called a **pixel** (picture element). A pixel stores a number representing brightness or color.

A chest X-ray might be 2,048 x 2,048 pixels — that is about 4 million numbers arranged in a grid. A bright pixel (high number) means dense tissue. A dark pixel (low number) means air or soft tissue.

CT and MRI scans are three-dimensional, so instead of pixels they have **voxels** (volume elements). A voxel is a tiny cube in 3D space, like one small block in a Minecraft world. A CT head scan might contain 512 x 512 x 300 voxels — about 78 million numbers. The computer processes every single one.

### What a Neural Network Is

A neural network is a type of computer program that learns patterns from examples. It is loosely inspired by how the brain works — layers of connected units that each do a small piece of math.

Think of it this way: imagine you are teaching a child to recognize dogs. You show them thousands of pictures labeled "dog" or "not dog." At first they guess randomly. But every time they get it wrong, you give them feedback, and they adjust. After seeing thousands of examples, they get very good at it. That is exactly how a neural network learns — except instead of a child, it is millions of math equations, and instead of "dog," it might be "brain bleed" or "lung nodule."

A neural network has layers:

- **Input layer:** Receives the image (all those pixel/voxel numbers).
- **Hidden layers:** Each layer extracts features — simple things like edges at first, then complex things like shapes and textures in deeper layers.
- **Output layer:** Produces the answer (a classification, a detection box, or a segmentation mask).

The more layers a network has, the "deeper" it is. That is why you hear the term **deep learning** — it just means a neural network with many layers.

### Training vs. Inference

These are two distinct phases:

**Training** is like studying for a test. You feed the network thousands of labeled examples (images where a human expert has marked the correct answer). The network adjusts its internal math to get better at predicting the right answer. Training is slow and computationally expensive — it can take hours or days on powerful hardware.

**Inference** is like taking the test. You give the trained network a new image it has never seen before, and it produces a prediction. Inference is fast — usually seconds per image.

The Imaging Intelligence Agent only does inference. The models were already trained by researchers (often using MONAI, NVIDIA's medical imaging AI framework). The agent loads those pre-trained models and runs them on new patient scans.

### Classification: Is There a Problem?

Classification answers a yes/no (or multi-category) question about the entire image. "Does this chest X-ray show pneumothorax?" The output is a **confidence score** — a number between 0 and 1 (or 0% to 100%) representing how sure the model is.

For example, the agent's CXR workflow runs a classifier called DenseNet-121. It looks at the whole chest X-ray and outputs confidence scores for five conditions: pneumothorax, consolidation, pleural effusion, cardiomegaly, and fracture.

### Detection: Where Is the Problem?

Detection goes further than classification. It not only says "yes, there is a problem" but also draws a **bounding box** — a rectangle — around where the problem is. Detection models are used for the lung nodule workflow: they scan the entire CT volume and draw boxes around every suspicious nodule.

Think of it like the difference between someone saying "there is a typo on this page" (classification) versus circling the exact word (detection).

### Segmentation: What Exact Shape Is the Problem?

Segmentation is the most detailed level. Instead of drawing a box, it colors in the exact shape of the abnormality — pixel by pixel (or voxel by voxel). The result is called a **segmentation mask** — an image the same size as the original, where each pixel is labeled as either "abnormality" or "normal."

The agent uses segmentation for hemorrhage volume measurement (coloring in the exact area of blood), nodule volumetrics (measuring the exact 3D size of a nodule), and MS lesion tracking (identifying each lesion precisely).

Models used for segmentation in this agent include 3D U-Net and SegResNet from the MONAI Model Zoo.

### Confidence Scores

Every prediction comes with a confidence score. A model might say "I am 95% confident this is a hemorrhage" or "I am 40% confident this is a nodule." The agent uses these scores to make decisions:

- If confidence is above a threshold (for example, 0.5), the finding is reported.
- If confidence is below the threshold, the finding is discarded or flagged as uncertain.

Confidence scores are not the same as accuracy. A model can be 95% confident and still be wrong. That is why a human radiologist always reviews the results.

### Why It Matters

Classification, detection, and segmentation are the three fundamental AI techniques the agent uses on every scan. Understanding them helps you explain exactly what the agent is doing at each step of each workflow — not just "AI magic" but specific, understandable operations.

### Test Yourself — Chapter 2

1. What is the difference between a pixel and a voxel?
2. Explain the difference between training and inference in your own words.
3. What does a segmentation mask look like?
4. Why are confidence scores not the same as accuracy?

---

## Chapter 3 — The DICOM Standard

### What DICOM Is

DICOM stands for Digital Imaging and Communications in Medicine. It is a universal standard for how medical images are stored, transmitted, and displayed. Every CT scanner, MRI machine, X-ray unit, and PACS system in the world uses DICOM.

Think of it this way: JPEG is the format your phone uses for photos. DICOM is the format hospitals use for medical images. But DICOM is much more than just an image file — it includes the image data AND a large amount of metadata (information about the image).

### DICOM Tags: Metadata Attached to Every Image

Every DICOM file has hundreds of **tags** — labeled pieces of information. Some important tags include:

| Tag | Name | Example |
|---|---|---|
| (0010,0010) | Patient Name | DOE^JOHN |
| (0010,0020) | Patient ID | MRN12345 |
| (0008,0020) | Study Date | 20260115 |
| (0008,0060) | Modality | CT, MR, CR, DX |
| (0020,000D) | Study Instance UID | 1.2.840.113619... |
| (0020,000E) | Series Instance UID | 1.2.840.113619... |
| (0008,0018) | SOP Instance UID | 1.2.840.113619... |
| (0018,0015) | Body Part Examined | CHEST, HEAD, BRAIN |

The numbers in parentheses (like 0010,0010) are the tag addresses. You do not need to memorize them, but knowing they exist helps you understand that every medical image carries its own "passport" of information.

### UIDs: Unique Identifiers

Three UIDs are critical:

- **Study Instance UID:** Identifies the entire study (one scan session). Every image from that session shares this UID.
- **Series Instance UID:** Identifies a series within the study. Different imaging protocols produce different series.
- **SOP Instance UID:** Identifies a single image file. Every individual slice has its own SOP Instance UID.

Think of it like an address: the Study UID is the city, the Series UID is the street, and the SOP Instance UID is the house number. Together they uniquely identify any image anywhere in the world.

The Python library **pydicom** (MIT license) is what the agent uses to read and write DICOM files.

### DICOM Servers: Where Images Live

A DICOM server is a specialized database for storing and retrieving medical images. The agent uses **Orthanc**, an open-source DICOM server (GPLv3 license). Orthanc runs as a container on port 8042 for its web interface and port 4242 for DICOM network communication.

When a scanner finishes a scan, it sends the images to Orthanc using a protocol called DIMSE C-STORE (the traditional way medical devices send images). Orthanc stores the images and notifies the agent that a new study has arrived.

### DICOMweb: The Modern Way

DICOMweb is a newer standard that lets you access DICOM images using normal web protocols (HTTP). Instead of the old networking protocol (DIMSE), you use web requests — the same kind your browser uses to load websites.

The three main DICOMweb services are:

- **STOW-RS (Store):** Upload images to the server. "Here is a new study — please store it."
- **WADO-RS (Retrieve):** Download images from the server. "Give me all the slices for this study."
- **QIDO-RS (Search):** Search for studies on the server. "Find all CT head studies from the last week."

The agent uses DICOMweb to send completed results (like DICOM Structured Reports) back to the hospital's PACS.

### DICOM Structured Reports (DICOM SR)

When the agent finishes analyzing a scan, it needs to send the results back in a format the hospital system understands. A DICOM SR is a structured document that contains findings, measurements, and references — all stored as DICOM, so it appears right alongside the original images in the PACS viewer.

The agent creates DICOM SRs using a Python library called **highdicom**. Each SR contains:
- What was found (e.g., "hemorrhage detected")
- Where it was found (references to specific image slices)
- Measurements (e.g., "volume: 25.3 mL")
- Confidence scores
- The model version that produced the finding

### Why It Matters

The agent must speak DICOM to work with any hospital. Without DICOM compliance, the agent's findings could not appear in the radiologist's viewer, could not be archived alongside the original images, and could not be searched in the hospital's systems. DICOM is the lingua franca of medical imaging.

### Test Yourself — Chapter 3

1. What is the difference between DICOM and JPEG?
2. Name the three types of UIDs and what level each identifies.
3. What is the difference between DIMSE C-STORE and DICOMweb STOW-RS?
4. Why does the agent output findings as DICOM SR instead of just a text file?

---

## Chapter 4 — What Is a GPU and Why Does It Matter?

### CPU vs. GPU

Every computer has a CPU (Central Processing Unit). The CPU is like one very smart worker — it can do almost any task, but it does them one at a time (or a few at a time with multiple cores). A modern CPU might have 8 to 64 cores.

A GPU (Graphics Processing Unit) is different. It has thousands of small, simple cores that all work at the same time. Each individual core is less capable than a CPU core, but the sheer number of them working in parallel makes the GPU incredibly fast for certain tasks.

Think of it this way: a CPU is like one expert chef who can make any dish. A GPU is like a kitchen with 1,000 line cooks — each one can only do simple tasks, but together they can prepare 1,000 plates simultaneously.

### Why GPUs Are Great for AI

Neural networks are built from millions of simple math operations — mostly multiplying numbers and adding them together. These operations are independent of each other, meaning they can all happen at the same time. That is exactly what a GPU is designed to do.

Processing a CT scan through a segmentation model requires billions of multiplications. On a CPU, this might take minutes. On a GPU, it takes seconds — because thousands of cores work on different parts of the calculation simultaneously.

This is why the Imaging Intelligence Agent runs on NVIDIA GPU hardware. Without a GPU, the agent would be too slow to be useful in a clinical setting.

### The DGX Spark

The DGX Spark is the specific computer this agent is designed to run on. Here are its key specs explained:

**ARM64 (aarch64) Architecture:** Most laptops and desktops use x86_64 processors (made by Intel or AMD). The DGX Spark uses an ARM processor — the same family of chips used in smartphones and newer Mac computers. This matters because all software (including Docker containers) must be compiled for ARM64. You cannot run x86_64 software on ARM64 without emulation.

**NVIDIA Grace CPU:** The specific ARM64 processor in the DGX Spark. It is named "Grace" after Grace Hopper, a computer science pioneer.

**NVIDIA Blackwell GB10 GPU:** The GPU chip. "Blackwell" is the architecture generation name (like how Intel has "Core i7" or "Core i9"). GB10 is the specific model designed for the DGX Spark.

**128 GB Unified LPDDR5x Memory:** This is the most important concept to understand. In most computers, the CPU has its own memory (RAM) and the GPU has its own separate memory (VRAM). Data must be copied from CPU memory to GPU memory before the GPU can work on it, and copied back when done. On the DGX Spark, the CPU and GPU share one single 128 GB pool. No copying back and forth. This is called **unified memory** and it makes everything faster and simpler.

**NVMe Storage:** NVMe is a type of very fast solid-state drive. Think of it as a USB flash drive but 100 times faster. The DGX Spark can have up to 4 TB of NVMe storage. This is where DICOM images and derived artifacts are stored locally.

**GPUDirect Storage:** A technology that lets data flow directly from the NVMe drive into GPU memory, skipping the CPU entirely. Normally data goes: Drive → CPU memory → GPU memory. With GPUDirect Storage it goes: Drive → GPU memory. This eliminates a bottleneck.

**Price: $4,699:** The DGX Spark is desktop-sized — it can sit on a desk in a hospital office. At $4,699, it is positioned as a proof-of-concept and development machine, not a full production server. There is no NVIDIA AI Enterprise (NVAIE) software license cost at the desktop class.

### DGX Compute Progression

The DGX Spark is just the starting point. As the agent scales from proof-of-concept to full hospital deployment, the hardware scales too:

| Phase | Hardware | Price | What Changes |
|---|---|---|---|
| Proof Build | DGX Spark (1 GPU) | $4,699 | 1-2 workflows on a desk |
| Departmental | 1-2x DGX B200 (8 GPUs each) | $500K-$1M | All workflows, integrated with PACS |
| Multi-Site | 4-8x DGX B200 + InfiniBand | $2M-$4M | Multiple hospitals, federated learning |
| AI Factory | DGX SuperPOD (thousands of GPUs) | $7M-$60M+ | Thousands of concurrent studies |

The software stays the same across all tiers. The containers, the database schema, the agent logic — all identical. Only the hardware scales.

### Why It Matters

The GPU is the engine that makes real-time medical image analysis possible. Without it, the AI models would run too slowly for clinical use. The DGX Spark's unified memory and GPUDirect Storage are specifically designed to eliminate the data transfer bottlenecks that slow down image processing. Understanding this hardware explains why the agent can process a brain CT in under 90 seconds.

### Test Yourself — Chapter 4

1. What is the fundamental difference between a CPU and a GPU?
2. What does "unified memory" mean on the DGX Spark, and why is it advantageous?
3. What is GPUDirect Storage and what bottleneck does it eliminate?
4. Why must all containers be ARM64-compatible on the DGX Spark?

---

## Chapter 5 — Containers and Docker

### What a Container Is

A container is a self-contained package that includes everything a piece of software needs to run — the code, the libraries, the configuration files, and even a minimal operating system. When you run a container, it behaves as if it is running on its own dedicated machine, even though it is actually sharing the computer with other containers.

Think of it this way: a container is like a lunchbox. Everything the program needs for its meal is packed inside. It does not matter what kitchen (computer) you open the lunchbox in — the meal is the same every time.

### The "Works on My Machine" Problem

Without containers, installing software can be a nightmare. A program might need a specific version of Python, a specific database driver, and a specific math library. If any of those pieces are missing or the wrong version, the program breaks. "It works on my machine" is one of the most common frustrations in software engineering.

Containers solve this by packaging everything together. If a container runs on your laptop, it runs the same way on a server in a hospital. No surprises.

### Docker: The Tool That Runs Containers

Docker is the most widely used tool for building and running containers. Here are the key concepts:

**Docker Image:** A read-only template — like a recipe. It lists every ingredient (file, library, configuration) needed. You build an image once.

**Docker Container:** A running instance of an image — like the actual meal cooked from the recipe. You can run many containers from the same image, and each one is independent.

**Dockerfile:** A text file with step-by-step instructions for building an image. It says things like "start with Ubuntu, install Python 3.11, copy my code in, run this command."

**Docker Hub / Container Registry:** A website where pre-built images are stored, like an app store for containers. NVIDIA hosts their images at nvcr.io (the NVIDIA Container Registry).

### docker-compose: Running Many Containers Together

The agent is not one container — it is eleven containers that work together. Managing eleven containers individually would be tedious. Docker Compose is a tool that lets you define all your containers in a single file called `docker-compose.yml` and start them all with one command: `docker compose up`.

Think of docker-compose as an orchestra conductor. Each musician (container) plays their own instrument, but the conductor ensures they all start at the right time and play together.

### Volumes: Persistent Storage

By default, when you stop a container, everything inside it disappears — like erasing a whiteboard. But the agent needs to keep data between restarts (database contents, DICOM images, model weights). **Volumes** are Docker's solution — they are storage areas outside the container that survive restarts.

Think of it like saving your game. The container is the game session, and the volume is the save file on your hard drive. You can close the game and reopen it later without losing progress.

### Health Checks

Each container in the agent has a **health check** — an automatic "are you alive?" test that Docker runs every few seconds. If a container stops responding, Docker can automatically restart it. For example, the PostgreSQL container's health check runs a simple database query. If the query fails, Docker knows the database is down.

### The Agent's 11 Containers

| Container | Port | Purpose |
|---|---|---|
| imaging-orthanc | 4242, 8042 | DICOM server — stores and serves medical images |
| imaging-postgres | 5432 | Database — stores findings, measurements, embeddings |
| imaging-nim-llm | 8520 | Large language model — generates clinical reports |
| imaging-embedding | 8521 | Embedding service — creates similarity vectors |
| imaging-dicom-listener | 8522 | Listens for new studies arriving in Orthanc |
| imaging-fhir-publisher | 8523 | Sends results to EHR systems via FHIR |
| imaging-agent | 8524 | AI reasoning engine — LangGraph clinical agent |
| imaging-portal | 8525 | Web dashboard — Streamlit user interface |
| imaging-dcgm | 9400 | GPU metrics exporter — monitors GPU health |
| imaging-prometheus | 9099 | Metrics collector — gathers system statistics |
| imaging-grafana | 3000 | Dashboard — visualizes metrics as charts |

Each container is a specialist. The DICOM listener only listens for new images. The database only stores data. The agent only reasons about findings. This separation means if one piece has a problem, the others keep running.

### ARM64 Compatibility

Because the DGX Spark uses an ARM64 processor, every container must be built for ARM64 architecture. This is usually handled by selecting the right base image (for example, `nvidia/cuda:12.x-runtime-ubuntu22.04` has ARM64 variants) and ensuring all Python packages have ARM64 wheels available.

If you try to run an x86_64 container on the DGX Spark, it will either fail to start or run extremely slowly through emulation.

### Why It Matters

Containers are how the agent is packaged and deployed. Every piece of the system — from the DICOM server to the AI models to the monitoring tools — runs in its own container. This makes the agent portable (runs on any machine with Docker), reliable (containers restart automatically), and maintainable (you can update one piece without touching the others).

### Test Yourself — Chapter 5

1. What is the difference between a Docker image and a Docker container?
2. Why does the agent use 11 separate containers instead of putting everything in one?
3. What happens to data inside a container when it stops, and how do volumes fix this?
4. Why is ARM64 compatibility important for the DGX Spark?

---

## Chapter 6 — Databases

### What a Database Is

A database is an organized collection of information that you can search, filter, and update efficiently. A spreadsheet is the simplest analogy — rows of data organized into columns. But a database is much more powerful than a spreadsheet: it can handle millions of rows, multiple users reading and writing at the same time, and complex queries that would crash Excel.

### Tables, Rows, and Columns

A database stores data in **tables**. Each table is like a spreadsheet with named columns and rows:

- **Column (field):** Defines what kind of information is stored. Example: "patient_id", "study_date", "finding_label".
- **Row (record):** One entry. Example: one specific finding from one specific study.

The agent's database has seven main tables. Here is a simplified view:

| Table | What It Stores | Example Row |
|---|---|---|
| studies | One row per imaging study | "CT Head, 2026-01-15, Patient MRN12345" |
| series | One row per image series | "Axial non-contrast, 300 slices, 0.5mm spacing" |
| findings | One row per AI-detected finding | "Hemorrhage, 95% confidence, CRITICAL" |
| measurements | One row per measurement | "Volume: 25.3 mL, Midline shift: 8.2 mm" |
| embeddings | One row per similarity vector | "384-dimensional vector for study-level search" |
| provenance | One row per processing event | "Model v2.1, processed at 14:32:05, took 42 seconds" |
| worklist_entries | One row per worklist item | "CRITICAL priority, assigned to on-call radiologist" |

### SQL: Asking the Database Questions

SQL (Structured Query Language) is the language you use to talk to a database. It reads almost like English:

- "Show me all critical findings from today" becomes:
  `SELECT * FROM findings WHERE severity = 'CRITICAL' AND created_at >= '2026-02-02';`
- "How many studies were processed this week?" becomes:
  `SELECT COUNT(*) FROM studies WHERE study_date >= '2026-01-27';`

You do not need to write SQL to use the agent — the agent writes SQL internally. But understanding that SQL exists helps you understand how the agent retrieves and stores information.

### PostgreSQL

The agent uses **PostgreSQL** (often called "Postgres") as its database. PostgreSQL is free, open-source, and trusted by organizations from banks to governments. It runs inside the `imaging-postgres` container on port 5432.

PostgreSQL was chosen because:
1. It is free (PostgreSQL License — very permissive).
2. It supports the **pgvector** extension for similarity search (explained in Chapter 9).
3. It is battle-tested and widely supported.
4. It runs on ARM64 (the DGX Spark's architecture).

### Views: Saved Queries That Look Like Tables

A **view** is a saved SQL query that you can treat as if it were a table. The agent has two views:

- **active_worklist:** Shows only the incomplete worklist entries, sorted by urgency. Instead of writing a complex query every time, the agent just queries `active_worklist`.
- **study_summary:** Joins multiple tables to show a complete picture of each study — its series, findings, and measurements — all in one view.

### Indexes: Making Searches Fast

An **index** is a data structure that makes database searches faster. Think of it like the index at the back of a textbook — instead of reading every page to find a topic, you look it up in the index and jump straight to the right page.

Without an index, searching 100,000 findings for "all CRITICAL hemorrhages" would require reading every single row. With an index on the severity column, the database jumps straight to the matching rows.

The agent uses multiple indexes, including a specialized **HNSW index** for vector similarity search (covered in Chapter 9).

### Why It Matters

The database is the agent's memory. Every finding, every measurement, every embedding, and every decision is stored here. Without the database, the agent would forget everything the moment a workflow finished. The database is also what enables longitudinal tracking — comparing today's scan to last month's scan — because all prior results are stored and queryable.

### Test Yourself — Chapter 6

1. Name the seven tables in the agent's database and what each stores.
2. What is SQL and what does it do?
3. Why was PostgreSQL chosen over other databases?
4. What is a database index and why does it matter for performance?

---

## Chapter 7 — The Four Clinical Workflows

Every workflow in the agent follows the same general pattern:

1. A new study arrives in Orthanc (the DICOM server).
2. The DICOM listener detects it and triggers processing.
3. GPU inference runs one or more AI models on the images.
4. Post-processing extracts measurements and classifications.
5. Results are stored in PostgreSQL.
6. Outputs are sent to PACS (as DICOM SR) and EHR (as FHIR).
7. The worklist is updated with the appropriate priority.

### Workflow 1: Chest X-Ray (CXR) Rapid Findings

**Target time:** Less than 30 seconds end-to-end.

**What it detects:** Five conditions visible on a chest X-ray:

| Condition | What It Means |
|---|---|
| Pneumothorax | Collapsed lung — air trapped between the lung and chest wall |
| Consolidation | Part of the lung is filled with fluid or pus (often pneumonia) |
| Pleural Effusion | Fluid buildup around the lungs |
| Cardiomegaly | Enlarged heart |
| Fracture | Broken rib or other bone |

**How it works:**

The model used is **DenseNet-121** — a classification neural network with 121 layers. It looks at the entire chest X-ray and outputs a confidence score (0 to 1) for each of the five conditions.

After classification, the agent generates a **GradCAM heatmap**. GradCAM (Gradient-weighted Class Activation Mapping) highlights the regions of the image that most influenced the model's decision. It produces a colored overlay — red/yellow for high influence, blue for low influence — like highlighting text with a marker. This helps the radiologist see exactly where the AI is "looking."

If any condition has a confidence score above the threshold (for example, 0.5), it is reported as a finding. Pneumothorax above 0.8 triggers an immediate CRITICAL alert.

**Why it matters:** Chest X-rays are the single most common imaging exam worldwide. Being able to screen every one in under 30 seconds and flag emergencies immediately can save lives.

### Workflow 2: CT Head Hemorrhage Triage

**Target time:** Less than 90 seconds end-to-end.

**What it detects:** Bleeding (hemorrhage) inside the skull.

**Why speed matters:** A brain hemorrhage is a medical emergency. Every minute that passes without treatment, more brain cells die. Getting the right scan to the top of the radiologist's worklist can be the difference between a good outcome and a devastating one.

**How it works — step by step:**

1. **Classification:** A DenseNet-121 model scans the CT and determines whether hemorrhage is present (yes/no with confidence).

2. **Segmentation:** If hemorrhage is detected, a 3D U-Net model outlines the exact shape of the bleed — voxel by voxel — producing a segmentation mask.

3. **Volume estimation:** The agent calculates the volume of the hemorrhage by counting the marked voxels and multiplying by the voxel size. The result is in milliliters (mL).

4. **Midline shift measurement:** The brain is normally symmetric — left and right halves are mirror images. A large hemorrhage pushes the brain to one side. The agent measures this **midline shift** in millimeters. A shift greater than 5 mm is a red flag.

5. **Urgency classification:** Based on the Brain Trauma Foundation guidelines:

| Volume | Midline Shift | Urgency |
|---|---|---|
| > 30 mL | Any | CRITICAL — potential surgical emergency |
| 5-30 mL | > 5 mm | URGENT — needs immediate attention |
| 5-30 mL | ≤ 5 mm | URGENT — close monitoring |
| < 5 mL | ≤ 5 mm | ROUTINE — standard review |

6. **Worklist prioritization:** CRITICAL studies jump to the top of the radiologist's worklist. An alert notification is sent to the on-call team.

**Why it matters:** This workflow demonstrates the agent's full pipeline — from raw DICOM images to an actionable triage decision in under 90 seconds. It directly impacts patient outcomes by ensuring the most urgent cases are seen first.

### Workflow 3: CT Chest Lung Nodule Tracking

**Target time:** Less than 5 minutes (multi-stage).

**What it detects:** Small spots in the lungs called nodules. Most nodules are harmless (benign), but some are early-stage lung cancer. The challenge is figuring out which ones need attention.

**How it works — step by step:**

1. **Detection:** A RetinaNet model scans the CT chest and draws bounding boxes around every suspicious nodule. Detection sensitivity target: > 90% for nodules 4 mm or larger.

2. **Segmentation:** For each detected nodule, a SegResNet model outlines the exact 3D shape.

3. **Volumetrics:** The agent calculates each nodule's volume (in mm³) and longest diameter (in mm) from the segmentation mask.

4. **Prior study retrieval:** The agent searches the database for previous CT chest studies of the same patient. If a prior exists, it compares the current nodule measurements to the previous ones.

5. **Volume Doubling Time (VDT):** If a prior study exists, the agent calculates how fast the nodule is growing:

   **VDT = (Δt × ln 2) / ln(V2 / V1)**

   Where:
   - Δt = time between scans (in days)
   - V1 = volume at the previous scan
   - V2 = volume at the current scan
   - ln = natural logarithm
   - ln 2 ≈ 0.693

   A short VDT means the nodule is growing fast (suspicious). A long VDT means slow growth (less concerning). VDT < 400 days upgrades the risk category.

6. **Lung-RADS classification:** Lung-RADS (Lung CT Screening Reporting and Data System) is a standardized scoring system from the American College of Radiology. The agent assigns a category based on nodule size and growth:

| Lung-RADS | What It Means | Solid Nodule Size |
|---|---|---|
| 1 | Negative — no nodules | N/A |
| 2 | Benign appearance | < 6 mm |
| 3 | Probably benign | 6-8 mm |
| 4A | Suspicious | 8-15 mm |
| 4B | Very suspicious | ≥ 15 mm |

   Ground-glass nodules (GGN) — nodules that appear hazy rather than solid — have different size thresholds because they behave differently.

7. **Genomics trigger:** If a nodule scores Lung-RADS 4B or higher, the agent automatically triggers the genomics pipeline (Parabricks) for tumor profiling. This cross-modal trigger connects imaging to molecular analysis.

**Why it matters:** Lung cancer is the leading cause of cancer death worldwide. Early detection through CT screening saves lives, but tracking hundreds of nodules across millions of patients manually is impossible. This workflow automates the tedious measurement and comparison work, letting radiologists focus on clinical decisions.

### Workflow 4: MRI Brain MS Lesion Tracking

**Target time:** Less than 5 minutes (multi-stage).

**What MS is:** Multiple sclerosis (MS) is a disease where the immune system mistakenly attacks the protective insulation (myelin) around nerve fibers in the brain and spinal cord. Each attack creates a damaged area called a **lesion**. Over time, lesions accumulate and cause neurological problems — vision issues, numbness, difficulty walking.

**What it detects:** Lesions on FLAIR MRI sequences.

**How it works — step by step:**

1. **3D U-Net segmentation:** The model processes the FLAIR MRI and segments every lesion — producing a 3D mask where each lesion is individually labeled.

2. **Lesion counting and measurement:** The agent counts the total number of lesions and measures each one's volume.

3. **Spatial registration to prior MRI:** If a prior MRI exists, the agent **registers** (aligns) the current scan to the previous one. Brain scans taken on different days are never perfectly aligned — the patient's head is in a slightly different position each time. Registration warps one scan to match the other so they can be compared voxel-by-voxel.

4. **New / enlarging lesion identification:** By comparing the registered scans, the agent identifies:
   - **New lesions** that were not present before
   - **Enlarging lesions** that have grown since the last scan
   - **Stable lesions** that have not changed

5. **Disease activity classification:**

| New / Enlarging Lesions | Classification |
|---|---|
| 0 | Stable |
| 1-2 | Active |
| 3 or more | Highly Active |

6. **Longitudinal trajectory tracking:** The agent tracks disease activity across every MRI over time, building a trajectory: "Stable for 2 years, then active at the most recent scan." This trajectory helps neurologists decide whether to change treatment.

**Why it matters:** MS patients get MRI scans every 6-12 months to monitor disease progression. Manually counting and measuring dozens of tiny lesions across two scans is tedious and error-prone. The agent automates this comparison, giving neurologists a clear picture of disease activity in minutes.

### Test Yourself — Chapter 7

1. What five conditions does the CXR workflow detect?
2. At what hemorrhage volume does the CT head workflow classify a finding as CRITICAL?
3. What is Volume Doubling Time and why does a short VDT indicate concern?
4. How does the MRI workflow determine whether MS is "Stable," "Active," or "Highly Active"?

---

## Chapter 8 — How the Agent Thinks (LangGraph and AI Reasoning)

### What Is an AI Agent?

In everyday language, "AI" usually means a program that can recognize patterns (like the image models from Chapter 2). An **AI agent** goes further — it is a program that can make decisions about what to do next based on what it finds.

A regular AI model works like a vending machine: you put in an image, you get out a prediction. An agent works more like a doctor: it receives information, decides what additional tests to run, interprets the results, and forms a conclusion.

The Imaging Intelligence Agent does not just detect hemorrhages — it detects them, decides whether they are critical, looks up prior scans for comparison, searches for similar cases, and generates a clinical report. Each step depends on the results of the previous step.

### LangGraph: The Framework

LangGraph is an open-source framework (MIT license) for building multi-step AI reasoning workflows. It is built on top of LangChain, which provides tools for connecting AI models to external data and services.

LangGraph uses a concept called a **StateGraph** — essentially a flowchart that the agent follows. The flowchart has:

- **Nodes:** Individual steps. Each node is a Python function that does one specific thing.
- **Edges:** Connections between steps. They define the order of operations.
- **Conditional edges:** Decision points where the agent takes different paths depending on the current state.

### The Agent's Nodes

The agent's StateGraph has four main nodes:

1. **triage_node:** Receives the study findings from the database. Determines severity (CRITICAL, URGENT, or ROUTINE). Decides which downstream analyses to run.

2. **longitudinal_node:** Retrieves prior studies for the same patient. Compares current findings to previous findings. Calculates changes (volume growth, new lesions, midline shift changes).

3. **population_node:** Searches the database for similar cases using embedding similarity (Chapter 9). Retrieves outcomes from similar patients. Provides population context.

4. **report_node:** Gathers all evidence from the previous nodes. Sends everything to the LLM (Chapter 11) with clinical guidelines. Generates a structured, evidence-grounded report.

### Conditional Routing

Not every study goes through every node. The agent uses conditional routing based on severity:

- **CRITICAL finding:** Runs all four nodes in sequence — full analysis pipeline. The patient gets longitudinal comparison, population context, and a comprehensive report.
- **ROUTINE finding:** May skip the longitudinal and population nodes and go straight to a brief report. No need for an extensive workup on a normal scan.

This saves compute time and focuses resources on the cases that need the most attention.

### State: The Agent's Working Memory

As the agent moves from node to node, it carries a **state** — a data structure that accumulates information. The state is defined as an `AgentState` TypedDict in Python:

- `study_id` — which study is being processed
- `findings` — list of findings from the AI models
- `prior_studies` — data from previous scans
- `similar_cases` — results from population search
- `severity` — current triage classification
- `report` — the final generated report
- `provenance` — tracking information for audit

Each node reads from the state, does its work, and writes its results back into the state. The next node picks up where the last one left off.

### Tools: What the Agent Can Do

The agent has access to **tools** — functions it can call to interact with external systems:

- **query_findings:** Search the PostgreSQL database for findings.
- **search_similar_studies:** Use embedding similarity to find matching cases.
- **calculate_vdt:** Compute volume doubling time for nodules.
- **retrieve_guidelines:** Look up clinical guidelines (Lung-RADS, Brain Trauma Foundation).
- **generate_report:** Call the LLM to write a clinical report.

These tools use MCP (Model Context Protocol) — a standardized way for AI agents to declare and use tools. MCP defines how a tool is described (its name, parameters, and expected output) so that any compatible agent can discover and use it.

### Why It Matters

The agent is not just a collection of AI models running independently. It is a reasoning system that orchestrates multiple models, databases, and services into a coherent clinical workflow. Understanding the graph structure helps you explain how the agent makes decisions — and why different cases trigger different levels of analysis.

### Test Yourself — Chapter 8

1. What is the difference between an AI model and an AI agent?
2. Name the four nodes in the agent's StateGraph and what each does.
3. Why does the agent skip some nodes for ROUTINE findings?
4. What is MCP and why is it useful?

---

## Chapter 9 — Embeddings and Similarity Search

### What Is an Embedding?

An embedding is a way to represent something complex (like a medical image) as a list of numbers. Not just pixel values — instead, numbers that capture the "essence" or "meaning" of the image.

Think of it this way: imagine you could describe a painting with GPS coordinates — not the location where it was painted, but a point in an abstract space where similar paintings are nearby. A landscape of a sunset and another landscape of a sunset would have coordinates close together. A landscape and a portrait of a person would be far apart.

An embedding does exactly this for medical images. The agent turns each study into a list of 384 numbers (a 384-dimensional vector). Studies with similar characteristics — same body part, similar findings, similar patient demographics — end up with similar vectors.

### BiomedCLIP: The Embedding Model

The specific model that creates embeddings is called **BiomedCLIP** (microsoft/BiomedCLIP-PubMedBERT_256-vit_base_patch16_224). It was trained on millions of medical images paired with text descriptions. Because it learned from medical data, it understands medical concepts — it knows that two CT scans showing hemorrhage are more similar to each other than a CT scan showing hemorrhage and one showing a normal brain.

BiomedCLIP runs in the `imaging-embedding` container on port 8521.

### Cosine Similarity: Measuring Closeness

When you have two embeddings (two lists of 384 numbers), you need a way to measure how similar they are. The agent uses **cosine similarity**.

Think of two arrows pointing from the same starting point. If they point in almost the same direction, they are similar (cosine similarity close to 1.0). If they point in opposite directions, they are dissimilar (cosine similarity close to -1.0). If they are at right angles, they are unrelated (cosine similarity 0).

The math is straightforward but not important to memorize. What matters is: cosine similarity close to 1.0 = very similar. Close to 0 = unrelated.

### pgvector: Similarity Search in the Database

The agent stores embeddings in the PostgreSQL database using the **pgvector** extension. pgvector adds vector data types and similarity search functions to PostgreSQL.

Without pgvector, finding the 10 most similar studies out of 100,000 would require comparing the query vector to every single stored vector — 100,000 comparisons. With pgvector's **HNSW index** (Hierarchical Navigable Small World), the database can find approximate nearest neighbors much faster by organizing vectors into a navigable graph structure.

HNSW works like this analogy: imagine you are in a large city trying to find the restaurant most similar to your favorite one. Instead of visiting every restaurant in the city, you start at a well-connected hub, then jump to a more specific neighborhood, then to a specific block — narrowing down quickly without checking every option.

### Three Levels of Embeddings

The agent creates embeddings at three levels:

1. **Study-level:** One vector per study. Used for "find studies similar to this one."
2. **Series-level:** One vector per series. Used for "find series with similar imaging characteristics."
3. **Lesion-level:** One vector per individual lesion/finding. Used for "find nodules that look like this one."

### Hybrid Queries: Combining Filters with Similarity

The real power comes from combining traditional database queries with similarity search. Examples:

- "Find the 10 CT chest studies most similar to this one that also had Lung-RADS 4A or higher."
- "Find all hemorrhage cases similar to this patient's scan where the patient was over 65."
- "Find growing nodules in patients with similar imaging phenotype who are APOE4 carriers."

This is called a **hybrid query** — it uses both SQL filtering (traditional) and vector similarity (AI-powered) in the same search.

### Why It Matters

Embeddings and similarity search enable "patients like this" queries. When a radiologist sees a complex case, the agent can instantly retrieve similar historical cases with known outcomes. This provides evidence for clinical decisions — "10 similar cases were found; 7 had benign outcomes, 3 required intervention." That context is vastly more useful than reading one scan in isolation.

### Test Yourself — Chapter 9

1. In your own words, what is an embedding?
2. What does cosine similarity measure?
3. What is the HNSW index and why is it needed?
4. Give an example of a hybrid query that combines filtering and similarity search.

---

## Chapter 10 — Talking to Hospital Systems (FHIR and Clinical Integration)

### What Is an EHR?

An EHR (Electronic Health Record) is the hospital's digital patient chart. It contains everything about a patient: demographics, lab results, medications, doctor's notes, imaging reports, and more. Epic, Cerner, and MEDITECH are common EHR systems.

The agent needs to send its findings to the EHR so that the patient's care team can see the results alongside all other clinical information.

### What Is FHIR?

FHIR (Fast Healthcare Interoperability Resources, pronounced "fire") is a modern standard for exchanging healthcare data. It defines how different healthcare systems can share information using web technology (HTTP, JSON).

Think of FHIR like a shared language between hospital systems. Without FHIR, every EHR, every lab system, and every imaging system would need a custom translation for every other system. With FHIR, they all speak the same language.

### FHIR Resources: Building Blocks

FHIR organizes healthcare data into **Resources** — standardized building blocks. The agent uses several:

| Resource | What It Represents |
|---|---|
| Patient | The person being imaged |
| ImagingStudy | Reference to the DICOM study |
| DiagnosticReport | The agent's findings and conclusions |
| Observation | Individual measurements (volume, diameter, shift) |

When the agent finishes analyzing a study, it creates a **DiagnosticReport** resource that links to the patient, the imaging study, and individual observations. This report is sent to the EHR via the `imaging-fhir-publisher` container (port 8523).

### Coding Systems: SNOMED CT and LOINC

Medical findings need to be coded in a standard way so that different systems interpret them identically.

**SNOMED CT** (Systematized Nomenclature of Medicine — Clinical Terms) is a coding system for clinical concepts. For example:
- "Intracranial hemorrhage" has a specific SNOMED code.
- "Pulmonary nodule" has a different specific code.

**LOINC** (Logical Observation Identifiers Names and Codes) is a coding system for lab tests and observations. For example:
- "Volume of hemorrhage" has a LOINC code.
- "Diameter of pulmonary nodule" has a different LOINC code.

The agent includes these codes in its FHIR output so that any EHR system can understand the findings without ambiguity.

### PACS: Where Radiologists View Images

PACS (Picture Archiving and Communication System) is the specialized system where radiologists view and interpret medical images. The agent sends two types of output to PACS:

1. **DICOM SR (Structured Report):** Contains the findings in a machine-readable format, viewable alongside the original images.
2. **GSPS (Grayscale Softcopy Presentation State):** Visual overlays that draw annotations, arrows, and labels directly on the images. When the radiologist opens the study, they see the AI's markings on the relevant slices.

These are pushed to PACS via DICOMweb STOW-RS.

### The Worklist: The Radiologist's To-Do List

The agent maintains a **worklist** — a prioritized queue of studies awaiting radiologist review. Each entry has a priority:

- **CRITICAL:** Jump to the top. Red flag. Immediate review needed.
- **URGENT:** High priority. Review within the hour.
- **ROUTINE:** Standard queue. Review when available.

The worklist is stored in the `worklist_entries` table and displayed through the Streamlit portal (port 8525). The agent's triage node assigns priority based on the clinical workflows' severity thresholds.

### Clinician-in-the-Loop

This is a fundamental design principle: the agent is decision **support**, not autonomous diagnosis. Every output is a recommendation for a human expert to review. The radiologist sees the agent's findings, evaluates them, and makes the final clinical decision.

The agent never says "this patient has a hemorrhage." It says "hemorrhage detected with 95% confidence, estimated volume 25.3 mL — recommend immediate review." The radiologist confirms or overrides.

### Why It Matters

An AI system that cannot communicate with existing hospital systems is useless in practice. FHIR, DICOM SR, GSPS, and worklist integration are what make the agent a practical clinical tool rather than a laboratory experiment. Understanding this integration layer is critical for explaining the agent's real-world value.

### Test Yourself — Chapter 10

1. What is the difference between an EHR and PACS?
2. What FHIR resource does the agent create to report its findings?
3. Why are coding systems like SNOMED CT important?
4. What does "clinician-in-the-loop" mean and why is it essential?

---

## Chapter 11 — RAG and Large Language Models

### What Is a Large Language Model (LLM)?

A Large Language Model is an AI that can read, understand, and generate text. ChatGPT is the most famous example. An LLM has been trained on enormous amounts of text and learned the patterns of language — grammar, facts, reasoning styles.

The agent uses an LLM to generate clinical reports. Instead of producing a rigid template, the LLM can write natural, context-aware text that incorporates all the evidence the agent has gathered.

### NVIDIA NIM: Running the LLM

The LLM runs inside the `imaging-nim-llm` container on port 8520. **NVIDIA NIM** (NVIDIA Inference Microservice) is the software that serves the model — it handles loading the model into GPU memory, processing requests, and returning responses.

The specific model used is **Meta Llama 3 8B Instruct** (container image: `nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark`). "8B" means it has 8 billion parameters — the numbers that define the model's learned patterns. "Instruct" means it has been fine-tuned to follow instructions.

On the DGX Spark, at the desktop class, there is no NVAIE license cost for running NIM.

### The Problem with LLMs: Hallucination

LLMs have a well-known weakness: they can make things up. This is called **hallucination**. An LLM might confidently state a measurement that does not exist, cite a guideline that was never published, or invent a patient history.

In a hospital setting, hallucination is dangerous. A report that says "hemorrhage volume 15 mL" when the actual volume is 35 mL could lead to a wrong treatment decision. The solution is RAG.

### What Is RAG?

RAG stands for **Retrieval-Augmented Generation**. Instead of asking the LLM to generate a report from its training memory alone, you first **retrieve** all the relevant facts, then give those facts to the LLM as context, and ask it to **generate** a report based only on those facts.

Think of it this way: instead of asking someone to write a book report from memory (where they might misremember details), you give them the book, their notes, and a highlighter, and ask them to write the report with the book open in front of them.

### The RAG Pipeline Step by Step

When the agent's report_node runs, here is what happens:

1. **Retrieve findings from the database:** "This study has 2 findings: hemorrhage with 95% confidence, midline shift of 8.2 mm."

2. **Retrieve prior measurements for comparison:** "Previous CT from 3 months ago showed hemorrhage volume of 18 mL. Current volume is 25.3 mL. Volume has increased by 40%."

3. **Retrieve similar cases and their outcomes:** "10 similar cases found via embedding search. 7 required surgical intervention. 3 were managed conservatively."

4. **Retrieve relevant clinical guidelines:** "Brain Trauma Foundation: hemorrhage > 30 mL with midline shift > 5 mm is a surgical indication."

5. **Construct the prompt:** All retrieved evidence is assembled into a structured prompt with instructions: "Based on the following evidence, generate a clinical report. Only include facts supported by the evidence. Do not add information not present in the evidence."

6. **Send to the LLM:** The prompt goes to the NIM-served Llama 3 model.

7. **Receive the grounded report:** The LLM generates a structured report that references the retrieved evidence. Every claim in the report can be traced back to a specific finding, measurement, or guideline.

### Evidence Grounding

The key concept is **evidence grounding** — every statement in the report is backed by data. The LLM is not generating new information; it is summarizing and contextualizing information that was already computed, measured, and retrieved. This dramatically reduces the risk of hallucination.

### Cross-Modal Context Enrichment

The RAG pipeline can pull evidence from multiple sources — not just imaging:

- **Imaging findings:** What the AI models detected.
- **Longitudinal data:** How findings have changed over time.
- **Population data:** What happened to similar patients.
- **Clinical guidelines:** What published standards recommend.
- **Genomic data** (when available): Molecular characteristics from Parabricks.

By combining evidence from all these sources, the LLM can generate a comprehensive report that no single data source could provide alone.

### Why It Matters

Reports with evidence are trustworthy. Hallucinated reports are dangerous. RAG is the mechanism that prevents the LLM from inventing facts. Understanding RAG is critical for explaining why the agent's reports are reliable — and what safeguards prevent it from producing misleading output.

### Test Yourself — Chapter 11

1. What is hallucination in the context of LLMs, and why is it dangerous in healthcare?
2. What does RAG stand for, and what problem does it solve?
3. List the five types of evidence the RAG pipeline retrieves before generating a report.
4. What is evidence grounding?

---

## Chapter 12 — Monitoring and Observability

### What Monitoring Means

Monitoring means continuously watching a system to make sure it is working correctly. For the Imaging Intelligence Agent, monitoring answers questions like:

- Is the system processing studies or is it stuck?
- How long is each workflow taking?
- Is the GPU running at full capacity or sitting idle?
- Are there errors that need attention?

### Metrics: Numbers About System Health

A **metric** is a single measurement tracked over time. Examples:

- **inference_latency_seconds:** How long GPU inference takes for each study.
- **pipeline_throughput:** How many studies are processed per hour.
- **error_rate:** How many processing attempts failed.
- **gpu_utilization_percent:** What percentage of the GPU's capacity is being used.
- **gpu_temperature_celsius:** How hot the GPU is running.

### Prometheus: The Metrics Collector

**Prometheus** is an open-source monitoring system that collects metrics at regular intervals (typically every 15-30 seconds). It runs in the `imaging-prometheus` container on port 9099.

Prometheus works by "scraping" — it sends HTTP requests to each container's metrics endpoint and collects the numbers. Each container exposes its metrics in a standard format that Prometheus understands.

Think of Prometheus as a hall monitor who walks around every 15 seconds, checks on every classroom, and writes down what they find.

### Grafana: The Dashboard

**Grafana** is an open-source visualization tool that turns Prometheus metrics into charts and dashboards. It runs in the `imaging-grafana` container on port 3000.

While Prometheus stores the raw numbers, Grafana makes them visual. You can see:

- A line chart of inference latency over the past 24 hours.
- A gauge showing current GPU utilization.
- A bar chart of studies processed per hour.
- A table of recent errors.

Grafana is what you would open on a screen in the hospital to see the agent's health at a glance.

### DCGM Exporter: GPU-Specific Metrics

**DCGM** (Data Center GPU Manager) is an NVIDIA tool specifically for monitoring GPU health. The `imaging-dcgm` container (port 9400) exports GPU metrics in Prometheus format:

- GPU utilization (compute percentage)
- GPU memory usage
- GPU temperature
- Power consumption
- Error counts

These metrics are critical because the GPU is the heart of the system. If the GPU is overheating, running out of memory, or encountering errors, the agent cannot process studies.

### Alerting

Grafana can be configured to send **alerts** — notifications that fire when a metric crosses a threshold. Examples:

- "GPU temperature exceeded 85°C" → alert to the IT team.
- "Pipeline error rate exceeded 5% in the last hour" → alert to the support team.
- "No studies processed in the last 30 minutes during business hours" → alert that the system may be down.

Alerts ensure that problems are caught quickly, even if no one is watching the dashboard.

### Provenance: The Complete Audit Trail

**Provenance** is a record of exactly how a result was produced. For every finding the agent generates, the provenance includes:

- **Model ID and version:** Which specific model produced this result.
- **Inference parameters:** What settings were used (thresholds, etc.).
- **Processing duration:** How long inference took.
- **Input data lineage:** Which specific DICOM images (by UID) were processed.
- **Timestamps:** When each step started and finished.
- **Operator approvals:** If a human reviewed and approved the result.
- **Predetermined change control plans:** Pre-approved rules for model updates (relevant for FDA compliance).

All provenance data is stored in the `provenance` table in PostgreSQL.

### Reproducibility

**Reproducibility** means the ability to re-run the exact same analysis and get the exact same result. The agent achieves this by:

1. Storing the original DICOM images as immutable files (never modified).
2. Recording every processing parameter in provenance.
3. Versioning every model.

If a question arises about a result from six months ago, you can reload the same images, use the same model version, apply the same parameters, and verify that the result is identical.

### Why It Matters

In healthcare, you cannot just deploy a system and hope it works. You must prove it is working correctly, continuously. Monitoring provides real-time visibility. Provenance provides a complete audit trail. Reproducibility provides verification. Together, they create the trust needed for clinical adoption.

### Test Yourself — Chapter 12

1. What is the difference between Prometheus and Grafana?
2. What does DCGM monitor and why is it important?
3. What information does a provenance record include?
4. Why is reproducibility essential in a medical AI system?

---

## Chapter 13 — The HCLS AI Factory (The Bigger Picture)

### What the HCLS AI Factory Is

The Imaging Intelligence Agent is not a standalone system. It is one piece of a larger platform called the **HCLS AI Factory** (Health and Life Sciences AI Factory). The AI Factory is a collection of specialized AI agents that work together on shared NVIDIA DGX hardware to provide comprehensive patient analysis.

Think of it like a hospital with multiple departments. The imaging department (this agent) handles scans. The genetics department handles DNA analysis. The pharmacy department handles drug matching. Each department is specialized, but they share patient records and collaborate on complex cases.

### Cross-Modal Integration

"Cross-modal" means combining information from different types of data. The HCLS AI Factory connects imaging to other domains:

**Imaging → Genomics (Parabricks):**
When the imaging agent finds a highly suspicious lung nodule (Lung-RADS 4B+), it automatically triggers the genomics pipeline. NVIDIA Parabricks performs whole-genome sequencing analysis — processing 30x WGS in approximately 10 minutes on DGX (compared to approximately 30 hours on CPU alone). This produces a molecular profile of the tumor: which mutations are present, which genes are affected, which treatments might work.

Up to 50% lower compute cost compared to CPU-only genomics processing.

**Imaging → Drug Discovery (BioNeMo):**
NVIDIA BioNeMo is used for drug discovery and molecular modeling. The imaging agent's quantitative measurements (tumor size, growth rate, response patterns) become endpoints in clinical trials. Drug candidates can be scored against the imaging phenotype — "does this drug work for tumors that look like this?"

BioNeMo has over 200 adopters including Eli Lilly, Astellas, Insilico Medicine, and Recursion.

**Imaging → Clinical Reasoning (NIM LLM):**
The RAG pipeline (Chapter 11) uses NVIDIA NIM to serve the LLM. Clinical reports are grounded in evidence from imaging findings, guidelines, and prior outcomes.

**Imaging → Longitudinal Care:**
Tracking patients across multiple time points — months or years of scans — to detect meaningful changes automatically. The agent's database retains all prior measurements, enabling trend analysis without manual chart review.

**Imaging → Outcomes:**
Linking imaging patterns to patient outcomes through cohort retrieval. "Of 50 patients with similar imaging findings, what were their outcomes?" This feeds back into clinical decision-making and care pathway optimization.

### The Biomarker Agent

Another agent in the AI Factory is the **Biomarker Intelligence Agent**. It combines genomic biomarkers with imaging biomarkers for a combined phenotype profile. For example: a patient has a specific genetic mutation (from Parabricks) AND a specific tumor morphology (from the imaging agent). Together, these might predict response to a specific therapy more accurately than either alone.

Both agents run on the same DGX Spark hardware and share data through the common database and data models.

### Scaling: From Desk to Data Center

| Phase | Hardware | Price | Capabilities |
|---|---|---|---|
| 1 — Proof Build | DGX Spark | $4,699 | 1-2 workflows, single user, local |
| 2 — Departmental | 1-2x DGX B200 | $500K-$1M | All workflows, multi-user, PACS-integrated |
| 3 — Multi-Site | 4-8x DGX B200 + InfiniBand | $2M-$4M | Multiple hospitals, federated learning |
| 4 — AI Factory | DGX SuperPOD | $7M-$60M+ | Thousands of concurrent studies, complete platform |

**NVAIE (NVIDIA AI Enterprise) licensing scales with GPUs:**
- Phase 1: $0 (desktop-class)
- Phase 2: $36K-$72K/year (8-16 GPUs at $4,500/GPU/year)
- Phase 3: $144K-$288K/year (32-64 GPUs)
- Phase 4: $576K-$1.15M/year (128-256 GPUs)

The software is identical at every tier. The same containers, the same database schema, the same agent code. Only the hardware changes.

### NVIDIA FLARE: Federated Learning

At Phase 3 and beyond, **NVIDIA FLARE** (Federated Learning Application Runtime Environment) enables multiple hospitals to collaboratively improve the AI models without sharing patient data.

How it works:
1. Each hospital trains (or fine-tunes) a model on their own local data.
2. Instead of sharing patient data, each hospital shares only the model's learned parameters (numbers, not patient information).
3. A central server aggregates the updates from all hospitals into an improved global model.
4. The improved model is sent back to each hospital.

This way, the model learns from the collective experience of many hospitals while patient data never leaves institutional control. NVIDIA FLARE is free (Apache 2.0 license).

### Why It Matters

The Imaging Intelligence Agent is designed from day one to be part of something bigger. Understanding the AI Factory context helps you explain:
- Why the agent's data model includes genomics triggers.
- Why the architecture uses standardized containers and data formats.
- Why the $4,699 DGX Spark is a starting point, not the destination.
- How imaging connects to the full spectrum of precision medicine.

### Test Yourself — Chapter 13

1. What happens when the imaging agent detects a Lung-RADS 4B+ nodule?
2. What is federated learning and how does NVIDIA FLARE protect patient privacy?
3. Name three cross-modal integrations the AI Factory enables.
4. Why is the DGX Spark described as a "proof build" rather than a production system?

---

## Chapter 14 — Trust, Safety, and Regulation

### Why Trust Matters

No hospital will adopt an AI system they do not trust. Trust in medical AI comes from three pillars: transparency (you can see what the system did), accountability (someone is responsible for the outcome), and reliability (the system performs consistently).

### FDA AI/ML SaMD Framework

The FDA (Food and Drug Administration) regulates medical AI as **Software as a Medical Device (SaMD)**. The FDA's AI/ML framework addresses a unique challenge: unlike a physical medical device that stays the same after manufacturing, AI models can be updated and improved.

Key concepts:

**Intended use:** What the software is designed to do. The imaging agent's intended use is decision support for radiologists — not autonomous diagnosis.

**Risk categorization:** Higher-risk applications (like triage decisions that affect patient care timing) face stricter requirements than lower-risk applications (like quality checks).

**Predetermined Change Control Plans (PCCP):** Pre-approved rules for when and how the AI can be updated without requiring a new regulatory submission each time. For example: "If the model's sensitivity on a validation dataset drops below 90%, it must be retrained and re-validated before deployment."

### Decision Support vs. Autonomous Diagnosis

This is a critical distinction:

- **Decision support:** The AI assists the clinician. It provides findings, measurements, and recommendations. The clinician makes the final diagnosis and treatment decision. The clinician is accountable.
- **Autonomous diagnosis:** The AI makes the diagnosis without human review. Currently, very few medical AI systems are authorized for fully autonomous operation.

The Imaging Intelligence Agent is designed as decision support. Every output includes the phrase "recommend review" or similar language. The radiologist always has the final say.

### Audit Trails: Traceability

Every output the agent produces can be traced back to:
- The exact model version (including a hash of the model weights).
- The exact input data (DICOM UIDs of every image processed).
- The exact configuration (thresholds, parameters, routing decisions).
- The exact timestamp of each processing step.
- Any human approvals or overrides.

This audit trail is stored in the provenance table (Chapter 12) and is immutable — once written, it cannot be changed or deleted. If a regulator or a legal proceeding asks "how did the system reach this conclusion on January 15?", the answer is fully reconstructable.

### Immutability: Preserving the Evidence

The original DICOM images are never modified. They are stored as immutable files — evidence that can always be re-examined. The agent's derived artifacts (segmentation masks, measurements, reports) are stored alongside the originals but never overwrite them.

Think of it like a crime scene: you take photos (original DICOM), you write analysis notes (derived artifacts), but you never alter the scene itself.

### Patient Data Security

The agent is designed with security principles:

- **Data stays local:** Patient data remains within the institution's control. It is not sent to cloud services or external servers (unless explicitly configured for multi-site deployment).
- **Least-privilege access:** Each container only has access to the data it needs. The monitoring containers cannot read patient data. The inference containers cannot modify the database directly.
- **Tenant isolation:** In multi-department deployments, different departments cannot see each other's data. Access is controlled at the database level.

### Controlled Rollouts

When a new model version is deployed, it does not replace the old version overnight. The agent supports controlled rollouts:

1. New model version is deployed alongside the existing version.
2. Both versions process the same studies in parallel.
3. Results are compared to verify the new version performs at least as well.
4. Only after validation is the new version promoted to production.
5. Previous version outputs are preserved — never deleted.

This approach ensures that model updates improve quality without introducing regressions.

### Why It Matters

Without trust, governance, and regulatory alignment, the most brilliant AI system will sit unused. Understanding these concepts lets you address the most common concern from hospital leadership: "Can we trust it?" The answer is not "yes, trust the AI" — the answer is "here is exactly how every result is produced, audited, and safeguarded."

### Test Yourself — Chapter 14

1. What is the difference between decision support and autonomous diagnosis?
2. What is a Predetermined Change Control Plan?
3. Why are original DICOM images stored as immutable files?
4. How does tenant isolation protect patient data in a multi-department deployment?

---

## Test Yourself — Comprehensive Review

These 30 questions cover all 14 chapters. Try answering them in your own words before checking the answer key.

### Questions

1. What imaging modality uses magnets and radio waves instead of radiation?
2. What is the difference between a pixel and a voxel?
3. Name the three fundamental AI operations used in the agent (yes/no, where, exact shape).
4. What does DICOM stand for?
5. What are the three types of UIDs in DICOM and what level does each identify?
6. What is DICOMweb and how does it differ from DIMSE?
7. How many cores does a GPU typically have compared to a CPU?
8. What does "unified memory" mean on the DGX Spark?
9. What is GPUDirect Storage?
10. What is a Docker container and how does it differ from a Docker image?
11. How many containers does the agent run, and why are they separated?
12. What is a Docker volume and why is it needed?
13. Name the seven tables in the agent's database.
14. What is SQL used for?
15. What is a database index and what real-world object is it analogous to?
16. What five conditions does the CXR workflow detect?
17. At what hemorrhage volume does the CT Head workflow classify a finding as CRITICAL?
18. What is Volume Doubling Time (VDT) and what does a short VDT indicate?
19. How does the MRI workflow determine disease activity classification?
20. What is a StateGraph and what are its main components?
21. Name the four nodes in the agent's reasoning graph.
22. What is an embedding?
23. What does cosine similarity measure?
24. What is a hybrid query?
25. What is FHIR and what problem does it solve?
26. What is the difference between an EHR and PACS?
27. What does RAG stand for and what problem does it solve?
28. What is the difference between Prometheus and Grafana?
29. What is federated learning and how does NVIDIA FLARE implement it?
30. What is the difference between decision support and autonomous diagnosis?

### Answer Key

1. **MRI** (Magnetic Resonance Imaging) uses magnets and radio waves. No radiation involved.

2. A **pixel** is a 2D picture element (one square in a flat image). A **voxel** is a 3D volume element (one cube in a 3D scan like CT or MRI).

3. **Classification** (is there a problem?), **detection** (where is it?), and **segmentation** (what exact shape is it?).

4. **Digital Imaging and Communications in Medicine.**

5. **Study Instance UID** (one scan session), **Series Instance UID** (one set of images within a study), **SOP Instance UID** (one individual image file).

6. **DICOMweb** uses standard web protocols (HTTP/REST) to access DICOM data. **DIMSE** uses a traditional networking protocol. DICOMweb is the modern approach; DIMSE is the legacy approach. Both are supported.

7. A GPU has **thousands** of simple cores. A CPU has **8-64** powerful cores. The GPU trades per-core capability for massive parallelism.

8. The 128 GB of memory is **shared between CPU and GPU** — one pool, no copying back and forth. In most computers, CPU and GPU have separate memory pools.

9. GPUDirect Storage lets data flow **directly from NVMe drive to GPU memory**, skipping the CPU. This eliminates a copy step and speeds up data loading.

10. A **Docker image** is a read-only template (like a recipe). A **Docker container** is a running instance of that image (like the cooked meal). You can run many containers from one image.

11. **11 containers.** They are separated so each is a specialist — if one fails, the others keep running. It also allows independent updates and scaling.

12. A Docker volume is **persistent storage** that survives container restarts. Without volumes, all data inside a container is lost when it stops.

13. **studies, series, findings, measurements, embeddings, provenance, worklist_entries.**

14. SQL is the language for **querying (asking questions of) and modifying data** in a relational database.

15. A database index is a data structure that **speeds up searches** — like the **index at the back of a textbook** that lets you jump to the right page instead of reading every page.

16. **Pneumothorax, consolidation, pleural effusion, cardiomegaly, and fracture.**

17. **> 30 mL** (any midline shift) = CRITICAL.

18. VDT measures **how fast a nodule's volume doubles**. A short VDT (< 400 days) means rapid growth, which is suspicious for malignancy.

19. **0 new/enlarging lesions = Stable, 1-2 = Active, 3+ = Highly Active.**

20. A StateGraph is a **flowchart** with **nodes** (individual processing steps), **edges** (connections between steps), and **conditional edges** (decision points).

21. **triage_node, longitudinal_node, population_node, report_node.**

22. An embedding is a **list of numbers** (vector) that captures the "essence" of an image. Similar images have similar embeddings.

23. Cosine similarity measures **how similar two vectors are** by comparing the angle between them. Close to 1.0 = very similar. Close to 0 = unrelated.

24. A hybrid query **combines traditional SQL filtering with vector similarity search** in the same query. Example: "CT chest studies with Lung-RADS 4A+ that are most similar to this scan."

25. FHIR (Fast Healthcare Interoperability Resources) is a **standard for exchanging healthcare data** using web technology. It solves the problem of different healthcare systems being unable to share information.

26. An **EHR** (Electronic Health Record) stores all patient clinical data. **PACS** (Picture Archiving and Communication System) specifically stores and displays medical images. Radiologists use PACS; the broader care team uses the EHR.

27. RAG = **Retrieval-Augmented Generation**. It solves LLM **hallucination** by first retrieving factual evidence, then asking the LLM to generate text based only on that evidence.

28. **Prometheus** collects and stores metrics (the numbers). **Grafana** visualizes those metrics as charts and dashboards. Prometheus is the data source; Grafana is the display.

29. Federated learning lets **multiple hospitals collaboratively train AI models without sharing patient data**. Each hospital trains locally and shares only model parameters. NVIDIA FLARE coordinates this process. A central server aggregates the updates.

30. **Decision support:** AI assists the clinician, who makes the final decision and is accountable. **Autonomous diagnosis:** AI makes the diagnosis without human review. The Imaging Intelligence Agent is decision support only.

---

## Glossary

| Term | Definition |
|---|---|
| Agent | An AI program that can make decisions about what to do next based on what it finds, rather than just producing a single output. |
| ARM64 / aarch64 | A processor architecture used by the DGX Spark. Different from the x86_64 architecture in most PCs. All software must be compiled for the correct architecture. |
| BiomedCLIP | A medical image embedding model that converts images into 384-dimensional vectors for similarity search. |
| BioNeMo | NVIDIA's platform for drug discovery and molecular modeling. Part of the HCLS AI Factory. |
| Bounding Box | A rectangle drawn around a detected object in an image (used in detection). |
| Classification | An AI task that answers "is there a problem?" with a yes/no (or multi-category) answer and a confidence score. |
| Clinician-in-the-Loop | A design principle where the AI assists but never replaces the human expert's judgment and accountability. |
| Confidence Score | A number between 0 and 1 representing how sure the AI model is about its prediction. |
| Container | A self-contained software package with everything needed to run — code, libraries, configuration, and a minimal OS. |
| Cosine Similarity | A mathematical measure of how similar two vectors are, based on the angle between them. 1.0 = identical direction. |
| CT (Computed Tomography) | An imaging technique that takes X-rays from many angles to create a 3D picture of the body. |
| DCGM | Data Center GPU Manager. An NVIDIA tool for monitoring GPU health metrics. |
| Deep Learning | Machine learning using neural networks with many layers (deep networks). |
| DenseNet-121 | A 121-layer classification neural network used for chest X-ray analysis and hemorrhage detection. |
| Detection | An AI task that answers "where is the problem?" by drawing a bounding box around it. |
| DGX Spark | NVIDIA's $4,699 desktop-class GPU computer with ARM64 Grace CPU and Blackwell GB10 GPU. |
| DICOM | Digital Imaging and Communications in Medicine. The universal standard for medical image storage and communication. |
| DICOM SR | DICOM Structured Report. A DICOM-formatted document containing findings and measurements. |
| DICOMweb | A modern web-based protocol for storing, retrieving, and searching DICOM images. |
| DIMSE | The traditional networking protocol for DICOM communication (C-STORE, C-FIND, etc.). |
| Docker | A tool for building and running containers. |
| docker-compose | A tool for defining and running multi-container applications from a single configuration file. |
| Dockerfile | A text file with instructions for building a Docker image. |
| Embedding | A list of numbers (vector) that represents the "essence" of a complex object like an image, enabling similarity comparison. |
| EHR | Electronic Health Record. The hospital's comprehensive digital patient chart. |
| Evidence Grounding | Ensuring every claim in an AI-generated report is supported by retrieved factual data. |
| FDA | Food and Drug Administration. The U.S. agency that regulates medical devices, including AI software. |
| Federated Learning | A technique where multiple institutions collaboratively train AI models without sharing raw data. |
| FHIR | Fast Healthcare Interoperability Resources. A standard for exchanging healthcare data via web protocols. |
| FLAIR | Fluid-Attenuated Inversion Recovery. An MRI sequence that highlights brain lesions by suppressing normal fluid signals. |
| FLARE | NVIDIA's Federated Learning Application Runtime Environment. Open-source (Apache 2.0). |
| GPU | Graphics Processing Unit. A processor with thousands of cores optimized for parallel computation. |
| GPUDirect Storage | A technology that moves data directly from NVMe storage to GPU memory, bypassing the CPU. |
| GradCAM | Gradient-weighted Class Activation Mapping. A technique that highlights which image regions influenced an AI's decision. |
| Grafana | An open-source visualization tool that displays metrics as charts and dashboards. |
| GSPS | Grayscale Softcopy Presentation State. DICOM overlays that draw annotations on images. |
| Hallucination | When an LLM generates text that is factually incorrect or fabricated. |
| Health Check | An automatic test Docker runs to verify a container is responding correctly. |
| highdicom | A Python library for creating DICOM Structured Reports and other DICOM objects. |
| HNSW | Hierarchical Navigable Small World. An index structure for fast approximate nearest-neighbor search in vector databases. |
| Hybrid Query | A database query that combines traditional SQL filtering with vector similarity search. |
| Immutable | Cannot be changed after creation. Original DICOM images are stored as immutable files. |
| Index (database) | A data structure that speeds up searches, like the index at the back of a textbook. |
| Inference | Running a trained AI model on new data to produce predictions. The "taking the test" phase. |
| LangChain | An open-source framework for building applications with LLMs. |
| LangGraph | An open-source framework (built on LangChain) for building multi-step AI reasoning workflows. |
| LLM | Large Language Model. An AI that reads, understands, and generates text. |
| LOINC | Logical Observation Identifiers Names and Codes. A coding system for lab tests and observations. |
| Lung-RADS | Lung CT Screening Reporting and Data System. A standardized scoring system for lung nodule risk. |
| MAP | MONAI Deploy Application Package. A containerized, portable inference pipeline. |
| MCP | Model Context Protocol. A standardized way for AI agents to declare and use tools. |
| Metric | A single numerical measurement tracked over time (e.g., inference latency, GPU temperature). |
| Midline Shift | The displacement of the brain's center line, measured in millimeters. Caused by swelling or hemorrhage. |
| MONAI | Medical Open Network for AI. NVIDIA's open-source framework for medical image analysis (Apache 2.0). |
| MRI | Magnetic Resonance Imaging. Uses magnets and radio waves (no radiation) to create detailed soft tissue images. |
| MS (Multiple Sclerosis) | A disease where the immune system attacks nerve fiber insulation (myelin) in the brain and spinal cord. |
| Neural Network | A computer program that learns patterns from examples, organized in layers of connected mathematical units. |
| NIM | NVIDIA Inference Microservice. Software for serving AI models in production. |
| Node (LangGraph) | One step in an agent's reasoning graph — a function that processes state and passes results forward. |
| NVMe | Non-Volatile Memory Express. A fast solid-state storage technology. |
| NVAIE | NVIDIA AI Enterprise. A software license ($4,500/GPU/year) required for production-scale NIM, Parabricks, and BioNeMo. |
| Orthanc | An open-source DICOM server (GPLv3) used to store and serve medical images. |
| PACS | Picture Archiving and Communication System. Where radiologists view and interpret medical images. |
| Parabricks | NVIDIA's GPU-accelerated genomics analysis platform. |
| PCCP | Predetermined Change Control Plan. Pre-approved rules for AI model updates under FDA regulation. |
| pgvector | A PostgreSQL extension that adds vector data types and similarity search. |
| Pixel | Picture element. One point in a 2D image, represented as a number. |
| PostgreSQL | An open-source relational database used by the agent. Runs on port 5432. |
| Prometheus | An open-source monitoring system that collects metrics at regular intervals. |
| Provenance | A complete record of how a result was produced — model version, input data, parameters, timestamps. |
| pydicom | A Python library (MIT license) for reading and writing DICOM files. |
| QIDO-RS | DICOMweb service for searching studies on a server. |
| RAG | Retrieval-Augmented Generation. Retrieving facts first, then generating text based on those facts to prevent hallucination. |
| Radiologist | A doctor who specializes in reading and interpreting medical images. |
| Reproducibility | The ability to re-run an analysis with the same inputs and get the same result. |
| RetinaNet | A neural network architecture used for object detection (finding and boxing nodules). |
| SaMD | Software as a Medical Device. The FDA's classification for medical AI software. |
| Segmentation | An AI task that labels every pixel/voxel as "abnormality" or "normal," producing an exact shape mask. |
| Segmentation Mask | An image the same size as the original where each pixel is labeled by category. |
| SegResNet | A neural network architecture used for medical image segmentation. |
| Series | A set of images within a study, acquired with the same settings. |
| SNOMED CT | Systematized Nomenclature of Medicine — Clinical Terms. A coding system for clinical concepts. |
| SQL | Structured Query Language. The language for querying and modifying data in a relational database. |
| StateGraph | A flowchart-like structure in LangGraph with nodes (steps), edges (connections), and conditional routing. |
| STOW-RS | DICOMweb service for uploading (storing) images to a server. |
| Study | The complete collection of images from one scan session. |
| SwinUNETR | A transformer-based neural network architecture for medical image segmentation. |
| Training | Teaching an AI model by showing it labeled examples. The "studying for the test" phase. |
| U-Net (3D) | A neural network architecture widely used for medical image segmentation, operating on 3D volumes. |
| UID | Unique Identifier. A globally unique number assigned to every DICOM study, series, and image. |
| Unified Memory | A memory architecture where CPU and GPU share the same physical memory pool (as on DGX Spark). |
| VDT | Volume Doubling Time. How long it takes a nodule to double in volume. VDT < 400 days is suspicious. |
| View (database) | A saved SQL query that can be treated as if it were a table. |
| Volume | (1) Docker volume: persistent storage. (2) Medical: a 3D measurement in milliliters or cubic millimeters. |
| Voxel | Volume element. One point in a 3D image, represented as a number. |
| WADO-RS | DICOMweb service for retrieving (downloading) images from a server. |
| Worklist | A prioritized queue of studies awaiting radiologist review, ordered by urgency. |
| X-Ray | An imaging technique that passes radiation through the body to create a 2D shadow image. |

---

> *This document was created for the HCLS AI Factory — Imaging Intelligence Agent.*
> *Apache 2.0 License | Author: Adam Jones | February 2026*
