#!/usr/bin/env nextflow

/*
 * ========================================================
 *  Healthcare & Life Sciences Pipeline Orchestrator
 * ========================================================
 *  Unified workflow: Genomics → RAG Chat → Drug Discovery
 *
 *  Author: HLS Pipeline Team
 *  Version: 1.1.0
 *
 *  NOTE (ORCH-2, 2026-06-29): Nextflow is NOT installed on the DGX Spark and this
 *  workflow does not run there (a prior attempt crashed at JVM/cgroup-v2 startup).
 *  The ACTIVE orchestrator is `run_pipeline.py` in this directory. This .nf is
 *  retained as the declarative reference / future-cluster target, not the live path.
 * ========================================================
 */

nextflow.enable.dsl = 2

// Pipeline metadata
def pipelineVersion = '1.1.0'
def pipelineName = 'HLS-Pipeline'

log.info """
╔═══════════════════════════════════════════════════════════════╗
║  Healthcare & Life Sciences Pipeline                          ║
║  ─────────────────────────────────────────────────────────── ║
║  Genomics → Target Discovery → Drug Design                    ║
╚═══════════════════════════════════════════════════════════════╝

  Pipeline:      ${pipelineName} v${pipelineVersion}
  Run ID:        ${workflow.runName}
  Session:       ${workflow.sessionId}

  Input:         ${params.input ?: 'Not specified'}
  Output:        ${params.outdir}
  Mode:          ${params.mode}

  Profile:       ${workflow.profile}
  Container:     ${workflow.containerEngine ?: 'None'}
──────────────────────────────────────────────────────────────────
"""

// ============================================================
//  Include modules
// ============================================================

include { GENOMICS_PIPELINE   } from './modules/genomics'
include { RAG_CHAT_PIPELINE   } from './modules/rag_chat'
include { DRUG_DISCOVERY      } from './modules/drug_discovery'
include { GENERATE_REPORT     } from './modules/reporting'

// ============================================================
//  Input validation
// ============================================================

def validateInputs() {
    if (params.mode == 'full' && !params.input) {
        error "Input FASTQ files required for full pipeline mode. Use --input"
    }
    if (params.mode == 'genomics_only' && !params.input) {
        error "Input FASTQ files required for genomics_only mode. Use --input"
    }
    if (params.mode == 'target' && !params.vcf) {
        error "VCF file required for target discovery mode. Use --vcf"
    }
    if (params.mode == 'drug' && !params.target) {
        error "Target hypothesis JSON required for drug discovery mode. Use --target"
    }
}

// ============================================================
//  Main workflow
// ============================================================

workflow {

    // Validate inputs based on mode
    validateInputs()

    // Channel definitions — pre-declare all to prevent undefined variable errors
    ch_versions  = Channel.empty()
    ch_molecules = Channel.empty()
    ch_targets   = Channel.empty()

    if (params.mode == 'full') {
        // ════════════════════════════════════════════════════
        //  FULL PIPELINE: Genomics → RAG Chat → Drug Discovery
        // ════════════════════════════════════════════════════

        // Create input channel from FASTQ files
        Channel
            .fromFilePairs(params.input, checkIfExists: true)
            .set { ch_reads }

        // Stage 1: Genomics Pipeline (FASTQ → VCF)
        GENOMICS_PIPELINE(ch_reads)
        ch_vcf = GENOMICS_PIPELINE.out.vcf
        ch_versions = ch_versions.mix(GENOMICS_PIPELINE.out.versions)

        // Stage 2: RAG Chat Pipeline (VCF → Target Hypothesis)
        RAG_CHAT_PIPELINE(ch_vcf)
        ch_targets = RAG_CHAT_PIPELINE.out.targets
        ch_versions = ch_versions.mix(RAG_CHAT_PIPELINE.out.versions)

        // Stage 3: Drug Discovery Pipeline (Target → Molecules)
        DRUG_DISCOVERY(ch_targets)
        ch_molecules = DRUG_DISCOVERY.out.molecules
        ch_versions = ch_versions.mix(DRUG_DISCOVERY.out.versions)

    } else if (params.mode == 'target') {
        // ════════════════════════════════════════════════════
        //  TARGET MODE: RAG Chat → Drug Discovery (skip genomics)
        // ════════════════════════════════════════════════════

        Channel
            .fromPath(params.vcf, checkIfExists: true)
            .set { ch_vcf }

        RAG_CHAT_PIPELINE(ch_vcf)
        ch_targets = RAG_CHAT_PIPELINE.out.targets

        DRUG_DISCOVERY(ch_targets)
        ch_molecules = DRUG_DISCOVERY.out.molecules

    } else if (params.mode == 'drug') {
        // ════════════════════════════════════════════════════
        //  DRUG MODE: Drug Discovery only (with provided target)
        // ════════════════════════════════════════════════════

        Channel
            .fromPath(params.target, checkIfExists: true)
            .set { ch_targets }

        DRUG_DISCOVERY(ch_targets)
        ch_molecules = DRUG_DISCOVERY.out.molecules

    } else if (params.mode == 'demo') {
        // ════════════════════════════════════════════════════
        //  DEMO MODE: VCP/FTD demonstration pipeline
        // ════════════════════════════════════════════════════

        // Use built-in VCP demo data
        Channel
            .fromPath("${projectDir}/data/demo/vcp_target.json")
            .set { ch_targets }

        DRUG_DISCOVERY(ch_targets)
        ch_molecules = DRUG_DISCOVERY.out.molecules

    } else if (params.mode == 'genomics_only') {
        // ════════════════════════════════════════════════════
        //  GENOMICS ONLY: FASTQ → VCF (no downstream analysis)
        // ════════════════════════════════════════════════════

        Channel
            .fromFilePairs(params.input, checkIfExists: true)
            .set { ch_reads }

        GENOMICS_PIPELINE(ch_reads)
        ch_versions = ch_versions.mix(GENOMICS_PIPELINE.out.versions)

        // Log completion with output location for genomics-only mode
        GENOMICS_PIPELINE.out.vcf
            .collect()
            .subscribe { vcf_files ->
                log.info """
  ────────────────────────────────────────────────────────
  Genomics-only mode complete.
  VCF outputs: ${params.outdir}/genomics/vcf
  Samples processed: ${vcf_files.size()}
  ────────────────────────────────────────────────────────
                """
            }
    }

    // Generate final report
    if (params.mode != 'genomics_only') {
        GENERATE_REPORT(
            ch_molecules,
            ch_targets
        )
    }
}

// ============================================================
//  Workflow completion handlers
// ============================================================

workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    def duration = workflow.duration
    def exitCode = workflow.exitStatus ?: 'N/A'

    log.info """
╔═══════════════════════════════════════════════════════════════╗
║  Pipeline Complete                                             ║
╚═══════════════════════════════════════════════════════════════╝

  Status:        ${status}
  Exit code:     ${exitCode}
  Duration:      ${duration}
  Pipeline:      ${pipelineName} v${pipelineVersion}
  Mode:          ${params.mode}
  Output:        ${params.outdir}

  Report:        ${params.outdir}/pipeline_report.html
  Timeline:      ${params.outdir}/timeline.html

  Nextflow:      ${workflow.nextflow.version}
  Run name:      ${workflow.runName}
  Session ID:    ${workflow.sessionId}

──────────────────────────────────────────────────────────────────
    """

    // Log failed process details when pipeline did not succeed
    if (!workflow.success) {
        log.warn """
  ┌─ Failure Details ──────────────────────────────────────────┐
  │  Error message: ${workflow.errorMessage ?: 'Unknown'}
  │  Error report:  ${workflow.errorReport  ?: 'None'}
  │
  │  Common causes:
  │    - Drug Discovery NIM service unavailable (check ports 8001/8002)
  │    - Insufficient GPU memory for MolMIM or DiffDock
  │    - Milvus not running for RAG Chat stage (check port 19530)
  │    - Missing reference files (BWA index, known sites VCF)
  └────────────────────────────────────────────────────────────┘
        """
    }

    // Send notification if configured
    if (params.email) {
        sendMail(
            to: params.email,
            subject: "[HLS-Pipeline] ${status}: ${workflow.runName}",
            body: "Pipeline completed with status: ${status}\nDuration: ${duration}\nVersion: ${pipelineVersion}"
        )
    }
}

workflow.onError {
    log.error """
  ╔═══════════════════════════════════════════════════════════════╗
  ║  PIPELINE ERROR                                                ║
  ╚═══════════════════════════════════════════════════════════════╝

  Error:         ${workflow.errorMessage}
  Pipeline:      ${pipelineName} v${pipelineVersion}
  Mode:          ${params.mode}

  Check the .nextflow.log file for full stack trace.
  """
}
