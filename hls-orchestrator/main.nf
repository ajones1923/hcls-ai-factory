#!/usr/bin/env nextflow

/*
 * ========================================================
 *  Healthcare & Life Sciences Pipeline Orchestrator
 * ========================================================
 *  Unified workflow: Genomics → RAG Chat → Drug Discovery
 *
 *  Author: Adam Jones
 *  Version: 1.0.0
 * ========================================================
 */

nextflow.enable.dsl = 2

// Pipeline metadata
def pipelineVersion = '1.0.0'
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

    // Channel definitions
    ch_versions = Channel.empty()

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

    log.info """
╔═══════════════════════════════════════════════════════════════╗
║  Pipeline Complete                                             ║
╚═══════════════════════════════════════════════════════════════╝

  Status:        ${status}
  Duration:      ${duration}
  Output:        ${params.outdir}

  Report:        ${params.outdir}/pipeline_report.html
  Timeline:      ${params.outdir}/timeline.html

──────────────────────────────────────────────────────────────────
    """

    // Send notification if configured
    if (params.email) {
        sendMail(
            to: params.email,
            subject: "[HLS-Pipeline] ${status}: ${workflow.runName}",
            body: "Pipeline completed with status: ${status}\nDuration: ${duration}"
        )
    }
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
