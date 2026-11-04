"""Regression tests for the Nextflow workflow's honesty defects.

These are source-level tests: Nextflow is not run here (no containers, no GPU). They pin the
exact shapes that produced fabricated or falsely-successful output.

compute#12 / C8 -- `RAG_IDENTIFY_TARGETS` was named RAG, commented "Identify drug targets using
RAG", wired in as "Stage 2 (RAG Chat)" and wrote `rag_chat: 1.0` into versions.yml, while
contacting no Milvus, no embedding model and no LLM: an 11-gene literal dict. `if gene in info`
was a substring test against the raw INFO column (VCP matched VCPIP1/VCPKMT, ALK matched
ALKBH5, GRN matched PGRN, KRAS matched KRASP1); the rsID clause sat inside the per-gene loop, so
one variant carrying rs188935092 emitted all eleven genes; every target was stamped confidence
'high', priority 5, status 'validated' regardless of evidence; and a VCF with no match at all
defaulted to a VCP/FTD target with the same stamps and a CB-5083 seed -- so a healthy control's
report read "Target Gene: VCP". `params.confidence_threshold` and `params.max_targets` were
never read.

compute#22 / C10 -- `params.mode` itself was never validated and the if/else-if chain had no
`else`, so `--mode genomics` (a typo for genomics_only) left both channels empty, executed
nothing, and printed "Status: SUCCESS" with paths to a report and timeline that were never
written.

compute#19 / C9 -- `GENERATE_REPORT` declared two unkeyed input channels, which Nextflow pairs
positionally in arrival order, so one patient's report could carry another patient's target and
PDB structures; `errorStrategy 'ignore'` let a failed report leave the workflow SUCCESS with no
report on disk; and the footer still read "Powered by NVIDIA BioNeMo".
"""
from pathlib import Path

import pytest

ORCH = Path(__file__).resolve().parent.parent


def _code(text: str) -> str:
    """Nextflow source with comments stripped.

    Each removed fabrication is quoted verbatim in the new comments so a reader can see what
    went, which means the comments must not count as occurrences. Groovy uses `/* */` and `//`;
    the Python process bodies use `#`.
    """
    out, in_block = [], False
    for line in text.splitlines():
        stripped = line.lstrip()
        if in_block:
            if "*/" in line:
                in_block = False
            continue
        if stripped.startswith("/*"):
            in_block = "*/" not in line
            continue
        if stripped.startswith(("//", "#", "*")):
            continue
        out.append(line)
    return "\n".join(out)


MAIN = (ORCH / "main.nf").read_text()
RAG = (ORCH / "modules" / "rag_chat.nf").read_text()
REPORTING = (ORCH / "modules" / "reporting.nf").read_text()
CONFIG = (ORCH / "nextflow.config").read_text()
MAIN_CODE = _code(MAIN)
RAG_CODE = _code(RAG)
REPORTING_CODE = _code(REPORTING)


class TestUnknownModeCannotSucceed:
    def test_the_mode_is_validated_against_an_allowed_set(self):
        assert "ALLOWED_MODES" in MAIN_CODE
        assert "ALLOWED_MODES.contains(params.mode)" in MAIN_CODE
        assert "Unknown --mode" in MAIN_CODE

    @pytest.mark.parametrize("mode", ["full", "target", "drug", "demo", "genomics_only"])
    def test_every_implemented_mode_is_allowed(self, mode):
        allowed = MAIN_CODE.split("def ALLOWED_MODES = [", 1)[1].split("]", 1)[0]
        assert f"'{mode}'" in allowed

    def test_the_demo_channel_checks_its_input_exists(self):
        line = next(ln for ln in MAIN_CODE.splitlines()
                    if "data/demo/vcp_target.json" in ln)
        assert "checkIfExists: true" in line


class TestTheReportIsKeyedToItsSample:
    def test_generate_report_takes_one_joined_channel(self):
        # Two unkeyed `tuple val(...)` inputs is the defect; one 3-element tuple is the fix.
        assert "tuple val(sample_id), path(ranked), path(targets)" in REPORTING_CODE
        assert "tuple val(sample_id2), path(targets)" not in REPORTING_CODE

    def test_main_joins_on_sample_id(self):
        assert "ch_molecules.join(ch_targets_keyed" in MAIN_CODE

    def test_a_failed_report_can_no_longer_be_ignored(self):
        assert "errorStrategy 'ignore'" not in REPORTING_CODE

    def test_the_bionemo_footer_is_gone(self):
        assert "Powered by NVIDIA BioNeMo" not in REPORTING_CODE


class TestTargetLookupDoesNotClaimRag:
    def test_the_process_is_not_called_rag(self):
        assert "process RAG_IDENTIFY_TARGETS" not in RAG_CODE
        assert "process LOOKUP_TARGETS_IN_GENE_TABLE" in RAG_CODE

    def test_versions_yml_does_not_claim_a_rag_component(self):
        assert "rag_chat: 1.0" not in RAG_CODE
        assert "gene_table_lookup: 1.0" in RAG_CODE

    def test_the_output_records_what_the_method_actually_was(self):
        assert "not RAG" in RAG_CODE


class TestGeneMatchingIsNotASubstringTest:
    def test_the_substring_test_is_gone(self):
        # Both executable spellings of the substring test. (`_code` strips Groovy and `#`
        # comments but not the Python docstrings inside the process body, which quote the old
        # expression deliberately.)
        assert "if gene in info or" not in RAG_CODE
        assert "if gene in info:" not in RAG_CODE

    def test_the_gene_symbol_is_parsed_from_the_structured_annotation(self):
        assert "ANN=" in RAG_CODE and "CSQ=" in RAG_CODE
        assert "def genes_in(variant)" in RAG_CODE

    def test_the_rsid_clause_is_no_longer_inside_the_gene_loop(self):
        # One variant carrying this rsID used to emit all eleven genes as targets.
        assert "'rs188935092' in str(variant.get('rsid'))" not in RAG_CODE
        # The only surviving mention is inside the opt-in demo default, after the lookup loop.
        assert "rs188935092" not in RAG_CODE.split("targets = sorted(")[0]

    def test_genes_in_matches_whole_symbols_only(self):
        """Executed against the same parser the process body uses."""
        ns = {}
        body = RAG.split("def genes_in(variant):", 1)[1].split("\n    # Identify target", 1)[0]
        src = "def genes_in(variant):" + body
        exec("\n".join(ln[4:] if ln.startswith("    ") else ln
                       for ln in src.splitlines()), ns)
        genes_in = ns["genes_in"]
        # snpEff ANN with a VCP-like neighbour: VCPIP1 must not read as VCP.
        assert genes_in({"info": "ANN=A|missense_variant|HIGH|VCPIP1|ENSG1|transcript"}) == {
            "VCPIP1"}
        assert genes_in({"info": "ANN=A|missense_variant|HIGH|VCP|ENSG1|transcript"}) == {"VCP"}
        assert genes_in({"info": "ANN=A|intron_variant|LOW|ALKBH5|x|y"}) == {"ALKBH5"}
        assert genes_in({"info": "DP=30;AF=0.5"}) == set()


class TestNoTargetIsStampedValidatedWithoutEvidence:
    def test_the_unconditional_stamps_are_gone(self):
        assert "'confidence': 'high',\n                    'priority': 5," not in RAG_CODE
        assert "'status': 'validated'," not in RAG_CODE

    def test_confidence_is_derived_and_thresholded(self):
        assert "CONFIDENCE_SCORE" in RAG_CODE
        assert "params.confidence_threshold" in RAG_CODE
        assert "params.max_targets" in RAG_CODE

    def test_the_status_says_what_was_actually_done(self):
        assert "'gene_table_match'" in RAG_CODE


class TestTheDemoTargetFallbackIsOptIn:
    def test_the_fallback_requires_an_explicit_flag(self):
        assert "params.demo_target_fallback" in RAG_CODE
        assert "demo_target_fallback = false" in CONFIG

    def test_the_fallback_labels_itself_as_a_default(self):
        block = RAG_CODE.split("params.demo_target_fallback", 1)[1]
        assert "'demo_default'" in block
        assert "DEMO DEFAULT" in block
        assert "'validated'" not in block.split("output = {")[0]
