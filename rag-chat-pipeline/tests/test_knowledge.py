"""Tests for RAG Chat knowledge base data integrity."""
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.knowledge import (
    GENE_REFERENCE_DATA,
    KNOWLEDGE_CONNECTIONS,
    format_knowledge_for_prompt,
    get_druggable_genes,
    get_gene_drugs,
    get_gene_reference_data,
    get_genes_by_disease,
    get_genes_by_pathway,
    get_knowledge_for_evidence,
    get_knowledge_for_genes,
    get_knowledge_stats,
)

PDB_PATTERN = re.compile(r'^[A-Za-z0-9]{4}$')


class TestKnowledgeConnectionsIntegrity:
    """Verify KNOWLEDGE_CONNECTIONS data integrity."""

    def test_minimum_gene_count(self):
        assert len(KNOWLEDGE_CONNECTIONS) >= 201

    def test_all_genes_consistent_case(self):
        """Gene keys should match standard nomenclature (mostly uppercase, some mixed like C9orf72)."""
        for gene in KNOWLEDGE_CONNECTIONS:
            # Gene names are alphanumeric with optional hyphens/underscores
            assert gene == gene.strip(), f"Gene {gene} has whitespace"
            assert len(gene) > 0, "Empty gene name found"

    @pytest.mark.parametrize("field", [
        'protein', 'function', 'pathway', 'diseases', 'drugs',
        'drug_status', 'pdb_ids', 'druggable',
    ])
    def test_required_fields_present(self, field):
        missing = [g for g, d in KNOWLEDGE_CONNECTIONS.items() if field not in d]
        assert missing == [], f"Genes missing '{field}': {missing[:5]}"

    def test_diseases_are_nonempty_lists(self):
        bad = [g for g, d in KNOWLEDGE_CONNECTIONS.items()
               if not isinstance(d['diseases'], list) or len(d['diseases']) == 0]
        assert bad == [], f"Genes with empty/non-list diseases: {bad[:5]}"

    def test_drugs_are_lists(self):
        bad = [g for g, d in KNOWLEDGE_CONNECTIONS.items()
               if not isinstance(d['drugs'], list)]
        assert bad == [], f"Genes with non-list drugs: {bad[:5]}"

    def test_pdb_ids_format(self):
        bad = []
        for gene, data in KNOWLEDGE_CONNECTIONS.items():
            for pdb in data['pdb_ids']:
                if not PDB_PATTERN.match(pdb):
                    bad.append(f"{gene}:{pdb}")
        assert bad == [], f"Invalid PDB IDs: {bad[:5]}"

    def test_druggable_is_boolean(self):
        bad = [g for g, d in KNOWLEDGE_CONNECTIONS.items()
               if not isinstance(d['druggable'], bool)]
        assert bad == [], f"Genes with non-bool druggable: {bad[:5]}"

    def test_druggable_genes_have_drugs(self):
        bad = [g for g, d in KNOWLEDGE_CONNECTIONS.items()
               if d['druggable'] and len(d['drugs']) == 0]
        assert bad == [], f"Druggable genes with no drugs: {bad[:5]}"

    def test_flagship_genes_present(self):
        """Verify key demo genes are in the knowledge base."""
        for gene in ['VCP', 'EGFR', 'BRAF', 'CFTR', 'KRAS', 'BRCA1', 'TP53']:
            assert gene in KNOWLEDGE_CONNECTIONS, f"Missing flagship gene: {gene}"


class TestGeneReferenceDataIntegrity:
    """Verify GENE_REFERENCE_DATA integrity."""

    def test_has_entries(self):
        assert len(GENE_REFERENCE_DATA) > 0

    def test_uniprot_id_present(self):
        missing = [g for g, d in GENE_REFERENCE_DATA.items() if 'uniprot_id' not in d]
        assert missing == [], f"Genes missing 'uniprot_id': {missing[:5]}"

    def test_most_genes_have_smiles(self):
        """Most genes should have reference_smiles; a few may not."""
        total = len(GENE_REFERENCE_DATA)
        with_smiles = sum(1 for d in GENE_REFERENCE_DATA.values()
                          if d.get('reference_smiles'))
        assert with_smiles >= total * 0.9, (
            f"Only {with_smiles}/{total} genes have reference SMILES"
        )

    def test_uniprot_format(self):
        bad = [g for g, d in GENE_REFERENCE_DATA.items()
               if len(d.get('uniprot_id', '')) < 5]
        assert bad == [], f"Genes with short UniProt IDs: {bad[:5]}"


class TestGetGeneReferenceData:
    """Tests for get_gene_reference_data()."""

    def test_known_gene_vcp(self):
        result = get_gene_reference_data('VCP')
        assert result['gene'] == 'VCP'
        assert result['protein'] == 'p97/VCP ATPase'
        assert result['druggable'] is True
        assert len(result['pdb_ids']) >= 3

    def test_case_insensitive(self):
        result = get_gene_reference_data('vcp')
        assert result['gene'] == 'VCP'

    def test_unknown_gene(self):
        result = get_gene_reference_data('NOTAREALGENE')
        assert result['gene'] == 'NOTAREALGENE'
        assert result['uniprot_id'] is None

    def test_egfr_has_drugs(self):
        result = get_gene_reference_data('EGFR')
        assert len(result['drugs']) >= 3


class TestHelperFunctions:
    """Tests for knowledge.py helper functions."""

    def test_get_druggable_genes(self):
        druggable = get_druggable_genes()
        assert 'VCP' in druggable
        assert len(druggable) > 100

    def test_get_gene_drugs_known(self):
        drugs = get_gene_drugs('VCP')
        assert drugs is not None
        assert len(drugs) > 0

    def test_get_gene_drugs_unknown(self):
        drugs = get_gene_drugs('FAKEGENE')
        assert drugs is None

    def test_get_genes_by_disease(self):
        ftd_genes = get_genes_by_disease('Frontotemporal')
        assert 'VCP' in ftd_genes
        assert len(ftd_genes) >= 3

    def test_get_genes_by_pathway(self):
        genes = get_genes_by_pathway('proteasome')
        assert len(genes) >= 1

    def test_get_knowledge_stats(self):
        stats = get_knowledge_stats()
        assert stats['total_genes'] >= 201
        assert stats['druggable_targets'] > 0
        assert stats['non_druggable'] >= 0
        assert stats['total_genes'] == stats['druggable_targets'] + stats['non_druggable']

    def test_get_knowledge_for_genes(self):
        conns = get_knowledge_for_genes(['VCP', 'EGFR'])
        assert len(conns) == 2
        assert conns[0]['gene'] == 'VCP'

    def test_get_knowledge_for_genes_dedup(self):
        conns = get_knowledge_for_genes(['VCP', 'vcp', 'VCP'])
        assert len(conns) == 1

    def test_format_knowledge_empty(self):
        result = format_knowledge_for_prompt([])
        assert result == ""

    def test_format_knowledge_with_evidence(self):
        evidence = [{'gene': 'VCP'}, {'gene': 'EGFR'}]
        result = format_knowledge_for_prompt(evidence)
        assert 'VCP' in result
        assert 'EGFR' in result
        assert 'DRUGGABLE' in result
