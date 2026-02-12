"""
RAG Engine - Retrieval-Augmented Generation for genomic evidence.
"""
import os
from typing import List, Dict, Any, Optional, Generator
from loguru import logger

from .embedder import EvidenceEmbedder
from .milvus_client import MilvusClient
from .llm_client import LLMClient, BaseLLMClient
from .knowledge import format_knowledge_for_prompt, get_knowledge_for_evidence


# System prompt for the RAG assistant
GENOMICS_RAG_SYSTEM_PROMPT = """You are a genomics expert assistant with access to variant evidence from whole genome sequencing data.

Your role is to:
1. Answer questions about genetic variants found in the sample (HG002)
2. Explain the potential significance of variants
3. Provide context about genes and their functions
4. Help interpret variant consequences and clinical significance
5. Provide pharmacogenomic insights for drug metabolism questions
6. Identify potential drug targets based on disease-associated variants

When responding:
- Base your answers on the provided evidence
- Cite specific variants when relevant (chromosome, position, gene, rsID)
- Explain technical terms in accessible language
- If the evidence doesn't contain relevant information, say so clearly
- Be precise about what is known vs. inferred
- For drug target questions, consider protein function and druggability

NEURODEGENERATION & FTD GUIDANCE:
When asked about dementia, neurodegeneration, or related conditions:
- VCP (Valosin-containing protein/p97): AAA+ ATPase critical for protein quality control
  * Mutations cause IBMPFD (Inclusion Body Myopathy with Paget's disease and Frontotemporal Dementia)
  * Also linked to ALS (amyotrophic lateral sclerosis)
  * VCP is a validated drug target with inhibitors in development (e.g., CB-5083)
  * Report ALL VCP variants as potentially relevant to neurodegeneration
- C9orf72: Most common genetic cause of FTD and ALS (repeat expansions)
- GRN (Progranulin): Haploinsufficiency causes FTD
- MAPT (Tau): Mutations cause frontotemporal dementia with parkinsonism
- TBK1, FUS, TARDBP: Additional FTD/ALS genes
- For VCP specifically, explain its role in:
  * Endoplasmic reticulum-associated degradation (ERAD)
  * Autophagy and mitophagy
  * DNA damage response
  * Why dysfunction leads to protein aggregation and neurodegeneration

PHARMACOGENOMICS GUIDANCE:
When asked about drug response or metabolism:
- CYP450 enzymes (CYP2C19, CYP2D6, CYP3A4, CYP2C9, etc.) are key drug metabolizers
- Even "benign" variants in CYP genes can affect drug metabolism rates
- Report ANY variants found in relevant CYP genes, regardless of clinical classification
- Common drug-gene associations:
  * CYP2C19: PPIs (omeprazole, pantoprazole), clopidogrel, antidepressants
  * CYP2D6: codeine, tramadol, tamoxifen, many antidepressants
  * CYP3A4: statins, calcium channel blockers, immunosuppressants
  * CYP2C9: warfarin, NSAIDs, sulfonylureas
  * VKORC1: warfarin dosing
  * SLCO1B1: statin myopathy risk
- If CYP variants are found, explain their potential impact on drug metabolism
- Distinguish between: poor metabolizers (reduced function), normal metabolizers, and ultrarapid metabolizers

Evidence format: Each piece of evidence represents a variant with its location,
gene (if known), consequence type, clinical significance, and quality score.
"""

# Neurodegeneration query expansion - maps disease keywords to relevant genes
NEURODEGENERATION_GENES = {
    # Frontotemporal dementia
    'frontotemporal': ['VCP', 'C9orf72', 'GRN', 'MAPT', 'TBK1', 'FUS', 'TARDBP'],
    'ftd': ['VCP', 'C9orf72', 'GRN', 'MAPT', 'TBK1', 'FUS', 'TARDBP'],
    'dementia': ['VCP', 'C9orf72', 'GRN', 'MAPT', 'PSEN1', 'PSEN2', 'APP'],
    # ALS
    'als': ['VCP', 'C9orf72', 'SOD1', 'FUS', 'TARDBP', 'TBK1'],
    'amyotrophic lateral sclerosis': ['VCP', 'C9orf72', 'SOD1', 'FUS', 'TARDBP', 'TBK1'],
    'motor neuron': ['VCP', 'C9orf72', 'SOD1', 'FUS', 'TARDBP'],
    # VCP-specific
    'vcp': ['VCP'],
    'p97': ['VCP'],
    'ibmpfd': ['VCP'],
    'inclusion body myopathy': ['VCP'],
    'paget': ['VCP'],
    # General neurodegeneration
    'neurodegeneration': ['VCP', 'C9orf72', 'GRN', 'MAPT', 'SOD1', 'FUS', 'TARDBP', 'TBK1'],
    'protein aggregation': ['VCP', 'MAPT', 'TARDBP', 'FUS'],
    'proteinopathy': ['VCP', 'MAPT', 'TARDBP'],
}

# Pharmacogenomic query expansion - maps drug keywords to relevant genes
PHARMACOGENOMIC_GENES = {
    # PPIs and acid reflux
    'ppi': ['CYP2C19', 'CYP3A4'],
    'proton pump': ['CYP2C19', 'CYP3A4'],
    'omeprazole': ['CYP2C19', 'CYP3A4'],
    'pantoprazole': ['CYP2C19', 'CYP3A4'],
    'esomeprazole': ['CYP2C19', 'CYP3A4'],
    'lansoprazole': ['CYP2C19', 'CYP3A4'],
    'acid reflux': ['CYP2C19', 'CYP3A4'],
    'gerd': ['CYP2C19', 'CYP3A4'],
    # Anticoagulants
    'warfarin': ['CYP2C9', 'VKORC1', 'CYP4F2'],
    'clopidogrel': ['CYP2C19'],
    'plavix': ['CYP2C19'],
    # Pain medications
    'codeine': ['CYP2D6'],
    'tramadol': ['CYP2D6', 'UGT2B7'],
    'oxycodone': ['CYP2D6', 'CYP3A4'],
    # Antidepressants
    'ssri': ['CYP2D6', 'CYP2C19'],
    'antidepressant': ['CYP2D6', 'CYP2C19'],
    'citalopram': ['CYP2C19', 'CYP2D6'],
    'sertraline': ['CYP2C19', 'CYP2D6'],
    # Statins
    'statin': ['SLCO1B1', 'CYP3A4'],
    'atorvastatin': ['CYP3A4', 'SLCO1B1'],
    'simvastatin': ['CYP3A4', 'SLCO1B1'],
    # Tamoxifen
    'tamoxifen': ['CYP2D6'],
    # General drug metabolism
    'drug metabolism': ['CYP2D6', 'CYP2C19', 'CYP3A4', 'CYP2C9'],
    'pharmacogenomic': ['CYP2D6', 'CYP2C19', 'CYP3A4', 'CYP2C9', 'VKORC1'],
}

# Metabolic/Endocrine query expansion
METABOLIC_GENES = {
    # Diabetes
    'diabetes': ['GLP1R', 'GIPR', 'GCGR', 'INSR', 'PPARG', 'DPP4', 'SGLT2', 'KCNJ11', 'ABCC8', 'GCK', 'HNF1A', 'HNF4A'],
    'type 2 diabetes': ['GLP1R', 'GIPR', 'GCGR', 'PPARG', 'DPP4', 'SGLT2'],
    'mody': ['GCK', 'HNF1A', 'HNF4A', 'KCNJ11', 'ABCC8'],
    'insulin': ['INSR', 'PPARG', 'KCNJ11', 'ABCC8'],
    'glp-1': ['GLP1R', 'GIPR', 'DPP4'],
    'semaglutide': ['GLP1R'],
    'tirzepatide': ['GLP1R', 'GIPR'],
    'ozempic': ['GLP1R'],
    'wegovy': ['GLP1R'],
    # Obesity
    'obesity': ['GLP1R', 'GIPR', 'GCGR', 'MC4R', 'LEPR', 'POMC', 'PCSK1'],
    'weight loss': ['GLP1R', 'GIPR', 'MC4R'],
    # Cholesterol/lipids
    'cholesterol': ['HMGCR', 'PCSK9', 'LDLR', 'APOB', 'NPC1L1', 'ANGPTL3', 'CETP'],
    'hypercholesterolemia': ['LDLR', 'PCSK9', 'APOB', 'ANGPTL3'],
    'statin': ['HMGCR', 'SLCO1B1'],
    'pcsk9': ['PCSK9'],
    # Thyroid
    'thyroid': ['TSHR'],
    'graves': ['TSHR'],
}

# Infectious Disease query expansion
INFECTIOUS_GENES = {
    # HIV
    'hiv': ['HIV1_RT', 'HIV1_PR', 'HIV1_IN', 'CCR5', 'CD4'],
    'aids': ['HIV1_RT', 'HIV1_PR', 'HIV1_IN', 'CCR5', 'CD4'],
    'antiretroviral': ['HIV1_RT', 'HIV1_PR', 'HIV1_IN'],
    # Hepatitis
    'hepatitis c': ['HCV_NS3', 'HCV_NS5A', 'HCV_NS5B'],
    'hcv': ['HCV_NS3', 'HCV_NS5A', 'HCV_NS5B'],
    'hepatitis b': ['HBV_RT'],
    'hbv': ['HBV_RT'],
    # COVID/SARS
    'covid': ['SARS2_MPRO', 'SARS2_RDRP', 'SARS2_SPIKE', 'ACE2', 'TMPRSS2'],
    'sars': ['SARS2_MPRO', 'SARS2_RDRP', 'SARS2_SPIKE', 'ACE2'],
    'coronavirus': ['SARS2_MPRO', 'SARS2_RDRP', 'SARS2_SPIKE'],
    'paxlovid': ['SARS2_MPRO'],
    # Bacterial
    'tuberculosis': ['RPOB', 'INHA'],
    'bacterial': ['DHFR', 'DHPS', 'GYRA', 'RPOB'],
    'antibiotic': ['DHFR', 'DHPS', 'GYRA', 'RPOB'],
    # Fungal
    'fungal': ['CYP51A1', 'FKS1'],
    'candida': ['CYP51A1', 'FKS1'],
    'antifungal': ['CYP51A1', 'FKS1'],
}

# Respiratory/Pulmonary query expansion
RESPIRATORY_GENES = {
    'asthma': ['ADRB2', 'CHRM3', 'IL4R', 'IL5', 'IL5RA', 'IL13', 'FCER1A', 'TSLP'],
    'copd': ['ADRB2', 'CHRM3', 'PDE4D', 'SERPINA1'],
    'alpha-1': ['SERPINA1'],
    'pulmonary hypertension': ['BMPR2', 'EDNRA', 'PDE5A', 'GUCY1A1'],
    'pah': ['BMPR2', 'EDNRA', 'PDE5A', 'GUCY1A1'],
    'bronchodilator': ['ADRB2', 'CHRM3'],
    'eosinophil': ['IL5', 'IL5RA', 'IL4R'],
    'allergy': ['FCER1A', 'TSLP', 'IL4R', 'IL13'],
}

# Ophthalmology query expansion
OPHTHALMOLOGY_GENES = {
    'macular degeneration': ['VEGFA', 'CFH', 'C3', 'C5'],
    'amd': ['VEGFA', 'CFH', 'C3', 'C5'],
    'diabetic retinopathy': ['VEGFA'],
    'glaucoma': ['ROCK1', 'CA2', 'OPTN'],
    'retinitis pigmentosa': ['RHO', 'RPE65', 'USH2A', 'ABCA4'],
    'leber': ['RPE65'],
    'stargardt': ['ABCA4'],
    'usher': ['USH2A'],
    'geographic atrophy': ['CFH', 'C3', 'C5'],
    'complement': ['CFH', 'C3', 'C5'],
}

# Pain/Migraine query expansion
PAIN_GENES = {
    'migraine': ['CGRP', 'CALCA'],
    'headache': ['CGRP', 'CALCA'],
    'pain': ['SCN9A', 'SCN10A', 'TRPV1'],
    'neuropathic pain': ['SCN9A', 'SCN10A', 'TRPV1'],
    'chronic pain': ['SCN9A', 'SCN10A', 'TRPV1'],
    'nav1.7': ['SCN9A'],
    'nav1.8': ['SCN10A'],
    'cgrp': ['CGRP', 'CALCA'],
}

# GI/Hepatology query expansion
GI_GENES = {
    'crohn': ['NOD2', 'ATG16L1', 'IL12B', 'ITGA4', 'ITGB7'],
    'ulcerative colitis': ['IL12B', 'ITGA4', 'ITGB7', 'S1PR1'],
    'ibd': ['NOD2', 'ATG16L1', 'IL12B', 'ITGA4', 'ITGB7', 'S1PR1'],
    'inflammatory bowel': ['NOD2', 'ATG16L1', 'IL12B', 'ITGA4', 'ITGB7'],
    'nash': ['PNPLA3', 'HSD17B13', 'FXR', 'THR_BETA', 'GLP1R'],
    'nafld': ['PNPLA3', 'HSD17B13', 'FXR', 'THR_BETA'],
    'fatty liver': ['PNPLA3', 'HSD17B13', 'FXR', 'THR_BETA'],
    'liver disease': ['PNPLA3', 'HSD17B13', 'FXR', 'THR_BETA', 'SERPINA1'],
    'gerd': ['ATP4A', 'CYP2C19'],
    'acid reflux': ['ATP4A', 'CYP2C19'],
    'ppi': ['ATP4A', 'CYP2C19'],
}

# Oncology query expansion
ONCOLOGY_GENES = {
    'breast cancer': ['BRCA1', 'BRCA2', 'ERBB2', 'PIK3CA', 'CDK4', 'CDK6'],
    'lung cancer': ['EGFR', 'ALK', 'KRAS', 'ROS1', 'MET', 'RET', 'BRAF'],
    'melanoma': ['BRAF', 'PDCD1', 'CD274', 'CTLA4', 'KIT'],
    'leukemia': ['BCR-ABL1', 'FLT3', 'IDH1', 'IDH2', 'BTK', 'BCL2'],
    'lymphoma': ['BTK', 'BCL2', 'CD274'],
    'colon cancer': ['KRAS', 'BRAF', 'EGFR'],
    'colorectal': ['KRAS', 'BRAF', 'EGFR'],
    'immunotherapy': ['PDCD1', 'CD274', 'CTLA4'],
    'checkpoint inhibitor': ['PDCD1', 'CD274', 'CTLA4'],
    'parp inhibitor': ['BRCA1', 'BRCA2'],
    'her2': ['ERBB2'],
    'egfr': ['EGFR'],
    'alk': ['ALK'],
    'kras': ['KRAS'],
    'braf': ['BRAF'],
}

# Hematology query expansion
HEMATOLOGY_GENES = {
    'thrombocytopenia': ['SYK', 'THPO', 'MPL'],
    'itp': ['SYK', 'THPO', 'MPL'],
    'anemia': ['EPOR', 'HIF2A', 'HBB'],
    'polycythemia': ['JAK2', 'EPOR'],
    'myelofibrosis': ['JAK2', 'CALR', 'MPL'],
    'thrombosis': ['F10', 'F2', 'SERPINC1'],
    'anticoagulant': ['F10', 'F2', 'VKORC1'],
    'ttp': ['ADAMTS13', 'VWF'],
    'von willebrand': ['VWF'],
    'hemophilia': ['F8', 'F9', 'SERPINC1'],
}


class RAGEngine:
    """
    RAG Engine for genomic evidence retrieval and question answering.
    """

    def __init__(
        self,
        milvus_client: MilvusClient,
        embedder: EvidenceEmbedder,
        llm_client: BaseLLMClient,
        top_k: int = 10,
        score_threshold: float = 0.5,
        system_prompt: Optional[str] = None,
    ):
        self.milvus = milvus_client
        self.embedder = embedder
        self.llm = llm_client
        self.top_k = top_k
        self.score_threshold = score_threshold
        self.system_prompt = system_prompt or GENOMICS_RAG_SYSTEM_PROMPT

        logger.info("RAG Engine initialized")

    def _get_expanded_genes(self, query: str) -> List[str]:
        """
        Extract relevant genes based on disease/drug keywords in the query.

        Checks all therapeutic area gene mappings for query expansion.
        Returns list of gene names to also search for.
        """
        query_lower = query.lower()
        relevant_genes = set()

        # All gene expansion dictionaries with their category names
        gene_mappings = [
            ('Neurodegeneration', NEURODEGENERATION_GENES),
            ('Pharmacogenomic', PHARMACOGENOMIC_GENES),
            ('Metabolic', METABOLIC_GENES),
            ('Infectious', INFECTIOUS_GENES),
            ('Respiratory', RESPIRATORY_GENES),
            ('Ophthalmology', OPHTHALMOLOGY_GENES),
            ('Pain', PAIN_GENES),
            ('GI', GI_GENES),
            ('Oncology', ONCOLOGY_GENES),
            ('Hematology', HEMATOLOGY_GENES),
        ]

        for category, mapping in gene_mappings:
            for keyword, genes in mapping.items():
                if keyword in query_lower:
                    relevant_genes.update(genes)
                    logger.info(f"{category} expansion: '{keyword}' -> {genes}")

        return list(relevant_genes)

    def _get_pharmacogenomic_genes(self, query: str) -> List[str]:
        """Deprecated: Use _get_expanded_genes instead. Kept for compatibility."""
        return self._get_expanded_genes(query)

    def _retrieve_by_genes(self, genes: List[str], limit_per_gene: int = 3) -> List[Dict[str, Any]]:
        """
        Retrieve variants for specific genes by direct query.
        """
        results = []
        collection = self.milvus.get_collection()

        for gene in genes:
            try:
                gene_results = collection.query(
                    expr=f'gene == "{gene}"',
                    output_fields=['id', 'gene', 'chrom', 'pos', 'ref', 'alt', 'rsid',
                                   'consequence', 'impact', 'clinical_significance',
                                   'text_summary', 'disease_associations'],
                    limit=limit_per_gene
                )
                for r in gene_results:
                    r['score'] = 0.85  # Assign high score for direct gene match
                    r['retrieval_method'] = 'pharmacogenomic_expansion'
                    results.append(r)
                logger.info(f"Found {len(gene_results)} variants for pharmacogenomic gene {gene}")
            except Exception as e:
                logger.warning(f"Error querying gene {gene}: {e}")

        return results

    def retrieve(
        self,
        query: str,
        top_k: Optional[int] = None,
        filter_expr: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Retrieve relevant evidence for a query.

        Args:
            query: User question or search query
            top_k: Number of results (default: self.top_k)
            filter_expr: Optional Milvus filter (e.g., "gene == 'EGFR'")

        Returns:
            List of relevant evidence items with scores
        """
        k = top_k or self.top_k

        # Embed the query
        query_embedding = self.embedder.embed_query(query)

        # Search Milvus with semantic similarity
        results = self.milvus.search(
            query_embedding=query_embedding,
            top_k=k,
            score_threshold=self.score_threshold,
            filter_expr=filter_expr,
        )

        logger.info(f"Semantic search retrieved {len(results)} evidence items")

        # Pharmacogenomic query expansion
        pgx_genes = self._get_expanded_genes(query)
        if pgx_genes:
            # Get variants from relevant pharmacogenomic genes
            pgx_results = self._retrieve_by_genes(pgx_genes, limit_per_gene=3)

            # Merge results, avoiding duplicates
            existing_ids = {r.get('id') for r in results}
            for pgx_result in pgx_results:
                if pgx_result.get('id') not in existing_ids:
                    results.append(pgx_result)
                    existing_ids.add(pgx_result.get('id'))

            logger.info(f"After pharmacogenomic expansion: {len(results)} total evidence items")

        # Sort by score and limit
        results = sorted(results, key=lambda x: x.get('score', 0), reverse=True)[:k + 5]

        logger.info(f"Retrieved {len(results)} evidence items for query")
        return results

    def _format_evidence_context(self, evidence_list: List[Dict[str, Any]]) -> str:
        """Format evidence for inclusion in LLM prompt."""
        if not evidence_list:
            return "No relevant evidence found in the database."

        context_parts = ["## Retrieved Evidence\n"]

        for i, ev in enumerate(evidence_list, 1):
            context_parts.append(f"### Evidence {i} (similarity: {ev['score']:.2f})")
            context_parts.append(ev['text_summary'])

            # Add structured details
            details = []
            if ev.get('gene'):
                details.append(f"Gene: {ev['gene']}")
            if ev.get('consequence'):
                details.append(f"Consequence: {ev['consequence']}")
            if ev.get('impact'):
                details.append(f"Impact: {ev['impact']}")
            if ev.get('clinical_significance'):
                details.append(f"Clinical: {ev['clinical_significance']}")

            if details:
                context_parts.append("Details: " + " | ".join(details))

            context_parts.append("")  # Empty line between evidence

        return "\n".join(context_parts)

    def _build_prompt(self, query: str, context: str, evidence: List[Dict[str, Any]] = None) -> str:
        """Build the full prompt for the LLM, including Clinker knowledge connections."""
        # Get Clinker knowledge context for genes in the evidence
        knowledge_context = ""
        if evidence:
            knowledge_context = format_knowledge_for_prompt(evidence)
            if knowledge_context:
                logger.info(f"Added Clinker knowledge for {len(get_knowledge_for_evidence(evidence))} genes")

        prompt_parts = [
            "Based on the following genomic evidence and therapeutic knowledge, please answer the user's question.",
            "",
            context,
        ]

        # Add Clinker knowledge if available
        if knowledge_context:
            prompt_parts.extend(["", knowledge_context])

        prompt_parts.extend([
            "",
            "## User Question",
            query,
            "",
            "## Instructions",
            "- Answer based on the evidence and therapeutic knowledge provided above",
            "- If genes have known drug targets or therapeutics in development, mention them",
            "- Highlight druggable targets and their development status when relevant",
            "- If the evidence doesn't contain relevant information, clearly state that",
            "- Explain any technical terms",
            "- Be specific about variant locations and genes when relevant",
        ])

        return "\n".join(prompt_parts)

    def query(
        self,
        question: str,
        filter_expr: Optional[str] = None,
        include_evidence: bool = True,
    ) -> Dict[str, Any]:
        """
        Answer a question using RAG.

        Args:
            question: User question
            filter_expr: Optional filter for evidence retrieval
            include_evidence: Whether to include evidence in response

        Returns:
            Dict with 'answer', 'evidence', and 'query' keys
        """
        # Retrieve relevant evidence
        evidence = self.retrieve(question, filter_expr=filter_expr)

        # Format context
        context = self._format_evidence_context(evidence)

        # Build prompt with Clinker knowledge
        prompt = self._build_prompt(question, context, evidence)

        # Generate answer
        answer = self.llm.generate(
            prompt=prompt,
            system_prompt=self.system_prompt,
        )

        result = {
            "query": question,
            "answer": answer,
        }

        if include_evidence:
            result["evidence"] = evidence
            result["evidence_count"] = len(evidence)

        return result

    def query_stream(
        self,
        question: str,
        filter_expr: Optional[str] = None,
    ) -> Generator[Dict[str, Any], None, None]:
        """
        Stream answer to a question using RAG.

        Yields:
            Dict with 'type' and 'content' keys:
            - {"type": "evidence", "content": evidence_list}
            - {"type": "token", "content": token_string}
            - {"type": "done", "content": full_answer}
        """
        # Retrieve relevant evidence
        evidence = self.retrieve(question, filter_expr=filter_expr)

        # Yield evidence first
        yield {"type": "evidence", "content": evidence}

        # Format context
        context = self._format_evidence_context(evidence)

        # Build prompt with Clinker knowledge
        prompt = self._build_prompt(question, context, evidence)

        # Stream answer
        full_answer = ""
        for token in self.llm.generate_stream(
            prompt=prompt,
            system_prompt=self.system_prompt,
        ):
            full_answer += token
            yield {"type": "token", "content": token}

        yield {"type": "done", "content": full_answer}

    def search_gene(self, gene: str) -> List[Dict[str, Any]]:
        """Search for all variants in a specific gene."""
        return self.milvus.search_by_gene(gene)

    def search_region(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> List[Dict[str, Any]]:
        """Search for variants in a genomic region."""
        return self.milvus.search_by_region(chrom, start, end)

    def summarize_gene(self, gene: str) -> str:
        """Generate a summary of variants in a gene."""
        variants = self.search_gene(gene)

        if not variants:
            return f"No variants found in gene {gene}."

        # Build summary prompt
        prompt = f"""Summarize the following {len(variants)} variants found in the {gene} gene:

"""
        for v in variants[:20]:  # Limit to first 20
            prompt += f"- {v.get('text_summary', 'No summary')}\n"

        prompt += """
Please provide:
1. Overview of variant types found
2. Most significant variants (if any)
3. Overall assessment of the gene's variant burden in this sample"""

        return self.llm.generate(
            prompt=prompt,
            system_prompt=self.system_prompt,
        )


def create_rag_engine(
    milvus_host: str = None,
    milvus_port: int = None,
    collection_name: str = "genomic_evidence",
    embedding_model: str = "BAAI/bge-small-en-v1.5",
    llm_provider: str = "anthropic",
    llm_model: Optional[str] = None,
    **kwargs
) -> RAGEngine:
    """
    Factory function to create a fully configured RAG engine.
    """
    # Initialize Milvus client
    milvus = MilvusClient(
        host=milvus_host,
        port=milvus_port,
        collection_name=collection_name,
    )
    milvus.connect()
    milvus.load()

    # Initialize embedder
    embedder = EvidenceEmbedder(model_name=embedding_model)

    # Initialize LLM client
    llm = LLMClient.create(
        provider=llm_provider,
        model=llm_model,
        **kwargs
    )

    # Create RAG engine
    return RAGEngine(
        milvus_client=milvus,
        embedder=embedder,
        llm_client=llm,
        top_k=kwargs.get("top_k", 10),
        score_threshold=kwargs.get("score_threshold", 0.5),
    )
