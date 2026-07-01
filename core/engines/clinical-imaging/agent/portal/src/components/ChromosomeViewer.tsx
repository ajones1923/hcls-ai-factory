import React, { useState } from 'react';
import { X, Dna, MapPin, AlertTriangle, BookOpen, ExternalLink } from 'lucide-react';

interface GeneInfo {
  symbol: string;
  name: string;
  chromosome: string;
  band: string;
  position_mb: number;       // Approximate position in megabases
  chr_length_mb: number;     // Chromosome length in Mb
  variant?: string;
  consequence?: string;
  clinical_significance?: string;
  condition?: string;
  inheritance?: string;
  drug_interaction?: string;
  omim?: string;
  description: string;
}

const GENE_DATABASE: Record<string, GeneInfo> = {
  LDLR: {
    symbol: 'LDLR',
    name: 'Low-Density Lipoprotein Receptor',
    chromosome: '19',
    band: 'p13.2',
    position_mb: 11.09,
    chr_length_mb: 58.6,
    variant: 'p.Gly592Glu (c.1775G>A)',
    consequence: 'Missense — disrupts LDL receptor folding',
    clinical_significance: 'Pathogenic',
    condition: 'Familial Hypercholesterolemia (FH)',
    inheritance: 'Autosomal Dominant',
    drug_interaction: 'Reduced response to statins alone; PCSK9 inhibitors (evolocumab, alirocumab) recommended as add-on therapy',
    omim: '606945',
    description: 'Encodes the LDL receptor responsible for clearing LDL cholesterol from the bloodstream. Pathogenic variants cause familial hypercholesterolemia, affecting 1 in 250 individuals. Each first-degree relative has a 50% chance of carrying the variant.',
  },
  PCSK9: {
    symbol: 'PCSK9',
    name: 'Proprotein Convertase Subtilisin/Kexin Type 9',
    chromosome: '1',
    band: 'p32.3',
    position_mb: 55.04,
    chr_length_mb: 248.9,
    variant: 'p.Asp374Tyr (c.1120G>T)',
    consequence: 'Gain-of-function — increased LDLR degradation',
    clinical_significance: 'Pathogenic',
    condition: 'Familial Hypercholesterolemia (FH)',
    inheritance: 'Autosomal Dominant',
    drug_interaction: 'PCSK9 inhibitors directly target this protein; enhanced therapeutic response expected',
    omim: '607786',
    description: 'Gain-of-function variants increase degradation of LDL receptors, raising circulating LDL cholesterol. Loss-of-function variants are cardioprotective — this biology led directly to the development of PCSK9 inhibitor drugs.',
  },
  APOB: {
    symbol: 'APOB',
    name: 'Apolipoprotein B',
    chromosome: '2',
    band: 'p24.1',
    position_mb: 21.0,
    chr_length_mb: 242.2,
    variant: 'p.Arg3527Gln (c.10580G>A)',
    consequence: 'Missense — impairs LDL receptor binding',
    clinical_significance: 'Pathogenic',
    condition: 'Familial Defective ApoB-100 (FDB)',
    inheritance: 'Autosomal Dominant',
    drug_interaction: 'Moderate statin response; may require combination therapy with ezetimibe',
    omim: '107730',
    description: 'ApoB-100 is the primary ligand for LDL receptor binding. The Arg3527Gln variant reduces receptor binding affinity by ~95%, causing LDL accumulation similar to FH but typically with a milder phenotype.',
  },
  LPA: {
    symbol: 'LPA',
    name: 'Lipoprotein(a)',
    chromosome: '6',
    band: 'q25.3',
    position_mb: 160.5,
    chr_length_mb: 170.8,
    variant: 'KIV-2 repeat expansion',
    consequence: 'Elevated Lp(a) levels (>50 mg/dL)',
    clinical_significance: 'Risk Factor',
    condition: 'Elevated Lipoprotein(a) — Independent CVD Risk',
    inheritance: 'Autosomal Dominant (>90% heritable)',
    drug_interaction: 'Not responsive to statins; antisense oligonucleotides (pelacarsen) in Phase 3 trials',
    omim: '152200',
    description: 'Lp(a) is a genetically determined, independent risk factor for atherosclerotic cardiovascular disease. The KIV-2 repeat number is inversely correlated with Lp(a) levels. Currently no approved targeted therapy, but pelacarsen and olpasiran are in late-stage clinical trials.',
  },
  '9p21': {
    symbol: '9p21.3',
    name: 'Chromosome 9p21.3 Risk Locus (CDKN2A/CDKN2B-AS1)',
    chromosome: '9',
    band: 'p21.3',
    position_mb: 22.0,
    chr_length_mb: 138.4,
    variant: 'rs10757274 / rs2383206 risk alleles',
    consequence: 'Non-coding — alters ANRIL lncRNA expression',
    clinical_significance: 'Risk Factor',
    condition: 'Coronary Artery Disease Susceptibility',
    inheritance: 'Complex (additive risk)',
    drug_interaction: 'No direct pharmacogenomic interaction; risk modifiable through aggressive risk factor management',
    omim: '611082',
    description: 'The 9p21.3 locus was the first genomic region identified in GWAS for coronary artery disease. The risk variants alter expression of ANRIL (CDKN2B-AS1), a long non-coding RNA involved in vascular smooth muscle cell proliferation and inflammation. Carriers of both risk alleles have ~1.6x increased CAD risk.',
  },
  SLCO1B1: {
    symbol: 'SLCO1B1',
    name: 'Solute Carrier Organic Anion Transporter 1B1',
    chromosome: '12',
    band: 'p12.2',
    position_mb: 21.3,
    chr_length_mb: 133.3,
    variant: 'rs4149056 (c.521T>C, Val174Ala)',
    consequence: 'Reduced OATP1B1 transporter function',
    clinical_significance: 'Pharmacogenomic',
    condition: 'Statin-Induced Myopathy Risk',
    inheritance: 'Autosomal Co-dominant',
    drug_interaction: 'TC heterozygous: reduce simvastatin (max 20mg), consider rosuvastatin or pravastatin. CC homozygous: avoid simvastatin entirely; use alternative statin at lower starting dose.',
    omim: '604843',
    description: 'OATP1B1 is the primary hepatic uptake transporter for statins. The Val174Ala variant reduces statin clearance from blood, increasing systemic exposure and myopathy risk by 4.5-fold (heterozygous) to 17-fold (homozygous). CPIC guideline provides evidence-based dosing recommendations.',
  },
  CYP2C19: {
    symbol: 'CYP2C19',
    name: 'Cytochrome P450 Family 2 Subfamily C Member 19',
    chromosome: '10',
    band: 'q23.33',
    position_mb: 94.8,
    chr_length_mb: 133.8,
    variant: '*2/*2 (c.681G>A, splicing defect)',
    consequence: 'No functional enzyme — poor metabolizer',
    clinical_significance: 'Pharmacogenomic',
    condition: 'Clopidogrel Resistance',
    inheritance: 'Autosomal Co-dominant',
    drug_interaction: 'Poor metabolizer: clopidogrel is ineffective. Use prasugrel 10mg or ticagrelor 90mg BID instead. Critical if patient undergoes PCI with stent placement.',
    omim: '124020',
    description: 'CYP2C19 converts the prodrug clopidogrel to its active antiplatelet metabolite. Poor metabolizers (*2/*2) have no enzyme function, rendering clopidogrel ineffective and increasing the risk of stent thrombosis by 3-4 fold. CPIC guideline recommends alternative antiplatelet agents.',
  },
  BRCA1: {
    symbol: 'BRCA1',
    name: 'BRCA1 DNA Repair Associated',
    chromosome: '17',
    band: 'q21.31',
    position_mb: 43.04,
    chr_length_mb: 83.3,
    variant: 'Pathogenic variant detected',
    consequence: 'Loss of DNA double-strand break repair',
    clinical_significance: 'Pathogenic',
    condition: 'Hereditary Breast and Ovarian Cancer Syndrome',
    inheritance: 'Autosomal Dominant',
    drug_interaction: 'PARP inhibitors (olaparib, talazoparib) show synthetic lethality in BRCA1-deficient tumors',
    omim: '113705',
    description: 'BRCA1 is a tumor suppressor involved in homologous recombination DNA repair. Pathogenic variants confer 60-80% lifetime risk of breast cancer and 40-60% risk of ovarian cancer. Cascade testing of at-risk relatives is strongly recommended.',
  },
  BRCA2: {
    symbol: 'BRCA2',
    name: 'BRCA2 DNA Repair Associated',
    chromosome: '13',
    band: 'q13.1',
    position_mb: 32.3,
    chr_length_mb: 114.4,
    variant: 'Pathogenic variant detected',
    consequence: 'Loss of DNA double-strand break repair',
    clinical_significance: 'Pathogenic',
    condition: 'Hereditary Breast and Ovarian Cancer Syndrome',
    inheritance: 'Autosomal Dominant',
    drug_interaction: 'PARP inhibitors effective; also associated with prostate and pancreatic cancer risk',
    omim: '600185',
    description: 'BRCA2 partners with RAD51 in homologous recombination repair. Pathogenic variants confer 45-65% lifetime breast cancer risk, 15-25% ovarian cancer risk, and elevated prostate cancer risk in males.',
  },
  EGFR: {
    symbol: 'EGFR',
    name: 'Epidermal Growth Factor Receptor',
    chromosome: '7',
    band: 'p11.2',
    position_mb: 55.1,
    chr_length_mb: 159.3,
    variant: 'Exon 19 deletion / L858R',
    consequence: 'Constitutive kinase activation',
    clinical_significance: 'Oncogenic Driver',
    condition: 'Non-Small Cell Lung Cancer (NSCLC)',
    inheritance: 'Somatic',
    drug_interaction: 'EGFR TKIs: osimertinib (3rd gen, preferred), erlotinib, gefitinib. T790M resistance mutation → osimertinib.',
    omim: '131550',
    description: 'EGFR activating mutations are found in 10-15% of NSCLC in Western populations and 40-50% in East Asian populations. These mutations predict dramatic response to EGFR tyrosine kinase inhibitors, with response rates >70%.',
  },
};

// Chromosome band colors for ideogram
const BAND_PATTERNS: Record<string, { bands: Array<{ start: number; end: number; color: string; label?: string }>; centromere: number }> = {
  '1':  { centromere: 0.49, bands: [{start:0,end:0.1,color:'#555'},{start:0.1,end:0.2,color:'#ccc'},{start:0.2,end:0.35,color:'#888'},{start:0.35,end:0.49,color:'#ccc'},{start:0.51,end:0.65,color:'#888'},{start:0.65,end:0.8,color:'#ccc'},{start:0.8,end:0.9,color:'#555'},{start:0.9,end:1,color:'#ccc'}] },
  '2':  { centromere: 0.39, bands: [{start:0,end:0.1,color:'#888'},{start:0.1,end:0.25,color:'#ccc'},{start:0.25,end:0.39,color:'#555'},{start:0.41,end:0.55,color:'#ccc'},{start:0.55,end:0.7,color:'#888'},{start:0.7,end:0.85,color:'#ccc'},{start:0.85,end:1,color:'#555'}] },
  '6':  { centromere: 0.39, bands: [{start:0,end:0.12,color:'#888'},{start:0.12,end:0.25,color:'#ccc'},{start:0.25,end:0.39,color:'#555'},{start:0.41,end:0.55,color:'#ccc'},{start:0.55,end:0.7,color:'#888'},{start:0.7,end:0.85,color:'#ccc'},{start:0.85,end:1,color:'#555'}] },
  '7':  { centromere: 0.38, bands: [{start:0,end:0.15,color:'#555'},{start:0.15,end:0.28,color:'#ccc'},{start:0.28,end:0.38,color:'#888'},{start:0.4,end:0.55,color:'#ccc'},{start:0.55,end:0.7,color:'#555'},{start:0.7,end:0.85,color:'#ccc'},{start:0.85,end:1,color:'#888'}] },
  '9':  { centromere: 0.36, bands: [{start:0,end:0.12,color:'#888'},{start:0.12,end:0.22,color:'#ccc'},{start:0.22,end:0.36,color:'#555'},{start:0.38,end:0.5,color:'#ccc'},{start:0.5,end:0.65,color:'#888'},{start:0.65,end:0.8,color:'#ccc'},{start:0.8,end:1,color:'#555'}] },
  '10': { centromere: 0.3,  bands: [{start:0,end:0.1,color:'#555'},{start:0.1,end:0.2,color:'#ccc'},{start:0.2,end:0.3,color:'#888'},{start:0.32,end:0.5,color:'#ccc'},{start:0.5,end:0.65,color:'#555'},{start:0.65,end:0.8,color:'#ccc'},{start:0.8,end:1,color:'#888'}] },
  '12': { centromere: 0.27, bands: [{start:0,end:0.12,color:'#888'},{start:0.12,end:0.27,color:'#ccc'},{start:0.29,end:0.45,color:'#555'},{start:0.45,end:0.6,color:'#ccc'},{start:0.6,end:0.75,color:'#888'},{start:0.75,end:0.9,color:'#ccc'},{start:0.9,end:1,color:'#555'}] },
  '13': { centromere: 0.15, bands: [{start:0,end:0.08,color:'#888'},{start:0.08,end:0.15,color:'#ccc'},{start:0.17,end:0.35,color:'#555'},{start:0.35,end:0.5,color:'#ccc'},{start:0.5,end:0.65,color:'#888'},{start:0.65,end:0.8,color:'#ccc'},{start:0.8,end:1,color:'#555'}] },
  '17': { centromere: 0.28, bands: [{start:0,end:0.1,color:'#555'},{start:0.1,end:0.2,color:'#ccc'},{start:0.2,end:0.28,color:'#888'},{start:0.3,end:0.45,color:'#ccc'},{start:0.45,end:0.6,color:'#555'},{start:0.6,end:0.8,color:'#ccc'},{start:0.8,end:1,color:'#888'}] },
  '19': { centromere: 0.45, bands: [{start:0,end:0.15,color:'#888'},{start:0.15,end:0.3,color:'#ccc'},{start:0.3,end:0.45,color:'#555'},{start:0.47,end:0.6,color:'#ccc'},{start:0.6,end:0.75,color:'#888'},{start:0.75,end:0.9,color:'#ccc'},{start:0.9,end:1,color:'#555'}] },
};

function getDefaultBands(centromere: number) {
  return {
    centromere,
    bands: [
      {start: 0, end: centromere * 0.4, color: '#888'},
      {start: centromere * 0.4, end: centromere * 0.75, color: '#ccc'},
      {start: centromere * 0.75, end: centromere, color: '#555'},
      {start: centromere + 0.02, end: centromere + 0.3, color: '#ccc'},
      {start: centromere + 0.3, end: centromere + 0.6, color: '#888'},
      {start: centromere + 0.6, end: 1, color: '#555'},
    ],
  };
}

const significanceColors: Record<string, { bg: string; text: string; border: string }> = {
  'Pathogenic': { bg: 'bg-red-500/15', text: 'text-red-400', border: 'border-red-500/30' },
  'Pharmacogenomic': { bg: 'bg-blue-500/15', text: 'text-blue-400', border: 'border-blue-500/30' },
  'Risk Factor': { bg: 'bg-amber-500/15', text: 'text-amber-400', border: 'border-amber-500/30' },
  'Oncogenic Driver': { bg: 'bg-purple-500/15', text: 'text-purple-400', border: 'border-purple-500/30' },
};

interface ChromosomeViewerProps {
  gene: string;
  onClose: () => void;
}

export default function ChromosomeViewer({ gene, onClose }: ChromosomeViewerProps) {
  const [activeTab, setActiveTab] = useState<'ideogram' | 'details'>('ideogram');
  const info = GENE_DATABASE[gene];

  if (!info) {
    return (
      <div className="fixed inset-0 z-50 bg-black/70 backdrop-blur-sm flex items-center justify-center p-4" onClick={onClose}>
        <div className="bg-[#1A1D23] rounded-xl border border-white/10 p-6 max-w-md" onClick={e => e.stopPropagation()}>
          <p className="text-white/60 text-sm">No genomic data available for {gene}</p>
          <button onClick={onClose} className="mt-4 text-xs text-[#76B900] hover:underline cursor-pointer">Close</button>
        </div>
      </div>
    );
  }

  const chrPattern = BAND_PATTERNS[info.chromosome] || getDefaultBands(0.35);
  const genePositionFrac = info.position_mb / info.chr_length_mb;
  const sigStyle = significanceColors[info.clinical_significance || ''] || significanceColors['Risk Factor'];

  // SVG dimensions
  const svgWidth = 600;
  const svgHeight = 220;
  const chrX = 50;
  const chrWidth = svgWidth - 100;
  const chrY = 80;
  const chrHeight = 28;
  const markerX = chrX + genePositionFrac * chrWidth;

  return (
    <div className="fixed inset-0 z-50 bg-black/70 backdrop-blur-sm flex items-center justify-center p-4" onClick={onClose}>
      <div
        className="bg-[#12151A] rounded-2xl border border-white/10 shadow-2xl shadow-black/50 w-full max-w-3xl max-h-[90vh] overflow-y-auto"
        onClick={e => e.stopPropagation()}
      >
        {/* Header */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-white/[0.08] bg-white/[0.02]">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 rounded-xl bg-purple-500/15 border border-purple-500/30 flex items-center justify-center">
              <Dna size={20} className="text-purple-400" />
            </div>
            <div>
              <h2 className="text-lg font-bold text-white">{info.symbol}</h2>
              <p className="text-xs text-[#9CA3AF]">{info.name}</p>
            </div>
          </div>
          <button
            onClick={onClose}
            className="w-8 h-8 rounded-lg bg-white/5 hover:bg-white/10 flex items-center justify-center transition-colors cursor-pointer"
          >
            <X size={16} className="text-[#9CA3AF]" />
          </button>
        </div>

        {/* Tabs */}
        <div className="flex items-center gap-0 px-6 pt-3 border-b border-white/[0.06]">
          <button
            onClick={() => setActiveTab('ideogram')}
            className={`px-4 py-2.5 text-sm font-medium border-b-2 transition-all cursor-pointer ${
              activeTab === 'ideogram'
                ? 'text-[#76B900] border-[#76B900]'
                : 'text-white/40 border-transparent hover:text-white/60'
            }`}
          >
            Chromosome View
          </button>
          <button
            onClick={() => setActiveTab('details')}
            className={`px-4 py-2.5 text-sm font-medium border-b-2 transition-all cursor-pointer ${
              activeTab === 'details'
                ? 'text-[#76B900] border-[#76B900]'
                : 'text-white/40 border-transparent hover:text-white/60'
            }`}
          >
            Clinical Details
          </button>
        </div>

        {/* Content */}
        <div className="p-6">
          {activeTab === 'ideogram' ? (
            <div className="space-y-5">
              {/* Chromosome ideogram */}
              <div className="bg-[#0E1117] rounded-xl border border-white/[0.06] p-4">
                <div className="flex items-center justify-between mb-3">
                  <span className="text-xs font-medium text-white/60">Chromosome {info.chromosome}</span>
                  <span className="text-[10px] text-[#9CA3AF] font-mono">{info.chr_length_mb.toFixed(1)} Mb</span>
                </div>
                <svg viewBox={`0 0 ${svgWidth} ${svgHeight}`} className="w-full" style={{ maxHeight: '200px' }}>
                  {/* Chromosome body with rounded ends */}
                  <defs>
                    <clipPath id="chr-clip">
                      <rect x={chrX} y={chrY} width={chrWidth} height={chrHeight} rx={chrHeight / 2} />
                    </clipPath>
                    <filter id="glow">
                      <feGaussianBlur stdDeviation="3" result="coloredBlur"/>
                      <feMerge>
                        <feMergeNode in="coloredBlur"/>
                        <feMergeNode in="SourceGraphic"/>
                      </feMerge>
                    </filter>
                    <filter id="marker-glow">
                      <feGaussianBlur stdDeviation="4" result="coloredBlur"/>
                      <feMerge>
                        <feMergeNode in="coloredBlur"/>
                        <feMergeNode in="SourceGraphic"/>
                      </feMerge>
                    </filter>
                  </defs>

                  {/* Chromosome outline */}
                  <rect x={chrX} y={chrY} width={chrWidth} height={chrHeight} rx={chrHeight / 2}
                    fill="#2a2d35" stroke="#444" strokeWidth="1" />

                  {/* Banding pattern */}
                  <g clipPath="url(#chr-clip)">
                    {chrPattern.bands.map((band, i) => (
                      <rect
                        key={i}
                        x={chrX + band.start * chrWidth}
                        y={chrY}
                        width={(band.end - band.start) * chrWidth}
                        height={chrHeight}
                        fill={band.color}
                        opacity={0.35}
                      />
                    ))}
                  </g>

                  {/* Centromere pinch */}
                  <g>
                    <polygon
                      points={`
                        ${chrX + chrPattern.centromere * chrWidth - 6},${chrY}
                        ${chrX + chrPattern.centromere * chrWidth},${chrY + chrHeight / 2}
                        ${chrX + chrPattern.centromere * chrWidth - 6},${chrY + chrHeight}
                      `}
                      fill="#12151A"
                    />
                    <polygon
                      points={`
                        ${chrX + chrPattern.centromere * chrWidth + 6},${chrY}
                        ${chrX + chrPattern.centromere * chrWidth},${chrY + chrHeight / 2}
                        ${chrX + chrPattern.centromere * chrWidth + 6},${chrY + chrHeight}
                      `}
                      fill="#12151A"
                    />
                  </g>

                  {/* p arm / q arm labels */}
                  <text x={chrX + chrPattern.centromere * chrWidth * 0.5} y={chrY - 8}
                    textAnchor="middle" fill="#666" fontSize="10" fontFamily="monospace">p arm</text>
                  <text x={chrX + chrPattern.centromere * chrWidth + (1 - chrPattern.centromere) * chrWidth * 0.5} y={chrY - 8}
                    textAnchor="middle" fill="#666" fontSize="10" fontFamily="monospace">q arm</text>

                  {/* Gene position marker */}
                  <line x1={markerX} y1={chrY - 2} x2={markerX} y2={chrY + chrHeight + 2}
                    stroke="#76B900" strokeWidth="3" filter="url(#marker-glow)" />

                  {/* Marker triangle above */}
                  <polygon
                    points={`${markerX - 6},${chrY - 18} ${markerX + 6},${chrY - 18} ${markerX},${chrY - 6}`}
                    fill="#76B900" filter="url(#marker-glow)"
                  />

                  {/* Gene label */}
                  <text x={markerX} y={chrY - 24} textAnchor="middle" fill="#76B900" fontSize="12" fontWeight="bold" fontFamily="monospace">
                    {info.symbol}
                  </text>

                  {/* Band label */}
                  <text x={markerX} y={chrY + chrHeight + 18} textAnchor="middle" fill="#9CA3AF" fontSize="10" fontFamily="monospace">
                    {info.chromosome}{info.band}
                  </text>

                  {/* Position label */}
                  <text x={markerX} y={chrY + chrHeight + 32} textAnchor="middle" fill="#666" fontSize="9" fontFamily="monospace">
                    ~{info.position_mb.toFixed(1)} Mb
                  </text>

                  {/* Scale markers */}
                  {[0, 0.25, 0.5, 0.75, 1].map(frac => (
                    <g key={frac}>
                      <line x1={chrX + frac * chrWidth} y1={chrY + chrHeight + 40} x2={chrX + frac * chrWidth} y2={chrY + chrHeight + 46}
                        stroke="#444" strokeWidth="1" />
                      <text x={chrX + frac * chrWidth} y={chrY + chrHeight + 56}
                        textAnchor="middle" fill="#555" fontSize="8" fontFamily="monospace">
                        {(frac * info.chr_length_mb).toFixed(0)}
                      </text>
                    </g>
                  ))}
                  <text x={chrX + chrWidth / 2} y={chrY + chrHeight + 68}
                    textAnchor="middle" fill="#555" fontSize="8" fontFamily="monospace">Position (Mb)</text>
                </svg>
              </div>

              {/* Quick info cards */}
              <div className="grid grid-cols-3 gap-3">
                <div className="bg-[#0E1117] rounded-lg border border-white/[0.06] p-3">
                  <div className="flex items-center gap-2 mb-1.5">
                    <MapPin size={12} className="text-[#76B900]" />
                    <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-medium">Location</span>
                  </div>
                  <p className="text-sm text-white font-mono font-semibold">Chr {info.chromosome}{info.band}</p>
                  <p className="text-[10px] text-[#9CA3AF] mt-0.5">~{info.position_mb} Mb</p>
                </div>
                <div className="bg-[#0E1117] rounded-lg border border-white/[0.06] p-3">
                  <div className="flex items-center gap-2 mb-1.5">
                    <AlertTriangle size={12} className={sigStyle.text} />
                    <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-medium">Significance</span>
                  </div>
                  <p className={`text-sm font-semibold ${sigStyle.text}`}>{info.clinical_significance}</p>
                  <p className="text-[10px] text-[#9CA3AF] mt-0.5">{info.inheritance}</p>
                </div>
                <div className="bg-[#0E1117] rounded-lg border border-white/[0.06] p-3">
                  <div className="flex items-center gap-2 mb-1.5">
                    <BookOpen size={12} className="text-purple-400" />
                    <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-medium">OMIM</span>
                  </div>
                  <p className="text-sm text-white font-mono font-semibold">{info.omim || 'N/A'}</p>
                  {info.omim && (
                    <a
                      href={`https://omim.org/entry/${info.omim}`}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="text-[10px] text-blue-400 hover:underline flex items-center gap-1 mt-0.5"
                    >
                      View <ExternalLink size={8} />
                    </a>
                  )}
                </div>
              </div>

              {/* Variant */}
              {info.variant && (
                <div className="bg-[#0E1117] rounded-lg border border-white/[0.06] p-4">
                  <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-medium mb-2">Variant Detected</h4>
                  <p className="text-sm text-white font-mono mb-1">{info.variant}</p>
                  <p className="text-xs text-[#9CA3AF]">{info.consequence}</p>
                </div>
              )}
            </div>
          ) : (
            <div className="space-y-4">
              {/* Gene description */}
              <div className="bg-[#0E1117] rounded-lg border border-white/[0.06] p-4">
                <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-medium mb-2">Gene Function</h4>
                <p className="text-sm text-[#E0E0E0] leading-relaxed">{info.description}</p>
              </div>

              {/* Condition */}
              <div className="bg-[#0E1117] rounded-lg border border-white/[0.06] p-4">
                <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-medium mb-2">Associated Condition</h4>
                <p className="text-sm text-white font-semibold mb-1">{info.condition}</p>
                <div className="flex items-center gap-2 mt-2">
                  <span className={`text-xs font-medium px-2.5 py-1 rounded-full border ${sigStyle.bg} ${sigStyle.text} ${sigStyle.border}`}>
                    {info.clinical_significance}
                  </span>
                  <span className="text-xs text-[#9CA3AF] bg-white/5 border border-white/10 px-2.5 py-1 rounded-full">
                    {info.inheritance}
                  </span>
                </div>
              </div>

              {/* Drug interaction */}
              {info.drug_interaction && (
                <div className="bg-[#0E1117] rounded-lg border border-amber-500/20 p-4">
                  <h4 className="text-[10px] uppercase tracking-wider text-amber-400 font-medium mb-2 flex items-center gap-1.5">
                    <AlertTriangle size={10} />
                    Pharmacogenomic Interaction
                  </h4>
                  <p className="text-sm text-[#E0E0E0] leading-relaxed">{info.drug_interaction}</p>
                </div>
              )}

              {/* Variant details */}
              {info.variant && (
                <div className="bg-[#0E1117] rounded-lg border border-white/[0.06] p-4">
                  <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-medium mb-2">Variant</h4>
                  <table className="w-full text-sm">
                    <tbody>
                      <tr className="border-b border-white/[0.04]">
                        <td className="py-2 text-[#9CA3AF]">Notation</td>
                        <td className="py-2 text-white font-mono text-right">{info.variant}</td>
                      </tr>
                      <tr className="border-b border-white/[0.04]">
                        <td className="py-2 text-[#9CA3AF]">Consequence</td>
                        <td className="py-2 text-white text-right">{info.consequence}</td>
                      </tr>
                      <tr className="border-b border-white/[0.04]">
                        <td className="py-2 text-[#9CA3AF]">Location</td>
                        <td className="py-2 text-white font-mono text-right">Chr {info.chromosome}{info.band} (~{info.position_mb} Mb)</td>
                      </tr>
                      {info.omim && (
                        <tr>
                          <td className="py-2 text-[#9CA3AF]">OMIM</td>
                          <td className="py-2 text-right">
                            <a href={`https://omim.org/entry/${info.omim}`} target="_blank" rel="noopener noreferrer"
                              className="text-blue-400 hover:underline font-mono">
                              #{info.omim}
                            </a>
                          </td>
                        </tr>
                      )}
                    </tbody>
                  </table>
                </div>
              )}
            </div>
          )}
        </div>

        {/* Footer */}
        <div className="px-6 py-3 border-t border-white/[0.06] bg-white/[0.02] flex items-center justify-between">
          <span className="text-[10px] text-[#9CA3AF]">
            Cross-modal genomic enrichment — Clinical Imaging Engine
          </span>
          <span className="text-[10px] text-[#9CA3AF] font-mono">
            {info.symbol} • Chr{info.chromosome}{info.band}
          </span>
        </div>
      </div>
    </div>
  );
}

export { GENE_DATABASE };
export type { GeneInfo };
