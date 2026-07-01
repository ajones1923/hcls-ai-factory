import React, { useEffect, useState } from 'react';
import {
  Play,
  Loader2,
  AlertTriangle,
  CheckCircle,
  Stethoscope,
  Heart,
  Brain,
  Eye,
  Scan,
  Activity,
  CircleDot,
  Search,
  Microscope,
  User,
  Dna,
  Zap,
  Download,
} from 'lucide-react';
import { fetchWorkflows, fetchDemoCases, runDemoCase, runWorkflow, generateReport } from '../lib/api';
import VolumeViewer3D from '../components/VolumeViewer3D';
import ChromosomeViewer from '../components/ChromosomeViewer';

interface Workflow {
  name: string;
  modality: string;
  body_region: string;
  target_latency_sec?: number;
  models_used?: string[];
  mock_mode?: boolean;
  scoring_system?: string;
}

interface DemoCase {
  case_id: string;
  title: string;
  clinical_scenario?: string;
  workflow_name?: string;
  expected_severity?: string;
  age?: number;
  sex?: string;
}

interface Finding {
  description: string;
  severity?: string;
  location?: string;
}

interface GenomicContext {
  genes?: string[];
  relevance?: string;
  cross_modal_queries?: string[];
}

interface WorkflowResult {
  classification?: string;
  severity?: string;
  findings?: Finding[];
  measurements?: Record<string, number>;
  inference_time_ms?: number;
  nim_services_used?: string[];
  is_mock?: boolean;
  [key: string]: unknown;
}

interface DemoCaseResult {
  case_id?: string;
  title?: string;
  workflow_name?: string;
  workflow_result?: WorkflowResult;
  genomic_context?: GenomicContext;
  talking_points?: string[];
}

const workflowIcons: Record<string, typeof Stethoscope> = {
  ct_head_hemorrhage: Brain,
  ct_chest_lung_nodule: Scan,
  ct_coronary_angiography: Heart,
  cxr_rapid_findings: Stethoscope,
  mri_brain_ms_lesion: Brain,
  mri_prostate_pirads: CircleDot,
  breast_birads: Eye,
  thyroid_tirads: Search,
  liver_lirads: Activity,
};

const workflowScoring: Record<string, string> = {
  ct_chest_lung_nodule: 'Lung-RADS',
  breast_birads: 'BI-RADS',
  mri_prostate_pirads: 'PI-RADS',
  thyroid_tirads: 'TI-RADS',
  liver_lirads: 'LI-RADS',
  ct_coronary_angiography: 'CAD-RADS',
  ct_head_hemorrhage: 'Marshall CT',
};

const modalityColors: Record<string, { bg: string; text: string; border: string }> = {
  ct: { bg: 'bg-blue-500/10', text: 'text-blue-400', border: 'border-blue-500/30' },
  mri: { bg: 'bg-purple-500/10', text: 'text-purple-400', border: 'border-purple-500/30' },
  cxr: { bg: 'bg-cyan-500/10', text: 'text-cyan-400', border: 'border-cyan-500/30' },
  mammography: { bg: 'bg-pink-500/10', text: 'text-pink-400', border: 'border-pink-500/30' },
  ultrasound: { bg: 'bg-teal-500/10', text: 'text-teal-400', border: 'border-teal-500/30' },
};

const severityStyles: Record<string, string> = {
  critical: 'bg-red-500/10 text-red-400 border border-red-500/30',
  urgent: 'bg-orange-500/10 text-orange-400 border border-orange-500/30',
  significant: 'bg-yellow-500/10 text-yellow-400 border border-yellow-500/30',
  routine: 'bg-green-500/10 text-green-400 border border-green-500/30',
  normal: 'bg-green-500/10 text-green-300 border border-green-500/30',
};

// Map workflows to sample images + AI segmentation overlays
const API_HOST = `${window.location.protocol}//${window.location.hostname}:8524`;
const workflowImages: Record<string, { primary: string; label: string; animation?: string; animLabel?: string }> = {
  ct_head_hemorrhage: {
    primary: `${API_HOST}/segmentation/highres_ct_head_overlay.png`,
    label: 'CT Head — AI Segmentation (Bone + Brain + CSF + Hemorrhage)',
    animation: `${API_HOST}/segmentation/highres_ct_head_segmented.gif`,
    animLabel: 'CT Head — AI Segmentation Scroll',
  },
  ct_chest_lung_nodule: {
    primary: `${API_HOST}/segmentation/highres_ct_chest_overlay.png`,
    label: 'CT Chest — AI Segmentation (Lung + Bone + Nodule)',
    animation: `${API_HOST}/segmentation/highres_ct_chest_segmented.gif`,
    animLabel: 'CT Chest — AI Segmentation Scroll',
  },
  ct_coronary_angiography: {
    primary: `${API_HOST}/segmentation/cardiac_ct_overlay.png`,
    label: 'Cardiac CT — AI Segmentation (Heart + Vessels + Calcification)',
    animation: `${API_HOST}/segmentation/cardiac_ct_segmented.gif`,
    animLabel: 'Cardiac CT — AI Segmentation Scroll (Real 512x512)',
  },
  cxr_rapid_findings: {
    primary: `${API_HOST}/images/annotated/cxr_pneumonia_bilateral_annotated.png`,
    label: 'Chest X-Ray — AI Pathology Detection',
  },
  mri_brain_ms_lesion: {
    primary: `${API_HOST}/segmentation/highres_brain_flair_overlay.png`,
    label: 'MRI FLAIR — AI Segmentation (Brain + CSF + MS Lesions)',
    animation: `${API_HOST}/segmentation/highres_brain_flair_segmented.gif`,
    animLabel: 'MRI Brain — AI Segmentation Scroll',
  },
  mri_prostate_pirads: { primary: `${API_HOST}/images/organ_ct_001.png`, label: 'MRI Prostate — AI PI-RADS Scoring' },
  breast_birads: { primary: `${API_HOST}/images/breast_us_000.png`, label: 'Breast Imaging — AI Mass Detection' },
  thyroid_tirads: { primary: `${API_HOST}/images/breast_us_002.png`, label: 'Thyroid Ultrasound — AI Nodule Analysis' },
  liver_lirads: { primary: `${API_HOST}/images/organ_ct_000.png`, label: 'Liver CT — AI HCC Screening' },
};

const severityDot: Record<string, string> = {
  critical: 'bg-red-500',
  urgent: 'bg-orange-500',
  significant: 'bg-yellow-500',
  routine: 'bg-green-500',
  normal: 'bg-green-400',
};

const severityAvatarColor: Record<string, string> = {
  critical: 'text-red-400 bg-red-500/15 border-red-500/40',
  urgent: 'text-orange-400 bg-orange-500/15 border-orange-500/40',
  significant: 'text-yellow-400 bg-yellow-500/15 border-yellow-500/40',
  routine: 'text-green-400 bg-green-500/15 border-green-500/40',
  normal: 'text-green-300 bg-green-500/15 border-green-500/40',
};

function getModalityStyle(modality: string) {
  const lower = modality.toLowerCase();
  for (const [key, style] of Object.entries(modalityColors)) {
    if (lower.includes(key)) return style;
  }
  return { bg: 'bg-white/5', text: 'text-white/60', border: 'border-white/10' };
}

function getWorkflowIcon(name: string) {
  for (const [key, Icon] of Object.entries(workflowIcons)) {
    if (name.toLowerCase().includes(key.toLowerCase())) return Icon;
  }
  // fallback by body region keyword
  const lower = name.toLowerCase();
  if (lower.includes('brain') || lower.includes('head') || lower.includes('neuro')) return Brain;
  if (lower.includes('heart') || lower.includes('cardiac') || lower.includes('coronary')) return Heart;
  if (lower.includes('lung') || lower.includes('chest')) return Scan;
  if (lower.includes('breast') || lower.includes('mammo')) return Eye;
  if (lower.includes('prostate')) return CircleDot;
  if (lower.includes('liver')) return Activity;
  if (lower.includes('thyroid')) return Search;
  return Stethoscope;
}

function getScoring(name: string): string | null {
  for (const [key, scoring] of Object.entries(workflowScoring)) {
    if (name.toLowerCase().includes(key.toLowerCase())) return scoring;
  }
  return null;
}

function getSeverityStyle(severity?: string): string {
  if (!severity) return 'bg-blue-500/10 text-blue-400 border border-blue-500/30';
  return severityStyles[severity.toLowerCase()] || 'bg-blue-500/10 text-blue-400 border border-blue-500/30';
}

function getTimelineSteps(workflowName: string, classification?: string, genomicGenes?: string[]) {
  const steps: { label: string; active: boolean; critical?: boolean }[] = [
    { label: 'DICOM Received', active: true },
    { label: 'AI Segmentation', active: true },
    { label: classification || 'Classification', active: true, critical: !!(classification && (classification.includes('4B') || classification.includes('critical') || classification.includes('LR-5') || classification.includes('TR5'))) },
    { label: 'Scoring Applied', active: true },
  ];

  if (genomicGenes && genomicGenes.length > 0) {
    steps.push({ label: 'Genomic Trigger', active: true, critical: false });
    steps.push({ label: genomicGenes.slice(0, 2).join(', '), active: true });
  }

  steps.push({ label: 'Report Generated', active: true });
  return steps;
}

function getSeverityDot(severity?: string): string {
  if (!severity) return 'bg-blue-500';
  return severityDot[severity.toLowerCase()] || 'bg-blue-500';
}

function getAvatarColor(severity?: string): string {
  if (!severity) return 'text-blue-400 bg-blue-500/15 border-blue-500/40';
  return severityAvatarColor[severity.toLowerCase()] || 'text-blue-400 bg-blue-500/15 border-blue-500/40';
}

export default function Workflows() {
  const [workflows, setWorkflows] = useState<Workflow[]>([]);
  const [demoCases, setDemoCases] = useState<DemoCase[]>([]);
  const [selectedWorkflow, setSelectedWorkflow] = useState<string | null>(null);
  const [result, setResult] = useState<DemoCaseResult | null>(null);
  const [loading, setLoading] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [imageMode, setImageMode] = useState<'ai' | 'raw'>('ai');
  const [cardiacVersion, setCardiacVersion] = useState<1 | 2>(2);
  const [selectedGene, setSelectedGene] = useState<string | null>(null);
  const [pipelineStage, setPipelineStage] = useState(0);

  useEffect(() => {
    fetchWorkflows()
      .then((data) => {
        const wfs = data.workflows || data;
        setWorkflows(Array.isArray(wfs) ? wfs : []);
      })
      .catch(() => {});
    fetchDemoCases()
      .then((data) => {
        const cases = Array.isArray(data) ? data : data.cases || [];
        setDemoCases(cases);
      })
      .catch(() => {});
  }, []);

  // Enhancement 4: Pipeline stage animation during loading
  useEffect(() => {
    if (loading) {
      setPipelineStage(0);
      const interval = setInterval(() => {
        setPipelineStage(prev => prev < 5 ? prev + 1 : prev);
      }, 400);
      return () => clearInterval(interval);
    }
  }, [loading]);

  // Enhancement 1: Reset image mode when new result loads
  useEffect(() => {
    if (result) setImageMode('ai');
  }, [result]);

  const handleRunDemo = async (id: string) => {
    setLoading(id);
    setError(null);
    setResult(null);
    try {
      const res = await runDemoCase(id);
      setResult(res);
      // Auto-scroll to results
      setTimeout(() => {
        document.getElementById('result-panel')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }, 100);
    } catch (err: unknown) {
      const msg = err instanceof Error ? err.message : 'Failed to run demo case';
      setError(msg);
    } finally {
      setLoading(null);
    }
  };

  const handleRunWorkflow = async (name: string) => {
    setLoading(name);
    setError(null);
    setResult(null);
    setSelectedWorkflow(name);
    try {
      const res = await runWorkflow(name);
      setResult({
        workflow_name: name,
        workflow_result: res,
        title: name.replace(/_/g, ' ').replace(/\b\w/g, (c: string) => c.toUpperCase()),
      });
      setTimeout(() => {
        document.getElementById('result-panel')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }, 100);
    } catch (err: unknown) {
      const msg = err instanceof Error ? err.message : 'Failed to run workflow';
      setError(msg);
    } finally {
      setLoading(null);
    }
  };

  const scoringSystems = new Set(workflows.map(wf => getScoring(wf.name)).filter(Boolean));

  return (
    <div className="space-y-8">
      {/* ── Header ─────────────────────────────────────────────── */}
      <div className="flex items-end justify-between">
        <div>
          <h1 className="text-2xl font-bold text-white tracking-tight flex items-center gap-3">
            <div className="p-2.5 bg-[#76B900]/10 rounded-xl border border-[#76B900]/20">
              <Microscope size={24} className="text-[#76B900]" />
            </div>
            Imaging Workflows
          </h1>
          <p className="text-sm text-[#9CA3AF] mt-1.5 ml-[52px]">
            {workflows.length > 0 ? `${workflows.length} workflows` : 'Loading workflows...'}
            {scoringSystems.size > 0 && (
              <span className="text-white/40 mx-2">|</span>
            )}
            {scoringSystems.size > 0 && `${scoringSystems.size} scoring systems`}
          </p>
        </div>
      </div>

      {/* ── Workflow Grid ──────────────────────────────────────── */}
      <div>
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          {workflows.map((wf) => {
            const Icon = getWorkflowIcon(wf.name);
            const isSelected = selectedWorkflow === wf.name;
            const mStyle = getModalityStyle(wf.modality);
            const scoring = wf.scoring_system || getScoring(wf.name);
            return (
              <button
                key={wf.name}
                onClick={() => setSelectedWorkflow(isSelected ? null : wf.name)}
                className={`text-left rounded-xl border p-5 transition-all duration-200 cursor-pointer group ${
                  isSelected
                    ? 'border-[#76B900]/50 bg-[#76B900]/5 shadow-[0_0_24px_rgba(118,185,0,0.12)]'
                    : `border-white/[0.08] bg-[#1A1D23] hover:border-[#76B900]/30 hover:shadow-lg hover:shadow-black/30`
                }`}
              >
                {/* Icon + Name */}
                <div className="flex items-start gap-4 mb-4">
                  <div className={`p-3 rounded-xl ${mStyle.bg} border ${mStyle.border} transition-colors group-hover:border-[#76B900]/30`}>
                    <Icon size={28} className={isSelected ? 'text-[#76B900]' : mStyle.text} />
                  </div>
                  <div className="flex-1 min-w-0 pt-0.5">
                    <h3 className="text-sm font-semibold text-white leading-tight">
                      {wf.name.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase())}
                    </h3>
                    {wf.models_used && wf.models_used.length > 0 && (
                      <p className="text-[11px] text-[#9CA3AF] mt-1 line-clamp-1">
                        {wf.models_used.join(' + ')}
                      </p>
                    )}
                  </div>
                </div>

                {/* Badges */}
                <div className="flex flex-wrap gap-1.5">
                  <span className={`text-[10px] font-medium px-2.5 py-1 rounded-full ${mStyle.bg} ${mStyle.text} border ${mStyle.border}`}>
                    {wf.modality.toUpperCase()}
                  </span>
                  <span className="text-[10px] font-medium px-2.5 py-1 rounded-full bg-white/5 text-white/60 border border-white/10">
                    {wf.body_region}
                  </span>
                  {scoring && (
                    <span className="text-[10px] font-medium px-2.5 py-1 rounded-full bg-[#76B900]/8 text-[#76B900] border border-[#76B900]/25">
                      {scoring}
                    </span>
                  )}
                  {wf.target_latency_sec != null && (
                    <span className="text-[10px] font-medium px-2.5 py-1 rounded-full bg-orange-500/8 text-orange-400 border border-orange-500/20">
                      {wf.target_latency_sec}s target
                    </span>
                  )}
                </div>

                {/* Run button */}
                <div className="mt-4 flex items-center gap-2">
                  <button
                    onClick={(e) => { e.stopPropagation(); handleRunWorkflow(wf.name); }}
                    disabled={loading === wf.name}
                    className="flex items-center gap-2 px-4 py-2 bg-[#76B900] hover:bg-[#5a9400] text-black text-xs font-semibold rounded-lg transition-all duration-200 cursor-pointer disabled:opacity-50"
                  >
                    {loading === wf.name ? (
                      <><Loader2 size={12} className="animate-spin" /> Running...</>
                    ) : (
                      <><Play size={12} /> Run Workflow</>
                    )}
                  </button>
                  {isSelected && result?.workflow_result && (
                    <span className="text-[10px] text-[#76B900]">
                      <CheckCircle size={12} className="inline mr-1" />Results below
                    </span>
                  )}
                </div>
              </button>
            );
          })}
          {workflows.length === 0 && (
            <div className="col-span-3 grid grid-cols-3 gap-4">
              {[...Array(6)].map((_, i) => (
                <div key={i} className="h-40 rounded-xl bg-[#1A1D23] border border-white/[0.06] animate-pulse" />
              ))}
            </div>
          )}
        </div>
      </div>

      {/* ── Demo Cases ─────────────────────────────────────────── */}
      <div>
        <div className="flex items-center gap-3 mb-4">
          <h2 className="text-lg font-semibold text-white">Clinical Demo Cases</h2>
          {demoCases.length > 0 && (
            <span className="text-[11px] font-medium px-2.5 py-1 rounded-full bg-[#76B900]/10 text-[#76B900] border border-[#76B900]/20">
              {demoCases.length} cases
            </span>
          )}
        </div>
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          {demoCases.map((dc) => {
            const avatarColor = getAvatarColor(dc.expected_severity);
            return (
              <div
                key={dc.case_id}
                className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 hover:border-[#76B900]/30 transition-all duration-200 shadow-lg shadow-black/20 flex flex-col"
              >
                {/* Header with avatar */}
                <div className="flex items-start gap-3 mb-3">
                  <div className={`p-2 rounded-lg border ${avatarColor} shrink-0`}>
                    <User size={18} />
                  </div>
                  <div className="flex-1 min-w-0">
                    <h3 className="text-sm font-semibold text-white leading-tight">{dc.title}</h3>
                    {(dc.age || dc.sex) && (
                      <p className="text-[11px] text-[#9CA3AF] mt-0.5">
                        {[dc.age && `${dc.age}y`, dc.sex].filter(Boolean).join(' / ')}
                      </p>
                    )}
                  </div>
                </div>

                {/* Clinical scenario */}
                {dc.clinical_scenario && (
                  <p className="text-xs text-[#9CA3AF] mb-3 line-clamp-2 leading-relaxed">{dc.clinical_scenario}</p>
                )}

                {/* Badges */}
                <div className="flex flex-wrap items-center gap-2 mb-4">
                  {dc.workflow_name && (
                    <span className="text-[10px] font-medium bg-blue-500/10 text-blue-400 border border-blue-500/25 px-2.5 py-1 rounded-full">
                      {dc.workflow_name.replace(/_/g, ' ')}
                    </span>
                  )}
                  {dc.expected_severity && (
                    <span className={`text-xs font-semibold px-3 py-1 rounded-full ${getSeverityStyle(dc.expected_severity)}`}>
                      {dc.expected_severity.charAt(0).toUpperCase() + dc.expected_severity.slice(1)}
                    </span>
                  )}
                </div>

                {/* Run button */}
                <div className="mt-auto">
                  <button
                    onClick={() => handleRunDemo(dc.case_id)}
                    disabled={loading === dc.case_id}
                    className="flex items-center justify-center gap-2 w-full bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-50 text-black text-sm font-semibold px-4 py-2.5 rounded-lg transition-all duration-200 cursor-pointer shadow-md shadow-[#76B900]/20"
                  >
                    {loading === dc.case_id ? (
                      <Loader2 size={16} className="animate-spin" />
                    ) : (
                      <Play size={16} />
                    )}
                    {loading === dc.case_id ? 'Running...' : 'Run Case'}
                  </button>
                </div>
              </div>
            );
          })}
          {demoCases.length === 0 && (
            <div className="col-span-3 grid grid-cols-3 gap-4">
              {[...Array(3)].map((_, i) => (
                <div key={i} className="h-52 rounded-xl bg-[#1A1D23] border border-white/[0.06] animate-pulse" />
              ))}
            </div>
          )}
        </div>
      </div>

      {/* ── Error ──────────────────────────────────────────────── */}
      {error && (
        <div className="bg-red-500/10 border border-red-500/30 rounded-xl p-4 flex items-center gap-3 shadow-lg shadow-black/20">
          <AlertTriangle size={18} className="text-red-400 shrink-0" />
          <p className="text-sm text-red-400">{error}</p>
        </div>
      )}

      {/* ── Enhancement 4: Pipeline Stage Animation ──────────── */}
      {loading && (
        <div className="bg-[#1A1D23] rounded-xl border border-[#76B900]/20 p-8 shadow-xl shadow-black/30">
          <div className="flex items-center justify-between mb-6">
            <h3 className="text-white font-semibold">Processing Pipeline</h3>
            <Loader2 size={16} className="animate-spin text-[#76B900]" />
          </div>
          <div className="flex items-center gap-2">
            {['DICOM Parse', 'Segment', 'Classify', 'Score', 'Cross-Modal', 'Report'].map((stage, i) => (
              <React.Fragment key={stage}>
                <div className={`flex-1 rounded-lg px-3 py-2 text-center text-[10px] font-medium transition-all duration-500 ${
                  i <= pipelineStage ? 'bg-[#76B900]/20 text-[#76B900] border border-[#76B900]/30' : 'bg-white/5 text-white/30 border border-white/10'
                }`}>
                  {stage}
                </div>
                {i < 5 && (
                  <div className={`w-4 h-0.5 transition-all duration-500 ${
                    i < pipelineStage ? 'bg-[#76B900]' : 'bg-white/10'
                  }`} />
                )}
              </React.Fragment>
            ))}
          </div>
        </div>
      )}

      {/* ── Result Panel ───────────────────────────────────────── */}
      {result && (
        <div id="result-panel" className={`bg-[#1A1D23] rounded-xl border border-[#76B900]/20 overflow-hidden shadow-xl shadow-black/30 ${result.workflow_result?.severity === 'critical' ? 'animate-critical-pulse' : ''}`}>
          {/* Result Header */}
          <div className="bg-[#76B900]/5 border-b border-[#76B900]/15 px-6 py-4 flex items-center justify-between">
            <div className="flex items-center gap-3">
              <CheckCircle size={20} className="text-[#76B900]" />
              <div>
                <h2 className="text-lg font-bold text-white">
                  {result.title || 'Analysis Complete'}
                </h2>
                {result.workflow_name && (
                  <p className="text-xs text-[#9CA3AF]">
                    {result.workflow_name.replace(/_/g, ' ')}
                  </p>
                )}
              </div>
            </div>
            <div className="flex items-center gap-2">
              {result.workflow_result?.inference_time_ms != null && (
                <span className="text-[11px] text-[#9CA3AF] bg-white/5 border border-white/10 px-3 py-1 rounded-full font-mono">
                  {Number(result.workflow_result.inference_time_ms).toFixed(1)}ms
                </span>
              )}
              <button
                onClick={async () => {
                  try {
                    const reportData = await generateReport({
                      format: 'markdown',
                      question: `Clinical imaging report for ${result.title || result.workflow_name}: ${result.workflow_result?.classification || 'Analysis'}`
                    });
                    const blob = new Blob([reportData.content || JSON.stringify(reportData, null, 2)], { type: 'text/markdown' });
                    const url = URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url;
                    a.download = `${result.workflow_name || 'imaging'}_report.md`;
                    a.click();
                    URL.revokeObjectURL(url);
                  } catch (err) {
                    console.error('Report download failed:', err);
                  }
                }}
                className="flex items-center gap-1.5 text-xs text-[#9CA3AF] hover:text-[#76B900] bg-white/5 hover:bg-[#76B900]/10 border border-white/10 hover:border-[#76B900]/30 px-3 py-1.5 rounded-lg transition-all duration-200 cursor-pointer"
              >
                <Download size={12} />
                <span>Report</span>
              </button>
            </div>
          </div>

          {/* Enhancement 2: Critical Alert Bar */}
          {result.workflow_result?.severity === 'critical' && (
            <div className="bg-red-500/10 border-l-4 border-red-500 px-4 py-2 flex items-center gap-3">
              <span className="relative flex h-3 w-3">
                <span className="animate-ping absolute inline-flex h-full w-full rounded-full bg-red-400 opacity-75"></span>
                <span className="relative inline-flex rounded-full h-3 w-3 bg-red-500"></span>
              </span>
              <span className="text-red-400 text-sm font-semibold tracking-wide">CRITICAL FINDING DETECTED</span>
            </div>
          )}

          {/* Patient Context (from demo cases) */}
          {(result as Record<string, unknown>).patient_demographics && (
            <div className="px-6 py-3 bg-white/[0.02] border-b border-white/[0.06] flex flex-wrap items-center gap-4">
              <div className="flex items-center gap-2">
                <User size={14} className="text-[#9CA3AF]" />
                <span className="text-sm text-white font-medium">
                  {((result as Record<string, unknown>).patient_demographics as Record<string, unknown>)?.age}
                  {((result as Record<string, unknown>).patient_demographics as Record<string, unknown>)?.sex === 'Male' ? 'M' : ((result as Record<string, unknown>).patient_demographics as Record<string, unknown>)?.sex === 'Female' ? 'F' : ''}
                </span>
              </div>
              {((result as Record<string, unknown>).patient_demographics as Record<string, unknown>)?.weight_kg && (
                <span className="text-xs text-[#9CA3AF]">
                  {String(((result as Record<string, unknown>).patient_demographics as Record<string, unknown>).weight_kg)}kg
                </span>
              )}
              {((result as Record<string, unknown>).patient_demographics as Record<string, unknown>)?.history && (
                <div className="flex flex-wrap gap-1">
                  {(((result as Record<string, unknown>).patient_demographics as Record<string, unknown>).history as string[]).map((h: string) => (
                    <span key={h} className="text-[10px] bg-amber-500/10 text-amber-400 border border-amber-500/20 px-2 py-0.5 rounded-full">
                      {h.replace(/_/g, ' ')}
                    </span>
                  ))}
                </div>
              )}
              {(result as Record<string, unknown>).clinical_scenario && (
                <p className="text-xs text-[#9CA3AF] flex-1 min-w-[200px]">
                  {String((result as Record<string, unknown>).clinical_scenario).slice(0, 150)}...
                </p>
              )}
            </div>
          )}

          {/* Enhancement 3: Patient Journey Timeline */}
          {(result as Record<string, unknown>).patient_demographics && (
            <div className="px-6 py-4 border-b border-white/[0.06]">
              <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] mb-3 font-medium">Patient Journey</h4>
              <div className="relative">
                {/* Timeline bar */}
                <div className="absolute top-3 left-0 right-0 h-0.5 bg-white/10" />
                <div className="absolute top-3 left-0 h-0.5 bg-[#76B900] transition-all duration-1000" style={{width: '100%'}} />

                {/* Timeline nodes */}
                <div className="flex justify-between relative">
                  {getTimelineSteps(
                    result.workflow_name || '',
                    result.workflow_result?.classification,
                    result.genomic_context?.genes
                  ).map((step, i) => (
                    <div key={i} className="flex flex-col items-center" style={{maxWidth: '120px'}}>
                      <div className={`w-6 h-6 rounded-full border-2 flex items-center justify-center text-[10px] font-bold z-10 ${
                        step.critical ? 'bg-red-500 border-red-500 text-white' :
                        step.active ? 'bg-[#76B900] border-[#76B900] text-black' :
                        'bg-[#1A1D23] border-white/20 text-white/40'
                      }`}>
                        {i + 1}
                      </div>
                      <span className={`text-[10px] mt-1.5 text-center leading-tight ${
                        step.critical ? 'text-red-400 font-semibold' : 'text-[#9CA3AF]'
                      }`}>
                        {step.label}
                      </span>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          )}

          {/* Result Body */}
          <div className="p-6">
            {/* Medical Image Display */}
            {(() => {
              const wfName = result.workflow_name || result.workflow_result?.workflow_name || '';
              const isCardiac = wfName === 'ct_coronary_angiography';

              // For cardiac workflow, select images based on version tab
              const cardiacV1Images = {
                primary: `${API_HOST}/images/annotated/ct_coronary_annotated.png`,
                label: 'CTA Coronary — AI Stenosis Analysis',
                animation: `${API_HOST}/segmentation/highres_ct_chest_segmented.gif`,
                animLabel: 'CT Chest — AI Segmentation Scroll (1024x1024)',
              };
              const cardiacV2Images = {
                primary: `${API_HOST}/segmentation/cardiac_ct_overlay.png`,
                label: 'Cardiac CT — AI Segmentation (Heart + Vessels + Calcification)',
                animation: `${API_HOST}/segmentation/cardiac_ct_segmented.gif`,
                animLabel: 'Cardiac CT — AI Segmentation Scroll (Real 512x512 Clinical CT)',
              };

              const imgData = isCardiac
                ? (cardiacVersion === 1 ? cardiacV1Images : cardiacV2Images)
                : workflowImages[wfName];
              if (!imgData) return null;
              const hasAnim = !!imgData.animation;
              const volumeWorkflows: Record<string, 'brain' | 'chest' | 'head'> = {
                ct_head_hemorrhage: 'head',
                mri_brain_ms_lesion: 'brain',
                ct_chest_lung_nodule: 'chest',
                ct_coronary_angiography: 'chest',
              };
              const volumeType = volumeWorkflows[wfName] || null;
              const colCount = hasAnim && volumeType ? 3 : hasAnim ? 2 : volumeType ? 2 : 1;
              return (
                <div className="mb-6">
                  {/* Version Tabs — only for cardiac workflow */}
                  {isCardiac && (
                    <div className="flex items-center gap-0 mb-0">
                      <button
                        onClick={() => setCardiacVersion(1)}
                        className={`flex items-center gap-2 px-5 py-2.5 text-sm font-medium rounded-t-xl border border-b-0 transition-all duration-200 cursor-pointer ${
                          cardiacVersion === 1
                            ? 'bg-[#1A1D23] text-white border-white/[0.08]'
                            : 'bg-[#0E1117] text-white/40 border-transparent hover:text-white/60'
                        }`}
                      >
                        <span className={`w-5 h-5 rounded-full text-[10px] font-bold flex items-center justify-center ${
                          cardiacVersion === 1 ? 'bg-[#76B900] text-black' : 'bg-white/10 text-white/40'
                        }`}>1</span>
                        Standard Demo
                      </button>
                      <button
                        onClick={() => setCardiacVersion(2)}
                        className={`flex items-center gap-2 px-5 py-2.5 text-sm font-medium rounded-t-xl border border-b-0 transition-all duration-200 cursor-pointer ${
                          cardiacVersion === 2
                            ? 'bg-[#1A1D23] text-white border-white/[0.08]'
                            : 'bg-[#0E1117] text-white/40 border-transparent hover:text-white/60'
                        }`}
                      >
                        <span className={`w-5 h-5 rounded-full text-[10px] font-bold flex items-center justify-center ${
                          cardiacVersion === 2 ? 'bg-[#76B900] text-black' : 'bg-white/10 text-white/40'
                        }`}>2</span>
                        Advanced Imagery
                        <span className="text-[9px] bg-blue-500/15 text-blue-400 border border-blue-500/25 px-1.5 py-0.5 rounded-full ml-1">1024px</span>
                      </button>
                    </div>
                  )}

                  <div className={`bg-[#0E1117] ${isCardiac ? 'rounded-b-xl rounded-tr-xl' : 'rounded-xl'} border border-white/[0.08] overflow-hidden`}>
                    <div className="flex items-center justify-between px-4 py-2 border-b border-white/[0.06] bg-white/[0.02]">
                      <div className="flex items-center gap-2">
                        <Scan size={14} className="text-[#76B900]" />
                        <span className="text-xs font-medium text-white">Medical Imaging — AI Analysis</span>
                        {isCardiac && (
                          <span className={`text-[10px] px-2 py-0.5 rounded-full font-medium ${
                            cardiacVersion === 1
                              ? 'bg-white/5 text-white/50 border border-white/10'
                              : 'bg-blue-500/10 text-blue-400 border border-blue-500/25'
                          }`}>
                            {cardiacVersion === 1 ? 'V1 — Demo-Matched' : 'V2 — Real Clinical CT'}
                          </span>
                        )}
                      </div>
                      <div className="flex items-center gap-3">
                        {/* Enhancement 1: Before/After Toggle */}
                        <div className="flex items-center gap-1 bg-[#0E1117] rounded-lg p-1">
                          <button
                            onClick={() => setImageMode('raw')}
                            className={`px-3 py-1.5 text-xs font-medium rounded-md transition-all duration-300 cursor-pointer ${
                              imageMode === 'raw' ? 'bg-white/10 text-white' : 'text-white/40 hover:text-white/60'
                            }`}
                          >
                            Raw Image
                          </button>
                          <button
                            onClick={() => setImageMode('ai')}
                            className={`px-3 py-1.5 text-xs font-medium rounded-md transition-all duration-300 cursor-pointer ${
                              imageMode === 'ai' ? 'bg-[#76B900]/20 text-[#76B900] border border-[#76B900]/30' : 'text-white/40 hover:text-white/60'
                            }`}
                          >
                            AI Analysis
                          </button>
                        </div>
                        <span className="text-[10px] text-[#9CA3AF] bg-white/5 border border-white/10 px-2 py-0.5 rounded-full">
                          VISTA-3D • NV-Segment-CT
                        </span>
                      </div>
                    </div>

                    {/* Version description banner for cardiac */}
                    {isCardiac && (
                      <div className={`px-4 py-2 border-b border-white/[0.06] ${
                        cardiacVersion === 1 ? 'bg-white/[0.02]' : 'bg-blue-500/[0.03]'
                      }`}>
                        <p className="text-[11px] text-[#9CA3AF] leading-relaxed">
                          {cardiacVersion === 1
                            ? 'Standard view using chest CT segmentation overlays. Matches the imagery used in the platform demo walkthrough for consistent presentation.'
                            : 'Advanced view using real 512x512 clinical cardiac CT from TCIA. AI segmentation highlights heart chambers, aorta, and coronary calcification at full diagnostic resolution.'}
                        </p>
                      </div>
                    )}

                    <div className={`grid ${colCount === 3 ? 'grid-cols-3' : colCount === 2 ? 'grid-cols-2' : 'grid-cols-1'} gap-0`}>
                      {/* Annotated AI Detection with Before/After toggle */}
                      <div className="flex flex-col items-center justify-center p-4 bg-black/30 border-r border-white/[0.06]">
                        <img
                          src={imgData.primary}
                          alt={imgData.label}
                          className={`max-h-[280px] rounded-lg shadow-2xl shadow-black/50 object-contain transition-all duration-500 ${
                            imageMode === 'raw' ? 'grayscale brightness-110' : ''
                          }`}
                          onError={(e) => { (e.target as HTMLImageElement).style.display = 'none'; }}
                        />
                        <span className="text-[10px] text-[#9CA3AF] mt-2">
                          {imageMode === 'raw' ? 'Original Grayscale' : imgData.label}
                        </span>
                      </div>
                      {/* Animated Volume Scroll */}
                      {hasAnim && (
                        <div className="flex flex-col items-center justify-center p-4 bg-black/40">
                          <img
                            src={imgData.animation}
                            alt={imgData.animLabel || 'Volume scroll'}
                            className="max-h-[280px] rounded-lg shadow-2xl shadow-black/50 object-contain"
                            onError={(e) => { (e.target as HTMLImageElement).style.display = 'none'; }}
                          />
                          <span className="text-[10px] text-[#9CA3AF] mt-2">{imgData.animLabel || 'Volume Scroll'}</span>
                        </div>
                      )}
                      {/* 3D AI Volume Rendering */}
                      {volumeType && (
                        <div className="flex flex-col items-center justify-center p-4 bg-black/40 border-l border-white/[0.06]">
                          <VolumeViewer3D volumeType={volumeType} height={280} />
                          <span className="text-[10px] text-[#9CA3AF] mt-2">3D AI Visualization</span>
                        </div>
                      )}
                    </div>
                  </div>
                </div>
              );
            })()}

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              {/* Left: Classification */}
              <div className="space-y-5">
                {/* Large Classification Badge */}
                {result.workflow_result?.classification && (
                  <div className="space-y-3">
                    <div className={`inline-flex items-center gap-3 px-6 py-3 rounded-xl text-2xl font-bold border-2 ${getSeverityStyle(result.workflow_result.severity)}`}>
                      {result.workflow_result.classification}
                    </div>
                    <div className="flex items-center gap-3">
                      {result.workflow_result.severity && (
                        <span className={`text-sm font-semibold px-3 py-1.5 rounded-lg ${getSeverityStyle(result.workflow_result.severity)}`}>
                          {result.workflow_result.severity.toUpperCase()}
                        </span>
                      )}
                      {result.workflow_result.is_mock && (
                        <span className="text-[10px] text-yellow-500 bg-yellow-500/10 border border-yellow-500/30 px-2 py-0.5 rounded-full font-medium">
                          SIMULATED
                        </span>
                      )}
                    </div>
                  </div>
                )}

                {/* NIM Services */}
                {result.workflow_result?.nim_services_used && result.workflow_result.nim_services_used.length > 0 && (
                  <div>
                    <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] mb-2">NIM Services Used</h4>
                    <div className="flex flex-wrap gap-1.5">
                      {result.workflow_result.nim_services_used.map((svc) => (
                        <span key={svc} className="text-[10px] font-mono bg-[#76B900]/8 text-[#76B900] border border-[#76B900]/20 px-2 py-0.5 rounded">
                          {svc}
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>

              {/* Right: Findings + Measurements */}
              <div className="space-y-5">
                {/* Findings */}
                {result.workflow_result?.findings && result.workflow_result.findings.length > 0 && (
                  <div>
                    <h3 className="text-sm font-semibold text-white mb-3 flex items-center gap-2">
                      <span>Findings</span>
                      <span className="text-[10px] bg-white/5 text-[#9CA3AF] border border-white/10 px-2 py-0.5 rounded-full">
                        {result.workflow_result.findings.length}
                      </span>
                    </h3>
                    <div className="space-y-2">
                      {result.workflow_result.findings.map((f, i) => (
                        <div
                          key={i}
                          className="bg-[#0E1117] rounded-lg p-3 border border-white/[0.06] hover:border-white/[0.12] transition-all duration-200 flex items-start gap-3"
                        >
                          <div className={`w-2 h-2 rounded-full mt-1.5 shrink-0 ${getSeverityDot(f.severity)}`} />
                          <div className="flex-1 min-w-0">
                            <p className="text-sm text-[#E0E0E0] leading-relaxed">{f.description}</p>
                            <div className="flex gap-2 mt-1.5">
                              {f.severity && (
                                <span className={`text-[10px] font-medium px-2 py-0.5 rounded-full ${getSeverityStyle(f.severity)}`}>
                                  {f.severity}
                                </span>
                              )}
                              {f.location && (
                                <span className="text-[10px] text-[#9CA3AF] bg-white/5 px-2 py-0.5 rounded-full">
                                  {f.location}
                                </span>
                              )}
                            </div>
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                )}

                {/* Measurements */}
                {result.workflow_result?.measurements && Object.keys(result.workflow_result.measurements).length > 0 && (
                  <div>
                    <h3 className="text-sm font-semibold text-white mb-3">Measurements</h3>
                    <div className="overflow-x-auto rounded-lg border border-white/[0.06]">
                      <table className="w-full text-sm">
                        <thead>
                          <tr className="border-b border-white/[0.08] bg-white/[0.03]">
                            <th className="text-left py-2.5 px-4 text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">
                              Metric
                            </th>
                            <th className="text-right py-2.5 px-4 text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">
                              Value
                            </th>
                          </tr>
                        </thead>
                        <tbody>
                          {Object.entries(result.workflow_result.measurements).map(([key, val]) => (
                            <tr key={key} className="border-b border-white/[0.04] hover:bg-white/[0.02] transition-colors">
                              <td className="py-2.5 px-4 text-[#E0E0E0]">
                                {key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase())}
                              </td>
                              <td className="py-2.5 px-4 text-[#76B900] font-semibold font-mono text-right">
                                {typeof val === 'number' ? val.toFixed(1) : String(val)}
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>
                )}
              </div>
            </div>

            {/* ── Enhancement 5: Cross-Modal Bridge Animation ──── */}
            {result.genomic_context && (
              <div className="mt-6 relative overflow-hidden rounded-xl border border-purple-500/20">
                {/* Animated bridge header */}
                <div className="bg-gradient-to-r from-blue-500/10 via-purple-500/10 to-[#76B900]/10 px-5 py-3 border-b border-purple-500/15">
                  <div className="flex items-center justify-between">
                    <div className="flex items-center gap-8">
                      {/* Imaging side */}
                      <div className="flex items-center gap-2">
                        <div className="w-8 h-8 rounded-lg bg-blue-500/20 border border-blue-500/30 flex items-center justify-center">
                          <Scan size={16} className="text-blue-400" />
                        </div>
                        <div>
                          <span className="text-[10px] uppercase tracking-wider text-blue-400 font-semibold">Engine 4</span>
                          <p className="text-[10px] text-white/50">Imaging</p>
                        </div>
                      </div>

                      {/* Animated bridge */}
                      <div className="flex items-center gap-1">
                        {[0,1,2,3,4].map(idx => (
                          <div key={idx} className="w-2 h-2 rounded-full bg-purple-400 animate-pulse" style={{animationDelay: `${idx * 200}ms`}} />
                        ))}
                        <Zap size={14} className="text-purple-400 mx-1" />
                        {[0,1,2,3,4].map(idx => (
                          <div key={idx} className="w-2 h-2 rounded-full bg-[#76B900] animate-pulse" style={{animationDelay: `${(idx+5) * 200}ms`}} />
                        ))}
                      </div>

                      {/* Genomics side */}
                      <div className="flex items-center gap-2">
                        <div className="w-8 h-8 rounded-lg bg-[#76B900]/20 border border-[#76B900]/30 flex items-center justify-center">
                          <Dna size={16} className="text-[#76B900]" />
                        </div>
                        <div>
                          <span className="text-[10px] uppercase tracking-wider text-[#76B900] font-semibold">Engine 2</span>
                          <p className="text-[10px] text-white/50">Genomics</p>
                        </div>
                      </div>
                    </div>

                    <span className="text-[10px] text-purple-300 bg-purple-500/10 border border-purple-500/20 px-2.5 py-1 rounded-full font-medium">
                      CROSS-MODAL TRIGGER
                    </span>
                  </div>
                </div>

                {/* Genomic context content */}
                <div className="p-5 bg-[#0E1117]/50">
                  {result.genomic_context.relevance && (
                    <p className="text-xs text-[#9CA3AF] mb-3 leading-relaxed">
                      {result.genomic_context.relevance}
                    </p>
                  )}
                  {result.genomic_context.genes && result.genomic_context.genes.length > 0 && (
                    <div className="flex flex-wrap gap-2 mb-3">
                      {result.genomic_context.genes.map((gene) => (
                        <button
                          key={gene}
                          onClick={() => setSelectedGene(gene)}
                          className="text-xs font-mono font-semibold bg-purple-500/15 text-purple-300 border border-purple-500/30 px-3 py-1 rounded-full hover:bg-purple-500/25 hover:border-purple-400/50 hover:text-purple-200 transition-all duration-200 cursor-pointer flex items-center gap-1.5"
                          title={`View ${gene} on chromosome`}
                        >
                          <Dna size={10} />
                          {gene}
                        </button>
                      ))}
                    </div>
                  )}
                  {result.genomic_context.cross_modal_queries && result.genomic_context.cross_modal_queries.length > 0 && (
                    <div>
                      <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] mb-2">Gene Queries</h4>
                      <div className="flex flex-wrap gap-2">
                        {result.genomic_context.cross_modal_queries.map((q, i) => (
                          <span
                            key={i}
                            className="text-[11px] bg-[#76B900]/10 text-[#76B900] border border-[#76B900]/25 px-3 py-1 rounded-full"
                          >
                            <Zap size={10} className="inline mr-1" />
                            {q}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}

            {/* ── Talking Points ───────────────────────────────── */}
            {result.talking_points && result.talking_points.length > 0 && (
              <div className="mt-5">
                <h3 className="text-sm font-semibold text-white mb-3">Key Points</h3>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-2">
                  {result.talking_points.map((tp, i) => (
                    <div key={i} className="flex items-start gap-2.5 bg-[#0E1117] rounded-lg p-3 border border-white/[0.06]">
                      <span className="text-[#76B900] font-bold text-xs mt-0.5 shrink-0">{i + 1}</span>
                      <p className="text-xs text-[#9CA3AF] leading-relaxed">{tp}</p>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* ── Chromosome Viewer Modal ───────────────────────── */}
      {selectedGene && (
        <ChromosomeViewer gene={selectedGene} onClose={() => setSelectedGene(null)} />
      )}
    </div>
  );
}
