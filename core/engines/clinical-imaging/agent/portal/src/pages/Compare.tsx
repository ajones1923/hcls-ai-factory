import { useEffect, useState } from 'react';
import {
  Play,
  Loader2,
  CheckCircle,
  GitCompare,
  Dna,
  AlertTriangle,
} from 'lucide-react';
import { fetchDemoCases, runDemoCase, describeApiError } from '../lib/api';

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

const severityStyles: Record<string, string> = {
  critical: 'bg-red-500/10 text-red-400 border border-red-500/30',
  urgent: 'bg-orange-500/10 text-orange-400 border border-orange-500/30',
  significant: 'bg-yellow-500/10 text-yellow-400 border border-yellow-500/30',
  routine: 'bg-green-500/10 text-green-400 border border-green-500/30',
  normal: 'bg-green-500/10 text-green-300 border border-green-500/30',
};

function getSeverityStyle(severity?: string): string {
  if (!severity) return 'bg-blue-500/10 text-blue-400 border border-blue-500/30';
  return severityStyles[severity.toLowerCase()] || 'bg-blue-500/10 text-blue-400 border border-blue-500/30';
}

const severityRank: Record<string, number> = {
  critical: 4,
  urgent: 3,
  significant: 2,
  routine: 1,
  normal: 0,
};

function PatientColumn({
  label,
  demoCases,
  selectedCase,
  setSelectedCase,
  result,
  loading,
  onRun,
}: {
  label: string;
  demoCases: DemoCase[];
  selectedCase: string;
  setSelectedCase: (v: string) => void;
  result: DemoCaseResult | null;
  loading: boolean;
  onRun: () => void;
}) {
  const wr = result?.workflow_result;
  return (
    <div className="flex-1 min-w-0 flex flex-col">
      {/* Column header */}
      <div className="px-5 py-3 border-b border-white/[0.06] bg-white/[0.02]">
        <span className="text-xs font-semibold text-[#76B900]">{label}</span>
      </div>

      {/* Selector + Run */}
      <div className="px-5 py-4 border-b border-white/[0.06] space-y-3">
        <select
          value={selectedCase}
          onChange={(e) => setSelectedCase(e.target.value)}
          className="w-full bg-[#0E1117] text-sm text-white border border-white/[0.12] rounded-lg px-3 py-2 focus:outline-none focus:border-[#76B900]/50 cursor-pointer"
        >
          <option value="">Select a case...</option>
          {demoCases.map((dc) => (
            <option key={dc.case_id} value={dc.case_id}>
              {dc.title}
            </option>
          ))}
        </select>
        <button
          onClick={onRun}
          disabled={loading || !selectedCase}
          className="flex items-center justify-center gap-2 w-full bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-40 text-black text-sm font-semibold px-4 py-2.5 rounded-lg transition-all duration-200 cursor-pointer"
        >
          {loading ? (
            <><Loader2 size={14} className="animate-spin" /> Running...</>
          ) : (
            <><Play size={14} /> Run Case</>
          )}
        </button>
      </div>

      {/* Results */}
      {result && (
        <div className="px-5 py-4 flex-1 space-y-4 overflow-y-auto">
          {/* Title */}
          <div className="flex items-center gap-2">
            <CheckCircle size={14} className="text-[#76B900] shrink-0" />
            <span className="text-sm font-semibold text-white truncate">{result.title}</span>
          </div>

          {/* Classification */}
          {wr?.classification && (
            <div className="space-y-2">
              <div className={`inline-flex items-center px-4 py-2 rounded-lg text-lg font-bold border ${getSeverityStyle(wr.severity)}`}>
                {wr.classification}
              </div>
              {wr.severity && (
                <div>
                  <span className={`text-xs font-semibold px-2.5 py-1 rounded-lg ${getSeverityStyle(wr.severity)}`}>
                    {wr.severity.toUpperCase()}
                  </span>
                </div>
              )}
            </div>
          )}

          {/* Findings */}
          {wr?.findings && wr.findings.length > 0 && (
            <div>
              <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] mb-2 font-medium">
                Findings ({wr.findings.length})
              </h4>
              <div className="space-y-1.5">
                {wr.findings.map((f, i) => (
                  <div key={i} className="bg-[#0E1117] rounded-lg p-2.5 border border-white/[0.06] text-xs text-[#E0E0E0] leading-relaxed">
                    {f.description}
                    {f.severity && (
                      <span className={`ml-2 text-[10px] font-medium px-1.5 py-0.5 rounded-full ${getSeverityStyle(f.severity)}`}>
                        {f.severity}
                      </span>
                    )}
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Genomic context */}
          {result.genomic_context && (
            <div>
              <h4 className="text-[10px] uppercase tracking-wider text-[#9CA3AF] mb-2 font-medium flex items-center gap-1">
                <Dna size={10} /> Genomic Context
              </h4>
              {result.genomic_context.genes && result.genomic_context.genes.length > 0 && (
                <div className="flex flex-wrap gap-1.5 mb-2">
                  {result.genomic_context.genes.map((gene) => (
                    <span key={gene} className="text-[11px] font-mono font-semibold bg-purple-500/15 text-purple-300 border border-purple-500/30 px-2.5 py-0.5 rounded-full">
                      {gene}
                    </span>
                  ))}
                </div>
              )}
              {result.genomic_context.relevance && (
                <p className="text-[11px] text-[#9CA3AF] leading-relaxed">{result.genomic_context.relevance}</p>
              )}
            </div>
          )}

          {/* Inference time */}
          {wr?.inference_time_ms != null && (
            <div className="text-[10px] text-[#9CA3AF] font-mono">
              Inference: {Number(wr.inference_time_ms).toFixed(1)}ms
            </div>
          )}
        </div>
      )}

      {!result && !loading && (
        <div className="flex-1 flex items-center justify-center p-8">
          <p className="text-sm text-white/20 text-center">Select a case and click Run</p>
        </div>
      )}

      {loading && (
        <div className="flex-1 flex items-center justify-center p-8">
          <Loader2 size={24} className="animate-spin text-[#76B900]" />
        </div>
      )}
    </div>
  );
}

export default function Compare() {
  const [demoCases, setDemoCases] = useState<DemoCase[]>([]);
  const [selectedA, setSelectedA] = useState('');
  const [selectedB, setSelectedB] = useState('');
  const [resultA, setResultA] = useState<DemoCaseResult | null>(null);
  const [resultB, setResultB] = useState<DemoCaseResult | null>(null);
  const [loadingA, setLoadingA] = useState(false);
  const [loadingB, setLoadingB] = useState(false);

  useEffect(() => {
    fetchDemoCases()
      .then((data) => {
        const cases = Array.isArray(data) ? data : data.cases || [];
        setDemoCases(cases);
      })
      .catch((e) => describeApiError('Loading demo cases', e));
  }, []);

  const runA = async () => {
    if (!selectedA) return;
    setLoadingA(true);
    setResultA(null);
    try {
      const res = await runDemoCase(selectedA);
      setResultA(res);
    } catch { /* ignore */ }
    finally { setLoadingA(false); }
  };

  const runB = async () => {
    if (!selectedB) return;
    setLoadingB(true);
    setResultB(null);
    try {
      const res = await runDemoCase(selectedB);
      setResultB(res);
    } catch { /* ignore */ }
    finally { setLoadingB(false); }
  };

  // Compare summary
  const bothDone = resultA && resultB;
  const sevA = resultA?.workflow_result?.severity?.toLowerCase() || '';
  const sevB = resultB?.workflow_result?.severity?.toLowerCase() || '';
  const genesA = resultA?.genomic_context?.genes || [];
  const genesB = resultB?.genomic_context?.genes || [];
  const commonGenes = genesA.filter((g) => genesB.includes(g));
  const wfA = resultA?.workflow_name || '';
  const wfB = resultB?.workflow_name || '';

  return (
    <div className="space-y-8">
      {/* Header */}
      <div>
        <h1 className="text-2xl font-bold text-white tracking-tight flex items-center gap-3">
          <div className="p-2.5 bg-[#76B900]/10 rounded-xl border border-[#76B900]/20">
            <GitCompare size={24} className="text-[#76B900]" />
          </div>
          Multi-Patient Comparison
        </h1>
        <p className="text-sm text-[#9CA3AF] mt-1.5 ml-[52px]">
          Side-by-side workflow comparison across demo cases
        </p>
      </div>

      {/* Two-column comparison */}
      <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] overflow-hidden shadow-xl shadow-black/30 flex flex-col md:flex-row min-h-[500px]">
        <PatientColumn
          label="Patient A"
          demoCases={demoCases}
          selectedCase={selectedA}
          setSelectedCase={setSelectedA}
          result={resultA}
          loading={loadingA}
          onRun={runA}
        />
        <div className="w-px bg-white/[0.08] hidden md:block" />
        <div className="h-px bg-white/[0.08] md:hidden" />
        <PatientColumn
          label="Patient B"
          demoCases={demoCases}
          selectedCase={selectedB}
          setSelectedCase={setSelectedB}
          result={resultB}
          loading={loadingB}
          onRun={runB}
        />
      </div>

      {/* Compare Summary */}
      {bothDone && (
        <div className="bg-[#1A1D23] rounded-xl border border-purple-500/20 overflow-hidden shadow-xl shadow-black/30">
          <div className="px-5 py-3 border-b border-purple-500/15 bg-purple-500/5 flex items-center gap-2">
            <GitCompare size={16} className="text-purple-400" />
            <span className="text-sm font-semibold text-white">Comparison Summary</span>
          </div>
          <div className="p-5 space-y-3">
            {/* Modality / workflow difference */}
            {wfA !== wfB && (
              <div className="flex items-start gap-2">
                <AlertTriangle size={12} className="text-yellow-400 mt-0.5 shrink-0" />
                <p className="text-xs text-[#E0E0E0]">
                  <span className="text-white font-medium">Different workflows: </span>
                  <span className="text-blue-400">{wfA.replace(/_/g, ' ')}</span>
                  {' vs '}
                  <span className="text-blue-400">{wfB.replace(/_/g, ' ')}</span>
                </p>
              </div>
            )}
            {wfA === wfB && wfA && (
              <div className="flex items-start gap-2">
                <CheckCircle size={12} className="text-[#76B900] mt-0.5 shrink-0" />
                <p className="text-xs text-[#E0E0E0]">
                  <span className="text-white font-medium">Same workflow: </span>
                  <span className="text-[#76B900]">{wfA.replace(/_/g, ' ')}</span>
                </p>
              </div>
            )}

            {/* Severity comparison */}
            {sevA && sevB && (
              <div className="flex items-start gap-2">
                {sevA !== sevB ? (
                  <AlertTriangle size={12} className="text-yellow-400 mt-0.5 shrink-0" />
                ) : (
                  <CheckCircle size={12} className="text-[#76B900] mt-0.5 shrink-0" />
                )}
                <p className="text-xs text-[#E0E0E0]">
                  <span className="text-white font-medium">Severity: </span>
                  <span className={`font-semibold ${getSeverityStyle(sevA).includes('red') ? 'text-red-400' : getSeverityStyle(sevA).includes('orange') ? 'text-orange-400' : 'text-green-400'}`}>
                    {sevA.toUpperCase()}
                  </span>
                  {' vs '}
                  <span className={`font-semibold ${getSeverityStyle(sevB).includes('red') ? 'text-red-400' : getSeverityStyle(sevB).includes('orange') ? 'text-orange-400' : 'text-green-400'}`}>
                    {sevB.toUpperCase()}
                  </span>
                  {sevA !== sevB && (
                    <span className="text-[#9CA3AF] ml-1">
                      ({(severityRank[sevA] || 0) > (severityRank[sevB] || 0) ? 'Patient A higher risk' : 'Patient B higher risk'})
                    </span>
                  )}
                </p>
              </div>
            )}

            {/* Genomic targets */}
            {(genesA.length > 0 || genesB.length > 0) && (
              <div className="flex items-start gap-2">
                <Dna size={12} className="text-purple-400 mt-0.5 shrink-0" />
                <div className="text-xs text-[#E0E0E0]">
                  <span className="text-white font-medium">Genomic targets: </span>
                  {genesA.length > 0 && (
                    <span>A: <span className="text-purple-300 font-mono">{genesA.join(', ')}</span> </span>
                  )}
                  {genesB.length > 0 && (
                    <span>B: <span className="text-purple-300 font-mono">{genesB.join(', ')}</span></span>
                  )}
                  {commonGenes.length > 0 && (
                    <div className="mt-1 flex items-center gap-1.5">
                      <CheckCircle size={10} className="text-[#76B900]" />
                      <span className="text-[#76B900]">Common: {commonGenes.join(', ')}</span>
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
}
