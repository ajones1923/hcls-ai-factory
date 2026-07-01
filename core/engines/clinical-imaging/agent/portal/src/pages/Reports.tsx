import { useState } from 'react';
import { FileText, Loader2, AlertTriangle, Download, Copy, Clock, BookOpen, FileJson, FileType, Activity, Database } from 'lucide-react';
import { generateReport } from '../lib/api';

interface ReportEvidence {
  collection?: string;
  id?: string;
  score?: number;
  text_snippet?: string;
  label?: string;
}

interface ReportResult {
  title?: string;
  generated_at?: string;
  question?: string;
  answer?: string;
  content?: string;
  format?: string;
  evidence?: ReportEvidence[];
  evidence_count?: number;
  generation_time_ms?: number;
  [key: string]: unknown;
}

const quickReports = [
  'Lung cancer screening guidelines and AI detection performance',
  'Radiation dose optimization strategies for pediatric CT',
  'BI-RADS classification overview and management recommendations',
  'Comparison of AI segmentation models for cardiac MRI',
  'ACR Appropriateness Criteria for acute headache workup',
  'DiffDock molecular docking performance benchmarks',
];

interface FormatOption {
  id: string;
  label: string;
  description: string;
  icon: typeof FileText;
}

const formatOptions: FormatOption[] = [
  { id: 'markdown', label: 'Markdown', description: 'Rich text with headers and lists', icon: FileText },
  { id: 'structured', label: 'JSON', description: 'Machine-readable structured data', icon: FileJson },
  { id: 'brief', label: 'Brief', description: 'Concise summary format', icon: FileType },
  { id: 'fhir_r4', label: 'FHIR R4', description: 'Interoperable clinical standard', icon: Activity },
  { id: 'dicom_sr', label: 'DICOM SR', description: 'Structured reporting standard', icon: Database },
];

export default function Reports() {
  const [formData, setFormData] = useState({
    question: '',
    format: 'markdown',
  });
  const [result, setResult] = useState<ReportResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [copied, setCopied] = useState(false);

  const handleGenerate = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    try {
      const res = await generateReport(formData);
      setResult(res);
    } catch {
      setError('Failed to generate report');
    } finally {
      setLoading(false);
    }
  };

  const handleQuickReport = (q: string) => {
    setFormData({ ...formData, question: q });
    setResult(null);
    setError(null);
  };

  const handleCopy = () => {
    const text = result?.answer || result?.content || '';
    navigator.clipboard.writeText(text);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  const handleDownload = () => {
    const text = result?.answer || result?.content || '';
    const ext = formData.format === 'structured' ? 'json' : 'md';
    const blob = new Blob([text], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `report_${Date.now()}.${ext}`;
    a.click();
    URL.revokeObjectURL(url);
  };

  const isJsonFormat = formData.format === 'structured';

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3">
        <div className="p-2.5 bg-[#76B900]/10 rounded-xl border border-[#76B900]/20">
          <FileText size={24} className="text-[#76B900]" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-white tracking-tight">Report Generation</h1>
          <p className="text-sm text-[#9CA3AF]">Generate structured clinical imaging reports with AI and RAG-retrieved evidence</p>
        </div>
      </div>

      {/* Quick Reports */}
      <div>
        <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2.5 block">Quick Reports</span>
        <div className="flex flex-wrap gap-2">
          {quickReports.map((q) => (
            <button
              key={q}
              onClick={() => handleQuickReport(q)}
              className={`text-xs px-3 py-2 rounded-full border transition-all duration-200 cursor-pointer ${
                formData.question === q
                  ? 'bg-[#76B900]/10 text-[#76B900] border-[#76B900]/30'
                  : 'bg-white/3 text-[#9CA3AF] border-white/10 hover:border-[#76B900]/20 hover:text-white'
              }`}
            >
              {q}
            </button>
          ))}
        </div>
      </div>

      {/* Two-column layout */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Left: Generation Form */}
        <form
          onSubmit={handleGenerate}
          className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 space-y-5 shadow-lg shadow-black/20"
        >
          {/* Question */}
          <div>
            <label className="block text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2">
              Report Topic / Question
            </label>
            <textarea
              value={formData.question}
              onChange={(e) => setFormData({ ...formData, question: e.target.value })}
              placeholder="e.g., Summarize lung cancer screening guidelines and AI-assisted detection performance"
              rows={4}
              className="w-full bg-[#0E1117] border border-white/[0.08] rounded-lg px-4 py-3 text-sm text-white placeholder-[#9CA3AF]/60 focus:outline-none focus:border-[#76B900]/50 focus:ring-1 focus:ring-[#76B900]/15 resize-none transition-all duration-200"
            />
          </div>

          {/* Format selector as visual cards */}
          <div>
            <label className="block text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2">
              Output Format
            </label>
            <div className="grid grid-cols-2 md:grid-cols-3 gap-2">
              {formatOptions.map(({ id, label, description, icon: Icon }) => (
                <button
                  key={id}
                  type="button"
                  onClick={() => setFormData({ ...formData, format: id })}
                  className={`text-left p-3 rounded-lg border transition-all duration-200 cursor-pointer ${
                    formData.format === id
                      ? 'bg-[#76B900]/8 border-[#76B900]/30 shadow-[0_0_12px_rgba(118,185,0,0.08)]'
                      : 'bg-[#0E1117] border-white/[0.06] hover:border-white/[0.15]'
                  }`}
                >
                  <Icon size={16} className={formData.format === id ? 'text-[#76B900]' : 'text-[#9CA3AF]'} />
                  <p className={`text-xs font-semibold mt-1.5 ${formData.format === id ? 'text-[#76B900]' : 'text-white'}`}>{label}</p>
                  <p className="text-[10px] text-[#9CA3AF] mt-0.5 line-clamp-1">{description}</p>
                </button>
              ))}
            </div>
          </div>

          {/* Submit */}
          <button
            type="submit"
            disabled={loading || !formData.question.trim()}
            className="flex items-center justify-center gap-2 w-full bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-40 text-black text-sm font-semibold px-5 py-3 rounded-lg transition-all duration-200 cursor-pointer shadow-md shadow-[#76B900]/20"
          >
            {loading ? <Loader2 size={18} className="animate-spin" /> : <FileText size={18} />}
            {loading ? 'Generating...' : 'Generate Report'}
          </button>
        </form>

        {/* Right: Result Panel */}
        <div>
          {error && (
            <div className="bg-red-500/10 border border-red-500/30 rounded-xl p-4 flex items-center gap-3 mb-4">
              <AlertTriangle size={18} className="text-red-400 shrink-0" />
              <p className="text-sm text-red-400">{error}</p>
            </div>
          )}

          {!result && !error && (
            <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-8 flex flex-col items-center justify-center text-center h-full min-h-[300px]">
              <div className="p-4 bg-white/5 rounded-2xl mb-4">
                <FileText size={40} className="text-white/20" />
              </div>
              <p className="text-sm text-[#9CA3AF]">Enter a topic and select a format to generate a clinical report</p>
            </div>
          )}

          {result && (
            <div className="bg-[#1A1D23] rounded-xl border border-[#76B900]/20 overflow-hidden shadow-lg shadow-black/20">
              {/* Result header */}
              <div className="bg-[#76B900]/5 border-b border-[#76B900]/15 px-5 py-3 flex items-center justify-between">
                <div>
                  <h3 className="text-sm font-semibold text-white">
                    {result.title || 'Generated Report'}
                  </h3>
                  <div className="flex items-center gap-3 mt-1">
                    {result.generated_at && (
                      <span className="text-[10px] text-[#9CA3AF] flex items-center gap-1">
                        <Clock size={10} />
                        {new Date(result.generated_at).toLocaleString()}
                      </span>
                    )}
                    {result.evidence_count != null && (
                      <span className="text-[10px] text-[#76B900] flex items-center gap-1 bg-[#76B900]/10 px-2 py-0.5 rounded-full">
                        <BookOpen size={10} />
                        {result.evidence_count} sources
                      </span>
                    )}
                    {result.generation_time_ms != null && (
                      <span className="text-[10px] text-[#9CA3AF]">
                        {result.generation_time_ms}ms
                      </span>
                    )}
                  </div>
                </div>
                <div className="flex items-center gap-2">
                  <button
                    onClick={handleCopy}
                    className="flex items-center gap-1 text-xs text-[#9CA3AF] hover:text-white px-2.5 py-1.5 rounded-lg bg-white/5 border border-white/[0.08] hover:border-white/[0.15] transition-all duration-200 cursor-pointer"
                  >
                    <Copy size={12} />
                    {copied ? 'Copied' : 'Copy'}
                  </button>
                  <button
                    onClick={handleDownload}
                    className="flex items-center gap-1 text-xs text-[#76B900] hover:text-white px-2.5 py-1.5 rounded-lg bg-[#76B900]/10 border border-[#76B900]/20 hover:bg-[#76B900]/20 transition-all duration-200 cursor-pointer"
                  >
                    <Download size={12} />
                    Download
                  </button>
                </div>
              </div>

              {/* Report content */}
              {(result.answer || result.content) && (
                <div className="p-5">
                  <div className={`bg-[#0E1117] rounded-lg p-5 border border-white/[0.06] max-h-[500px] overflow-auto ${
                    isJsonFormat ? 'font-mono text-xs' : 'text-sm'
                  }`}>
                    <pre className={`text-[#E0E0E0] whitespace-pre-wrap ${isJsonFormat ? 'font-mono' : 'font-sans'} leading-relaxed`}>
                      {result.answer || result.content}
                    </pre>
                  </div>
                </div>
              )}

              {/* Evidence sources */}
              {result.evidence && result.evidence.length > 0 && (
                <div className="px-5 pb-5">
                  <div className="border-t border-white/[0.06] pt-4">
                    <h4 className="text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-3 flex items-center gap-2">
                      <BookOpen size={12} className="text-[#76B900]" />
                      Evidence Sources
                    </h4>
                    <div className="grid grid-cols-1 gap-2">
                      {result.evidence.slice(0, 5).map((e, i) => (
                        <div key={i} className="bg-[#0E1117] rounded-lg p-3 border border-white/[0.06] hover:border-white/[0.12] transition-all duration-200">
                          <div className="flex items-center justify-between mb-1">
                            <span className="text-xs text-white font-medium">{e.label || e.collection || 'Source'}</span>
                            {e.score != null && (
                              <div className="flex items-center gap-2">
                                <div className="w-16 h-1.5 bg-white/5 rounded-full overflow-hidden">
                                  <div
                                    className="h-full bg-[#76B900] rounded-full"
                                    style={{ width: `${e.score * 100}%` }}
                                  />
                                </div>
                                <span className="text-[10px] text-[#76B900] font-mono shrink-0">
                                  {(e.score * 100).toFixed(0)}%
                                </span>
                              </div>
                            )}
                          </div>
                          {e.text_snippet && (
                            <p className="text-[11px] text-[#9CA3AF] line-clamp-2 leading-relaxed">{e.text_snippet}</p>
                          )}
                        </div>
                      ))}
                    </div>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
