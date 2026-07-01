import { useCallback, useEffect, useState } from 'react';
import {
  Upload,
  Zap,
  Loader2,
  AlertTriangle,
  CheckCircle,
  Cpu,
  FileImage,
  Play,
  Info,
} from 'lucide-react';
import axios from 'axios';

const API_BASE =
  import.meta.env.VITE_API_URL ||
  `${window.location.protocol}//${window.location.hostname}:8524`;
const api = axios.create({ baseURL: API_BASE });

interface Finding {
  category?: string;
  description?: string;
  severity?: string;
  confidence?: number;
}

interface AnalysisResult {
  workflow_name?: string;
  status?: string;
  classification?: string;
  severity?: string;
  findings?: Finding[];
  measurements?: Record<string, number>;
  inference_time_ms?: number;
  is_mock?: boolean;
  nim_services_used?: string[];
  dicom_metadata?: Record<string, string | number>;
  all_pathology_scores?: Record<string, number>;
  sample_name?: string;
  error?: string;
  analysis_note?: string;
}

interface AnalysisStatus {
  gpu_available?: boolean;
  gpu_name?: string;
  gpu_memory_mb?: number;
  device?: string;
  models_loaded?: string[];
  supported_modalities?: string[];
  live_analysis_enabled?: boolean;
  available_samples?: Record<
    string,
    { path: string; exists: boolean; is_series: boolean }
  >;
}

const severityStyles: Record<string, string> = {
  critical: 'bg-red-500/10 text-red-400 border border-red-500/30',
  urgent: 'bg-orange-500/10 text-orange-400 border border-orange-500/30',
  significant: 'bg-yellow-500/10 text-yellow-400 border border-yellow-500/30',
  routine: 'bg-green-500/10 text-green-400 border border-green-500/30',
  normal: 'bg-green-500/10 text-green-300 border border-green-500/30',
};

const severityDot: Record<string, string> = {
  critical: 'bg-red-500',
  urgent: 'bg-orange-500',
  significant: 'bg-yellow-500',
  routine: 'bg-green-500',
  normal: 'bg-green-400',
};

export default function LiveAnalysis() {
  const [status, setStatus] = useState<AnalysisStatus | null>(null);
  const [result, setResult] = useState<AnalysisResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [dragOver, setDragOver] = useState(false);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [dicomMeta, setDicomMeta] = useState<Record<string, string | number> | null>(null);

  // Fetch analyzer status on mount
  useEffect(() => {
    api
      .get('/analyze/status')
      .then((r) => setStatus(r.data))
      .catch(() => {});
  }, []);

  // ── File Upload Handlers ──

  const handleFiles = useCallback((files: FileList | null) => {
    if (!files || files.length === 0) return;
    const file = files[0];
    setSelectedFile(file);
    setResult(null);
    setError(null);
    setDicomMeta(null);
  }, []);

  const handleDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      setDragOver(false);
      handleFiles(e.dataTransfer.files);
    },
    [handleFiles]
  );

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setDragOver(true);
  }, []);

  const handleDragLeave = useCallback(() => {
    setDragOver(false);
  }, []);

  // ── Analysis ──

  const runAnalysis = async (file?: File, sampleName?: string) => {
    setLoading(true);
    setError(null);
    setResult(null);

    try {
      let res;
      if (sampleName) {
        res = await api.post(`/analyze/analyze-sample/${sampleName}`);
      } else if (file) {
        const formData = new FormData();
        formData.append('file', file);
        // Don't set Content-Type — let axios set it with the correct boundary
        res = await api.post('/analyze/upload', formData);
      } else {
        throw new Error('No file or sample selected');
      }

      setResult(res.data);
      if (res.data.dicom_metadata) {
        setDicomMeta(res.data.dicom_metadata);
      }
    } catch (err: unknown) {
      const msg =
        err instanceof Error
          ? err.message
          : axios.isAxiosError(err) && err.response?.data?.detail
            ? err.response.data.detail
            : 'Analysis failed';
      setError(msg);
    } finally {
      setLoading(false);
    }
  };

  // ── Render Helpers ──

  const getSeverityStyle = (severity?: string) => {
    if (!severity) return 'bg-blue-500/10 text-blue-400 border border-blue-500/30';
    return severityStyles[severity.toLowerCase()] || 'bg-blue-500/10 text-blue-400 border border-blue-500/30';
  };

  const getSeverityDot = (severity?: string) => {
    if (!severity) return 'bg-blue-500';
    return severityDot[severity.toLowerCase()] || 'bg-blue-500';
  };

  const availableSamples = status?.available_samples
    ? Object.entries(status.available_samples).filter(([, v]) => v.exists)
    : [];

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold text-white flex items-center gap-3">
            <Zap className="text-[#76B900]" size={28} />
            Live DICOM Analysis
          </h1>
          <p className="text-sm text-[#9CA3AF] mt-1">
            Upload real DICOM images for AI-powered analysis with DenseNet-121 inference
          </p>
        </div>
        <div className="flex items-center gap-2">
          <span
            className={`w-2.5 h-2.5 rounded-full ${
              status?.gpu_available ? 'bg-[#76B900]' : 'bg-yellow-500'
            }`}
          />
          <span className="text-xs text-[#9CA3AF]">
            {status?.gpu_name || (status?.gpu_available ? 'GPU Ready' : 'CPU Mode')}
          </span>
        </div>
      </div>

      {/* GPU Status Bar */}
      {status && (
        <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-4">
          <div className="flex items-center gap-6 text-xs">
            <div className="flex items-center gap-2">
              <Cpu size={14} className="text-[#76B900]" />
              <span className="text-[#9CA3AF]">Device:</span>
              <span className="text-white font-medium">{status.device || 'cpu'}</span>
            </div>
            {status.gpu_memory_mb && (
              <div className="flex items-center gap-2">
                <span className="text-[#9CA3AF]">VRAM:</span>
                <span className="text-white font-medium">
                  {(status.gpu_memory_mb / 1024).toFixed(0)} GB
                </span>
              </div>
            )}
            <div className="flex items-center gap-2">
              <span className="text-[#9CA3AF]">Models loaded:</span>
              <span className="text-white font-medium">
                {status.models_loaded?.length || 0}
              </span>
            </div>
            <div className="flex items-center gap-2">
              <span className="text-[#9CA3AF]">Modalities:</span>
              <span className="text-white font-medium">
                {status.supported_modalities?.join(', ') || 'None'}
              </span>
            </div>
          </div>
        </div>
      )}

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Left Column: Upload + Samples */}
        <div className="space-y-6">
          {/* Upload Drop Zone */}
          <div
            onDrop={handleDrop}
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            className={`bg-[#1A1D23] border-2 border-dashed rounded-lg p-8 text-center transition-all cursor-pointer ${
              dragOver
                ? 'border-[#76B900] bg-[#76B900]/[0.05]'
                : 'border-white/[0.15] hover:border-white/[0.3]'
            }`}
            onClick={() => document.getElementById('dicom-file-input')?.click()}
          >
            <input
              id="dicom-file-input"
              type="file"
              accept=".dcm,.dicom,.DCM,application/dicom,*/*"
              className="hidden"
              onChange={(e) => handleFiles(e.target.files)}
            />
            <Upload
              size={40}
              className={`mx-auto mb-3 ${
                dragOver ? 'text-[#76B900]' : 'text-[#9CA3AF]'
              }`}
            />
            <p className="text-white font-medium mb-1">
              {selectedFile
                ? selectedFile.name
                : 'Drop DICOM file here or click to browse'}
            </p>
            <p className="text-xs text-[#9CA3AF]">
              {selectedFile
                ? `${(selectedFile.size / 1024).toFixed(0)} KB`
                : 'Supports .dcm files -- CXR, CT, MRI'}
            </p>

            {selectedFile && (
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  runAnalysis(selectedFile);
                }}
                disabled={loading}
                className="mt-4 px-6 py-2.5 bg-[#76B900] text-black font-semibold rounded-lg hover:bg-[#8DD100] transition-colors disabled:opacity-50 flex items-center gap-2 mx-auto"
              >
                {loading ? (
                  <Loader2 size={16} className="animate-spin" />
                ) : (
                  <Zap size={16} />
                )}
                {loading ? 'Analyzing...' : 'Run Live Analysis'}
              </button>
            )}
          </div>

          {/* Sample DICOMs */}
          {availableSamples.length > 0 && (
            <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-5">
              <h3 className="text-sm font-semibold text-white mb-3 flex items-center gap-2">
                <FileImage size={16} className="text-[#76B900]" />
                Pre-loaded Samples
              </h3>
              <div className="space-y-2">
                {availableSamples.map(([name, info]) => (
                  <button
                    key={name}
                    onClick={() => runAnalysis(undefined, name)}
                    disabled={loading}
                    className="w-full flex items-center justify-between px-4 py-3 bg-white/[0.03] hover:bg-white/[0.06] border border-white/[0.08] rounded-lg transition-colors disabled:opacity-50 cursor-pointer"
                  >
                    <div className="flex items-center gap-3">
                      <Play size={14} className="text-[#76B900]" />
                      <div className="text-left">
                        <p className="text-sm text-white font-medium">
                          {name.replace(/_/g, ' ')}
                        </p>
                        <p className="text-xs text-[#9CA3AF]">
                          {info.is_series ? 'CT Series' : 'Single DICOM'} --{' '}
                          {info.path}
                        </p>
                      </div>
                    </div>
                    {loading ? (
                      <Loader2 size={14} className="animate-spin text-[#9CA3AF]" />
                    ) : (
                      <span className="text-xs text-[#76B900] font-medium">Analyze</span>
                    )}
                  </button>
                ))}
              </div>
            </div>
          )}

          {/* DICOM Metadata */}
          {dicomMeta && (
            <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-5">
              <h3 className="text-sm font-semibold text-white mb-3 flex items-center gap-2">
                <Info size={16} className="text-blue-400" />
                DICOM Metadata
              </h3>
              <div className="grid grid-cols-2 gap-2 text-xs">
                {Object.entries(dicomMeta)
                  .filter(([, v]) => v !== '' && v !== 0)
                  .map(([key, val]) => (
                    <div key={key} className="flex justify-between gap-2">
                      <span className="text-[#9CA3AF] truncate">
                        {key.replace(/_/g, ' ')}
                      </span>
                      <span className="text-white font-medium truncate max-w-[140px]">
                        {String(val)}
                      </span>
                    </div>
                  ))}
              </div>
            </div>
          )}
        </div>

        {/* Right Column: Results */}
        <div className="space-y-6">
          {/* Error */}
          {error && (
            <div className="bg-red-500/10 border border-red-500/30 rounded-lg p-4 flex items-start gap-3">
              <AlertTriangle size={18} className="text-red-400 mt-0.5 shrink-0" />
              <div>
                <p className="text-sm font-medium text-red-400">Analysis Failed</p>
                <p className="text-xs text-red-400/80 mt-1">{error}</p>
              </div>
            </div>
          )}

          {/* Loading */}
          {loading && !result && (
            <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-8 text-center">
              <Loader2 size={32} className="animate-spin text-[#76B900] mx-auto mb-3" />
              <p className="text-white font-medium">Running AI Inference...</p>
              <p className="text-xs text-[#9CA3AF] mt-1">
                DenseNet-121 model loading and running on GPU
              </p>
            </div>
          )}

          {/* Results */}
          {result && (
            <div className="space-y-4">
              {/* Header Card */}
              <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-5">
                <div className="flex items-center justify-between mb-3">
                  <div className="flex items-center gap-3">
                    <CheckCircle size={20} className="text-[#76B900]" />
                    <h3 className="text-base font-semibold text-white">Analysis Complete</h3>
                  </div>
                  <div className="flex items-center gap-2">
                    {result.is_mock === false && (
                      <span className="px-2.5 py-0.5 text-[10px] font-bold uppercase tracking-wider bg-[#76B900]/15 text-[#76B900] border border-[#76B900]/40 rounded-full">
                        LIVE INFERENCE
                      </span>
                    )}
                    {result.is_mock === true && (
                      <span className="px-2.5 py-0.5 text-[10px] font-bold uppercase tracking-wider bg-yellow-500/15 text-yellow-400 border border-yellow-500/40 rounded-full">
                        SIMULATED
                      </span>
                    )}
                  </div>
                </div>

                <div className="grid grid-cols-3 gap-4 text-xs">
                  <div>
                    <span className="text-[#9CA3AF]">Workflow</span>
                    <p className="text-white font-medium mt-0.5">
                      {result.workflow_name?.replace(/_/g, ' ') || 'Unknown'}
                    </p>
                  </div>
                  <div>
                    <span className="text-[#9CA3AF]">Inference Time</span>
                    <p className="text-white font-medium mt-0.5">
                      {result.inference_time_ms
                        ? `${result.inference_time_ms.toFixed(0)} ms`
                        : '--'}
                    </p>
                  </div>
                  <div>
                    <span className="text-[#9CA3AF]">Models Used</span>
                    <p className="text-white font-medium mt-0.5 truncate" title={result.nim_services_used?.join(', ')}>
                      {result.nim_services_used?.[0] || '--'}
                    </p>
                  </div>
                </div>
              </div>

              {/* Classification + Severity */}
              <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-5">
                <div className="flex items-center gap-4 mb-4">
                  <div className="flex items-center gap-2">
                    <span
                      className={`w-3 h-3 rounded-full ${getSeverityDot(
                        result.severity
                      )}`}
                    />
                    <span
                      className={`px-3 py-1 rounded-full text-xs font-semibold ${getSeverityStyle(
                        result.severity
                      )}`}
                    >
                      {result.severity?.toUpperCase() || 'PENDING'}
                    </span>
                  </div>
                  <span className="text-sm text-white font-medium">
                    {result.classification || 'No classification'}
                  </span>
                </div>

                {/* Findings */}
                <h4 className="text-xs font-semibold text-[#9CA3AF] uppercase tracking-wider mb-2">
                  Findings
                </h4>
                <div className="space-y-2">
                  {result.findings?.map((finding, i) => (
                    <div
                      key={i}
                      className="flex items-start gap-3 px-3 py-2.5 bg-white/[0.02] rounded-lg border border-white/[0.05]"
                    >
                      <span
                        className={`w-2 h-2 rounded-full mt-1.5 shrink-0 ${getSeverityDot(
                          finding.severity
                        )}`}
                      />
                      <div className="min-w-0">
                        <div className="flex items-center gap-2 mb-0.5">
                          <span className="text-xs font-semibold text-white">
                            {finding.category?.replace(/_/g, ' ') || 'Finding'}
                          </span>
                          {finding.confidence !== undefined && (
                            <span className="text-[10px] text-[#9CA3AF]">
                              {(finding.confidence * 100).toFixed(1)}%
                            </span>
                          )}
                          <span
                            className={`px-1.5 py-0.5 rounded text-[9px] font-medium ${getSeverityStyle(
                              finding.severity
                            )}`}
                          >
                            {finding.severity}
                          </span>
                        </div>
                        <p className="text-xs text-[#9CA3AF] leading-relaxed">
                          {finding.description}
                        </p>
                      </div>
                    </div>
                  ))}
                </div>
              </div>

              {/* Pathology Scores (for CXR) */}
              {result.all_pathology_scores && (
                <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-5">
                  <h4 className="text-xs font-semibold text-[#9CA3AF] uppercase tracking-wider mb-3">
                    All Pathology Scores (18 labels)
                  </h4>
                  <div className="space-y-1.5">
                    {Object.entries(result.all_pathology_scores)
                      .sort(([, a], [, b]) => b - a)
                      .map(([label, score]) => (
                        <div key={label} className="flex items-center gap-3">
                          <span className="text-xs text-[#9CA3AF] w-44 truncate">
                            {label}
                          </span>
                          <div className="flex-1 h-2 bg-white/[0.05] rounded-full overflow-hidden">
                            <div
                              className={`h-full rounded-full transition-all ${
                                score > 0.5
                                  ? 'bg-red-500'
                                  : score > 0.3
                                    ? 'bg-yellow-500'
                                    : 'bg-[#76B900]/60'
                              }`}
                              style={{ width: `${Math.max(score * 100, 1)}%` }}
                            />
                          </div>
                          <span
                            className={`text-xs font-mono w-14 text-right ${
                              score > 0.5 ? 'text-red-400 font-bold' : 'text-[#9CA3AF]'
                            }`}
                          >
                            {(score * 100).toFixed(1)}%
                          </span>
                        </div>
                      ))}
                  </div>
                </div>
              )}

              {/* Measurements */}
              {result.measurements &&
                Object.keys(result.measurements).length > 0 &&
                !result.all_pathology_scores && (
                  <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-5">
                    <h4 className="text-xs font-semibold text-[#9CA3AF] uppercase tracking-wider mb-3">
                      Measurements
                    </h4>
                    <div className="grid grid-cols-2 gap-2 text-xs">
                      {Object.entries(result.measurements).map(([key, val]) => (
                        <div key={key} className="flex justify-between gap-2">
                          <span className="text-[#9CA3AF] truncate">
                            {key.replace(/_/g, ' ')}
                          </span>
                          <span className="text-white font-medium">
                            {typeof val === 'number' ? val.toFixed(4) : String(val)}
                          </span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}

              {/* Analysis Note */}
              {result.analysis_note && (
                <div className="bg-blue-500/5 border border-blue-500/20 rounded-lg p-3 flex items-start gap-2">
                  <Info size={14} className="text-blue-400 mt-0.5 shrink-0" />
                  <p className="text-xs text-blue-400/80">{result.analysis_note}</p>
                </div>
              )}
            </div>
          )}

          {/* Empty state */}
          {!result && !loading && !error && (
            <div className="bg-[#1A1D23] border border-white/[0.08] rounded-lg p-12 text-center">
              <Zap size={40} className="text-[#9CA3AF]/30 mx-auto mb-3" />
              <p className="text-[#9CA3AF] text-sm">
                Upload a DICOM file or select a sample to run live AI analysis
              </p>
              <p className="text-[#9CA3AF]/60 text-xs mt-2">
                CXR analysis uses DenseNet-121 with real torchxrayvision inference
              </p>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
