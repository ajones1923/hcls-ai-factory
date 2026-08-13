import { useState, useEffect } from 'react';
import { FileCheck, Loader2, AlertTriangle, Shield, ChevronRight, Gauge, Zap } from 'lucide-react';
import { recommendProtocol, fetchIndications, describeApiError } from '../lib/api';

interface Alternative {
  modality?: string;
  protocol?: string;
  rating?: number;
  note?: string;
}

interface ProtocolResult {
  indication?: string;
  recommended_protocol?: string;
  recommended_modality?: string;
  acr_appropriateness_rating?: number;
  parameters?: Record<string, unknown>;
  contrast?: string | null;
  dose_estimate_msv?: number | null;
  rationale?: string;
  warnings?: string[];
  alternatives?: Alternative[];
  [key: string]: unknown;
}

function ACRGauge({ rating }: { rating: number }) {
  const segments = Array.from({ length: 9 }, (_, i) => i + 1);
  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">ACR Appropriateness</span>
        <span className="text-lg font-bold text-white">{rating}<span className="text-sm text-[#9CA3AF] font-normal">/9</span></span>
      </div>
      <div className="flex gap-1">
        {segments.map((seg) => {
          let color = 'bg-white/10';
          if (seg <= rating) {
            if (seg <= 3) color = 'bg-red-500';
            else if (seg <= 6) color = 'bg-yellow-500';
            else color = 'bg-green-500';
          }
          return <div key={seg} className={`h-3 flex-1 rounded-sm ${color} transition-all duration-300`} />;
        })}
      </div>
      <div className="flex justify-between text-[9px] text-[#9CA3AF]">
        <span>Usually not appropriate</span>
        <span>May be appropriate</span>
        <span>Usually appropriate</span>
      </div>
    </div>
  );
}

export default function Protocol() {
  const [indications, setIndications] = useState<string[]>([]);
  const [formData, setFormData] = useState({
    indication: '',
    body_region: '',
    age: '',
    sex: '',
    weight: '',
    egfr: '',
    clinical_context: '',
    pregnancy: false,
    contrast_allergy: false,
    pediatric: false,
    allergy_severity: 'mild',
  });
  const [result, setResult] = useState<ProtocolResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetchIndications()
      .then((data) => {
        const list = data.indications || data;
        setIndications(Array.isArray(list) ? list : []);
      })
      .catch((e) => describeApiError('Loading protocol data', e));
  }, []);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    try {
      const res = await recommendProtocol({
        ...formData,
        age: formData.age ? parseInt(formData.age) : undefined,
        weight: formData.weight ? parseFloat(formData.weight) : undefined,
        egfr: formData.egfr ? parseFloat(formData.egfr) : undefined,
      });
      setResult(res);
    } catch {
      setError('Failed to get protocol recommendation');
    } finally {
      setLoading(false);
    }
  };

  const handleIndicationClick = (ind: string) => {
    setFormData({ ...formData, indication: ind });
    setResult(null);
    setError(null);
  };

  const inputClass = "w-full bg-[#0E1117] border border-white/[0.08] rounded-lg px-3 py-2.5 text-sm text-white placeholder-[#9CA3AF]/60 focus:outline-none focus:border-[#76B900]/50 focus:ring-1 focus:ring-[#76B900]/15 transition-all duration-200";

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3">
        <div className="p-2.5 bg-[#76B900]/10 rounded-xl border border-[#76B900]/20">
          <FileCheck size={24} className="text-[#76B900]" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-white tracking-tight">Protocol Recommendation</h1>
          <p className="text-sm text-[#9CA3AF]">AI-assisted imaging protocol selection based on ACR Appropriateness Criteria</p>
        </div>
      </div>

      {/* Two-column layout */}
      <div className="grid grid-cols-1 lg:grid-cols-5 gap-6">
        {/* Left column: Form (3/5 width) */}
        <div className="lg:col-span-3">
          <form
            onSubmit={handleSubmit}
            className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 space-y-5 shadow-lg shadow-black/20"
          >
            {/* Clinical Indication */}
            <div>
              <label className="block text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2">
                Clinical Indication
              </label>
              {indications.length > 0 ? (
                <select
                  value={formData.indication}
                  onChange={(e) => setFormData({ ...formData, indication: e.target.value })}
                  className={inputClass}
                >
                  <option value="">Select indication...</option>
                  {indications.map((ind) => (
                    <option key={ind} value={ind}>{ind}</option>
                  ))}
                </select>
              ) : (
                <input
                  type="text"
                  value={formData.indication}
                  onChange={(e) => setFormData({ ...formData, indication: e.target.value })}
                  placeholder="e.g., Lung cancer screening"
                  className={inputClass}
                />
              )}
            </div>

            {/* Patient Factors - compact grid */}
            <div>
              <label className="block text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2">
                Patient Factors
              </label>
              <div className="grid grid-cols-2 md:grid-cols-3 gap-3">
                <div>
                  <label className="block text-[11px] text-[#9CA3AF] mb-1">Age</label>
                  <input
                    type="number"
                    value={formData.age}
                    onChange={(e) => setFormData({ ...formData, age: e.target.value })}
                    placeholder="55"
                    className={inputClass}
                  />
                </div>
                <div>
                  <label className="block text-[11px] text-[#9CA3AF] mb-1">Sex</label>
                  <select
                    value={formData.sex}
                    onChange={(e) => setFormData({ ...formData, sex: e.target.value })}
                    className={inputClass}
                  >
                    <option value="">Select...</option>
                    <option value="M">Male</option>
                    <option value="F">Female</option>
                  </select>
                </div>
                <div>
                  <label className="block text-[11px] text-[#9CA3AF] mb-1">Weight (kg)</label>
                  <input
                    type="number"
                    value={formData.weight}
                    onChange={(e) => setFormData({ ...formData, weight: e.target.value })}
                    placeholder="70"
                    className={inputClass}
                  />
                </div>
                <div>
                  <label className="block text-[11px] text-[#9CA3AF] mb-1">eGFR</label>
                  <input
                    type="number"
                    value={formData.egfr}
                    onChange={(e) => setFormData({ ...formData, egfr: e.target.value })}
                    placeholder="90"
                    className={inputClass}
                  />
                </div>
                <div>
                  <label className="block text-[11px] text-[#9CA3AF] mb-1">Body Region</label>
                  <input
                    type="text"
                    value={formData.body_region}
                    onChange={(e) => setFormData({ ...formData, body_region: e.target.value })}
                    placeholder="Chest"
                    className={inputClass}
                  />
                </div>
              </div>
            </div>

            {/* Checkboxes */}
            <div className="flex flex-wrap gap-4">
              {[
                { key: 'pregnancy', label: 'Pregnancy' },
                { key: 'contrast_allergy', label: 'Contrast Allergy' },
                { key: 'pediatric', label: 'Pediatric' },
              ].map(({ key, label }) => (
                <label key={key} className="flex items-center gap-2 cursor-pointer group">
                  <div className={`w-4 h-4 rounded border transition-all duration-200 flex items-center justify-center ${
                    formData[key as keyof typeof formData]
                      ? 'bg-[#76B900] border-[#76B900]'
                      : 'border-white/20 group-hover:border-white/40'
                  }`}>
                    {formData[key as keyof typeof formData] && (
                      <svg width="10" height="8" viewBox="0 0 10 8" fill="none">
                        <path d="M1 4L3.5 6.5L9 1" stroke="black" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
                      </svg>
                    )}
                  </div>
                  <span className="text-sm text-[#E0E0E0]">{label}</span>
                </label>
              ))}
            </div>

            {/* Contrast allergy severity */}
            {formData.contrast_allergy && (
              <div className="bg-amber-500/5 border border-amber-500/20 rounded-lg p-3">
                <label className="block text-[11px] text-amber-400 mb-1.5 font-medium">Allergy Severity</label>
                <select
                  value={formData.allergy_severity}
                  onChange={(e) => setFormData({ ...formData, allergy_severity: e.target.value })}
                  className={inputClass}
                >
                  <option value="mild">Mild</option>
                  <option value="moderate">Moderate</option>
                  <option value="severe">Severe</option>
                </select>
              </div>
            )}

            {/* Clinical Context */}
            <div>
              <label className="block text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2">
                Clinical Context
              </label>
              <textarea
                value={formData.clinical_context}
                onChange={(e) => setFormData({ ...formData, clinical_context: e.target.value })}
                placeholder="Additional clinical details, history, prior imaging..."
                rows={2}
                className={`${inputClass} resize-none`}
              />
            </div>

            {/* Submit */}
            <button
              type="submit"
              disabled={loading || !formData.indication}
              className="flex items-center justify-center gap-2 w-full bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-40 text-black text-sm font-semibold px-5 py-3 rounded-lg transition-all duration-200 cursor-pointer shadow-md shadow-[#76B900]/20"
            >
              {loading ? (
                <Loader2 size={18} className="animate-spin" />
              ) : (
                <FileCheck size={18} />
              )}
              {loading ? 'Analyzing...' : 'Get Recommendation'}
            </button>
          </form>
        </div>

        {/* Right column: Results (2/5 width) */}
        <div className="lg:col-span-2">
          {error && (
            <div className="bg-red-500/10 border border-red-500/30 rounded-xl p-4 flex items-center gap-3 mb-4">
              <AlertTriangle size={18} className="text-red-400 shrink-0" />
              <p className="text-sm text-red-400">{error}</p>
            </div>
          )}

          {!result && !error && (
            <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-8 flex flex-col items-center justify-center text-center h-full min-h-[300px]">
              <div className="p-4 bg-white/5 rounded-2xl mb-4">
                <Gauge size={40} className="text-white/20" />
              </div>
              <p className="text-sm text-[#9CA3AF]">Select an indication and submit to receive an AI-assisted protocol recommendation</p>
            </div>
          )}

          {result && (
            <div className="bg-[#1A1D23] rounded-xl border border-[#76B900]/20 overflow-hidden shadow-lg shadow-black/20 space-y-0">
              {/* Protocol Name */}
              {result.recommended_protocol && (
                <div className="bg-[#76B900]/5 border-b border-[#76B900]/15 px-5 py-4">
                  <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Recommended Protocol</span>
                  <h3 className="text-xl font-bold text-[#76B900] mt-1">
                    {result.recommended_protocol}
                  </h3>
                  {result.recommended_modality && (
                    <span className="inline-block mt-2 text-[11px] font-medium bg-blue-500/10 text-blue-400 border border-blue-500/25 px-2.5 py-1 rounded-full">
                      {result.recommended_modality}
                    </span>
                  )}
                </div>
              )}

              <div className="p-5 space-y-5">
                {/* ACR Rating Gauge */}
                {result.acr_appropriateness_rating != null && (
                  <ACRGauge rating={result.acr_appropriateness_rating} />
                )}

                {/* Dose Estimate */}
                {result.dose_estimate_msv != null && (
                  <div className="bg-[#0E1117] rounded-lg p-4 border border-white/[0.06]">
                    <div className="flex items-center justify-between">
                      <div className="flex items-center gap-2">
                        <Zap size={14} className="text-yellow-400" />
                        <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Estimated Dose</span>
                      </div>
                      <span className="text-lg font-bold text-white">{result.dose_estimate_msv} <span className="text-sm text-[#9CA3AF] font-normal">mSv</span></span>
                    </div>
                    <div className="mt-2 h-2 bg-white/5 rounded-full overflow-hidden">
                      <div
                        className={`h-full rounded-full transition-all duration-500 ${
                          result.dose_estimate_msv <= 3 ? 'bg-green-500' :
                          result.dose_estimate_msv <= 10 ? 'bg-yellow-500' :
                          result.dose_estimate_msv <= 20 ? 'bg-orange-500' : 'bg-red-500'
                        }`}
                        style={{ width: `${Math.min(100, (result.dose_estimate_msv / 30) * 100)}%` }}
                      />
                    </div>
                  </div>
                )}

                {/* Warnings */}
                {result.warnings && result.warnings.length > 0 && (
                  <div className="space-y-2">
                    {result.warnings.map((w, i) => (
                      <div key={i} className="bg-amber-500/8 border border-amber-500/25 rounded-lg p-3 flex items-start gap-2.5">
                        <AlertTriangle size={14} className="text-amber-400 shrink-0 mt-0.5" />
                        <p className="text-xs text-amber-300 leading-relaxed">{w}</p>
                      </div>
                    ))}
                  </div>
                )}

                {/* Rationale */}
                {result.rationale && (
                  <div className="border-l-2 border-[#76B900]/30 pl-4">
                    <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Rationale</span>
                    <p className="text-sm text-[#E0E0E0] mt-1 leading-relaxed">{result.rationale}</p>
                  </div>
                )}

                {/* Parameters */}
                {result.parameters && Object.keys(result.parameters).length > 0 && (
                  <div>
                    <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Parameters</span>
                    <div className="grid grid-cols-2 gap-2 mt-2">
                      {Object.entries(result.parameters).map(([k, v]) => (
                        <div key={k} className="bg-[#0E1117] rounded-lg p-2.5 border border-white/[0.06]">
                          <span className="text-[10px] text-[#9CA3AF] block">{k}</span>
                          <span className="text-sm text-white font-medium">
                            {Array.isArray(v) ? v.join(', ') : String(v)}
                          </span>
                        </div>
                      ))}
                    </div>
                  </div>
                )}

                {/* Alternatives */}
                {result.alternatives && result.alternatives.length > 0 && (
                  <div>
                    <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Alternatives</span>
                    <div className="space-y-2 mt-2">
                      {result.alternatives.map((alt, i) => (
                        <div key={i} className="bg-[#0E1117] rounded-lg p-3 border border-white/[0.06] flex items-center justify-between hover:border-white/[0.12] transition-all duration-200">
                          <div className="flex items-center gap-3 min-w-0">
                            <span className="text-xs font-bold text-[#9CA3AF] shrink-0">#{i + 1}</span>
                            <div className="min-w-0">
                              <span className="text-sm text-white font-medium block truncate">{alt.protocol}</span>
                              {alt.modality && (
                                <span className="text-[10px] text-[#9CA3AF]">{alt.modality}</span>
                              )}
                              {alt.note && (
                                <p className="text-[10px] text-[#9CA3AF] mt-0.5 line-clamp-1">{alt.note}</p>
                              )}
                            </div>
                          </div>
                          {alt.rating != null && (
                            <div className="flex items-center gap-1 shrink-0 ml-2">
                              <div className={`w-2 h-2 rounded-full ${
                                alt.rating >= 7 ? 'bg-green-500' : alt.rating >= 4 ? 'bg-yellow-500' : 'bg-red-500'
                              }`} />
                              <span className="text-xs font-mono text-white">{alt.rating}/9</span>
                            </div>
                          )}
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Supported Indications */}
      {indications.length > 0 && (
        <div>
          <div className="flex items-center gap-2 mb-3">
            <Shield size={14} className="text-[#9CA3AF]" />
            <span className="text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold">
              Supported Indications
            </span>
            <span className="text-[10px] bg-white/5 text-[#9CA3AF] border border-white/10 px-2 py-0.5 rounded-full">
              {indications.length}
            </span>
          </div>
          <div className="flex flex-wrap gap-2">
            {indications.map((ind) => (
              <button
                key={ind}
                onClick={() => handleIndicationClick(ind)}
                className={`text-xs px-3 py-1.5 rounded-full border transition-all duration-200 cursor-pointer ${
                  formData.indication === ind
                    ? 'bg-[#76B900]/10 text-[#76B900] border-[#76B900]/30'
                    : 'bg-white/3 text-[#9CA3AF] border-white/10 hover:border-[#76B900]/20 hover:text-white'
                }`}
              >
                {ind}
                <ChevronRight size={10} className="inline ml-1 opacity-50" />
              </button>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
