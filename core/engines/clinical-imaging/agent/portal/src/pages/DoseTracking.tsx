import { useState } from 'react';
import { Radiation, Loader2, AlertTriangle, Search, Activity, CheckCircle, BarChart3, Scale } from 'lucide-react';
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, Cell } from 'recharts';
import { recordDose, fetchCumulativeDose, compareDRL } from '../lib/api';

interface CumulativeResult {
  patient_id?: string;
  total_effective_dose_msv?: number;
  study_count?: number;
  date_range?: Record<string, string>;
  by_modality?: Record<string, number>;
  by_body_region?: Record<string, number>;
  alert_level?: string;
  alert_message?: string | null;
  [key: string]: unknown;
}

interface DRLResult {
  protocol?: string;
  patient_dose_msv?: number;
  drl_msv?: number;
  achievable_dose_msv?: number;
  ratio?: number;
  status?: string;
  optimization_suggestions?: string[];
}

const tabs = [
  { id: 'lookup', label: 'Patient Lookup', icon: Search },
  { id: 'record', label: 'Record Dose', icon: Radiation },
  { id: 'drl', label: 'DRL Comparison', icon: Scale },
] as const;

type TabId = typeof tabs[number]['id'];

const alertStyles: Record<string, { bg: string; text: string; border: string; label: string }> = {
  high: { bg: 'bg-red-500/10', text: 'text-red-400', border: 'border-red-500/30', label: 'HIGH RISK' },
  warning: { bg: 'bg-red-500/10', text: 'text-red-400', border: 'border-red-500/30', label: 'WARNING' },
  moderate: { bg: 'bg-orange-500/10', text: 'text-orange-400', border: 'border-orange-500/30', label: 'MODERATE' },
  caution: { bg: 'bg-yellow-500/10', text: 'text-yellow-400', border: 'border-yellow-500/30', label: 'CAUTION' },
  normal: { bg: 'bg-green-500/10', text: 'text-green-400', border: 'border-green-500/30', label: 'NORMAL' },
  low: { bg: 'bg-green-500/10', text: 'text-green-400', border: 'border-green-500/30', label: 'LOW RISK' },
};

function getAlertStyle(level?: string) {
  if (!level) return alertStyles.normal;
  return alertStyles[level.toLowerCase()] || alertStyles.normal;
}

function DoseGauge({ dose }: { dose: number }) {
  // Thresholds: <20=green, 20-50=yellow, 50-100=orange, >100=red
  const maxDisplay = 120;
  const pct = Math.min(100, (dose / maxDisplay) * 100);
  const color = dose < 20 ? '#22c55e' : dose < 50 ? '#eab308' : dose < 100 ? '#f97316' : '#ef4444';

  return (
    <div className="space-y-2">
      <div className="h-4 bg-white/5 rounded-full overflow-hidden relative">
        <div
          className="h-full rounded-full transition-all duration-700 ease-out"
          style={{ width: `${pct}%`, backgroundColor: color }}
        />
        {/* Threshold markers */}
        {[20, 50, 100].map(t => (
          <div
            key={t}
            className="absolute top-0 bottom-0 w-px bg-white/20"
            style={{ left: `${(t / maxDisplay) * 100}%` }}
          />
        ))}
      </div>
      <div className="flex justify-between text-[9px] text-[#9CA3AF]">
        <span>0</span>
        <span>20 mSv</span>
        <span>50 mSv</span>
        <span>100 mSv</span>
        <span>120+</span>
      </div>
    </div>
  );
}

const chartColors = ['#76B900', '#22d3ee', '#a855f7', '#f97316', '#ec4899', '#6366f1', '#14b8a6'];

export default function DoseTracking() {
  const [activeTab, setActiveTab] = useState<TabId>('lookup');
  const [patientId, setPatientId] = useState('');
  const [cumulative, setCumulative] = useState<CumulativeResult | null>(null);
  const [loading, setLoading] = useState(false);
  const [recordLoading, setRecordLoading] = useState(false);
  const [drlLoading, setDrlLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [recordForm, setRecordForm] = useState({
    patient_id: '',
    modality: 'CT',
    body_region: 'Chest',
    dose_mgy: '',
    dlp_mgycm: '',
  });
  const [recordSuccess, setRecordSuccess] = useState(false);
  const [drlForm, setDrlForm] = useState({
    protocol: 'CT Chest',
    effective_dose_msv: '',
    body_region: 'chest',
    modality: 'CT',
    pediatric: false,
  });
  const [drlResult, setDrlResult] = useState<DRLResult | null>(null);

  const handleLookup = async () => {
    if (!patientId.trim()) return;
    setLoading(true);
    setError(null);
    try {
      const data = await fetchCumulativeDose(patientId);
      setCumulative(data);
    } catch {
      setError('Failed to fetch cumulative dose');
    } finally {
      setLoading(false);
    }
  };

  const handleRecord = async (e: React.FormEvent) => {
    e.preventDefault();
    setRecordLoading(true);
    setRecordSuccess(false);
    setError(null);
    try {
      await recordDose({
        ...recordForm,
        dose_mgy: recordForm.dose_mgy ? parseFloat(recordForm.dose_mgy) : undefined,
        dlp_mgycm: recordForm.dlp_mgycm ? parseFloat(recordForm.dlp_mgycm) : undefined,
      });
      setRecordSuccess(true);
    } catch {
      setError('Failed to record dose');
    } finally {
      setRecordLoading(false);
    }
  };

  const handleDRLCompare = async (e: React.FormEvent) => {
    e.preventDefault();
    setDrlLoading(true);
    setError(null);
    setDrlResult(null);
    try {
      const data = await compareDRL({
        ...drlForm,
        effective_dose_msv: parseFloat(drlForm.effective_dose_msv),
      });
      setDrlResult(data);
    } catch {
      setError('Failed to compare to DRL');
    } finally {
      setDrlLoading(false);
    }
  };

  const inputClass = "w-full bg-[#0E1117] border border-white/[0.08] rounded-lg px-3 py-2.5 text-sm text-white placeholder-[#9CA3AF]/60 focus:outline-none focus:border-[#76B900]/50 focus:ring-1 focus:ring-[#76B900]/15 transition-all duration-200";

  const modalityData = cumulative?.by_modality
    ? Object.entries(cumulative.by_modality).map(([name, value]) => ({ name, value: Number(value) }))
    : [];
  const regionData = cumulative?.by_body_region
    ? Object.entries(cumulative.by_body_region).map(([name, value]) => ({ name, value: Number(value) }))
    : [];

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3">
        <div className="p-2.5 bg-[#76B900]/10 rounded-xl border border-[#76B900]/20">
          <Activity size={24} className="text-[#76B900]" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-white tracking-tight">Dose Tracking</h1>
          <p className="text-sm text-[#9CA3AF]">Patient radiation dose monitoring, recording, and DRL benchmarking</p>
        </div>
      </div>

      {/* Tab Bar */}
      <div className="flex gap-1 bg-[#1A1D23] p-1 rounded-xl border border-white/[0.08] w-fit">
        {tabs.map(({ id, label, icon: Icon }) => (
          <button
            key={id}
            onClick={() => { setActiveTab(id); setError(null); }}
            className={`flex items-center gap-2 px-4 py-2.5 rounded-lg text-sm font-medium transition-all duration-200 cursor-pointer ${
              activeTab === id
                ? 'bg-[#76B900]/15 text-[#76B900] border border-[#76B900]/25'
                : 'text-[#9CA3AF] hover:text-white border border-transparent'
            }`}
          >
            <Icon size={16} />
            {label}
          </button>
        ))}
      </div>

      {/* Error */}
      {error && (
        <div className="bg-red-500/10 border border-red-500/30 rounded-xl p-4 flex items-center gap-3">
          <AlertTriangle size={18} className="text-red-400 shrink-0" />
          <p className="text-sm text-red-400">{error}</p>
        </div>
      )}

      {/* ── Patient Lookup Tab ─────────────────────────────────── */}
      {activeTab === 'lookup' && (
        <div className="space-y-6">
          {/* Search bar */}
          <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 shadow-lg shadow-black/20">
            <div className="flex gap-3">
              <input
                type="text"
                value={patientId}
                onChange={(e) => setPatientId(e.target.value)}
                onKeyDown={(e) => e.key === 'Enter' && handleLookup()}
                placeholder="Enter Patient ID..."
                className={`flex-1 ${inputClass}`}
              />
              <button
                onClick={handleLookup}
                disabled={loading || !patientId.trim()}
                className="flex items-center gap-2 bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-40 text-black font-semibold px-6 py-2.5 rounded-lg transition-all duration-200 cursor-pointer shadow-md shadow-[#76B900]/20"
              >
                {loading ? <Loader2 size={16} className="animate-spin" /> : <Search size={16} />}
                Lookup
              </button>
            </div>
          </div>

          {/* Results */}
          {cumulative && (
            <div className="space-y-4">
              {/* Top stats row */}
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                {/* Cumulative Dose */}
                <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 shadow-lg shadow-black/20">
                  <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Cumulative Dose</span>
                  <div className="mt-2 flex items-baseline gap-1">
                    <span className="text-4xl font-bold text-white tracking-tight">
                      {cumulative.total_effective_dose_msv?.toFixed(1) ?? '--'}
                    </span>
                    <span className="text-lg text-[#9CA3AF]">mSv</span>
                  </div>
                </div>

                {/* Alert Level */}
                <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 shadow-lg shadow-black/20">
                  <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Alert Level</span>
                  <div className="mt-2">
                    {(() => {
                      const style = getAlertStyle(cumulative.alert_level);
                      return (
                        <span className={`inline-flex items-center gap-2 text-lg font-bold px-4 py-1.5 rounded-lg ${style.bg} ${style.text} border ${style.border}`}>
                          {style.label}
                        </span>
                      );
                    })()}
                  </div>
                  {cumulative.alert_message && (
                    <p className="text-xs text-yellow-400 mt-2">{cumulative.alert_message}</p>
                  )}
                </div>

                {/* Study Count & Range */}
                <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 shadow-lg shadow-black/20">
                  <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Studies</span>
                  <div className="mt-2 flex items-baseline gap-1">
                    <span className="text-4xl font-bold text-white tracking-tight">
                      {cumulative.study_count ?? '--'}
                    </span>
                    <span className="text-sm text-[#9CA3AF]">studies</span>
                  </div>
                  {cumulative.date_range && (
                    <p className="text-[11px] text-[#9CA3AF] mt-2">
                      {cumulative.date_range.earliest || '?'} - {cumulative.date_range.latest || '?'}
                    </p>
                  )}
                </div>
              </div>

              {/* Dose Gauge */}
              {cumulative.total_effective_dose_msv != null && (
                <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 shadow-lg shadow-black/20">
                  <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold mb-3 block">Dose Threshold Position</span>
                  <DoseGauge dose={cumulative.total_effective_dose_msv} />
                </div>
              )}

              {/* Breakdown Charts */}
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {modalityData.length > 0 && (
                  <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 shadow-lg shadow-black/20">
                    <div className="flex items-center gap-2 mb-4">
                      <BarChart3 size={14} className="text-[#76B900]" />
                      <span className="text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold">By Modality (mSv)</span>
                    </div>
                    <ResponsiveContainer width="100%" height={180}>
                      <BarChart data={modalityData} layout="vertical" margin={{ left: 60, right: 20, top: 0, bottom: 0 }}>
                        <XAxis type="number" tick={{ fill: '#9CA3AF', fontSize: 11 }} axisLine={false} tickLine={false} />
                        <YAxis type="category" dataKey="name" tick={{ fill: '#E0E0E0', fontSize: 11 }} axisLine={false} tickLine={false} width={55} />
                        <Tooltip
                          contentStyle={{ backgroundColor: '#1A1D23', border: '1px solid rgba(255,255,255,0.08)', borderRadius: '8px', fontSize: '12px' }}
                          labelStyle={{ color: '#E0E0E0' }}
                          formatter={(value) => [`${Number(value).toFixed(1)} mSv`, 'Dose']}
                        />
                        <Bar dataKey="value" radius={[0, 4, 4, 0]} maxBarSize={20}>
                          {modalityData.map((_, i) => (
                            <Cell key={i} fill={chartColors[i % chartColors.length]} />
                          ))}
                        </Bar>
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                )}

                {regionData.length > 0 && (
                  <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-5 shadow-lg shadow-black/20">
                    <div className="flex items-center gap-2 mb-4">
                      <BarChart3 size={14} className="text-cyan-400" />
                      <span className="text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold">By Body Region (mSv)</span>
                    </div>
                    <ResponsiveContainer width="100%" height={180}>
                      <BarChart data={regionData} layout="vertical" margin={{ left: 60, right: 20, top: 0, bottom: 0 }}>
                        <XAxis type="number" tick={{ fill: '#9CA3AF', fontSize: 11 }} axisLine={false} tickLine={false} />
                        <YAxis type="category" dataKey="name" tick={{ fill: '#E0E0E0', fontSize: 11 }} axisLine={false} tickLine={false} width={55} />
                        <Tooltip
                          contentStyle={{ backgroundColor: '#1A1D23', border: '1px solid rgba(255,255,255,0.08)', borderRadius: '8px', fontSize: '12px' }}
                          labelStyle={{ color: '#E0E0E0' }}
                          formatter={(value) => [`${Number(value).toFixed(1)} mSv`, 'Dose']}
                        />
                        <Bar dataKey="value" radius={[0, 4, 4, 0]} maxBarSize={20}>
                          {regionData.map((_, i) => (
                            <Cell key={i} fill={chartColors[(i + 3) % chartColors.length]} />
                          ))}
                        </Bar>
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                )}
              </div>
            </div>
          )}
        </div>
      )}

      {/* ── Record Dose Tab ────────────────────────────────────── */}
      {activeTab === 'record' && (
        <div className="max-w-xl">
          <form
            onSubmit={handleRecord}
            className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 shadow-lg shadow-black/20 space-y-4"
          >
            <div>
              <label className="block text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2">Patient ID</label>
              <input
                type="text"
                value={recordForm.patient_id}
                onChange={(e) => setRecordForm({ ...recordForm, patient_id: e.target.value })}
                placeholder="Enter patient identifier"
                className={inputClass}
              />
            </div>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="block text-[11px] text-[#9CA3AF] mb-1">Modality</label>
                <select
                  value={recordForm.modality}
                  onChange={(e) => setRecordForm({ ...recordForm, modality: e.target.value })}
                  className={inputClass}
                >
                  <option>CT</option>
                  <option>PET/CT</option>
                  <option>Fluoroscopy</option>
                  <option>X-Ray</option>
                </select>
              </div>
              <div>
                <label className="block text-[11px] text-[#9CA3AF] mb-1">Body Region</label>
                <select
                  value={recordForm.body_region}
                  onChange={(e) => setRecordForm({ ...recordForm, body_region: e.target.value })}
                  className={inputClass}
                >
                  <option>Chest</option>
                  <option>Abdomen</option>
                  <option>Head</option>
                  <option>Pelvis</option>
                  <option>Spine</option>
                </select>
              </div>
            </div>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="block text-[11px] text-[#9CA3AF] mb-1">Dose (mGy)</label>
                <input
                  type="number"
                  step="0.1"
                  value={recordForm.dose_mgy}
                  onChange={(e) => setRecordForm({ ...recordForm, dose_mgy: e.target.value })}
                  placeholder="0.0"
                  className={inputClass}
                />
              </div>
              <div>
                <label className="block text-[11px] text-[#9CA3AF] mb-1">DLP (mGy*cm)</label>
                <input
                  type="number"
                  step="0.1"
                  value={recordForm.dlp_mgycm}
                  onChange={(e) => setRecordForm({ ...recordForm, dlp_mgycm: e.target.value })}
                  placeholder="0.0"
                  className={inputClass}
                />
              </div>
            </div>
            <button
              type="submit"
              disabled={recordLoading || !recordForm.patient_id}
              className="flex items-center justify-center gap-2 w-full bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-40 text-black font-semibold px-5 py-3 rounded-lg transition-all duration-200 cursor-pointer shadow-md shadow-[#76B900]/20"
            >
              {recordLoading ? <Loader2 size={16} className="animate-spin" /> : <Radiation size={16} />}
              {recordLoading ? 'Recording...' : 'Record Dose'}
            </button>
            {recordSuccess && (
              <div className="bg-green-500/10 border border-green-500/30 rounded-lg p-3 flex items-center gap-2">
                <CheckCircle size={16} className="text-green-400" />
                <p className="text-sm text-green-400">Dose recorded successfully.</p>
              </div>
            )}
          </form>
        </div>
      )}

      {/* ── DRL Comparison Tab ─────────────────────────────────── */}
      {activeTab === 'drl' && (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Form */}
          <form
            onSubmit={handleDRLCompare}
            className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 shadow-lg shadow-black/20 space-y-4"
          >
            <div>
              <label className="block text-xs uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2">Protocol</label>
              <select
                value={drlForm.protocol}
                onChange={(e) => setDrlForm({ ...drlForm, protocol: e.target.value })}
                className={inputClass}
              >
                {['CT Chest', 'CT Abdomen-Pelvis', 'CT Head', 'CT Coronary Angiography', 'CT Spine', 'CT Chest-Abdomen-Pelvis'].map(p => (
                  <option key={p} value={p}>{p}</option>
                ))}
              </select>
            </div>
            <div>
              <label className="block text-[11px] text-[#9CA3AF] mb-1">Patient Dose (mSv)</label>
              <input
                type="number"
                step="0.1"
                value={drlForm.effective_dose_msv}
                onChange={(e) => setDrlForm({ ...drlForm, effective_dose_msv: e.target.value })}
                placeholder="e.g., 8.5"
                className={inputClass}
              />
            </div>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="block text-[11px] text-[#9CA3AF] mb-1">Modality</label>
                <select
                  value={drlForm.modality}
                  onChange={(e) => setDrlForm({ ...drlForm, modality: e.target.value })}
                  className={inputClass}
                >
                  <option>CT</option>
                  <option>PET/CT</option>
                  <option>Fluoroscopy</option>
                </select>
              </div>
              <div>
                <label className="block text-[11px] text-[#9CA3AF] mb-1">Body Region</label>
                <select
                  value={drlForm.body_region}
                  onChange={(e) => setDrlForm({ ...drlForm, body_region: e.target.value })}
                  className={inputClass}
                >
                  <option value="chest">Chest</option>
                  <option value="abdomen">Abdomen</option>
                  <option value="head">Head</option>
                  <option value="pelvis">Pelvis</option>
                  <option value="spine">Spine</option>
                </select>
              </div>
            </div>
            <label className="flex items-center gap-2 cursor-pointer group">
              <div className={`w-4 h-4 rounded border transition-all duration-200 flex items-center justify-center ${
                drlForm.pediatric ? 'bg-[#76B900] border-[#76B900]' : 'border-white/20 group-hover:border-white/40'
              }`}>
                {drlForm.pediatric && (
                  <svg width="10" height="8" viewBox="0 0 10 8" fill="none">
                    <path d="M1 4L3.5 6.5L9 1" stroke="black" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round"/>
                  </svg>
                )}
              </div>
              <span className="text-sm text-[#E0E0E0]">Pediatric patient</span>
            </label>
            <button
              type="submit"
              disabled={drlLoading || !drlForm.effective_dose_msv}
              className="flex items-center justify-center gap-2 w-full bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-40 text-black font-semibold px-5 py-3 rounded-lg transition-all duration-200 cursor-pointer shadow-md shadow-[#76B900]/20"
            >
              {drlLoading ? <Loader2 size={16} className="animate-spin" /> : <Scale size={16} />}
              {drlLoading ? 'Comparing...' : 'Compare to DRL'}
            </button>
          </form>

          {/* DRL Results */}
          {drlResult ? (
            <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 shadow-lg shadow-black/20 space-y-5">
              {/* Protocol */}
              <div>
                <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Protocol</span>
                <h3 className="text-lg font-bold text-white mt-1">{drlResult.protocol}</h3>
              </div>

              {/* Visual comparison bar */}
              {drlResult.patient_dose_msv != null && drlResult.drl_msv != null && drlResult.achievable_dose_msv != null && (
                <div className="space-y-3">
                  <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Dose Comparison</span>
                  <div className="relative h-10 bg-white/5 rounded-lg overflow-visible">
                    {/* Achievable marker */}
                    {(() => {
                      const maxVal = Math.max(drlResult.patient_dose_msv!, drlResult.drl_msv!) * 1.3;
                      const achievPct = Math.min(100, (drlResult.achievable_dose_msv! / maxVal) * 100);
                      const drlPct = Math.min(100, (drlResult.drl_msv! / maxVal) * 100);
                      const patientPct = Math.min(100, (drlResult.patient_dose_msv! / maxVal) * 100);
                      const patientColor = drlResult.patient_dose_msv! <= drlResult.achievable_dose_msv! ? '#22c55e' :
                        drlResult.patient_dose_msv! <= drlResult.drl_msv! ? '#eab308' : '#ef4444';
                      return (
                        <>
                          {/* Achievable zone */}
                          <div className="absolute top-0 bottom-0 left-0 bg-green-500/10 rounded-l-lg" style={{ width: `${achievPct}%` }} />
                          {/* DRL line */}
                          <div className="absolute top-0 bottom-0 w-0.5 bg-yellow-500" style={{ left: `${drlPct}%` }}>
                            <span className="absolute -top-5 left-1/2 -translate-x-1/2 text-[9px] text-yellow-400 whitespace-nowrap">DRL {drlResult.drl_msv} mSv</span>
                          </div>
                          {/* Achievable line */}
                          <div className="absolute top-0 bottom-0 w-0.5 bg-green-500" style={{ left: `${achievPct}%` }}>
                            <span className="absolute -bottom-5 left-1/2 -translate-x-1/2 text-[9px] text-green-400 whitespace-nowrap">Achievable {drlResult.achievable_dose_msv} mSv</span>
                          </div>
                          {/* Patient dose marker */}
                          <div
                            className="absolute top-1/2 -translate-y-1/2 w-4 h-4 rounded-full border-2 shadow-lg"
                            style={{ left: `${patientPct}%`, backgroundColor: patientColor, borderColor: patientColor, transform: `translate(-50%, -50%)` }}
                          />
                        </>
                      );
                    })()}
                  </div>
                </div>
              )}

              {/* Ratio */}
              {drlResult.ratio != null && (
                <div className="bg-[#0E1117] rounded-lg p-4 border border-white/[0.06] flex items-center justify-between">
                  <div>
                    <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold">Dose-to-DRL Ratio</span>
                    <p className="text-2xl font-bold text-white mt-1">{drlResult.ratio.toFixed(2)}x</p>
                  </div>
                  <span className={`text-sm font-semibold px-3 py-1.5 rounded-lg ${
                    drlResult.status === 'below_achievable' ? 'bg-green-500/10 text-green-400 border border-green-500/30' :
                    drlResult.status === 'below_drl' ? 'bg-yellow-500/10 text-yellow-400 border border-yellow-500/30' :
                    drlResult.status === 'above_drl' ? 'bg-orange-500/10 text-orange-400 border border-orange-500/30' :
                    'bg-red-500/10 text-red-400 border border-red-500/30'
                  }`}>
                    {drlResult.status?.replace(/_/g, ' ').toUpperCase()}
                  </span>
                </div>
              )}

              {/* Optimization suggestions */}
              {drlResult.optimization_suggestions && drlResult.optimization_suggestions.length > 0 && (
                <div>
                  <span className="text-[10px] uppercase tracking-wider text-[#9CA3AF] font-semibold mb-2 block">Optimization Suggestions</span>
                  <div className="space-y-2">
                    {drlResult.optimization_suggestions.map((s, i) => (
                      <div key={i} className="bg-[#0E1117] rounded-lg p-3 border border-white/[0.06] flex items-start gap-3">
                        <span className="text-[#76B900] font-bold text-xs mt-0.5 shrink-0">{i + 1}</span>
                        <p className="text-sm text-[#E0E0E0] leading-relaxed">{s}</p>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          ) : (
            <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-8 flex flex-col items-center justify-center text-center">
              <div className="p-4 bg-white/5 rounded-2xl mb-4">
                <Scale size={40} className="text-white/20" />
              </div>
              <p className="text-sm text-[#9CA3AF]">Enter a protocol and dose to compare against national Diagnostic Reference Levels</p>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
