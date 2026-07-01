import { useEffect, useState } from 'react';
import { BarChart3, Loader2, RefreshCw } from 'lucide-react';
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  PieChart,
  Pie,
  Cell,
} from 'recharts';
import { fetchAnalyticsPopulation, generateDemoData, fetchTrends } from '../lib/api';

const COLORS = ['#76B900', '#4FC3F7', '#FF9800', '#EF5350', '#AB47BC', '#26A69A'];

interface PopulationData {
  total_studies?: number;
  modality_distribution?: Record<string, number>;
  severity_distribution?: Record<string, number>;
  body_region_distribution?: Record<string, number>;
  [key: string]: unknown;
}

interface TrendPoint {
  period: string;
  count: number;
  rate?: number;
}

export default function Analytics() {
  const [population, setPopulation] = useState<PopulationData | null>(null);
  const [trends, setTrends] = useState<TrendPoint[]>([]);
  const [loading, setLoading] = useState(true);
  const [generating, setGenerating] = useState(false);
  const [selectedMetric, setSelectedMetric] = useState('studies');
  const [refreshKey, setRefreshKey] = useState(0);

  useEffect(() => {
    let cancelled = false;
    const load = async () => {
      setLoading(true);
      try {
        const [pop, tr] = await Promise.all([
          fetchAnalyticsPopulation().catch(() => null),
          fetchTrends(selectedMetric).catch(() => ({ data_points: [] })),
        ]);
        if (!cancelled) {
          setPopulation(pop);
          setTrends(tr?.data_points || tr?.data || tr?.trends || []);
        }
      } finally {
        if (!cancelled) setLoading(false);
      }
    };
    load();
    return () => { cancelled = true; };
  }, [selectedMetric, refreshKey]);

  const handleGenerate = async () => {
    setGenerating(true);
    try {
      await generateDemoData(50);
      setRefreshKey((k) => k + 1);
    } finally {
      setGenerating(false);
    }
  };

  const toChartData = (dist?: Record<string, number>) =>
    dist
      ? Object.entries(dist).map(([name, value]) => ({ name, value }))
      : [];

  if (loading) {
    return (
      <div className="flex items-center justify-center h-64">
        <Loader2 size={24} className="animate-spin text-[#76B900]" />
      </div>
    );
  }

  const modalityData = toChartData(population?.modality_distribution);
  const severityData = toChartData(population?.severity_distribution);

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-lg font-semibold text-white">Population Analytics</h2>
          <p className="text-sm text-[#9CA3AF]">
            {population?.total_studies?.toLocaleString() ?? 0} total studies analyzed
          </p>
        </div>
        <div className="flex gap-3">
          <button
            onClick={handleGenerate}
            disabled={generating}
            className="flex items-center gap-2 bg-[#1A1D23] border border-white/[0.08] hover:border-[#76B900]/30 text-white text-sm px-4 py-2 rounded-lg transition-all duration-200 cursor-pointer"
          >
            {generating ? (
              <Loader2 size={14} className="animate-spin" />
            ) : (
              <RefreshCw size={14} />
            )}
            Generate Demo Data
          </button>
          <select
            value={selectedMetric}
            onChange={(e) => setSelectedMetric(e.target.value)}
            className="bg-[#1A1D23] border border-white/[0.08] rounded-lg px-3 py-2 text-sm text-white focus:outline-none focus:border-[#76B900]/50"
          >
            <option value="studies">Studies</option>
            <option value="dose">Dose</option>
            <option value="findings">Findings</option>
          </select>
        </div>
      </div>

      {/* Charts */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Modality Distribution */}
        <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 hover:border-[#76B900]/30 transition-all duration-200 shadow-lg shadow-black/20">
          <h3 className="text-sm font-medium text-white mb-4 flex items-center gap-2">
            <BarChart3 size={16} className="text-[#76B900]" />
            Modality Distribution
          </h3>
          {modalityData.length > 0 ? (
            <ResponsiveContainer width="100%" height={250}>
              <BarChart data={modalityData}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.06)" />
                <XAxis
                  dataKey="name"
                  tick={{ fill: '#9CA3AF', fontSize: 11 }}
                  axisLine={{ stroke: 'rgba(255,255,255,0.08)' }}
                />
                <YAxis
                  tick={{ fill: '#9CA3AF', fontSize: 11 }}
                  axisLine={{ stroke: 'rgba(255,255,255,0.08)' }}
                />
                <Tooltip
                  contentStyle={{
                    backgroundColor: '#1E2229',
                    border: '1px solid rgba(255,255,255,0.08)',
                    borderRadius: 8,
                    color: '#E0E0E0',
                    fontSize: 12,
                  }}
                />
                <Bar dataKey="value" fill="#76B900" radius={[4, 4, 0, 0]} />
              </BarChart>
            </ResponsiveContainer>
          ) : (
            <p className="text-sm text-[#9CA3AF] text-center py-10">
              No data available. Generate demo data to see charts.
            </p>
          )}
        </div>

        {/* Severity Distribution */}
        <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 hover:border-[#76B900]/30 transition-all duration-200 shadow-lg shadow-black/20">
          <h3 className="text-sm font-medium text-white mb-4">
            Severity Distribution
          </h3>
          {severityData.length > 0 ? (
            <ResponsiveContainer width="100%" height={250}>
              <PieChart>
                <Pie
                  data={severityData}
                  cx="50%"
                  cy="50%"
                  innerRadius={60}
                  outerRadius={100}
                  paddingAngle={2}
                  dataKey="value"
                  label={({ name, percent }: { name?: string; percent?: number }) =>
                    `${name ?? ''} ${((percent ?? 0) * 100).toFixed(0)}%`
                  }
                >
                  {severityData.map((_, i) => (
                    <Cell key={i} fill={COLORS[i % COLORS.length]} />
                  ))}
                </Pie>
                <Tooltip
                  contentStyle={{
                    backgroundColor: '#1E2229',
                    border: '1px solid rgba(255,255,255,0.08)',
                    borderRadius: 8,
                    color: '#E0E0E0',
                    fontSize: 12,
                  }}
                />
              </PieChart>
            </ResponsiveContainer>
          ) : (
            <p className="text-sm text-[#9CA3AF] text-center py-10">
              No data available.
            </p>
          )}
        </div>
      </div>

      {/* Trends */}
      {trends.length > 0 && (
        <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 hover:border-[#76B900]/30 transition-all duration-200 shadow-lg shadow-black/20">
          <h3 className="text-sm font-medium text-white mb-4">
            Trend: {selectedMetric}
          </h3>
          <ResponsiveContainer width="100%" height={200}>
            <BarChart data={trends}>
              <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.06)" />
              <XAxis
                dataKey="period"
                tick={{ fill: '#9CA3AF', fontSize: 10 }}
                axisLine={{ stroke: 'rgba(255,255,255,0.08)' }}
              />
              <YAxis
                tick={{ fill: '#9CA3AF', fontSize: 10 }}
                axisLine={{ stroke: 'rgba(255,255,255,0.08)' }}
              />
              <Tooltip
                contentStyle={{
                  backgroundColor: '#1E2229',
                  border: '1px solid rgba(255,255,255,0.08)',
                  borderRadius: 8,
                  color: '#E0E0E0',
                  fontSize: 12,
                }}
              />
              <Bar dataKey="count" fill="#4FC3F7" radius={[4, 4, 0, 0]} />
            </BarChart>
          </ResponsiveContainer>
        </div>
      )}
    </div>
  );
}
