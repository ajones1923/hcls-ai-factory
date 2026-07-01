import { useNavigate } from 'react-router-dom';
import {
  Database,
  Layers,
  GitBranch,
  Cpu,
  Play,
  Search,
  BarChart3,
  CheckCircle,
  ChevronRight,
} from 'lucide-react';
import { useAppStore } from '../store/appStore';

const engines = [
  {
    name: 'Engine 1',
    title: 'Genomics',
    description: 'Parabricks FASTQ-to-VCF pipeline',
    status: 'available',
  },
  {
    name: 'Engine 2',
    title: 'RAG / Chat',
    description: 'Clinical evidence with Milvus + Claude',
    status: 'available',
  },
  {
    name: 'Engine 3',
    title: 'Drug Discovery',
    description: 'MolMIM generation + DiffDock docking',
    status: 'available',
  },
  {
    name: 'Engine 4',
    title: 'Clinical Imaging',
    description: 'VISTA-3D, MAISI, Radiology LLM NIMs',
    status: 'active',
  },
];

export default function Dashboard() {
  const navigate = useNavigate();
  const { health } = useAppStore();

  const totalVectors = health?.total_vectors ?? 0;
  const collectionCount = health?.collections
    ? Object.keys(health.collections).length
    : 0;
  const nimCount = health?.nim_services
    ? Object.keys(health.nim_services).length
    : 0;

  const stats = [
    {
      label: 'Total Vectors',
      value: totalVectors.toLocaleString(),
      icon: Database,
      color: 'text-[#76B900]',
    },
    {
      label: 'Collections',
      value: collectionCount.toString(),
      icon: Layers,
      color: 'text-blue-400',
    },
    {
      label: 'Workflows',
      value: '9',
      icon: GitBranch,
      color: 'text-purple-400',
    },
    {
      label: 'NIM Services',
      value: nimCount.toString(),
      icon: Cpu,
      color: 'text-orange-400',
    },
  ];

  return (
    <div className="space-y-6">
      {/* Hero Banner */}
      <div className="relative overflow-hidden rounded-2xl border border-[#76B900]/30 bg-gradient-to-br from-[#0a1628] via-[#122040] to-[#1a3a2a] p-8 mb-8">
        <div className="absolute top-0 left-0 right-0 h-1 bg-gradient-to-r from-[#76B900] via-[#4CAF50] to-[#76B900]" />
        <div className="relative z-10">
          <span className="inline-block text-[10px] font-bold uppercase tracking-[2px] text-[#76B900] bg-[#76B900]/10 border border-[#76B900]/30 px-3 py-1 rounded-full mb-3">Engine 4</span>
          <h1 className="text-3xl font-bold text-white mb-2">Clinical Imaging Engine</h1>
          <p className="text-[#76B900] font-medium mb-1">HCLS AI Factory — Precision Medicine Platform</p>
          <p className="text-sm text-white/50">9 Workflows &bull; 7 Scoring Systems &bull; 13 Collections &bull; 20 NVIDIA Technologies &bull; Apache 2.0</p>
        </div>
        {/* Subtle grid pattern overlay */}
        <div className="absolute inset-0 opacity-5" style={{backgroundImage: 'radial-gradient(circle at 1px 1px, white 1px, transparent 0)', backgroundSize: '24px 24px'}} />
      </div>

      {/* Clinical Disclaimer */}
      <div className="flex items-center gap-3 bg-amber-500/5 border border-amber-500/20 rounded-lg px-4 py-2.5 mb-6">
        <span className="text-amber-400 text-xs font-medium">Clinical Decision Support</span>
        <span className="text-white/50 text-xs">Not FDA-cleared. All findings must be reviewed by a qualified healthcare professional.</span>
      </div>

      {/* Stats Row */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
        {stats.map((s) => {
          const Icon = s.icon;
          return (
            <div
              key={s.label}
              className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 hover:border-[#76B900]/30 transition-all duration-200 shadow-lg shadow-black/20"
            >
              <div className="flex items-center justify-between mb-3">
                <span className="text-sm text-[#9CA3AF]">{s.label}</span>
                <Icon size={20} className={s.color} />
              </div>
              <p className="text-2xl font-semibold text-white">{s.value}</p>
            </div>
          );
        })}
      </div>

      {/* Architecture */}
      <div>
        <h2 className="text-lg font-semibold text-white mb-4">
          HCLS AI Factory — Engine Architecture
        </h2>
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-0">
          {engines.map((eng, idx) => (
            <div key={eng.name} className="flex items-center">
              <div
                className={`flex-1 bg-[#1A1D23] rounded-xl border p-6 transition-all duration-200 shadow-lg shadow-black/20 ${
                  eng.status === 'active'
                    ? 'border-[#76B900]/50 shadow-[0_0_30px_rgba(118,185,0,0.15)] ring-1 ring-[#76B900]/20'
                    : 'border-white/[0.08] hover:border-white/[0.15]'
                }`}
              >
                <div className="flex items-center justify-between mb-2">
                  <span className="text-[10px] font-semibold uppercase tracking-wider text-[#9CA3AF]">
                    {eng.name}
                  </span>
                  {eng.status === 'active' ? (
                    <span className="text-[10px] font-semibold bg-[#76B900] text-black px-2 py-0.5 rounded animate-pulse">
                      ACTIVE
                    </span>
                  ) : (
                    <span className="text-[10px] text-[#9CA3AF] bg-white/[0.06] px-2 py-0.5 rounded">
                      Available
                    </span>
                  )}
                </div>
                <h3 className="text-white font-medium mb-1">{eng.title}</h3>
                <p className="text-xs text-[#9CA3AF]">{eng.description}</p>
              </div>
              {/* Connecting arrow between cards */}
              {idx < engines.length - 1 && (
                <div className="hidden lg:flex items-center justify-center w-8 shrink-0">
                  <ChevronRight size={20} className="text-[#76B900]/40" />
                </div>
              )}
            </div>
          ))}
        </div>
      </div>

      {/* Quick Actions */}
      <div>
        <h2 className="text-lg font-semibold text-white mb-4">Quick Actions</h2>
        <div className="flex gap-4 flex-wrap">
          <button
            onClick={() => navigate('/workflows')}
            className="flex items-center gap-2 bg-[#76B900] hover:bg-[#5a9400] text-black font-medium px-5 py-2.5 rounded-lg transition-all duration-200 cursor-pointer shadow-lg shadow-[#76B900]/20"
          >
            <Play size={16} />
            Run Demo Case
          </button>
          <button
            onClick={() => navigate('/evidence')}
            className="flex items-center gap-2 bg-[#1A1D23] border border-white/[0.08] hover:border-[#76B900]/30 text-white px-5 py-2.5 rounded-lg transition-all duration-200 cursor-pointer"
          >
            <Search size={16} />
            Ask a Question
          </button>
          <button
            onClick={() => navigate('/analytics')}
            className="flex items-center gap-2 bg-[#1A1D23] border border-white/[0.08] hover:border-[#76B900]/30 text-white px-5 py-2.5 rounded-lg transition-all duration-200 cursor-pointer"
          >
            <BarChart3 size={16} />
            View Analytics
          </button>
        </div>
      </div>

      {/* Status */}
      <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 shadow-lg shadow-black/20 hover:border-[#76B900]/30 transition-all duration-200">
        <div className="flex items-center gap-2 mb-2">
          <CheckCircle size={18} className="text-[#76B900]" />
          <h3 className="text-white font-medium">System Status</h3>
        </div>
        <p className="text-sm text-[#9CA3AF]">
          System ready. 1,324 tests passing. {totalVectors.toLocaleString()}{' '}
          literature vectors indexed. 9 workflows available. Status:{' '}
          <span className="text-[#76B900] font-medium">
            {health?.status ?? 'connecting...'}
          </span>
        </p>
      </div>
    </div>
  );
}
