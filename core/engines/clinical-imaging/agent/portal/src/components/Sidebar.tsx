import { useEffect, useState } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import {
  LayoutDashboard,
  Stethoscope,
  Search,
  FileCheck,
  Radiation,
  BarChart3,
  FileText,
  Award,
  GitCompare,
  Zap,
} from 'lucide-react';
import { fetchHealth, describeApiError } from '../lib/api';
import { useAppStore } from '../store/appStore';

const navItems = [
  { id: 'dashboard', label: 'Dashboard', icon: LayoutDashboard, path: '/' },
  { id: 'workflows', label: 'Workflows', icon: Stethoscope, path: '/workflows' },
  { id: 'live', label: 'Live Analysis', icon: Zap, path: '/live-analysis' },
  { id: 'evidence', label: 'Evidence', icon: Search, path: '/evidence' },
  { id: 'protocol', label: 'Protocol', icon: FileCheck, path: '/protocol' },
  { id: 'dose', label: 'Dose Tracking', icon: Radiation, path: '/dose' },
  { id: 'analytics', label: 'Analytics', icon: BarChart3, path: '/analytics' },
  { id: 'reports', label: 'Reports', icon: FileText, path: '/reports' },
  { id: 'benchmarks', label: 'Benchmarks', icon: Award, path: '/benchmarks' },
  { id: 'compare', label: 'Compare', icon: GitCompare, path: '/compare' },
];

const nimServiceLabels: Record<string, string> = {
  vista3d: 'VISTA-3D',
  maisi: 'MAISI',
  vila_m3: 'VILA-M3',
  llm: 'Llama-3 / Claude',
};

export default function Sidebar() {
  const location = useLocation();
  const navigate = useNavigate();
  const { health, setHealth } = useAppStore();
  const [nimStatus, setNimStatus] = useState<Record<string, string>>({});

  useEffect(() => {
    const load = () => {
      fetchHealth()
        .then((data) => {
          setHealth(data);
          setNimStatus(data.nim_services || {});
        })
        .catch((e) => describeApiError('Loading engine status', e));
    };
    load();
    const interval = setInterval(load, 15000);
    return () => clearInterval(interval);
  }, [setHealth]);

  const statusDot = (status: string) => {
    if (status === 'available' || status === 'running')
      return 'bg-[#76B900]';
    if (status === 'cloud')
      return 'bg-blue-400';
    if (status === 'mock')
      return 'bg-yellow-500';
    if (status === 'degraded' || status === 'warming')
      return 'bg-yellow-500';
    return 'bg-red-500';
  };

  const statusLabel = (status: string) => {
    if (status === 'available' || status === 'running') return '';
    if (status === 'cloud') return '(cloud)';
    if (status === 'mock') return '(mock)';
    return '';
  };

  return (
    <aside className="w-60 min-h-screen bg-[#1A1D23] border-r border-white/[0.08] flex flex-col">
      {/* Logo */}
      <div className="p-5 border-b border-white/[0.08]">
        <div className="flex items-center gap-2 mb-1">
          <h1 className="text-base font-semibold text-white">Clinical Imaging Engine</h1>
        </div>
        <div className="flex items-center gap-2">
          <span className="text-[10px] font-semibold bg-[#76B900] text-black px-2 py-0.5 rounded">
            ENGINE 4
          </span>
          <span className="text-xs text-[#9CA3AF]">HCLS AI Factory</span>
        </div>
      </div>

      {/* Navigation */}
      <nav className="flex-1 py-3">
        {navItems.map((item) => {
          const isActive =
            item.path === '/'
              ? location.pathname === '/'
              : location.pathname.startsWith(item.path);
          const Icon = item.icon;
          return (
            <button
              key={item.id}
              onClick={() => navigate(item.path)}
              className={`w-full flex items-center gap-3 px-5 py-2.5 text-sm transition-all duration-200 cursor-pointer ${
                isActive
                  ? 'text-[#76B900] bg-[#76B900]/[0.08] border-l-2 border-[#76B900]'
                  : 'text-[#9CA3AF] hover:text-white hover:bg-white/[0.04] border-l-2 border-transparent'
              }`}
            >
              <Icon size={18} />
              <span>{item.label}</span>
            </button>
          );
        })}
      </nav>

      {/* NIM Status */}
      <div className="px-5 py-3 border-t border-white/[0.08]">
        <p className="text-[10px] uppercase tracking-wider text-[#9CA3AF] mb-2 font-medium">
          NIM Services
        </p>
        <div className="space-y-1.5">
          {Object.entries(nimServiceLabels).map(([key, label]) => {
            const status = nimStatus[key] || 'offline';
            return (
              <div key={key} className="flex items-center gap-2">
                <span className={`w-2 h-2 rounded-full ${statusDot(status)}`} />
                <span className="text-xs text-[#9CA3AF]">{label}</span>
                {statusLabel(status) && (
                  <span className="text-[10px] text-[#9CA3AF]/60">{statusLabel(status)}</span>
                )}
              </div>
            );
          })}
        </div>
      </div>

      {/* Collection Stats */}
      {health?.collections && (
        <div className="px-5 py-3 border-t border-white/[0.08]">
          <p className="text-[10px] uppercase tracking-wider text-[#9CA3AF] mb-2 font-medium">
            Collections
          </p>
          <div className="space-y-1">
            {Object.entries(health.collections).map(([name, count]) => {
              const shortName = name.replace('imaging_', '').replace('genomic_evidence', 'genomic');
              return (
                <div key={name} className="flex items-center justify-between">
                  <span className="text-xs text-[#9CA3AF] truncate" title={name}>{shortName}</span>
                  <span className="text-xs text-white font-medium">
                    {(count as number).toLocaleString()}
                  </span>
                </div>
              );
            })}
          </div>
        </div>
      )}

      {/* Footer */}
      <div className="px-5 py-3 border-t border-white/[0.08]">
        <p className="text-[10px] text-[#9CA3AF]">
          Apache 2.0 | 1,365 tests passing
        </p>
      </div>
    </aside>
  );
}
