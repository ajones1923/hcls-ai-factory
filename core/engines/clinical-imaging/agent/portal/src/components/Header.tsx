import { useLocation } from 'react-router-dom';

const pageTitles: Record<string, string> = {
  '/': 'Dashboard',
  '/workflows': 'Workflows',
  '/evidence': 'Evidence Explorer',
  '/protocol': 'Protocol Recommendation',
  '/dose': 'Dose Tracking',
  '/analytics': 'Analytics',
  '/reports': 'Reports',
  '/benchmarks': 'Benchmarks',
};

export default function Header() {
  const location = useLocation();
  const title = pageTitles[location.pathname] || 'Dashboard';

  return (
    <header className="h-12 border-b border-white/[0.06] flex items-center justify-between px-6 bg-[#0E1117]/80 backdrop-blur-sm">
      <div className="flex items-center gap-2 text-sm text-[#9CA3AF]">
        <span className="text-white font-medium">Clinical Imaging Engine</span>
        <span className="text-white/20">/</span>
        <span>{title}</span>
      </div>
      <span className="text-[10px] text-amber-400/70 bg-amber-400/5 border border-amber-400/10 px-2 py-0.5 rounded">
        CDS — Not FDA-cleared
      </span>
    </header>
  );
}
