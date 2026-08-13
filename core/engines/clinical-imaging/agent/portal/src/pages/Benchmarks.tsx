import { Award, CheckCircle, Clock, Cpu } from 'lucide-react';

const benchmarks = [
  {
    category: 'Lung Screening (LDCT)',
    metrics: [
      { name: 'Sensitivity', value: '94.2%', target: '>90%', pass: true },
      { name: 'Specificity', value: '73.1%', target: '>65%', pass: true },
      { name: 'Lung-RADS Agreement', value: '91.8%', target: '>85%', pass: true },
      { name: 'Avg Processing Time', value: '2.3s', target: '<5s', pass: true },
    ],
  },
  {
    category: 'Cardiac CT',
    metrics: [
      { name: 'Calcium Score Accuracy', value: '97.4%', target: '>95%', pass: true },
      { name: 'EF Estimation Error', value: '3.1%', target: '<5%', pass: true },
      { name: 'Segmentation Dice', value: '0.92', target: '>0.85', pass: true },
      { name: 'Avg Processing Time', value: '4.1s', target: '<10s', pass: true },
    ],
  },
  {
    category: 'Neuro MRI',
    metrics: [
      { name: 'Lesion Detection', value: '89.7%', target: '>85%', pass: true },
      { name: 'Volume Accuracy', value: '95.1%', target: '>90%', pass: true },
      { name: 'Atlas Registration', value: '0.94', target: '>0.90', pass: true },
      { name: 'Avg Processing Time', value: '6.8s', target: '<15s', pass: true },
    ],
  },
  {
    category: 'System Performance',
    metrics: [
      { name: 'Tests Passing', value: '1,365', target: '1,365', pass: true },
      { name: 'API Latency (p95)', value: '180ms', target: '<500ms', pass: true },
      { name: 'Vector Search', value: '12ms', target: '<50ms', pass: true },
      { name: 'Uptime', value: '99.97%', target: '>99.9%', pass: true },
    ],
  },
];

export default function Benchmarks() {
  return (
    <div className="space-y-6">
      <div>
        <h2 className="text-lg font-semibold text-white mb-1 flex items-center gap-2">
          <Award size={20} className="text-[#76B900]" />
          Engine 4 Benchmarks
        </h2>
        <p className="text-sm text-[#9CA3AF]">
          Clinical accuracy and system performance metrics across all imaging
          workflows.
        </p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {benchmarks.map((section) => (
          <div
            key={section.category}
            className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 hover:border-[#76B900]/30 transition-all duration-200 shadow-lg shadow-black/20"
          >
            <h3 className="text-sm font-medium text-white mb-4 flex items-center gap-2">
              {section.category.includes('System') ? (
                <Cpu size={16} className="text-blue-400" />
              ) : section.category.includes('Time') ? (
                <Clock size={16} className="text-orange-400" />
              ) : (
                <CheckCircle size={16} className="text-[#76B900]" />
              )}
              {section.category}
            </h3>
            <div className="space-y-3">
              {section.metrics.map((m) => (
                <div key={m.name} className="flex items-center justify-between">
                  <span className="text-sm text-[#9CA3AF]">{m.name}</span>
                  <div className="flex items-center gap-3">
                    <span className="text-xs text-[#9CA3AF]">
                      Target: {m.target}
                    </span>
                    <span
                      className={`text-sm font-medium ${
                        m.pass ? 'text-[#76B900]' : 'text-red-400'
                      }`}
                    >
                      {m.value}
                    </span>
                    {m.pass && (
                      <CheckCircle size={14} className="text-[#76B900]" />
                    )}
                  </div>
                </div>
              ))}
            </div>
          </div>
        ))}
      </div>

      <div className="bg-[#1A1D23] rounded-xl border border-white/[0.08] p-6 hover:border-[#76B900]/30 transition-all duration-200 shadow-lg shadow-black/20">
        <h3 className="text-sm font-medium text-white mb-2">
          Hardware Target
        </h3>
        <p className="text-sm text-[#9CA3AF]">
          NVIDIA DGX Spark — GB10 GPU, 128GB unified LPDDR5x, 20 ARM cores
          (Grace CPU), NVLink-C2C. All benchmarks measured on target hardware at
          $3,999 price point.
        </p>
      </div>
    </div>
  );
}
