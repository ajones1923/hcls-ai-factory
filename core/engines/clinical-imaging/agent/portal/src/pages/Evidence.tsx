import { useState, useRef, useEffect } from 'react';
import { Send, Loader2, BookOpen, Clock, Database, Search, Zap, MessageSquare, Sparkles } from 'lucide-react';
import { askQuestion } from '../lib/api';

interface Citation {
  title?: string;
  source?: string;
  pubmed_id?: string;
  relevance_score?: number;
  collection?: string;
  snippet?: string;
  id?: string;
  score?: number;
  text_snippet?: string;
  metadata?: {
    title?: string;
    source_type?: string;
    year?: number;
    modality?: string;
    body_region?: string;
    journal?: string;
    [key: string]: unknown;
  };
  label?: string;
}

interface Message {
  role: 'user' | 'assistant';
  content: string;
  citations?: Citation[];
  stats?: {
    evidence_count?: number;
    collections_searched?: number;
    search_time_ms?: number;
    total_time_ms?: number;
  };
  follow_up_questions?: string[];
}

const exampleQueries = [
  'What is the Lung-RADS classification system?',
  'How does VISTA-3D segment cardiac structures?',
  'What are the ACR TI-RADS criteria for thyroid nodules?',
  'Explain BI-RADS 4 breast lesion management',
  'What radiation dose limits apply to pediatric CT?',
  'How does DiffDock assist in drug-target docking?',
];

const suggestedFollowUps: Record<string, string[]> = {
  'lung': ['What are Lung-RADS 4A management recommendations?', 'How does AI improve lung nodule detection sensitivity?'],
  'cardiac': ['What is the CAD-RADS scoring system?', 'Compare CT vs MRI for cardiac assessment'],
  'breast': ['When should breast MRI supplement mammography?', 'BI-RADS 3 vs 4 management differences'],
  'thyroid': ['TI-RADS 4 vs 5 distinction criteria', 'FNA biopsy indications for thyroid nodules'],
  'dose': ['ALARA principle in pediatric imaging', 'Iterative reconstruction dose reduction strategies'],
};

function getFollowUps(question: string): string[] {
  const lower = question.toLowerCase();
  for (const [key, questions] of Object.entries(suggestedFollowUps)) {
    if (lower.includes(key)) return questions;
  }
  return ['What are the clinical implications?', 'How does AI improve diagnostic accuracy here?'];
}

export default function Evidence() {
  const [messages, setMessages] = useState<Message[]>([]);
  const [input, setInput] = useState('');
  const [loading, setLoading] = useState(false);
  const [selectedCitations, setSelectedCitations] = useState<Citation[]>([]);
  const messagesEndRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  const handleAsk = async (question?: string) => {
    const q = question || input.trim();
    if (!q || loading) return;

    const userMsg: Message = { role: 'user', content: q };
    setMessages((prev) => [...prev, userMsg]);
    setInput('');
    setLoading(true);

    try {
      const data = await askQuestion(q);
      const rawSources = data.sources || data.citations || data.evidence || [];
      const citations: Citation[] = rawSources.map((s: Citation) => ({
        ...s,
        title: s.metadata?.title || s.title || s.label || 'Source',
        snippet: s.text_snippet || s.snippet || '',
        relevance_score: s.score ?? s.relevance_score,
        collection: s.collection || s.metadata?.collection,
        pubmed_id: s.id || s.pubmed_id,
      }));
      const assistantMsg: Message = {
        role: 'assistant',
        content: data.answer || data.response || JSON.stringify(data),
        citations,
        stats: {
          evidence_count: data.evidence_count ?? citations.length,
          collections_searched: data.collections_searched ?? 0,
          search_time_ms: data.search_time_ms ?? 0,
          total_time_ms: data.total_time_ms ?? 0,
        },
        follow_up_questions: data.follow_up_questions || getFollowUps(q),
      };
      setMessages((prev) => [...prev, assistantMsg]);
      setSelectedCitations(assistantMsg.citations || []);
    } catch {
      setMessages((prev) => [
        ...prev,
        {
          role: 'assistant',
          content:
            'Unable to reach the evidence engine. Ensure the FastAPI backend is running on port 8524.',
        },
      ]);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="flex gap-6 h-[calc(100vh-7rem)]">
      {/* Chat Area */}
      <div className="flex-1 flex flex-col min-w-0">
        {/* Messages */}
        <div className="flex-1 overflow-auto space-y-4 pb-4">
          {messages.length === 0 && (
            <div className="flex flex-col items-center justify-center h-full text-center">
              <div className="p-5 bg-[#76B900]/10 rounded-2xl border border-[#76B900]/15 mb-5">
                <Search size={44} className="text-[#76B900]/60" />
              </div>
              <h2 className="text-2xl font-bold text-white mb-2 tracking-tight">
                Evidence Explorer
              </h2>
              <p className="text-sm text-[#9CA3AF] mb-8 max-w-md leading-relaxed">
                Ask clinical imaging questions backed by RAG-retrieved literature
                and guideline evidence.
              </p>
              <div className="flex flex-wrap gap-2.5 max-w-xl justify-center">
                {exampleQueries.map((q) => (
                  <button
                    key={q}
                    onClick={() => handleAsk(q)}
                    className="text-xs text-[#9CA3AF] bg-[#1A1D23] border border-white/[0.08] hover:border-[#76B900]/30 hover:text-white hover:bg-[#76B900]/5 px-4 py-2.5 rounded-full transition-all duration-200 cursor-pointer group"
                  >
                    <Sparkles size={10} className="inline mr-1.5 text-[#76B900]/40 group-hover:text-[#76B900]" />
                    {q}
                  </button>
                ))}
              </div>
            </div>
          )}

          {messages.map((msg, i) => (
            <div
              key={i}
              className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}
            >
              {msg.role === 'user' ? (
                <div className="max-w-[80%] bg-[#76B900]/10 border border-[#76B900]/20 rounded-2xl rounded-tr-sm px-4 py-3">
                  <p className="text-sm text-white whitespace-pre-wrap">{msg.content}</p>
                </div>
              ) : (
                <div className="max-w-[85%] space-y-0">
                  <div className="bg-[#1A1D23] border border-white/[0.08] rounded-2xl rounded-tl-sm overflow-hidden">
                    {/* Answer content */}
                    <div className="px-5 py-4">
                      <p className="text-sm text-[#E0E0E0] whitespace-pre-wrap leading-relaxed">{msg.content}</p>
                    </div>

                    {/* Stats bar - visually separated */}
                    {msg.stats && (
                      <div className="flex items-center gap-4 px-5 py-3 border-t border-white/[0.06] bg-white/[0.01] text-xs text-[#9CA3AF]">
                        {msg.stats.evidence_count != null && (
                          <span className="flex items-center gap-1.5">
                            <BookOpen size={12} className="text-[#76B900]" />
                            {msg.stats.evidence_count} sources
                          </span>
                        )}
                        {msg.stats.collections_searched != null && msg.stats.collections_searched > 0 && (
                          <span className="flex items-center gap-1.5">
                            <Database size={12} className="text-blue-400" />
                            {msg.stats.collections_searched} collections
                          </span>
                        )}
                        {msg.stats.search_time_ms != null && msg.stats.search_time_ms > 0 && (
                          <span className="flex items-center gap-1.5">
                            <Zap size={12} className="text-orange-400" />
                            {msg.stats.search_time_ms}ms search
                          </span>
                        )}
                        {msg.stats.total_time_ms != null && (
                          <span className="flex items-center gap-1.5">
                            <Clock size={12} className="text-purple-400" />
                            {msg.stats.total_time_ms}ms total
                          </span>
                        )}
                      </div>
                    )}

                    {/* Citation chips */}
                    {msg.citations && msg.citations.length > 0 && (
                      <div className="flex flex-wrap gap-1.5 px-5 py-3 border-t border-white/[0.06]">
                        {msg.citations.map((c, j) => (
                          <button
                            key={j}
                            onClick={() => setSelectedCitations(msg.citations || [])}
                            className="text-[10px] text-blue-400 bg-blue-500/8 border border-blue-500/20 px-2.5 py-1 rounded-full hover:bg-blue-500/15 hover:border-blue-500/30 transition-all duration-200 cursor-pointer"
                          >
                            [{j + 1}] {c.title?.slice(0, 40) || c.source || 'Source'}
                            {c.title && c.title.length > 40 ? '...' : ''}
                          </button>
                        ))}
                      </div>
                    )}
                  </div>

                  {/* Follow-up questions */}
                  {msg.follow_up_questions && msg.follow_up_questions.length > 0 && (
                    <div className="flex flex-wrap gap-2 mt-3 ml-2">
                      {msg.follow_up_questions.map((fq, j) => (
                        <button
                          key={j}
                          onClick={() => handleAsk(fq)}
                          className="text-[11px] text-[#76B900] bg-[#76B900]/8 border border-[#76B900]/20 px-3 py-1.5 rounded-full hover:bg-[#76B900]/15 hover:border-[#76B900]/30 transition-all duration-200 cursor-pointer flex items-center gap-1.5"
                        >
                          <MessageSquare size={10} />
                          {fq}
                        </button>
                      ))}
                    </div>
                  )}
                </div>
              )}
            </div>
          ))}

          {loading && (
            <div className="flex justify-start">
              <div className="bg-[#1A1D23] border border-white/[0.08] rounded-2xl rounded-tl-sm p-4 flex items-center gap-3">
                <Loader2 size={16} className="animate-spin text-[#76B900]" />
                <span className="text-sm text-[#9CA3AF]">Searching evidence...</span>
                <div className="flex gap-1">
                  <div className="w-1.5 h-1.5 rounded-full bg-[#76B900]/40 animate-pulse" />
                  <div className="w-1.5 h-1.5 rounded-full bg-[#76B900]/40 animate-pulse" style={{ animationDelay: '0.2s' }} />
                  <div className="w-1.5 h-1.5 rounded-full bg-[#76B900]/40 animate-pulse" style={{ animationDelay: '0.4s' }} />
                </div>
              </div>
            </div>
          )}
          <div ref={messagesEndRef} />
        </div>

        {/* Input */}
        <div className="border-t border-white/[0.08] pt-4">
          <form
            onSubmit={(e) => {
              e.preventDefault();
              handleAsk();
            }}
            className="flex gap-3"
          >
            <input
              type="text"
              value={input}
              onChange={(e) => setInput(e.target.value)}
              placeholder="Ask a clinical imaging question..."
              className="flex-1 bg-[#1A1D23] border border-white/[0.08] rounded-xl px-4 py-3 text-sm text-white placeholder-[#9CA3AF]/60 focus:outline-none focus:border-[#76B900]/50 focus:ring-1 focus:ring-[#76B900]/20 transition-all duration-200"
            />
            <button
              type="submit"
              disabled={loading || !input.trim()}
              className="bg-[#76B900] hover:bg-[#5a9400] disabled:opacity-40 text-black p-3 rounded-xl transition-all duration-200 cursor-pointer shadow-md shadow-[#76B900]/20"
            >
              <Send size={18} />
            </button>
          </form>
        </div>
      </div>

      {/* Citations Sidebar */}
      <div className="w-80 bg-[#1A1D23] rounded-xl border border-white/[0.08] overflow-auto hidden lg:flex lg:flex-col shadow-lg shadow-black/20">
        <div className="p-4 border-b border-white/[0.06] bg-white/[0.01]">
          <h3 className="text-sm font-semibold text-white flex items-center gap-2">
            <BookOpen size={16} className="text-[#76B900]" />
            Evidence Sources
          </h3>
          {selectedCitations.length > 0 && (
            <p className="text-[10px] text-[#9CA3AF] mt-1">{selectedCitations.length} citations found</p>
          )}
        </div>
        <div className="p-4 flex-1 overflow-auto">
          {selectedCitations.length === 0 ? (
            <div className="flex flex-col items-center justify-center h-full text-center py-12">
              <BookOpen size={28} className="text-white/10 mb-3" />
              <p className="text-xs text-[#9CA3AF] leading-relaxed">
                Ask a question to see supporting evidence and citations.
              </p>
            </div>
          ) : (
            <div className="space-y-3">
              {selectedCitations.map((c, i) => (
                <div
                  key={i}
                  className="bg-[#0E1117] rounded-xl p-3.5 border border-white/[0.06] hover:border-[#76B900]/20 transition-all duration-200"
                >
                  <div className="flex items-start justify-between gap-2 mb-2">
                    <h4 className="text-xs font-semibold text-white line-clamp-2 leading-relaxed">
                      {c.title || 'Untitled source'}
                    </h4>
                  </div>

                  {/* Relevance bar */}
                  {c.relevance_score != null && (
                    <div className="flex items-center gap-2 mb-2">
                      <div className="flex-1 h-1.5 bg-white/5 rounded-full overflow-hidden">
                        <div
                          className="h-full bg-[#76B900] rounded-full transition-all duration-500"
                          style={{ width: `${c.relevance_score * 100}%` }}
                        />
                      </div>
                      <span className="text-[10px] text-[#76B900] font-mono shrink-0 font-semibold">
                        {(c.relevance_score * 100).toFixed(0)}%
                      </span>
                    </div>
                  )}

                  {c.snippet && (
                    <p className="text-[11px] text-[#9CA3AF] line-clamp-3 mb-2 leading-relaxed">
                      {c.snippet}
                    </p>
                  )}
                  <div className="flex items-center gap-2">
                    {c.collection && (
                      <span className="text-[10px] bg-purple-500/10 text-purple-400 border border-purple-500/20 px-2 py-0.5 rounded-full">
                        {c.collection}
                      </span>
                    )}
                    {c.metadata?.year && (
                      <span className="text-[10px] text-[#9CA3AF]">
                        {c.metadata.year}
                      </span>
                    )}
                    {c.pubmed_id && (
                      <a
                        href={`https://pubmed.ncbi.nlm.nih.gov/${c.pubmed_id}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-[10px] text-blue-400 hover:underline"
                      >
                        PubMed
                      </a>
                    )}
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
