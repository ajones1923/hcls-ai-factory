import { create } from 'zustand';

interface HealthData {
  status: string;
  total_vectors: number;
  collections: Record<string, number>;
  nim_services: Record<string, string>;
}

interface AppState {
  health: HealthData | null;
  setHealth: (h: HealthData) => void;
  activeTab: string;
  setActiveTab: (t: string) => void;
  lastWorkflowResult: Record<string, unknown> | null;
  setLastWorkflowResult: (r: Record<string, unknown> | null) => void;
  lastQueryResult: Record<string, unknown> | null;
  setLastQueryResult: (r: Record<string, unknown> | null) => void;
}

export const useAppStore = create<AppState>((set) => ({
  health: null,
  setHealth: (h) => set({ health: h }),
  activeTab: 'dashboard',
  setActiveTab: (t) => set({ activeTab: t }),
  lastWorkflowResult: null,
  setLastWorkflowResult: (r) => set({ lastWorkflowResult: r }),
  lastQueryResult: null,
  setLastQueryResult: (r) => set({ lastQueryResult: r }),
}));
