import axios from 'axios';

const API_BASE = import.meta.env.VITE_API_URL || `${window.location.protocol}//${window.location.hostname}:8524`;
const api = axios.create({ baseURL: API_BASE });

/**
 * Turn a failed API call into something a human can act on.
 *
 * Every loader in this portal used to end in `.catch(() => {})`. That made a CORS block, a stopped
 * API and an empty database look identical: the page simply sat on "Loading..." with nothing in the
 * console. A one-line origin mismatch (vite fell back from :3001 to :3002, which was not in
 * CORS_ORIGINS) cost real debugging time for exactly that reason.
 *
 * Returns a short message suitable for display, and always logs the detail.
 */
export function describeApiError(what: string, err: unknown): string {
  const e = err as { message?: string; code?: string; response?: { status?: number } };
  const status = e?.response?.status;
  let msg: string;
  if (status) {
    msg = `${what} failed: HTTP ${status} from ${API_BASE}`;
  } else if (e?.code === 'ERR_NETWORK' || /Network Error/i.test(e?.message ?? '')) {
    // No status at all means the browser never got a usable response: the API is down, or it
    // answered without an Access-Control-Allow-Origin header for this exact origin.
    msg = `${what} failed: cannot reach ${API_BASE} — API down, or this origin `
        + `(${window.location.origin}) is not in the API's CORS_ORIGINS.`;
  } else {
    msg = `${what} failed: ${e?.message ?? 'unknown error'}`;
  }
  console.error(`[imaging-portal] ${msg}`, err);
  return msg;
}


export const fetchHealth = () => api.get('/health').then(r => r.data);
export const fetchWorkflows = () => api.get('/workflows').then(r => r.data);
export const runWorkflow = (name: string) => api.post(`/workflow/${name}/run`, {}).then(r => r.data);
export const askQuestion = (question: string) => api.post('/api/ask', { question }).then(r => r.data);
export const fetchNIMStatus = () => api.get('/nim/status').then(r => r.data);
export const fetchDemoCases = () => api.get('/demo-cases').then(r => r.data);
export const runDemoCase = (id: string) => api.post(`/demo-cases/${id}/run`).then(r => r.data);
export const fetchCollections = () => api.get('/collections').then(r => r.data);
export const fetchKnowledgeStats = () => api.get('/knowledge/stats').then(r => r.data);
export const recommendProtocol = (data: Record<string, unknown>) => api.post('/protocol/recommend', data).then(r => r.data);
export const fetchIndications = () => api.get('/protocol/indications').then(r => r.data);
export const recordDose = (data: Record<string, unknown>) => api.post('/dose/record', data).then(r => r.data);
export const fetchCumulativeDose = (patientId: string) => api.get(`/dose/cumulative/${patientId}`).then(r => r.data);
export const compareDRL = (data: Record<string, unknown>) => api.post('/dose/compare-drl', data).then(r => r.data);
export const generateReport = (data: Record<string, unknown>) => api.post('/reports/generate', data).then(r => r.data);
export const fetchAnalyticsPopulation = () => api.get('/api/analytics/population').then(r => r.data);
export const generateDemoData = (n: number) => api.post(`/api/analytics/generate-demo-data?n_studies=${n}`).then(r => r.data);
export const fetchTrends = (metric: string) => api.get(`/api/analytics/trends/${metric}`).then(r => r.data);
