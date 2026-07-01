import axios from 'axios';

const API_BASE = import.meta.env.VITE_API_URL || `${window.location.protocol}//${window.location.hostname}:8524`;
const api = axios.create({ baseURL: API_BASE });

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
