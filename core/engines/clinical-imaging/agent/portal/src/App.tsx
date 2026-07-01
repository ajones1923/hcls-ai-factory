import { BrowserRouter, Routes, Route } from 'react-router-dom';
import Layout from './components/Layout';
import PageTransition from './components/PageTransition';
import Dashboard from './pages/Dashboard';
import Workflows from './pages/Workflows';
import Evidence from './pages/Evidence';
import Protocol from './pages/Protocol';
import DoseTracking from './pages/DoseTracking';
import Analytics from './pages/Analytics';
import Reports from './pages/Reports';
import Benchmarks from './pages/Benchmarks';
import Compare from './pages/Compare';
import LiveAnalysis from './pages/LiveAnalysis';

function App() {
  return (
    <BrowserRouter>
      <Layout>
        <PageTransition>
          <Routes>
            <Route path="/" element={<Dashboard />} />
            <Route path="/workflows" element={<Workflows />} />
            <Route path="/live-analysis" element={<LiveAnalysis />} />
            <Route path="/evidence" element={<Evidence />} />
            <Route path="/protocol" element={<Protocol />} />
            <Route path="/dose" element={<DoseTracking />} />
            <Route path="/analytics" element={<Analytics />} />
            <Route path="/reports" element={<Reports />} />
            <Route path="/benchmarks" element={<Benchmarks />} />
            <Route path="/compare" element={<Compare />} />
          </Routes>
        </PageTransition>
      </Layout>
    </BrowserRouter>
  );
}

export default App;
