import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'

export default defineConfig({
  plugins: [react(), tailwindcss()],
  server: {
    host: true,
    port: 3001,
    proxy: {
      '/api': 'http://localhost:8524',
      '/health': 'http://localhost:8524',
      '/workflows': 'http://localhost:8524',
      '/workflow': 'http://localhost:8524',
      '/nim': 'http://localhost:8524',
      '/protocol': 'http://localhost:8524',
      '/dose': 'http://localhost:8524',
      '/demo-cases': 'http://localhost:8524',
      '/collections': 'http://localhost:8524',
      '/knowledge': 'http://localhost:8524',
      '/reports': 'http://localhost:8524',
      '/metrics': 'http://localhost:8524',
      '/events': 'http://localhost:8524',
    }
  }
})
