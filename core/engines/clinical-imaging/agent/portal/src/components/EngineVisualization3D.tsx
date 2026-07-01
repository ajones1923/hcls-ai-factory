import { useRef, useEffect, useCallback } from 'react';
import * as THREE from 'three';

interface EngineConfig {
  id: number;
  name: string;
  title: string;
  color: string;
  position: [number, number, number];
  active: boolean;
}

const ENGINES: EngineConfig[] = [
  { id: 1, name: 'Engine 1', title: 'Genomics',         color: '#00BCD4', position: [-2.0,  1.0,  0.0], active: false },
  { id: 2, name: 'Engine 2', title: 'Intelligence',     color: '#2196F3', position: [-2.0, -1.0,  0.0], active: false },
  { id: 3, name: 'Engine 3', title: 'Drug Discovery',   color: '#9C27B0', position: [ 2.0, -1.0,  0.0], active: false },
  { id: 4, name: 'Engine 4', title: 'Clinical Imaging',  color: '#76B900', position: [ 0.0, -2.2,  0.3], active: true },
  { id: 5, name: 'Engine 5', title: 'Orchestration',    color: '#FF9800', position: [ 0.0,  2.2,  0.0], active: false },
  { id: 6, name: 'Engine 6', title: 'Monitoring',       color: '#607D8B', position: [ 2.0,  1.0,  0.0], active: false },
];

interface StreamConfig {
  source: number;
  target: number;
  color: string;
  particleCount: number;
  speed: number;
  alwaysOn: boolean;
}

const STREAMS: StreamConfig[] = [
  { source: 4, target: 2, color: '#9C27B0', particleCount: 70, speed: 2.5, alwaysOn: true },
  { source: 2, target: 3, color: '#2196F3', particleCount: 50, speed: 3.0, alwaysOn: false },
  { source: 1, target: 2, color: '#00BCD4', particleCount: 45, speed: 2.8, alwaysOn: false },
  { source: 5, target: 4, color: '#FF9800', particleCount: 35, speed: 2.0, alwaysOn: true },
  { source: 4, target: 6, color: '#607D8B', particleCount: 25, speed: 4.0, alwaysOn: true },
  { source: 3, target: 4, color: '#76B900', particleCount: 40, speed: 3.0, alwaysOn: false },
  { source: 1, target: 4, color: '#00BCD4', particleCount: 30, speed: 3.2, alwaysOn: false },
];

function createCircleTexture(): THREE.Texture {
  const size = 128;
  const canvas = document.createElement('canvas');
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext('2d')!;
  const center = size / 2;
  const gradient = ctx.createRadialGradient(center, center, 0, center, center, center);
  gradient.addColorStop(0, 'rgba(255,255,255,1)');
  gradient.addColorStop(0.2, 'rgba(255,255,255,0.85)');
  gradient.addColorStop(0.5, 'rgba(255,255,255,0.4)');
  gradient.addColorStop(0.8, 'rgba(255,255,255,0.1)');
  gradient.addColorStop(1, 'rgba(255,255,255,0)');
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  return new THREE.CanvasTexture(canvas);
}

// Glow shader
const glowVertexShader = `
  varying vec3 vNormal;
  void main() {
    vNormal = normalize(normalMatrix * normal);
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
  }
`;

const glowFragmentShader = `
  uniform vec3 glowColor;
  uniform float intensity;
  varying vec3 vNormal;
  void main() {
    float glow = pow(0.6 - dot(vNormal, vec3(0.0, 0.0, 1.0)), 2.0);
    gl_FragColor = vec4(glowColor, glow * intensity);
  }
`;

function computeBezierControl(
  source: THREE.Vector3,
  target: THREE.Vector3,
  index: number
): THREE.Vector3 {
  const mid = new THREE.Vector3().addVectors(source, target).multiplyScalar(0.5);
  const dir = new THREE.Vector3().subVectors(target, source).normalize();
  const perp = new THREE.Vector3(-dir.y, dir.x, 0);
  const sign = index % 2 === 0 ? 1 : -1;
  mid.add(perp.multiplyScalar(0.4 * sign));
  mid.z += 0.6;
  return mid;
}

function sampleBezier(p0: THREE.Vector3, p1: THREE.Vector3, p2: THREE.Vector3, t: number): THREE.Vector3 {
  const omt = 1 - t;
  return new THREE.Vector3(
    omt * omt * p0.x + 2 * omt * t * p1.x + t * t * p2.x,
    omt * omt * p0.y + 2 * omt * t * p1.y + t * t * p2.y,
    omt * omt * p0.z + 2 * omt * t * p1.z + t * t * p2.z,
  );
}

interface EngineVisualization3DProps {
  height?: number;
  className?: string;
}

export default function EngineVisualization3D({ height = 420, className }: EngineVisualization3DProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const frameRef = useRef<number>(0);
  const clockRef = useRef(new THREE.Clock());

  const setupScene = useCallback(() => {
    const container = containerRef.current;
    if (!container) return;

    const width = container.clientWidth || 800;
    const h = height;

    // Scene
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0a0d14);

    // Camera
    const camera = new THREE.PerspectiveCamera(50, width / h, 0.1, 100);
    camera.position.set(0, -0.5, 7);
    camera.lookAt(0, 0, 0);

    // Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, h);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Lighting
    scene.add(new THREE.AmbientLight(0x76B900, 0.2));
    const keyLight = new THREE.DirectionalLight(0xffffff, 0.6);
    keyLight.position.set(3, 4, 5);
    scene.add(keyLight);

    // Green accent light on Engine 4
    const accentLight = new THREE.PointLight(0x76B900, 0.5, 8);
    accentLight.position.set(0, -2.2, 1.5);
    scene.add(accentLight);

    // Subtle grid
    const grid = new THREE.GridHelper(16, 32, 0x76B900, 0x76B900);
    grid.rotation.x = Math.PI / 2;
    grid.position.z = -1.2;
    (grid.material as THREE.Material).transparent = true;
    (grid.material as THREE.Material).opacity = 0.03;
    scene.add(grid);

    // Main group for rotation
    const mainGroup = new THREE.Group();
    scene.add(mainGroup);

    // Sprite texture for particles
    const sprite = createCircleTexture();

    // ── Engine Nodes ──────────────────────────────────────────
    const engineGroups: THREE.Group[] = [];
    const glowMaterials: THREE.ShaderMaterial[] = [];

    ENGINES.forEach((eng) => {
      const group = new THREE.Group();
      group.position.set(...eng.position);

      // Glow sphere (outer)
      const glowGeo = new THREE.SphereGeometry(eng.active ? 0.7 : 0.55, 32, 32);
      const glowMat = new THREE.ShaderMaterial({
        uniforms: {
          glowColor: { value: new THREE.Color(eng.color) },
          intensity: { value: eng.active ? 0.7 : 0.35 },
        },
        vertexShader: glowVertexShader,
        fragmentShader: glowFragmentShader,
        side: THREE.BackSide,
        blending: THREE.AdditiveBlending,
        transparent: true,
        depthWrite: false,
      });
      glowMaterials.push(glowMat);
      group.add(new THREE.Mesh(glowGeo, glowMat));

      // Core sphere
      const coreGeo = new THREE.SphereGeometry(eng.active ? 0.38 : 0.3, 32, 32);
      const coreMat = new THREE.MeshPhongMaterial({
        color: new THREE.Color(eng.color),
        emissive: new THREE.Color(eng.color),
        emissiveIntensity: eng.active ? 0.5 : 0.15,
        shininess: 90,
        transparent: true,
        opacity: 0.92,
      });
      group.add(new THREE.Mesh(coreGeo, coreMat));

      // Orbiting ring for active engine
      if (eng.active) {
        const ringGeo = new THREE.RingGeometry(0.52, 0.56, 64);
        const ringMat = new THREE.MeshBasicMaterial({
          color: new THREE.Color(eng.color),
          transparent: true,
          opacity: 0.4,
          side: THREE.DoubleSide,
        });
        const ring = new THREE.Mesh(ringGeo, ringMat);
        ring.rotation.x = Math.PI / 3;
        ring.name = 'orbit-ring';
        group.add(ring);
      }

      mainGroup.add(group);
      engineGroups.push(group);
    });

    // ── Connection lines (faint) ──────────────────────────────
    STREAMS.forEach((stream, idx) => {
      const srcEng = ENGINES.find(e => e.id === stream.source)!;
      const tgtEng = ENGINES.find(e => e.id === stream.target)!;
      const src = new THREE.Vector3(...srcEng.position);
      const tgt = new THREE.Vector3(...tgtEng.position);
      const ctrl = computeBezierControl(src, tgt, idx);

      const curvePoints: THREE.Vector3[] = [];
      for (let t = 0; t <= 1; t += 0.02) {
        curvePoints.push(sampleBezier(src, ctrl, tgt, t));
      }

      const lineGeo = new THREE.BufferGeometry().setFromPoints(curvePoints);
      const lineMat = new THREE.LineBasicMaterial({
        color: new THREE.Color(stream.color),
        transparent: true,
        opacity: 0.08,
      });
      mainGroup.add(new THREE.Line(lineGeo, lineMat));
    });

    // ── Particle Streams ──────────────────────────────────────
    interface StreamState {
      tValues: Float32Array;
      speeds: Float32Array;
      positions: Float32Array;
      geometry: THREE.BufferGeometry;
      source: THREE.Vector3;
      control: THREE.Vector3;
      target: THREE.Vector3;
      baseSpeed: number;
      alwaysOn: boolean;
    }

    const streamStates: StreamState[] = [];

    STREAMS.forEach((stream, idx) => {
      const srcEng = ENGINES.find(e => e.id === stream.source)!;
      const tgtEng = ENGINES.find(e => e.id === stream.target)!;
      const src = new THREE.Vector3(...srcEng.position);
      const tgt = new THREE.Vector3(...tgtEng.position);
      const ctrl = computeBezierControl(src, tgt, idx);

      const count = stream.particleCount;
      const positions = new Float32Array(count * 3);
      const tValues = new Float32Array(count);
      const speeds = new Float32Array(count);

      for (let i = 0; i < count; i++) {
        tValues[i] = Math.random();
        speeds[i] = 0.75 + Math.random() * 0.5;
        const pos = sampleBezier(src, ctrl, tgt, tValues[i]);
        positions[i * 3] = pos.x;
        positions[i * 3 + 1] = pos.y;
        positions[i * 3 + 2] = pos.z;
      }

      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));

      const material = new THREE.PointsMaterial({
        size: 0.05,
        color: new THREE.Color(stream.color),
        transparent: true,
        opacity: stream.alwaysOn ? 0.7 : 0.5,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        sizeAttenuation: true,
        map: sprite,
      });

      const points = new THREE.Points(geometry, material);
      mainGroup.add(points);

      streamStates.push({
        tValues, speeds, positions, geometry,
        source: src, control: ctrl, target: tgt,
        baseSpeed: stream.speed,
        alwaysOn: stream.alwaysOn,
      });
    });

    // ── HTML Labels (via CSS) ─────────────────────────────────
    // We'll render engine labels as absolutely positioned divs
    const labelContainer = document.createElement('div');
    labelContainer.style.cssText = 'position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none;overflow:hidden;';
    container.style.position = 'relative';
    container.appendChild(labelContainer);

    const labels: HTMLDivElement[] = [];
    ENGINES.forEach((eng) => {
      const label = document.createElement('div');
      label.style.cssText = `
        position: absolute;
        text-align: center;
        pointer-events: none;
        user-select: none;
        transform: translate(-50%, -50%);
        white-space: nowrap;
        transition: opacity 0.3s;
      `;
      label.innerHTML = `
        <div style="
          background: rgba(14,17,23,0.85);
          border: 1px solid rgba(255,255,255,0.08);
          border-radius: 8px;
          padding: 5px 10px;
          backdrop-filter: blur(6px);
        ">
          <div style="font-size:11px;font-weight:700;color:${eng.color};letter-spacing:0.3px;">${eng.title}</div>
          <div style="font-size:9px;color:rgba(255,255,255,0.4);margin-top:1px;">${eng.name}</div>
          ${eng.active ? '<div style="font-size:8px;font-weight:700;color:#76B900;margin-top:2px;letter-spacing:1px;">● ACTIVE</div>' : ''}
        </div>
      `;
      labelContainer.appendChild(label);
      labels.push(label);
    });

    // ── Mouse Interaction ─────────────────────────────────────
    let isDragging = false;
    let prevX = 0;
    let prevY = 0;
    let lastInteraction = Date.now();

    const onMouseDown = (e: MouseEvent) => {
      isDragging = true;
      prevX = e.clientX;
      prevY = e.clientY;
      lastInteraction = Date.now();
    };
    const onMouseMove = (e: MouseEvent) => {
      if (!isDragging) return;
      const dx = e.clientX - prevX;
      const dy = e.clientY - prevY;
      mainGroup.rotation.y += dx * 0.008;
      mainGroup.rotation.x += dy * 0.005;
      mainGroup.rotation.x = Math.max(-0.5, Math.min(0.5, mainGroup.rotation.x));
      prevX = e.clientX;
      prevY = e.clientY;
      lastInteraction = Date.now();
    };
    const onMouseUp = () => { isDragging = false; };

    const onWheel = (e: WheelEvent) => {
      e.preventDefault();
      camera.position.z = Math.max(4, Math.min(12, camera.position.z + e.deltaY * 0.005));
      lastInteraction = Date.now();
    };

    renderer.domElement.addEventListener('mousedown', onMouseDown);
    window.addEventListener('mousemove', onMouseMove);
    window.addEventListener('mouseup', onMouseUp);
    renderer.domElement.addEventListener('wheel', onWheel, { passive: false });

    // ── Animation Loop ────────────────────────────────────────
    let elapsedTotal = 0;

    const animate = () => {
      frameRef.current = requestAnimationFrame(animate);
      const delta = clockRef.current.getDelta();
      elapsedTotal += delta;

      // Auto-rotate when idle
      const idleTime = (Date.now() - lastInteraction) / 1000;
      if (!isDragging && idleTime > 3) {
        const rotSpeed = Math.min(0.0015, 0.0003 * idleTime);
        mainGroup.rotation.y += rotSpeed;
      }

      // Pulse Engine 4 glow
      const engine4Glow = glowMaterials[3]; // Engine 4 is index 3
      engine4Glow.uniforms.intensity.value = 0.55 + Math.sin(elapsedTotal * 2.5) * 0.2;

      // Rotate orbit ring on Engine 4
      const ring = engineGroups[3].getObjectByName('orbit-ring');
      if (ring) {
        ring.rotation.z = elapsedTotal * 0.8;
      }

      // Subtle float for all engines
      ENGINES.forEach((eng, i) => {
        const group = engineGroups[i];
        group.position.y = eng.position[1] + Math.sin(elapsedTotal * 0.5 + i * 1.2) * 0.04;
      });

      // Update particle streams
      streamStates.forEach((state) => {
        const count = state.tValues.length;
        const speedFactor = delta / state.baseSpeed;
        const activityMult = state.alwaysOn ? 1.0 : 0.6 + Math.sin(elapsedTotal * 0.3) * 0.4;

        for (let i = 0; i < count; i++) {
          state.tValues[i] += speedFactor * state.speeds[i] * activityMult;
          if (state.tValues[i] > 1.0) state.tValues[i] -= 1.0;

          const pos = sampleBezier(state.source, state.control, state.target, state.tValues[i]);
          state.positions[i * 3] = pos.x;
          state.positions[i * 3 + 1] = pos.y;
          state.positions[i * 3 + 2] = pos.z;
        }
        state.geometry.attributes.position.needsUpdate = true;
      });

      // Update label positions (project 3D to 2D)
      ENGINES.forEach((eng, i) => {
        const group = engineGroups[i];
        const worldPos = new THREE.Vector3();
        group.getWorldPosition(worldPos);
        worldPos.project(camera);

        const x = (worldPos.x * 0.5 + 0.5) * width;
        const y = (-worldPos.y * 0.5 + 0.5) * h;
        const label = labels[i];
        label.style.left = `${x}px`;
        label.style.top = `${y + 45}px`;

        // Fade labels when behind camera
        label.style.opacity = worldPos.z < 1 ? '1' : '0';
      });

      renderer.render(scene, camera);
    };

    animate();

    // Resize handler
    const onResize = () => {
      const w = container.clientWidth || 800;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    window.addEventListener('resize', onResize);

    // Cleanup
    return () => {
      cancelAnimationFrame(frameRef.current);
      renderer.domElement.removeEventListener('mousedown', onMouseDown);
      window.removeEventListener('mousemove', onMouseMove);
      window.removeEventListener('mouseup', onMouseUp);
      renderer.domElement.removeEventListener('wheel', onWheel);
      window.removeEventListener('resize', onResize);

      // Dispose
      renderer.dispose();
      sprite.dispose();
      engineGroups.forEach(g => {
        g.traverse(child => {
          if ((child as THREE.Mesh).geometry) (child as THREE.Mesh).geometry.dispose();
          if ((child as THREE.Mesh).material) {
            const mat = (child as THREE.Mesh).material;
            if (Array.isArray(mat)) mat.forEach(m => m.dispose());
            else (mat as THREE.Material).dispose();
          }
        });
      });
      streamStates.forEach(s => s.geometry.dispose());

      if (container.contains(renderer.domElement)) container.removeChild(renderer.domElement);
      if (container.contains(labelContainer)) container.removeChild(labelContainer);
    };
  }, [height]);

  useEffect(() => {
    const cleanup = setupScene();
    return cleanup;
  }, [setupScene]);

  return (
    <div ref={containerRef} className={`relative rounded-xl overflow-hidden border border-white/[0.06] ${className || ''}`} style={{ height }}>
      {/* Legend overlay */}
      <div className="absolute bottom-3 left-3 flex flex-wrap gap-2 pointer-events-none z-10">
        {ENGINES.map(eng => (
          <div key={eng.id} className="flex items-center gap-1.5 bg-[#0E1117]/80 backdrop-blur-sm rounded-full px-2.5 py-1 border border-white/[0.06]">
            <div className="w-2 h-2 rounded-full" style={{ backgroundColor: eng.color, boxShadow: `0 0 4px ${eng.color}` }} />
            <span className="text-[9px] text-white/50 font-medium">{eng.title}</span>
          </div>
        ))}
      </div>
      {/* Interaction hint */}
      <div className="absolute top-3 right-3 pointer-events-none z-10">
        <span className="text-[9px] text-white/20 bg-white/[0.03] px-2 py-1 rounded-full">
          Drag to rotate • Scroll to zoom
        </span>
      </div>
    </div>
  );
}
