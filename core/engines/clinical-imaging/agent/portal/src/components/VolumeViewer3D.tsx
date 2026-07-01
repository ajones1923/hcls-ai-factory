import { useRef, useEffect } from 'react';
import * as THREE from 'three';

interface VolumeViewer3DProps {
  volumeType: 'brain' | 'chest' | 'head';
  className?: string;
  autoRotate?: boolean;
  height?: number;
}

function createCircleTexture(): THREE.Texture {
  const size = 128;
  const canvas = document.createElement('canvas');
  canvas.width = size;
  canvas.height = size;
  const ctx = canvas.getContext('2d')!;
  const center = size / 2;
  const gradient = ctx.createRadialGradient(center, center, 0, center, center, center);
  gradient.addColorStop(0, 'rgba(255,255,255,1)');
  gradient.addColorStop(0.2, 'rgba(255,255,255,0.9)');
  gradient.addColorStop(0.5, 'rgba(255,255,255,0.5)');
  gradient.addColorStop(0.8, 'rgba(255,255,255,0.15)');
  gradient.addColorStop(1, 'rgba(255,255,255,0)');
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, size, size);
  const tex = new THREE.CanvasTexture(canvas);
  return tex;
}

function randomEllipsoidPoint(rx: number, ry: number, rz: number): [number, number, number] {
  const u = Math.random();
  const v = Math.random();
  const theta = 2 * Math.PI * u;
  const phi = Math.acos(2 * v - 1);
  const r = Math.cbrt(Math.random());
  const x = r * rx * Math.sin(phi) * Math.cos(theta);
  const y = r * ry * Math.sin(phi) * Math.sin(theta);
  const z = r * rz * Math.cos(phi);
  return [x, y, z];
}

function randomShellPoint(rx: number, ry: number, rz: number, minR: number): [number, number, number] {
  const u = Math.random();
  const v = Math.random();
  const theta = 2 * Math.PI * u;
  const phi = Math.acos(2 * v - 1);
  const r = minR + (1 - minR) * Math.cbrt(Math.random());
  return [r * rx * Math.sin(phi) * Math.cos(theta), r * ry * Math.sin(phi) * Math.sin(theta), r * rz * Math.cos(phi)];
}

function buildBrainGeometry(): { positions: Float32Array; colors: Float32Array } {
  const positions: number[] = [];
  const colors: number[] = [];

  // Gray matter shell (blue-purple) — denser
  for (let i = 0; i < 2400; i++) {
    const [x, y, z] = randomShellPoint(1.2, 0.9, 1.0, 0.75);
    positions.push(x, y, z);
    colors.push(0.45 + Math.random() * 0.15, 0.2 + Math.random() * 0.1, 0.8 + Math.random() * 0.15);
  }

  // White matter core (yellow) — denser
  for (let i = 0; i < 1000; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.7, 0.5, 0.6);
    positions.push(x, y, z);
    colors.push(0.9 + Math.random() * 0.1, 0.85 + Math.random() * 0.1, 0.3 + Math.random() * 0.15);
  }

  // Hemorrhage cluster (red) — brighter
  for (let i = 0; i < 160; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.25, 0.2, 0.2);
    positions.push(x + 0.5, y + 0.3, z);
    colors.push(1.0, 0.08 + Math.random() * 0.08, 0.05);
  }

  // CSF ventricles (cyan) — denser
  for (let i = 0; i < 400; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.15, 0.3, 0.5);
    positions.push(x - 0.1, y, z);
    colors.push(0.2, 0.85 + Math.random() * 0.15, 0.95 + Math.random() * 0.05);
  }

  return { positions: new Float32Array(positions), colors: new Float32Array(colors) };
}

function buildChestGeometry(): { positions: Float32Array; colors: Float32Array } {
  const positions: number[] = [];
  const colors: number[] = [];

  // Spine (vertical column of bright points)
  for (let i = 0; i < 120; i++) {
    const y = (i / 120) * 1.6 - 0.8;
    const jitter = () => (Math.random() - 0.5) * 0.06;
    positions.push(jitter(), y, -0.85 + jitter());
    colors.push(0.95, 0.93, 0.87);
  }

  // Rib cage (cream/white arcs) — denser, thicker
  for (let rib = 0; rib < 12; rib++) {
    const yOff = (rib - 5.5) * 0.14;
    for (let j = 0; j < 60; j++) {
      const angle = (j / 60) * Math.PI - Math.PI / 2;
      const rx = 1.1 + Math.random() * 0.03;
      const ry = 0.8 + Math.random() * 0.03;
      positions.push(rx * Math.cos(angle), yOff + (Math.random() - 0.5) * 0.02, ry * Math.sin(angle));
      colors.push(0.95, 0.92, 0.85);
      // Add thickness
      if (j % 3 === 0) {
        positions.push(rx * Math.cos(angle) * 0.97, yOff + (Math.random() - 0.5) * 0.02, ry * Math.sin(angle) * 0.97);
        colors.push(0.90, 0.87, 0.80);
      }
    }
  }

  // Left lung (green) — much denser
  for (let i = 0; i < 900; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.4, 0.7, 0.5);
    positions.push(x - 0.55, y, z);
    const intensity = 0.55 + Math.random() * 0.25;
    colors.push(0.15 + Math.random() * 0.1, intensity, 0.15 + Math.random() * 0.1);
  }

  // Right lung (green) — much denser
  for (let i = 0; i < 900; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.4, 0.7, 0.5);
    positions.push(x + 0.55, y, z);
    const intensity = 0.55 + Math.random() * 0.25;
    colors.push(0.15 + Math.random() * 0.1, intensity, 0.15 + Math.random() * 0.1);
  }

  // Lung surface definition — add shell points for visible boundary
  for (let i = 0; i < 300; i++) {
    const [x, y, z] = randomShellPoint(0.4, 0.7, 0.5, 0.88);
    positions.push(x - 0.55, y, z);
    colors.push(0.25, 0.75 + Math.random() * 0.1, 0.3);
  }
  for (let i = 0; i < 300; i++) {
    const [x, y, z] = randomShellPoint(0.4, 0.7, 0.5, 0.88);
    positions.push(x + 0.55, y, z);
    colors.push(0.25, 0.75 + Math.random() * 0.1, 0.3);
  }

  // Heart (red) — denser, more defined
  for (let i = 0; i < 500; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.25, 0.3, 0.25);
    positions.push(x - 0.05, y - 0.05, z);
    colors.push(0.85 + Math.random() * 0.15, 0.1 + Math.random() * 0.08, 0.1 + Math.random() * 0.08);
  }
  // Heart surface shell
  for (let i = 0; i < 200; i++) {
    const [x, y, z] = randomShellPoint(0.25, 0.3, 0.25, 0.85);
    positions.push(x - 0.05, y - 0.05, z);
    colors.push(1.0, 0.2, 0.15);
  }

  // Coronary arteries (bright yellow-orange lines on heart surface)
  // LAD — runs down the front of the heart
  for (let i = 0; i < 40; i++) {
    const t = i / 40;
    const x = -0.05 + Math.sin(t * 0.5) * 0.05;
    const y = -0.05 + 0.28 - t * 0.45;
    const z = 0.22 + Math.sin(t * 2) * 0.03;
    positions.push(x + (Math.random() - 0.5) * 0.02, y, z + (Math.random() - 0.5) * 0.02);
    colors.push(1.0, 0.75 + Math.random() * 0.15, 0.1);
  }
  // LCx — wraps around left side
  for (let i = 0; i < 30; i++) {
    const t = i / 30;
    const angle = t * Math.PI * 0.6;
    const x = -0.05 - 0.22 * Math.sin(angle);
    const y = -0.05 + 0.1 - t * 0.15;
    const z = 0.22 * Math.cos(angle);
    positions.push(x + (Math.random() - 0.5) * 0.02, y, z + (Math.random() - 0.5) * 0.02);
    colors.push(1.0, 0.7 + Math.random() * 0.15, 0.1);
  }
  // RCA — right side
  for (let i = 0; i < 30; i++) {
    const t = i / 30;
    const angle = t * Math.PI * 0.5;
    const x = -0.05 + 0.23 * Math.sin(angle);
    const y = -0.05 + 0.1 - t * 0.2;
    const z = 0.2 * Math.cos(angle);
    positions.push(x + (Math.random() - 0.5) * 0.02, y, z + (Math.random() - 0.5) * 0.02);
    colors.push(1.0, 0.7 + Math.random() * 0.15, 0.1);
  }

  // Aorta (bright red-pink arc from heart upward)
  for (let i = 0; i < 80; i++) {
    const t = i / 80;
    const angle = t * Math.PI * 0.7;
    const x = -0.05 + 0.3 * Math.sin(angle) - 0.15;
    const y = 0.25 + 0.5 * Math.sin(angle);
    const z = -0.1 + t * 0.05;
    positions.push(x + (Math.random() - 0.5) * 0.04, y, z + (Math.random() - 0.5) * 0.04);
    colors.push(0.9 + Math.random() * 0.1, 0.25 + Math.random() * 0.1, 0.3 + Math.random() * 0.1);
  }

  // Nodule in upper right lung (bright orange, pulsing cluster) — denser
  for (let i = 0; i < 80; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.1, 0.08, 0.08);
    positions.push(x + 0.6, y + 0.4, z + 0.1);
    colors.push(1.0, 0.55 + Math.random() * 0.25, 0.05);
  }
  // Nodule highlight ring
  for (let i = 0; i < 30; i++) {
    const [x, y, z] = randomShellPoint(0.12, 0.1, 0.1, 0.9);
    positions.push(x + 0.6, y + 0.4, z + 0.1);
    colors.push(1.0, 0.8, 0.2);
  }

  // Sternum (front midline)
  for (let i = 0; i < 40; i++) {
    const y = (i / 40) * 1.0 - 0.3;
    positions.push((Math.random() - 0.5) * 0.04, y, 0.78 + (Math.random() - 0.5) * 0.03);
    colors.push(0.93, 0.90, 0.83);
  }

  return { positions: new Float32Array(positions), colors: new Float32Array(colors) };
}

function buildHeadGeometry(): { positions: Float32Array; colors: Float32Array } {
  const positions: number[] = [];
  const colors: number[] = [];

  // Skull outline (cream shell) — denser
  for (let i = 0; i < 1600; i++) {
    const [x, y, z] = randomShellPoint(1.0, 1.2, 1.0, 0.92);
    positions.push(x, y, z);
    colors.push(0.93, 0.9, 0.82);
  }

  // Brain tissue (blue/purple inside) — denser
  for (let i = 0; i < 2000; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.85, 1.0, 0.85);
    positions.push(x, y + 0.1, z);
    colors.push(0.4 + Math.random() * 0.15, 0.2 + Math.random() * 0.1, 0.7 + Math.random() * 0.2);
  }

  // Hemorrhage in right hemisphere (bright red) — denser
  for (let i = 0; i < 200; i++) {
    const [x, y, z] = randomEllipsoidPoint(0.2, 0.25, 0.2);
    positions.push(x + 0.4, y + 0.2, z);
    colors.push(1.0, 0.05 + Math.random() * 0.05, 0.05);
  }
  // Hemorrhage highlight ring
  for (let i = 0; i < 60; i++) {
    const [x, y, z] = randomShellPoint(0.22, 0.27, 0.22, 0.9);
    positions.push(x + 0.4, y + 0.2, z);
    colors.push(1.0, 0.3, 0.1);
  }

  return { positions: new Float32Array(positions), colors: new Float32Array(colors) };
}

export default function VolumeViewer3D({ volumeType, className, autoRotate = true, height = 300 }: VolumeViewer3DProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const frameRef = useRef<number>(0);

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const width = container.clientWidth || 280;
    const h = height;

    // Scene setup
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0a0d12);

    const camera = new THREE.PerspectiveCamera(45, width / h, 0.1, 100);
    camera.position.set(0, 0, 3.5);

    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, h);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // Lighting — brighter
    const ambient = new THREE.AmbientLight(0x76b900, 0.5);
    scene.add(ambient);
    const directional = new THREE.DirectionalLight(0xffffff, 0.8);
    directional.position.set(2, 3, 4);
    scene.add(directional);

    // Build geometry
    let data: { positions: Float32Array; colors: Float32Array };
    if (volumeType === 'brain') data = buildBrainGeometry();
    else if (volumeType === 'chest') data = buildChestGeometry();
    else data = buildHeadGeometry();

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(data.positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(data.colors, 3));

    const sprite = createCircleTexture();
    const material = new THREE.PointsMaterial({
      size: 0.055,
      map: sprite,
      vertexColors: true,
      transparent: true,
      opacity: 0.88,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
      sizeAttenuation: true,
    });

    const points = new THREE.Points(geometry, material);
    scene.add(points);

    // Mouse interaction
    let isDragging = false;
    let prevX = 0;
    let prevY = 0;

    const onMouseDown = (e: MouseEvent) => {
      isDragging = true;
      prevX = e.clientX;
      prevY = e.clientY;
    };
    const onMouseMove = (e: MouseEvent) => {
      if (!isDragging) return;
      const dx = e.clientX - prevX;
      const dy = e.clientY - prevY;
      points.rotation.y += dx * 0.01;
      points.rotation.x += dy * 0.01;
      prevX = e.clientX;
      prevY = e.clientY;
    };
    const onMouseUp = () => { isDragging = false; };

    renderer.domElement.addEventListener('mousedown', onMouseDown);
    window.addEventListener('mousemove', onMouseMove);
    window.addEventListener('mouseup', onMouseUp);

    // Animation loop
    const animate = () => {
      frameRef.current = requestAnimationFrame(animate);
      if (autoRotate && !isDragging) {
        points.rotation.y += 0.004;
      }
      renderer.render(scene, camera);
    };
    animate();

    // Resize handler
    const onResize = () => {
      const w = container.clientWidth || 280;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    window.addEventListener('resize', onResize);

    return () => {
      cancelAnimationFrame(frameRef.current);
      renderer.domElement.removeEventListener('mousedown', onMouseDown);
      window.removeEventListener('mousemove', onMouseMove);
      window.removeEventListener('mouseup', onMouseUp);
      window.removeEventListener('resize', onResize);
      renderer.dispose();
      geometry.dispose();
      material.dispose();
      sprite.dispose();
      if (container.contains(renderer.domElement)) {
        container.removeChild(renderer.domElement);
      }
    };
  }, [volumeType, autoRotate, height]);

  return (
    <div ref={containerRef} className={`relative ${className || ''}`} style={{ height }}>
      <span className="absolute bottom-2 left-2 text-[9px] text-[#76B900]/60 font-mono pointer-events-none select-none">
        3D AI Visualization
      </span>
    </div>
  );
}
