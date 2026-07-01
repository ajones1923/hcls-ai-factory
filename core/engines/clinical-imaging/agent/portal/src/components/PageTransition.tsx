import { useLocation } from 'react-router-dom';
import { useEffect, useState } from 'react';

export default function PageTransition({ children }: { children: React.ReactNode }) {
  const location = useLocation();
  const [show, setShow] = useState(false);

  useEffect(() => {
    setShow(false);
    const timer = requestAnimationFrame(() => {
      requestAnimationFrame(() => setShow(true));
    });
    return () => cancelAnimationFrame(timer);
  }, [location.pathname]);

  return (
    <div className={`transition-all duration-300 ease-out ${
      show ? 'opacity-100 translate-y-0' : 'opacity-0 translate-y-2'
    }`}>
      {children}
    </div>
  );
}
