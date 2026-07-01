import os
import sys
from pathlib import Path

# Tests NEVER make real API calls (deterministic + free), even with a key in .env.
os.environ.setdefault("TSC_OFFLINE", "1")

# engine root on path so `config` and `src` import as top-level packages
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
