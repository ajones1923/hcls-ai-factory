#!/usr/bin/env python3
"""
Start the RAG Chat Streamlit interface.

Usage:
    python scripts/run_chat.py [--port PORT]
"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from config.settings import settings


def main():
    parser = argparse.ArgumentParser(description="Start RAG Chat UI")
    parser.add_argument(
        "--port",
        type=int,
        default=settings.STREAMLIT_PORT,
        help="Port for Streamlit server"
    )
    args = parser.parse_args()

    app_path = Path(__file__).parent.parent / "app" / "chat_ui.py"

    print(f"Starting RAG Chat on port {args.port}...")
    print(f"Open http://localhost:{args.port} in your browser")
    print()

    subprocess.run([
        sys.executable, "-m", "streamlit", "run",
        str(app_path),
        "--server.port", str(args.port),
        "--server.headless", "true",
    ])


if __name__ == "__main__":
    main()
