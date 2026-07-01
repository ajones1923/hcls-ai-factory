#!/usr/bin/env python3
"""Test NIM service connectivity for the Imaging Intelligence Agent.

Instantiates the NIMServiceManager and checks the health of all four
NIM services: VISTA-3D, MAISI, VILA-M3, and Llama3 LLM.

Usage:
    python scripts/test_nim_services.py
    python scripts/test_nim_services.py --mock

Options:
    --mock    Force mock mode for all NIM services

Author: Adam Jones
Date: February 2026
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

from config.settings import ImagingSettings, settings
from src.nim.service_manager import NIMServiceManager


def main():
    parser = argparse.ArgumentParser(
        description="Test NIM service connectivity"
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Force mock mode for all NIM services",
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  Imaging Intelligence Agent — NIM Service Test")
    print("=" * 65)
    print(f"  Mock mode: {args.mock}")
    print()

    # Optionally override settings for mock mode
    if args.mock:
        test_settings = ImagingSettings(
            NIM_MODE="mock",
            NIM_ALLOW_MOCK_FALLBACK=True,
        )
    else:
        test_settings = settings

    logger.info("Initializing NIM Service Manager...")
    nim_manager = NIMServiceManager(test_settings)

    # Check all services
    logger.info("Checking all NIM services...")
    status = nim_manager.check_all_services()

    # Report results
    print()
    print("  NIM Service Status:")
    print("  " + "-" * 40)

    status_symbols = {
        "available": "[OK]",
        "mock": "[MOCK]",
        "unavailable": "[FAIL]",
    }

    service_names = {
        "vista3d": "VISTA-3D (Segmentation)",
        "maisi": "MAISI (Synthetic CT)",
        "vila_m3": "VILA-M3 (Vision LLM)",
        "llm": "Llama3 / LLM",
    }

    all_ok = True
    for service_key, service_status in status.items():
        symbol = status_symbols.get(service_status, "[?]")
        name = service_names.get(service_key, service_key)
        print(f"  {symbol:8s} {name}: {service_status}")
        if service_status == "unavailable" and not args.mock:
            all_ok = False

    print("  " + "-" * 40)

    # Summary
    available = sum(1 for s in status.values() if s == "available")
    mock = sum(1 for s in status.values() if s == "mock")
    unavailable = sum(1 for s in status.values() if s == "unavailable")

    print()
    print(f"  Available: {available}  |  Mock: {mock}  |  Unavailable: {unavailable}")

    if all_ok or args.mock:
        print("  Overall: PASS")
    else:
        print("  Overall: PARTIAL (some services unavailable)")

    # List available services
    available_services = nim_manager.get_available_services()
    if available_services:
        logger.info(f"Available services: {available_services}")

    print()
    print("=" * 65)


if __name__ == "__main__":
    main()
