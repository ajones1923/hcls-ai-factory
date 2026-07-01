"""Demo case manager for pre-loaded clinical scenarios."""
import json
from pathlib import Path
from typing import List, Dict, Any, Optional
from src.models import DemoCaseDefinition

DEMO_CASES_PATH = Path(__file__).parent.parent / "data" / "reference" / "demo_cases.json"


def load_demo_cases() -> List[DemoCaseDefinition]:
    """Load all demo case definitions from JSON."""
    with open(DEMO_CASES_PATH) as f:
        raw = json.load(f)
    return [DemoCaseDefinition(**case) for case in raw]


def get_demo_case(case_id: str) -> Optional[DemoCaseDefinition]:
    """Get a specific demo case by ID."""
    for case in load_demo_cases():
        if case.case_id == case_id:
            return case
    return None


def run_demo_case(case_id: str, workflows: dict) -> Dict[str, Any]:
    """Run a demo case through its workflow with mock overrides.

    Args:
        case_id: Demo case identifier (DEMO-001, DEMO-002, DEMO-003)
        workflows: Dict of workflow_name -> workflow_instance (from WORKFLOW_REGISTRY)

    Returns:
        Dict with workflow_result, genomic_context, case_definition, and narrative
    """
    case = get_demo_case(case_id)
    if not case:
        raise ValueError(f"Demo case {case_id} not found")

    workflow_cls = workflows.get(case.workflow_name)
    if not workflow_cls:
        raise ValueError(f"Workflow {case.workflow_name} not registered")

    # Create workflow instance with mock_overrides
    workflow = workflow_cls(mock_overrides=case.workflow_overrides)
    result = workflow.run()

    case_data = case.model_dump() if hasattr(case, 'model_dump') else case.dict()
    result_data = result.model_dump() if hasattr(result, 'model_dump') else result
    return {
        "case": case_data,
        "workflow_result": result_data,
        "genomic_context": case.genomic_context,
        "talking_points": case_data.get("talking_points", []),
    }
