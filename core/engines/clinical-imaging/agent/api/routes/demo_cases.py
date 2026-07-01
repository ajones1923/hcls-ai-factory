"""Demo case endpoints for pre-loaded clinical scenarios."""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional

router = APIRouter(prefix="/demo-cases", tags=["Demo Cases"])

class DemoCaseListItem(BaseModel):
    case_id: str
    title: str
    clinical_scenario: str
    workflow_name: str
    expected_severity: str

class DemoCaseRunResponse(BaseModel):
    case_id: str
    title: str
    clinical_scenario: str
    patient_demographics: Dict[str, Any]
    workflow_name: str
    workflow_result: Dict[str, Any]
    genomic_context: Optional[Dict[str, Any]] = None
    talking_points: List[str] = []

@router.get("", response_model=List[DemoCaseListItem])
def list_demo_cases():
    """List all available demo cases."""
    from src.demo_cases import load_demo_cases
    cases = load_demo_cases()
    return [
        DemoCaseListItem(
            case_id=c.case_id,
            title=c.title,
            clinical_scenario=c.clinical_scenario,
            workflow_name=c.workflow_name,
            expected_severity=c.expected_severity,
        )
        for c in cases
    ]

@router.post("/{case_id}/run", response_model=DemoCaseRunResponse)
def run_demo_case(case_id: str):
    """Run a demo case through its workflow with mock data."""
    from src.demo_cases import get_demo_case

    case = get_demo_case(case_id)
    if not case:
        raise HTTPException(status_code=404, detail=f"Demo case {case_id} not found")

    # Import workflows - use mock mode
    from src.workflows import WORKFLOW_REGISTRY

    workflow_cls = WORKFLOW_REGISTRY.get(case.workflow_name)
    if not workflow_cls:
        raise HTTPException(status_code=500, detail=f"Workflow {case.workflow_name} not registered")

    try:
        workflow = workflow_cls(mock_mode=True, mock_overrides=case.workflow_overrides)
        result = workflow.run()
        result_dict = result.dict() if hasattr(result, 'dict') else (result.model_dump() if hasattr(result, 'model_dump') else result)
    except Exception as e:
        result_dict = {"error": str(e), "status": "mock_fallback", "findings": case.workflow_overrides}

    return DemoCaseRunResponse(
        case_id=case.case_id,
        title=case.title,
        clinical_scenario=case.clinical_scenario,
        patient_demographics=case.patient_demographics,
        workflow_name=case.workflow_name,
        workflow_result=result_dict,
        genomic_context=case.genomic_context,
        talking_points=case.dict().get("talking_points", []) if hasattr(case, 'dict') else [],
    )
