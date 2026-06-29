"""PDF rendering engine for synthetic clinical documents.

Generates realistic clinical PDFs using ReportLab:
- Progress notes (PCP, specialist)
- Laboratory reports
- Imaging reports
- ER visit notes
- Pathology reports
- Genetic testing reports
- Referral letters
- Medication reconciliation reports

Author: Adam Jones
Date: March 2026
"""

from typing import Dict, List

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    HRFlowable,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

# =====================================================================
# Clinic / Provider Database
# =====================================================================

CLINICS = {
    "pcp": {
        "name": "Greenfield Family Medicine",
        "address": "4200 Wellness Blvd, Suite 101, Austin, TX 78745",
        "phone": "(512) 555-0142",
        "fax": "(512) 555-0143",
    },
    "rheumatology": {
        "name": "Central Texas Rheumatology Associates",
        "address": "1800 Medical Parkway, Suite 300, Austin, TX 78756",
        "phone": "(512) 555-0287",
        "fax": "(512) 555-0288",
    },
    "dermatology": {
        "name": "Austin Dermatology & Skin Center",
        "address": "5500 Balcones Drive, Suite 210, Austin, TX 78731",
        "phone": "(512) 555-0319",
        "fax": "(512) 555-0320",
    },
    "nephrology": {
        "name": "Kidney Care Specialists of Texas",
        "address": "3200 Red River St, Suite 400, Austin, TX 78705",
        "phone": "(512) 555-0456",
        "fax": "(512) 555-0457",
    },
    "cardiology": {
        "name": "Heart of Texas Cardiology",
        "address": "7100 Woodrow Ave, Suite 500, Austin, TX 78757",
        "phone": "(512) 555-0512",
        "fax": "(512) 555-0513",
    },
    "neurology": {
        "name": "Austin Neuroscience Associates",
        "address": "2400 Seton Parkway, Suite 200, Austin, TX 78728",
        "phone": "(512) 555-0634",
        "fax": "(512) 555-0635",
    },
    "ophthalmology": {
        "name": "Lone Star Eye Physicians",
        "address": "6300 Bee Caves Rd, Suite 110, Austin, TX 78746",
        "phone": "(512) 555-0721",
        "fax": "(512) 555-0722",
    },
    "gi": {
        "name": "Capitol Gastroenterology",
        "address": "4500 S Lamar Blvd, Suite 350, Austin, TX 78745",
        "phone": "(512) 555-0845",
        "fax": "(512) 555-0846",
    },
    "orthopedics": {
        "name": "Texas Orthopedic & Sports Medicine",
        "address": "8800 N MoPac Expy, Suite 220, Austin, TX 78759",
        "phone": "(512) 555-0933",
        "fax": "(512) 555-0934",
    },
    "er": {
        "name": "St. David's Medical Center -- Emergency Department",
        "address": "919 E 32nd St, Austin, TX 78705",
        "phone": "(512) 555-0100",
        "fax": "(512) 555-0101",
    },
    "pathology": {
        "name": "Austin Pathology Associates",
        "address": "919 E 32nd St, Austin, TX 78705",
        "phone": "(512) 555-0199",
        "fax": "(512) 555-0200",
    },
    "genetics": {
        "name": "Genomic Health Diagnostics",
        "address": "9500 Arboretum Blvd, Suite 100, Austin, TX 78759",
        "phone": "(512) 555-1050",
        "fax": "(512) 555-1051",
    },
    "lab": {
        "name": "Quest Diagnostics -- Austin Regional Lab",
        "address": "3100 Industrial Terrace, Austin, TX 78758",
        "phone": "(512) 555-0777",
        "fax": "(512) 555-0778",
    },
    "allergy": {
        "name": "Austin Allergy & Immunology Clinic",
        "address": "3600 N Lamar Blvd, Suite 180, Austin, TX 78756",
        "phone": "(512) 555-1122",
        "fax": "(512) 555-1123",
    },
    "dentistry": {
        "name": "Austin Oral Health Partners",
        "address": "2200 S Congress Ave, Suite 105, Austin, TX 78704",
        "phone": "(512) 555-1200",
        "fax": "(512) 555-1201",
    },
    "psychiatry": {
        "name": "Austin Behavioral Health Associates",
        "address": "3800 N Lamar Blvd, Suite 400, Austin, TX 78756",
        "phone": "(512) 555-1300",
        "fax": "(512) 555-1301",
    },
    "pulmonology": {
        "name": "Capital Pulmonary & Critical Care",
        "address": "5600 Bee Caves Rd, Suite 320, Austin, TX 78746",
        "phone": "(512) 555-1410",
        "fax": "(512) 555-1411",
    },
    "urgent_care": {
        "name": "MinuteClinic -- Austin South",
        "address": "4801 S Congress Ave, Austin, TX 78745",
        "phone": "(512) 555-1500",
        "fax": "(512) 555-1501",
    },
    "sleep_medicine": {
        "name": "Austin Sleep & Respiratory Center",
        "address": "7200 N MoPac Expy, Suite 150, Austin, TX 78731",
        "phone": "(512) 555-1600",
        "fax": "(512) 555-1601",
    },
    "pharmacy": {
        "name": "Greenfield Family Medicine -- Pharmacy Services",
        "address": "4200 Wellness Blvd, Suite 101, Austin, TX 78745",
        "phone": "(512) 555-0142",
        "fax": "(512) 555-0143",
    },
    "endocrinology": {
        "name": "Austin Endocrine & Diabetes Specialists",
        "address": "4700 Seton Center Pkwy, Suite 250, Austin, TX 78759",
        "phone": "(512) 555-1700",
        "fax": "(512) 555-1701",
    },
    "nuclear_medicine": {
        "name": "St. David's Nuclear Medicine & Molecular Imaging",
        "address": "919 E 32nd St, Austin, TX 78705",
        "phone": "(512) 555-1800",
        "fax": "(512) 555-1801",
    },
    "dietitian": {
        "name": "Austin Nutrition & Wellness Center",
        "address": "3200 Bee Caves Rd, Suite 140, Austin, TX 78746",
        "phone": "(512) 555-1900",
        "fax": "(512) 555-1901",
    },
}

PROVIDERS = {
    "pcp": "Jennifer Martinez, MD",
    "rheumatology": "Robert Chen, MD, FACR",
    "dermatology": "Amanda Foster, MD, FAAD",
    "nephrology": "Priya Sharma, MD",
    "cardiology": "Michael Thompson, MD, FACC",
    "neurology": "David Kim, MD, PhD",
    "ophthalmology": "Catherine Wells, MD",
    "gi": "James Okafor, MD",
    "orthopedics": "Brian Mitchell, MD",
    "er": "Lisa Nguyen, MD, FACEP",
    "pathology": "Harold Greene, MD",
    "genetics": "Sarah Goldman, MD, FACMG",
    "lab": "Quest Diagnostics",
    "allergy": "Karen Whitfield, MD",
    "dentistry": "Paul Rivera, DDS",
    "psychiatry": "Diane Crawford, MD, FAPA",
    "pulmonology": "Richard Osei, MD, FCCP",
    "urgent_care": "Nurse Practitioner on Duty",
    "sleep_medicine": "Angela Wu, MD, DABSM",
    "pharmacy": "Jennifer Martinez, MD",
    "endocrinology": "Patricia Huang, MD, FACE",
    "nuclear_medicine": "Samuel Beck, MD, ABNM",
    "dietitian": "Emily Carson, RD, LD, CDE",
}


# =====================================================================
# Styles
# =====================================================================

def get_styles():
    """Return custom paragraph styles for clinical documents."""
    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(
        "ClinicName", parent=styles["Normal"],
        fontSize=14, leading=17, textColor=colors.HexColor("#1a3c5e"),
        fontName="Helvetica-Bold", spaceAfter=2,
    ))
    styles.add(ParagraphStyle(
        "ClinicAddr", parent=styles["Normal"],
        fontSize=8, leading=10, textColor=colors.HexColor("#555555"),
        spaceAfter=1,
    ))
    styles.add(ParagraphStyle(
        "DocTitle", parent=styles["Normal"],
        fontSize=13, leading=16, fontName="Helvetica-Bold",
        spaceAfter=8, spaceBefore=6,
        textColor=colors.HexColor("#222222"),
    ))
    styles.add(ParagraphStyle(
        "SectionHead", parent=styles["Normal"],
        fontSize=10, leading=13, fontName="Helvetica-Bold",
        spaceBefore=10, spaceAfter=3,
        textColor=colors.HexColor("#1a3c5e"),
    ))
    styles.add(ParagraphStyle(
        "ClinicalText", parent=styles["Normal"],
        fontSize=10, leading=13, spaceAfter=4,
    ))
    styles.add(ParagraphStyle(
        "SmallText", parent=styles["Normal"],
        fontSize=8, leading=10, textColor=colors.HexColor("#777777"),
    ))
    styles.add(ParagraphStyle(
        "LabValue", parent=styles["Normal"],
        fontSize=9, leading=11,
    ))
    styles.add(ParagraphStyle(
        "AlertText", parent=styles["Normal"],
        fontSize=10, leading=13, fontName="Helvetica-Bold",
        textColor=colors.HexColor("#cc0000"),
    ))
    return styles


# =====================================================================
# Header / Footer Builders
# =====================================================================

def build_header(story, clinic_key: str, styles):
    """Add clinic letterhead to the story."""
    clinic = CLINICS[clinic_key]
    story.append(Paragraph(clinic["name"], styles["ClinicName"]))
    story.append(Paragraph(clinic["address"], styles["ClinicAddr"]))
    story.append(Paragraph(
        f'Phone: {clinic["phone"]}  |  Fax: {clinic["fax"]}',
        styles["ClinicAddr"],
    ))
    story.append(HRFlowable(width="100%", thickness=1.5, color=colors.HexColor("#1a3c5e")))
    story.append(Spacer(1, 8))


def build_patient_banner(story, patient: Dict, date: str, styles):
    """Add patient demographics banner."""
    data = [[
        f"Patient: {patient['name']}",
        f"DOB: {patient['dob']}",
        f"MRN: {patient['mrn']}",
        f"Date: {date}",
    ]]
    t = Table(data, colWidths=[2.0 * inch, 1.5 * inch, 1.5 * inch, 2.3 * inch])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#f0f4f8")),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("FONTNAME", (0, 0), (-1, -1), "Helvetica"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#cccccc")),
    ]))
    story.append(t)
    story.append(Spacer(1, 8))


def build_signature(story, provider_key: str, styles):
    """Add provider signature block."""
    story.append(Spacer(1, 20))
    story.append(HRFlowable(width="40%", thickness=0.5, color=colors.grey))
    story.append(Paragraph(PROVIDERS[provider_key], styles["ClinicalText"]))
    story.append(Paragraph("Electronically signed", styles["SmallText"]))


def build_disclaimer(story, styles):
    """Add standard disclaimer footer."""
    story.append(Spacer(1, 16))
    story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#cccccc")))
    story.append(Paragraph(
        "CONFIDENTIAL: This document contains protected health information. "
        "Unauthorized disclosure is prohibited under HIPAA (45 CFR Parts 160-164). "
        "SYNTHETIC DATA -- FOR DEMONSTRATION PURPOSES ONLY.",
        styles["SmallText"],
    ))


# =====================================================================
# Document Generators
# =====================================================================

def generate_progress_note(
    outpath: str,
    patient: Dict,
    date: str,
    provider_key: str,
    doc_title: str,
    sections: Dict[str, str],
):
    """Generate a clinical progress note PDF.

    sections: dict of section_name -> text content
    Common keys: chief_complaint, hpi, ros, exam, assessment, plan, vitals
    """
    doc = SimpleDocTemplate(outpath, pagesize=letter,
                            topMargin=36, bottomMargin=36,
                            leftMargin=40, rightMargin=40)
    styles = get_styles()
    story = []

    build_header(story, provider_key, styles)
    build_patient_banner(story, patient, date, styles)
    story.append(Paragraph(doc_title, styles["DocTitle"]))

    section_order = [
        "chief_complaint", "hpi", "vitals", "ros", "exam",
        "labs_reviewed", "imaging_reviewed", "assessment", "plan",
        "follow_up", "referrals", "medications",
    ]
    labels = {
        "chief_complaint": "CHIEF COMPLAINT",
        "hpi": "HISTORY OF PRESENT ILLNESS",
        "vitals": "VITAL SIGNS",
        "ros": "REVIEW OF SYSTEMS",
        "exam": "PHYSICAL EXAMINATION",
        "labs_reviewed": "LABORATORY DATA REVIEWED",
        "imaging_reviewed": "IMAGING REVIEWED",
        "assessment": "ASSESSMENT",
        "plan": "PLAN",
        "follow_up": "FOLLOW-UP",
        "referrals": "REFERRALS",
        "medications": "MEDICATIONS",
    }

    for key in section_order:
        if key in sections:
            story.append(Paragraph(labels.get(key, key.upper()), styles["SectionHead"]))
            for line in sections[key].split("\n"):
                line = line.strip()
                if line:
                    story.append(Paragraph(line, styles["ClinicalText"]))

    # Any extra sections not in standard order
    for key, text in sections.items():
        if key not in section_order:
            story.append(Paragraph(key.upper().replace("_", " "), styles["SectionHead"]))
            for line in text.split("\n"):
                if line.strip():
                    story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    build_signature(story, provider_key, styles)
    build_disclaimer(story, styles)
    doc.build(story)


def generate_lab_report(
    outpath: str,
    patient: Dict,
    date: str,
    order_provider: str,
    panels: List[Dict],
):
    """Generate a laboratory report PDF.

    panels: list of dicts, each with:
        panel_name: str
        results: list of {test, value, unit, ref_range, flag}
    """
    doc = SimpleDocTemplate(outpath, pagesize=letter,
                            topMargin=36, bottomMargin=36,
                            leftMargin=40, rightMargin=40)
    styles = get_styles()
    story = []

    build_header(story, "lab", styles)
    build_patient_banner(story, patient, date, styles)
    story.append(Paragraph("LABORATORY REPORT", styles["DocTitle"]))
    story.append(Paragraph(f"Ordering Provider: {order_provider}", styles["ClinicalText"]))
    story.append(Paragraph(f"Collected: {date}  |  Reported: {date}", styles["SmallText"]))
    story.append(Spacer(1, 8))

    for panel in panels:
        story.append(Paragraph(panel["panel_name"], styles["SectionHead"]))

        header = ["Test", "Result", "Units", "Reference Range", "Flag"]
        rows = [header]
        for r in panel["results"]:
            flag = r.get("flag", "")
            rows.append([
                r["test"], str(r["value"]), r.get("unit", ""),
                r.get("ref_range", ""), flag,
            ])

        t = Table(rows, colWidths=[2.2 * inch, 1.0 * inch, 0.7 * inch, 1.6 * inch, 0.6 * inch])
        style_cmds = [
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3c5e")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cccccc")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f8f8f8")]),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ]
        # Highlight abnormal flags
        for row_idx, r in enumerate(panel["results"], 1):
            flag = r.get("flag", "")
            if flag in ("H", "HH", "A", "L", "LL"):
                style_cmds.append(
                    ("TEXTCOLOR", (4, row_idx), (4, row_idx), colors.HexColor("#cc0000"))
                )
                style_cmds.append(
                    ("FONTNAME", (4, row_idx), (4, row_idx), "Helvetica-Bold")
                )
                style_cmds.append(
                    ("TEXTCOLOR", (1, row_idx), (1, row_idx), colors.HexColor("#cc0000"))
                )

        t.setStyle(TableStyle(style_cmds))
        story.append(t)
        story.append(Spacer(1, 10))

    build_disclaimer(story, styles)
    doc.build(story)


def generate_imaging_report(
    outpath: str,
    patient: Dict,
    date: str,
    provider_key: str,
    modality: str,
    body_part: str,
    indication: str,
    technique: str,
    findings: str,
    impression: str,
    reading_radiologist: str = "Thomas Wright, MD, FACR",
):
    """Generate a radiology/imaging report PDF."""
    doc = SimpleDocTemplate(outpath, pagesize=letter,
                            topMargin=36, bottomMargin=36,
                            leftMargin=40, rightMargin=40)
    styles = get_styles()
    story = []

    build_header(story, provider_key if provider_key in CLINICS else "er", styles)
    build_patient_banner(story, patient, date, styles)
    story.append(Paragraph(f"RADIOLOGY REPORT -- {modality.upper()}", styles["DocTitle"]))

    story.append(Paragraph("EXAM", styles["SectionHead"]))
    story.append(Paragraph(f"{modality} {body_part}", styles["ClinicalText"]))

    story.append(Paragraph("CLINICAL INDICATION", styles["SectionHead"]))
    story.append(Paragraph(indication, styles["ClinicalText"]))

    story.append(Paragraph("TECHNIQUE", styles["SectionHead"]))
    story.append(Paragraph(technique, styles["ClinicalText"]))

    story.append(Paragraph("FINDINGS", styles["SectionHead"]))
    for line in findings.split("\n"):
        if line.strip():
            story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    story.append(Paragraph("IMPRESSION", styles["SectionHead"]))
    for line in impression.split("\n"):
        if line.strip():
            story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    story.append(Spacer(1, 20))
    story.append(HRFlowable(width="40%", thickness=0.5, color=colors.grey))
    story.append(Paragraph(reading_radiologist, styles["ClinicalText"]))
    story.append(Paragraph("Board Certified Radiologist  |  Electronically signed", styles["SmallText"]))
    build_disclaimer(story, styles)
    doc.build(story)


def generate_pathology_report(
    outpath: str,
    patient: Dict,
    date: str,
    specimen: str,
    clinical_history: str,
    gross_description: str,
    microscopic: str,
    diagnosis: str,
    comment: str = "",
    pathologist: str = "Harold Greene, MD",
):
    """Generate a pathology report PDF."""
    doc = SimpleDocTemplate(outpath, pagesize=letter,
                            topMargin=36, bottomMargin=36,
                            leftMargin=40, rightMargin=40)
    styles = get_styles()
    story = []

    build_header(story, "pathology", styles)
    build_patient_banner(story, patient, date, styles)
    story.append(Paragraph("SURGICAL PATHOLOGY REPORT", styles["DocTitle"]))

    for label, text in [
        ("SPECIMEN", specimen),
        ("CLINICAL HISTORY", clinical_history),
        ("GROSS DESCRIPTION", gross_description),
        ("MICROSCOPIC DESCRIPTION", microscopic),
        ("DIAGNOSIS", diagnosis),
    ]:
        story.append(Paragraph(label, styles["SectionHead"]))
        for line in text.split("\n"):
            if line.strip():
                story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    if comment:
        story.append(Paragraph("COMMENT", styles["SectionHead"]))
        for line in comment.split("\n"):
            if line.strip():
                story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    story.append(Spacer(1, 20))
    story.append(HRFlowable(width="40%", thickness=0.5, color=colors.grey))
    story.append(Paragraph(pathologist, styles["ClinicalText"]))
    story.append(Paragraph("Board Certified Pathologist  |  Electronically signed", styles["SmallText"]))
    build_disclaimer(story, styles)
    doc.build(story)


def generate_genetic_report(
    outpath: str,
    patient: Dict,
    date: str,
    test_name: str,
    indication: str,
    methodology: str,
    results: List[Dict],
    interpretation: str,
    recommendations: str = "",
):
    """Generate a genetic/HLA testing report PDF."""
    doc = SimpleDocTemplate(outpath, pagesize=letter,
                            topMargin=36, bottomMargin=36,
                            leftMargin=40, rightMargin=40)
    styles = get_styles()
    story = []

    build_header(story, "genetics", styles)
    build_patient_banner(story, patient, date, styles)
    story.append(Paragraph(f"GENETIC TESTING REPORT -- {test_name}", styles["DocTitle"]))

    story.append(Paragraph("INDICATION", styles["SectionHead"]))
    story.append(Paragraph(indication, styles["ClinicalText"]))

    story.append(Paragraph("METHODOLOGY", styles["SectionHead"]))
    story.append(Paragraph(methodology, styles["ClinicalText"]))

    story.append(Paragraph("RESULTS", styles["SectionHead"]))

    if results:
        header = list(results[0].keys())
        rows = [header]
        for r in results:
            rows.append([str(v) for v in r.values()])

        t = Table(rows, colWidths=[1.5 * inch] * len(header))
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3c5e")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 9),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cccccc")),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ]))
        story.append(t)
        story.append(Spacer(1, 8))

    story.append(Paragraph("INTERPRETATION", styles["SectionHead"]))
    for line in interpretation.split("\n"):
        if line.strip():
            story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    if recommendations:
        story.append(Paragraph("RECOMMENDATIONS", styles["SectionHead"]))
        for line in recommendations.split("\n"):
            if line.strip():
                story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    story.append(Spacer(1, 20))
    story.append(HRFlowable(width="40%", thickness=0.5, color=colors.grey))
    story.append(Paragraph(PROVIDERS["genetics"], styles["ClinicalText"]))
    story.append(Paragraph("ABMGG Board Certified  |  Electronically signed", styles["SmallText"]))
    build_disclaimer(story, styles)
    doc.build(story)


def generate_referral_letter(
    outpath: str,
    patient: Dict,
    date: str,
    from_provider_key: str,
    to_provider_key: str,
    to_provider_name: str,
    reason: str,
    clinical_summary: str,
    specific_questions: str = "",
    urgency: str = "Routine",
):
    """Generate a physician-to-physician referral letter PDF."""
    doc = SimpleDocTemplate(outpath, pagesize=letter,
                            topMargin=36, bottomMargin=36,
                            leftMargin=40, rightMargin=40)
    styles = get_styles()
    story = []

    build_header(story, from_provider_key, styles)
    story.append(Spacer(1, 8))
    story.append(Paragraph(f"Date: {date}", styles["ClinicalText"]))
    story.append(Spacer(1, 4))

    to_clinic = CLINICS.get(to_provider_key, CLINICS["pcp"])
    story.append(Paragraph(f"TO: {to_provider_name}", styles["ClinicalText"]))
    story.append(Paragraph(to_clinic["name"], styles["ClinicalText"]))
    story.append(Paragraph(f"Fax: {to_clinic['fax']}", styles["SmallText"]))
    story.append(Spacer(1, 4))

    story.append(Paragraph(
        f"RE: {patient['name']}  |  DOB: {patient['dob']}  |  MRN: {patient['mrn']}",
        styles["ClinicalText"],
    ))
    story.append(Spacer(1, 4))

    # Urgency badge
    urgency_color = "#cc0000" if urgency.upper() in ("URGENT", "STAT") else "#1a3c5e"
    story.append(Paragraph(
        f"<b>REFERRAL -- {urgency.upper()}</b>",
        ParagraphStyle("Urgency", parent=styles["Normal"],
                       fontSize=12, leading=15, fontName="Helvetica-Bold",
                       textColor=colors.HexColor(urgency_color), spaceAfter=8),
    ))

    story.append(Paragraph("REASON FOR REFERRAL", styles["SectionHead"]))
    for line in reason.split("\n"):
        if line.strip():
            story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    story.append(Paragraph("CLINICAL SUMMARY", styles["SectionHead"]))
    for line in clinical_summary.split("\n"):
        if line.strip():
            story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    if specific_questions:
        story.append(Paragraph("SPECIFIC QUESTIONS", styles["SectionHead"]))
        for line in specific_questions.split("\n"):
            if line.strip():
                story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    story.append(Spacer(1, 12))
    story.append(Paragraph(
        f"Thank you for seeing this patient. Please send consultation report to our office "
        f"at fax {CLINICS[from_provider_key]['fax']}.",
        styles["ClinicalText"],
    ))

    build_signature(story, from_provider_key, styles)
    build_disclaimer(story, styles)
    doc.build(story)


def generate_medication_list(
    outpath: str,
    patient: Dict,
    date: str,
    provider_key: str,
    medications: List[Dict],
    allergies: str = "NKDA",
    pharmacy_notes: str = "",
):
    """Generate a medication reconciliation report PDF.

    medications: list of dicts with:
        name, dose, frequency, prescriber, start_date, indication
    """
    doc = SimpleDocTemplate(outpath, pagesize=letter,
                            topMargin=36, bottomMargin=36,
                            leftMargin=40, rightMargin=40)
    styles = get_styles()
    story = []

    build_header(story, provider_key, styles)
    build_patient_banner(story, patient, date, styles)
    story.append(Paragraph("MEDICATION RECONCILIATION", styles["DocTitle"]))

    # Allergies banner
    allergy_data = [[f"ALLERGIES: {allergies}"]]
    allergy_color = colors.HexColor("#fff3f3") if allergies != "NKDA" else colors.HexColor("#f0f8f0")
    at = Table(allergy_data, colWidths=[7.3 * inch])
    at.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), allergy_color),
        ("FONTSIZE", (0, 0), (-1, -1), 10),
        ("FONTNAME", (0, 0), (-1, -1), "Helvetica-Bold"),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#cccccc")),
    ]))
    story.append(at)
    story.append(Spacer(1, 10))

    # Medication table
    header = ["Medication", "Dose / Route", "Frequency", "Prescriber", "Start Date", "Indication"]
    rows = [header]
    for m in medications:
        rows.append([
            m["name"], m["dose"], m["frequency"],
            m.get("prescriber", ""), m.get("start_date", ""),
            m.get("indication", ""),
        ])

    t = Table(rows, colWidths=[1.4 * inch, 1.0 * inch, 0.8 * inch, 1.0 * inch, 0.8 * inch, 1.3 * inch])
    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1a3c5e")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cccccc")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f8f8f8")]),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]
    t.setStyle(TableStyle(style_cmds))
    story.append(t)
    story.append(Spacer(1, 8))

    story.append(Paragraph(
        f"Total active medications: {len(medications)}",
        styles["ClinicalText"],
    ))

    if pharmacy_notes:
        story.append(Paragraph("PHARMACY NOTES", styles["SectionHead"]))
        for line in pharmacy_notes.split("\n"):
            if line.strip():
                story.append(Paragraph(line.strip(), styles["ClinicalText"]))

    build_signature(story, provider_key, styles)
    build_disclaimer(story, styles)
    doc.build(story)
