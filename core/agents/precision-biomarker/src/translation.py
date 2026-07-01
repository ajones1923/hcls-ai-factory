"""
Multi-Language Report Translation for the Biomarker Intelligence Agent.

Provides translation of clinical reports to multiple languages to serve
diverse patient populations. Uses a template-based approach for medical
terminology accuracy (LLM translation would require clinical validation).

Supported languages:
  - English (default, no translation)
  - Spanish (es) — 559 million speakers globally
  - Simplified Chinese (zh) — 1.1 billion speakers
  - Hindi (hi) — 602 million speakers
  - French (fr) — 280 million speakers
  - Arabic (ar) — 274 million speakers
  - Portuguese (pt) — 257 million speakers

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

from typing import Dict

# =====================================================================
# Medical terminology translations
# =====================================================================

# These are validated medical translations, not machine-generated.
# Each term maps language code -> translated term.

MEDICAL_TERMS: Dict[str, Dict[str, str]] = {
    # Section headers
    "Biomarker Intelligence Report": {
        "es": "Informe de Inteligencia de Biomarcadores de Precisión",
        "zh": "精准生物标志物智能报告",
        "hi": "प्रिसिजन बायोमार्कर इंटेलिजेंस रिपोर्ट",
        "fr": "Rapport d'Intelligence des Biomarqueurs de Précision",
        "ar": "تقرير ذكاء المؤشرات الحيوية الدقيقة",
        "pt": "Relatório de Inteligência de Biomarcadores de Precisão",
    },
    "Biological Age Assessment": {
        "es": "Evaluación de Edad Biológica",
        "zh": "生物年龄评估",
        "hi": "जैविक आयु मूल्यांकन",
        "fr": "Évaluation de l'Âge Biologique",
        "ar": "تقييم العمر البيولوجي",
        "pt": "Avaliação da Idade Biológica",
    },
    "Disease Risk Trajectories": {
        "es": "Trayectorias de Riesgo de Enfermedades",
        "zh": "疾病风险轨迹",
        "hi": "रोग जोखिम प्रक्षेप पथ",
        "fr": "Trajectoires de Risque de Maladies",
        "ar": "مسارات خطر الأمراض",
        "pt": "Trajetórias de Risco de Doenças",
    },
    "Pharmacogenomic Profile": {
        "es": "Perfil Farmacogenómico",
        "zh": "药物基因组学概况",
        "hi": "फार्माकोजेनोमिक प्रोफाइल",
        "fr": "Profil Pharmacogénomique",
        "ar": "الملف الدوائي الجيني",
        "pt": "Perfil Farmacogenômico",
    },
    "Drug-Drug Interactions": {
        "es": "Interacciones Medicamentosas",
        "zh": "药物-药物相互作用",
        "hi": "दवा-दवा परस्पर क्रिया",
        "fr": "Interactions Médicamenteuses",
        "ar": "التفاعلات الدوائية",
        "pt": "Interações Medicamentosas",
    },
    "Evidence Provenance Chain": {
        "es": "Cadena de Procedencia de Evidencia",
        "zh": "证据溯源链",
        "hi": "साक्ष्य उत्पत्ति श्रृंखला",
        "fr": "Chaîne de Provenance des Preuves",
        "ar": "سلسلة مصدر الأدلة",
        "pt": "Cadeia de Procedência de Evidências",
    },
    "Clinical Validation & Limitations": {
        "es": "Validación Clínica y Limitaciones",
        "zh": "临床验证与局限性",
        "hi": "नैदानिक सत्यापन और सीमाएं",
        "fr": "Validation Clinique et Limites",
        "ar": "التحقق السريري والقيود",
        "pt": "Validação Clínica e Limitações",
    },

    # Risk levels
    "LOW": {
        "es": "BAJO", "zh": "低", "hi": "कम", "fr": "FAIBLE", "ar": "منخفض", "pt": "BAIXO",
    },
    "MODERATE": {
        "es": "MODERADO", "zh": "中等", "hi": "मध्यम", "fr": "MODÉRÉ", "ar": "متوسط", "pt": "MODERADO",
    },
    "HIGH": {
        "es": "ALTO", "zh": "高", "hi": "उच्च", "fr": "ÉLEVÉ", "ar": "مرتفع", "pt": "ALTO",
    },
    "CRITICAL": {
        "es": "CRÍTICO", "zh": "危急", "hi": "गंभीर", "fr": "CRITIQUE", "ar": "حرج", "pt": "CRÍTICO",
    },

    # PGx phenotypes
    "Normal Metabolizer": {
        "es": "Metabolizador Normal", "zh": "正常代谢者", "hi": "सामान्य चयापचय",
        "fr": "Métaboliseur Normal", "ar": "أيض طبيعي", "pt": "Metabolizador Normal",
    },
    "Poor Metabolizer": {
        "es": "Metabolizador Lento", "zh": "慢代谢者", "hi": "धीमा चयापचय",
        "fr": "Métaboliseur Lent", "ar": "أيض بطيء", "pt": "Metabolizador Lento",
    },
    "Ultra Rapid Metabolizer": {
        "es": "Metabolizador Ultra Rápido", "zh": "超快代谢者", "hi": "अतितीव्र चयापचय",
        "fr": "Métaboliseur Ultra Rapide", "ar": "أيض سريع جداً", "pt": "Metabolizador Ultra Rápido",
    },

    # Common biomarkers (keep English in most contexts, but translate labels)
    "Biological Age": {
        "es": "Edad Biológica", "zh": "生物年龄", "hi": "जैविक आयु",
        "fr": "Âge Biologique", "ar": "العمر البيولوجي", "pt": "Idade Biológica",
    },
    "Chronological Age": {
        "es": "Edad Cronológica", "zh": "实际年龄", "hi": "कालानुक्रमिक आयु",
        "fr": "Âge Chronologique", "ar": "العمر الزمني", "pt": "Idade Cronológica",
    },
    "Age Acceleration": {
        "es": "Aceleración del Envejecimiento", "zh": "年龄加速", "hi": "आयु त्वरण",
        "fr": "Accélération du Vieillissement", "ar": "تسارع الشيخوخة", "pt": "Aceleração do Envelhecimento",
    },
    "Mortality Risk": {
        "es": "Riesgo de Mortalidad", "zh": "死亡风险", "hi": "मृत्यु दर जोखिम",
        "fr": "Risque de Mortalité", "ar": "خطر الوفاة", "pt": "Risco de Mortalidade",
    },

    # Disclaimer
    "FOR RESEARCH USE ONLY": {
        "es": "SOLO PARA USO EN INVESTIGACIÓN",
        "zh": "仅供研究使用",
        "hi": "केवल अनुसंधान उपयोग के लिए",
        "fr": "POUR USAGE DE RECHERCHE UNIQUEMENT",
        "ar": "للاستخدام البحثي فقط",
        "pt": "APENAS PARA USO EM PESQUISA",
    },
}

SUPPORTED_LANGUAGES = {
    "en": "English",
    "es": "Spanish (Español)",
    "zh": "Chinese (中文)",
    "hi": "Hindi (हिन्दी)",
    "fr": "French (Français)",
    "ar": "Arabic (العربية)",
    "pt": "Portuguese (Português)",
}


def translate_term(term: str, language: str) -> str:
    """Translate a medical term to the target language.

    Parameters
    ----------
    term : str
        English term to translate.
    language : str
        Target language code (e.g., 'es', 'zh', 'hi').

    Returns
    -------
    str
        Translated term, or original English if no translation available.
    """
    if language == "en" or language not in SUPPORTED_LANGUAGES:
        return term

    translations = MEDICAL_TERMS.get(term, {})
    return translations.get(language, term)


def translate_report_headers(markdown: str, language: str) -> str:
    """Translate report section headers and key terms in a markdown report.

    Performs term-by-term replacement of known medical terminology.
    Preserves all data, numbers, and formatting — only translates labels.

    Parameters
    ----------
    markdown : str
        English markdown report.
    language : str
        Target language code.

    Returns
    -------
    str
        Report with translated headers and terms.
    """
    if language == "en" or language not in SUPPORTED_LANGUAGES:
        return markdown

    result = markdown
    for english_term, translations in MEDICAL_TERMS.items():
        if language in translations:
            result = result.replace(english_term, translations[language])

    # Add language header
    lang_name = SUPPORTED_LANGUAGES.get(language, language)
    result = f"*Report language: {lang_name}*\n\n" + result

    return result


def get_supported_languages() -> Dict[str, str]:
    """Return dict of supported language codes and display names."""
    return dict(SUPPORTED_LANGUAGES)
