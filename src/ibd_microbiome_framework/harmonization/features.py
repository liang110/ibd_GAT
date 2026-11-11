"""Deterministic feature-name handling used before any fitted transform."""

from __future__ import annotations

import re
from typing import Optional


FORBIDDEN_CLINICAL_FEATURES = {
    "age",
    "age years",
    "age_years",
    "age at sampling",
    "fecal.calprotectin",
    "fecal_calprotectin",
    "fecal calprotectin",
    "fcp",
    "fecal calprotectin concentration",
    "diagnosis",
    "group",
    "label",
    "target",
    "participant id",
    "sample_id",
    "sample id",
    "participant id",
    "participant",
    "antibiotic",
    "antibiotics",
    "steroid",
    "steroids",
    "medication",
    "drug",
    "immunosuppressant",
    "mesalamine",
    "sra metagenome name",
}


def normalize_text(value: object) -> str:
    """Normalize punctuation and whitespace without assigning a biological ID."""
    text = str(value).strip().lower()
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text)
    return text


def canonical_taxon(value: object) -> Optional[str]:
    """Return a conservative species label or None for non-taxon/ambiguous rows."""
    # Strip MetaPhlAn rank prefixes before converting underscores to spaces.
    text = str(value).strip().lower()
    text = re.sub(r"(^|\|)[kpcofgs]__", r"\1", text)
    text = normalize_text(text)
    text = text.replace("unclassified", "").strip()
    if not text or text in {"feature", "species", "metabolite"}:
        return None
    if normalize_text(text) in FORBIDDEN_CLINICAL_FEATURES:
        return None
    return text


def canonical_metabolite(value: object) -> Optional[str]:
    """Create an exact-name key; chemical crosswalks are intentionally separate."""
    text = normalize_text(value)
    if not text or text in {"feature", "metabolite"}:
        return None
    if text in FORBIDDEN_CLINICAL_FEATURES:
        return None
    return text


def is_forbidden_feature(value: object) -> bool:
    return normalize_text(value) in FORBIDDEN_CLINICAL_FEATURES
