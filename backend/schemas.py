"""Pydantic schemas for backend API."""

from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel, Field


class GeneSetInput(BaseModel):
    name: str
    description: Optional[str] = ""
    genes: List[str]


class ScoreRequest(BaseModel):
    expression: Dict[str, Dict[str, float]]  # sample -> gene -> value
    gene_sets: Dict[str, List[str]] | List[GeneSetInput]  # pathway -> genes or list
    method: str = Field(default="both", description="ssgsea, zscore, or both")
    alpha: float = Field(default=0.25, description="ssGSEA alpha weighting")

    def normalized_gene_sets(self) -> Dict[str, List[str]]:
        if isinstance(self.gene_sets, dict):
            return self.gene_sets
        return {gs.name: gs.genes for gs in self.gene_sets}


class ScoreResponse(BaseModel):
    samples: List[str]
    pathways: List[str]
    scores: Dict[str, Dict[str, float]]  # method -> pathway -> sample? Actually sample->pathway->score nested
    # We'll return method -> sample -> pathway -> score


class DifferentialRequest(BaseModel):
    scores: Dict[str, Dict[str, float]]  # sample -> pathway -> score
    groups: Dict[str, str]  # sample -> group
    group_a: Optional[str] = None
    group_b: Optional[str] = None


class HealthResponse(BaseModel):
    status: str
    gate: dict
    version: str
