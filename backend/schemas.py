"""Pydantic schemas for backend API."""

from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


class GeneSetInput(BaseModel):
    name: str = Field(min_length=1, description="Pathway name, non-empty")
    description: Optional[str] = ""
    genes: List[str] = Field(min_length=1, description="At least one gene required")


class ScoreRequest(BaseModel):
    expression: Dict[str, Dict[str, float]]  # sample -> gene -> value
    gene_sets: Dict[str, List[str]] | List[GeneSetInput]  # pathway -> genes or list
    method: str = Field(default="both", description="ssgsea, zscore, or both")
    alpha: float = Field(default=0.25, gt=0, le=5, description="ssGSEA alpha weighting, (0,5]")

    def normalized_gene_sets(self) -> Dict[str, List[str]]:
        if isinstance(self.gene_sets, dict):
            return self.gene_sets
        return {gs.name: gs.genes for gs in self.gene_sets}

    @field_validator("expression")
    @classmethod
    def validate_expression(cls, v: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
        import math
        if not isinstance(v, dict) or len(v) == 0:
            raise ValueError("expression must be non-empty dict")
        if len(v) > 5000:
            raise ValueError("expression too large: max 5000 samples")
        for sample, genes in v.items():
            if not isinstance(sample, str) or not sample.strip():
                raise ValueError("Sample name must be non-empty string")
            if not isinstance(genes, dict) or len(genes) == 0:
                raise ValueError(f"Expression for sample '{sample}' must be non-empty dict")
            if len(genes) > 30000:
                raise ValueError(f"Too many genes for sample '{sample}': max 30000")
            for gene, val in genes.items():
                if not isinstance(gene, str) or not gene.strip():
                    raise ValueError(f"Gene name empty in sample '{sample}'")
                if not isinstance(val, (int, float)):
                    raise ValueError(f"Expression value for {sample}/{gene} must be numeric")
                if not math.isfinite(val):
                    raise ValueError(f"Expression value for {sample}/{gene} must be finite, got {val}")
        return v

    @field_validator("gene_sets")
    @classmethod
    def validate_gene_sets(cls, v: Dict[str, List[str]] | List[GeneSetInput]) -> Dict[str, List[str]] | List[GeneSetInput]:
        if isinstance(v, dict):
            if len(v) == 0:
                raise ValueError("gene_sets must be non-empty")
            if len(v) > 1000:
                raise ValueError("Too many pathways: max 1000")
            for name, genes in v.items():
                if not isinstance(name, str) or not name.strip():
                    raise ValueError("Pathway name must be non-empty")
                if not isinstance(genes, list) or len(genes) == 0:
                    raise ValueError(f"Gene set '{name}' must have at least one gene")
                if len(genes) > 10000:
                    raise ValueError(f"Gene set '{name}' too large: max 10000 genes")
                for g in genes:
                    if not isinstance(g, str) or not g.strip():
                        raise ValueError(f"Gene name in pathway '{name}' must be non-empty string")
        else:
            if len(v) == 0:
                raise ValueError("gene_sets list must be non-empty")
        return v


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

    @field_validator("scores")
    @classmethod
    def validate_scores(cls, v: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
        import math
        if not isinstance(v, dict) or len(v) == 0:
            raise ValueError("scores must be non-empty")
        for sample, pathways in v.items():
            if not isinstance(sample, str) or not sample.strip():
                raise ValueError("Sample name must be non-empty")
            if not isinstance(pathways, dict) or len(pathways) == 0:
                raise ValueError(f"Scores for sample '{sample}' must be non-empty")
            for pw, val in pathways.items():
                if not isinstance(pw, str) or not pw.strip():
                    raise ValueError("Pathway name must be non-empty")
                if not isinstance(val, (int, float)):
                    raise ValueError(f"Score for {sample}/{pw} must be numeric")
                if not math.isfinite(val) and not (isinstance(val, float) and val != val):  # allow NaN? but finite otherwise
                    # Allow NaN? No, require finite or NaN for missing; but reject inf
                    if val in (float("inf"), float("-inf")):
                        raise ValueError(f"Score for {sample}/{pw} must be finite, got {val}")
        return v

    @field_validator("groups")
    @classmethod
    def validate_groups(cls, v: Dict[str, str]) -> Dict[str, str]:
        if not isinstance(v, dict) or len(v) == 0:
            raise ValueError("groups must be non-empty")
        for sample, group in v.items():
            if not isinstance(sample, str) or not sample.strip():
                raise ValueError("Sample name in groups must be non-empty")
            if not isinstance(group, str) or not group.strip():
                raise ValueError(f"Group for sample '{sample}' must be non-empty string")
        return v


class HealthResponse(BaseModel):
    status: str
    gate: dict
    version: str
