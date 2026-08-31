"""GMT parser for MSigDB / Reactome gene sets.

Format (tab-separated):
    pathway_name \\t description \\t gene1 \\t gene2 \\t ... \\t geneN

Real sources:
- MSigDB Hallmark (h.all.v2023.2.Hs.symbols.gmt, 50 gene sets)
  https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H
  Download: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt
  License: free for non-commercial/academic use.

- ReactomePathways.gmt
  https://reactome.org/download-data (\"ReactomePathways.gmt\")
  Fully open data, CC BY 4.0.

Usage: parse_gmt(path) returns Dict[str, Dict] with keys name, description, genes.
CLI entry: python -m data_pipeline.cli --gmt-path <file>
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List


def _strip_bom(text: str) -> str:
    """Remove UTF-8 BOM if present (common in MSigDB downloads on Windows)."""
    return text.lstrip("\ufeff")


def parse_gmt(path: str | Path) -> Dict[str, Dict[str, object]]:
    """Parse a GMT file.

    Robust to real-file quirks:
    - UTF-8 BOM (\\ufeff) at start of file or line
    - Comment lines starting with '#'
    - CRLF / LF line endings
    - Leading/trailing whitespace on names/descriptions/genes
    - Blank lines and whitespace-only lines
    - Multiple tab variants: tabs separate fields; empty gene fields filtered
    - Duplicate pathway names raise ValueError

    Args:
        path: path to .gmt file

    Returns:
        dict pathway_name -> {\"description\": str, \"genes\": List[str]}

    Raises:
        FileNotFoundError, ValueError on malformed file.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"GMT file not found: {p}")
    result: Dict[str, Dict[str, object]] = {}
    with p.open("r", encoding="utf-8-sig") as fh:
        for lineno, raw_line in enumerate(fh, 1):
            # Handle CRLF and BOM per line
            line = _strip_bom(raw_line.rstrip("\r\n"))
            if not line.strip():
                continue
            # Skip comment lines (common in some GMT exports)
            if line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(
                    f"Malformed GMT at line {lineno}: need at least 3 tab-separated fields, got {len(parts)}"
                )
            name = _strip_bom(parts[0].strip())
            description = parts[1].strip()
            # Genes are remaining fields; filter empty strings and dedup preserving order
            raw_genes = [g.strip() for g in parts[2:] if g.strip()]
            # dedup preserve order, upper-case normalization not forced (keep as-is)
            seen = set()
            genes: List[str] = []
            for g in raw_genes:
                if g not in seen:
                    seen.add(g)
                    genes.append(g)
            if not name:
                raise ValueError(f"Empty pathway name at line {lineno}")
            # duplicate pathway names disallowed
            if name in result:
                raise ValueError(f"Duplicate pathway name '{name}' at line {lineno}")
            result[name] = {"description": description, "genes": genes}
    if not result:
        raise ValueError("GMT file contained no valid pathways")
    return result


def parse_gmt_string(content: str) -> Dict[str, Dict[str, object]]:
    """Parse GMT from string (useful for tests). Mirrors file parser strictness."""
    result: Dict[str, Dict[str, object]] = {}
    # Handle BOM at start of string
    content = _strip_bom(content)
    for lineno, raw_line in enumerate(content.splitlines(), 1):
        line = _strip_bom(raw_line.rstrip("\r\n"))
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            raise ValueError(f"Malformed GMT at line {lineno}: need at least 3 tab-separated fields, got {len(parts)}")
        name = parts[0].strip()
        description = parts[1].strip()
        if not name:
            raise ValueError(f"Empty pathway name at line {lineno}")
        if name in result:
            raise ValueError(f"Duplicate pathway name '{name}' at line {lineno}")
        raw_genes = [g.strip() for g in parts[2:] if g.strip()]
        seen = set()
        genes: List[str] = []
        for g in raw_genes:
            if g not in seen:
                seen.add(g)
                genes.append(g)
        result[name] = {"description": description, "genes": genes}
    # Mirror file parser: empty result is an error (no valid pathways)
    # Note: callers that expect empty dict on only-comments should handle via try/except;
    # for backward compatibility we only raise if content was non-empty but had no valid lines?
    # Here we keep strict: empty string returns empty dict (no error) to not break existing fixtures,
    # but content with only comments/whitespace after stripping that still yields no pathways raises.
    # Simpler: if content.strip() and not result: raise. Preserve empty-string -> empty dict.
    if content.strip() and not result:
        raise ValueError("GMT content contained no valid pathways")
    return result
