"""Homology search of a candidate protein set with MMseqs2."""

from __future__ import annotations
import gzip
import os
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from protflow.types import Table, Provenance

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SWISSPROT_URL = "https://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
COLUMNS = ["query", "target", "pident", "alnlen", "evalue", "bits", "qcov", "tcov", "theader"]


def find_mmseqs() -> list[str]:
    """MMSEQS env var, PATH, or the bundled Windows build in tools/mmseqs (its .bat unpacks busybox on first run)."""
    exe = os.environ.get("MMSEQS") or shutil.which("mmseqs")
    if exe:
        return [exe]
    for cand in (PROJECT_ROOT / "tools" / "mmseqs" / "mmseqs.bat", PROJECT_ROOT / "tools" / "mmseqs2" / "mmseqs.bat"):
        if cand.exists():
            return [str(cand)]
    raise FileNotFoundError("mmseqs bulunamadı: PATH'e ekleyin, MMSEQS ortam değişkenini ayarlayın veya tools/mmseqs altına kurun.")


def resolve_target(target: str) -> Path:
    """'swissprot' downloads UniProt Swiss-Prot once (~90 MB) into tools/db; anything else is a FASTA path."""
    if target != "swissprot":
        p = Path(target)
        if not p.is_file():
            raise FileNotFoundError(f"Hedef FASTA yok: {target}")
        return p
    db = PROJECT_ROOT / "tools" / "db"
    db.mkdir(parents=True, exist_ok=True)
    fa = db / "uniprot_sprot.fasta"
    if not fa.exists():
        tmp_gz = db / "uniprot_sprot.fasta.gz.tmp"
        tmp_fa = db / "uniprot_sprot.fasta.tmp"
        try:
            urllib.request.urlretrieve(SWISSPROT_URL, tmp_gz)
            with gzip.open(tmp_gz, "rb") as src, open(tmp_fa, "wb") as dst:
                shutil.copyfileobj(src, dst)
            tmp_fa.replace(fa)
        finally:
            if tmp_gz.exists():
                tmp_gz.unlink()
            if tmp_fa.exists():
                tmp_fa.unlink()
    return fa


def parse_mmseqs_homology(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Searches every protein in the candidate set against a target FASTA (default Swiss-Prot).
    Output: one row per hit. Queries without hits are counted in provenance (no hit = no known homolog above the e-value).
    Parameters: target_fasta ('swissprot' or a FASTA path), evalue, max_hits (per query), sensitivity.
    """
    coll = context["inputs"]["proteins"]
    params = context.get("parameters", {})
    target = resolve_target(str(params.get("target_fasta") or "swissprot"))
    evalue = float(params.get("evalue", 1e-3))
    max_hits = int(params.get("max_hits", 5))
    sens = float(params.get("sensitivity", 5.7))

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        q = tmp / "query.fasta"
        q.write_text("".join(f">{s.protein_id}\n{s.sequence}\n" for s in coll.sequences.values()), encoding="utf-8")
        out = tmp / "hits.m8"
        cmd = find_mmseqs() + ["easy-search", str(q), str(target), str(out), str(tmp / "work"),
                               "--format-output", ",".join(COLUMNS), "-e", str(evalue),
                               "--max-seqs", str(max_hits * 20), "-s", str(sens), "-v", "1"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"mmseqs hata verdi ({res.returncode}): {(res.stderr or res.stdout)[-500:]}")
        df = pd.read_csv(out, sep="\t", names=COLUMNS) if out.stat().st_size else pd.DataFrame(columns=COLUMNS)

    df = df.sort_values(["query", "evalue"]).groupby("query", sort=False).head(max_hits).reset_index(drop=True)
    table = Table(protein_id=coll.protein_id, df=df)
    table.provenance = Provenance(
        tool_name="mmseqs2",
        parameters={"input_count": len(coll.sequences), "queries_with_hits": int(df["query"].nunique()),
                    "no_hit": len(coll.sequences) - int(df["query"].nunique()),
                    "target": target.name, "evalue": evalue},
    )
    return {"homology": table}
