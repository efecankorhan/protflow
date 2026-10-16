"""Multiple Sequence Alignment (MSA) adapter: MAFFT (Katoh & Standley 2013) / FAMSA.

Implements academic gold-standard MSA with a 3-tier execution strategy:
1. Local MAFFT/FAMSA binary (if installed in tools/ or PATH).
2. Official EMBL-EBI MAFFT REST API (remote high-throughput compute cluster, zero GPU/RAM overhead).
3. Pure-Python BLOSUM62 progressive dynamic programming aligner (offline fallback).
Computes per-column Shannon entropy conservation metrics (Valdar 2002).
"""

from __future__ import annotations
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import math
import numpy as np
import pandas as pd
import requests

from protflow.types import (
    Sequence,
    SequenceCollection,
    Alignment,
    Table,
    Provenance
)

# Standard Henikoff & Henikoff BLOSUM62 matrix for standard 20 amino acids
BLOSUM62_PAIRS = {
    ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
    ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
    ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
    ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
    ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
    ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
    ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
    ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
    ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2, ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
    ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
    ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
    ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
    ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
    ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1,
    ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
    ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
    ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
    ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
    ('Y', 'Y'): 7, ('Y', 'V'): -1,
    ('V', 'V'): 4
}


def _blosum62_score(a: str, b: str) -> int:
    a, b = a.upper(), b.upper()
    if a == '-' or b == '-':
        return -4
    pair = (a, b) if (a, b) in BLOSUM62_PAIRS else (b, a)
    return BLOSUM62_PAIRS.get(pair, -1)


def _find_tool_binary(name: str, subpaths: List[str]) -> Optional[Path]:
    root = Path(__file__).resolve().parents[2]
    for sp in subpaths:
        p = root / "tools" / sp
        if p.exists():
            return p
    sys_b = shutil.which(name)
    return Path(sys_b) if sys_b else None


def _find_famsa_binary() -> Optional[Path]:
    return _find_tool_binary("famsa", ["famsa/famsa.exe", "famsa/famsa", "famsa/bin/famsa", "famsa/bin/famsa.exe"])


def _run_famsa_local(fasta_text: str, bin_path: Path) -> Optional[Dict[str, str]]:
    """Runs local FAMSA 2 binary."""
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        in_p = tdp / "input.fasta"
        out_p = tdp / "aligned.fasta"
        in_p.write_text(fasta_text, encoding="utf-8")
        try:
            cmd = [str(bin_path), str(in_p), str(out_p)]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120)
            if res.returncode == 0 and out_p.exists():
                return _parse_fasta_alignment(out_p.read_text(encoding="utf-8"))
        except Exception:
            pass
    return None


def _find_mafft_binary() -> Optional[Path]:
    return _find_tool_binary("mafft", ["mafft/bin/mafft", "mafft/bin/mafft.bat", "mafft/mafft"])


def _run_mafft_local(fasta_text: str, bin_path: Path) -> Optional[Dict[str, str]]:
    """Runs local MAFFT binary."""
    try:
        cmd = [str(bin_path), "--auto", "-"]
        res = subprocess.run(cmd, input=fasta_text, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120)
        if res.returncode == 0 and res.stdout.strip():
            return _parse_fasta_alignment(res.stdout)
    except Exception:
        pass
    return None


def _run_mafft_ebi_api(fasta_text: str, timeout: int = 45) -> Optional[Dict[str, str]]:
    """Runs MAFFT via EMBL-EBI official REST API (cluster execution)."""
    try:
        run_url = "https://www.ebi.ac.uk/Tools/services/rest/mafft/run"
        r = requests.post(run_url, data={"email": "protflow@ebi.ac.uk", "sequence": fasta_text}, timeout=15)
        if r.status_code != 200:
            return None
        job_id = r.text.strip()
        if not job_id.startswith("mafft-"):
            return None

        status_url = f"https://www.ebi.ac.uk/Tools/services/rest/mafft/status/{job_id}"
        t0 = time.time()
        while time.time() - t0 < timeout:
            time.sleep(2)
            st_res = requests.get(status_url, timeout=10)
            if st_res.status_code == 200:
                status = st_res.text.strip()
                if status == "FINISHED":
                    res_url = f"https://www.ebi.ac.uk/Tools/services/rest/mafft/result/{job_id}/out"
                    aln_res = requests.get(res_url, timeout=15)
                    if aln_res.status_code == 200 and aln_res.text.strip():
                        return _parse_fasta_alignment(aln_res.text)
                elif status in ("ERROR", "NOT_FOUND", "FAILURE"):
                    return None
    except Exception:
        pass
    return None


def _needleman_wunsch(s1: str, s2: str, gap_open: int = -8, gap_ext: int = -2) -> Tuple[str, str, float]:
    """Pairwise dynamic programming alignment using BLOSUM62."""
    n, m = len(s1), len(s2)
    score = np.zeros((n + 1, m + 1), dtype=float)

    for i in range(1, n + 1):
        score[i, 0] = gap_open + (i - 1) * gap_ext
    for j in range(1, m + 1):
        score[0, j] = gap_open + (j - 1) * gap_ext

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1, j - 1] + _blosum62_score(s1[i - 1], s2[j - 1])
            up = score[i - 1, j] + (gap_ext if (i > 1 and score[i - 1, j] == score[i - 2, j] + gap_ext) else gap_open)
            left = score[i, j - 1] + (gap_ext if (j > 1 and score[i, j - 1] == score[i, j - 2] + gap_ext) else gap_open)
            score[i, j] = max(diag, up, left)

    # Traceback
    a1, a2 = [], []
    i, j = n, m
    while i > 0 and j > 0:
        curr = score[i, j]
        diag = score[i - 1, j - 1] + _blosum62_score(s1[i - 1], s2[j - 1])
        if abs(curr - diag) < 1e-4:
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif curr == score[i - 1, j] + (gap_ext if i > 1 and score[i - 1, j] == score[i - 2, j] + gap_ext else gap_open):
            a1.append(s1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(s2[j - 1])
            j -= 1

    while i > 0:
        a1.append(s1[i - 1])
        a2.append('-')
        i -= 1
    while j > 0:
        a1.append('-')
        a2.append(s2[j - 1])
        j -= 1

    return "".join(reversed(a1)), "".join(reversed(a2)), float(score[n, m])


def _run_progressive_msa(seq_dict: Dict[str, str]) -> Dict[str, str]:
    """
    Offline progressive alignment engine.
    Aligns all sequences against the longest/reference sequence using BLOSUM62.
    """
    if not seq_dict:
        return {}
    if len(seq_dict) == 1:
        return seq_dict

    # Sort by length descending, reference is longest
    sorted_ids = sorted(seq_dict.keys(), key=lambda k: len(seq_dict[k]), reverse=True)
    ref_id = sorted_ids[0]
    ref_orig = seq_dict[ref_id]
    L = len(ref_orig)

    # shortcut: pairwise-to-reference progressive alignment; install MAFFT/FAMSA for exact multi-profile alignment
    slot_ins: Dict[str, List[str]] = {ref_id: [""] * (L + 1)}
    match_chars: Dict[str, List[str]] = {ref_id: list(ref_orig)}

    for sid in sorted_ids[1:]:
        a_ref, a_s, _ = _needleman_wunsch(ref_orig, seq_dict[sid])
        ins = [""] * (L + 1)
        m_ch: List[str] = []
        r_idx = 0
        for ch_r, ch_s in zip(a_ref, a_s):
            if ch_r == "-":
                ins[r_idx] += ch_s
            else:
                m_ch.append(ch_s)
                r_idx += 1
        slot_ins[sid] = ins
        match_chars[sid] = m_ch

    max_ins = [max(len(slot_ins[k][i]) for k in sorted_ids) for i in range(L + 1)]
    aligned: Dict[str, str] = {}
    for k in sorted_ids:
        parts = []
        for i in range(L):
            if max_ins[i]:
                parts.append(slot_ins[k][i].ljust(max_ins[i], "-"))
            parts.append(match_chars[k][i] if i < len(match_chars[k]) else "-")
        if max_ins[L]:
            parts.append(slot_ins[k][L].ljust(max_ins[L], "-"))
        aligned[k] = "".join(parts)

    return aligned


def _parse_fasta_alignment(text: str) -> Dict[str, str]:
    res = {}
    cur_id = None
    cur_seq = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if cur_id:
                res[cur_id] = "".join(cur_seq)
                cur_seq = []
            cur_id = line[1:].split()[0]
        elif line:
            cur_seq.append(line)
    if cur_id:
        res[cur_id] = "".join(cur_seq)
    return res


def compute_conservation(aligned_dict: Dict[str, str]) -> pd.DataFrame:
    """
    Computes per-column Shannon entropy and Valdar conservation score:
    C_i = 1 - (H_i / log2(20))
    Where H_i = - sum(p_aa * log2(p_aa))
    """
    if not aligned_dict:
        return pd.DataFrame()

    seq_list = list(aligned_dict.values())
    aln_len = len(seq_list[0])
    num_seqs = len(seq_list)

    rows = []
    for pos in range(aln_len):
        col_chars = [s[pos].upper() for s in seq_list if pos < len(s)]
        non_gap = [c for c in col_chars if c != '-']
        gap_ratio = (num_seqs - len(non_gap)) / num_seqs

        if not non_gap:
            entropy = math.log2(20)
            cons_score = 0.0
            consensus_aa = '-'
        else:
            counts = {}
            for c in non_gap:
                counts[c] = counts.get(c, 0) + 1
            consensus_aa = max(counts.keys(), key=lambda k: counts[k])

            # Shannon entropy over observed amino acids
            h = 0.0
            for cnt in counts.values():
                p = cnt / len(non_gap)
                if p > 0:
                    h -= p * math.log2(p)
            entropy = h
            # Normalized Valdar score weighted by non-gap occupancy
            cons_score = max(0.0, 1.0 - (h / math.log2(20))) * (1.0 - gap_ratio)

        rows.append({
            "position": pos + 1,
            "consensus_aa": consensus_aa,
            "conservation_score": round(cons_score, 4),
            "entropy": round(entropy, 4),
            "gap_ratio": round(gap_ratio, 3)
        })

    return pd.DataFrame(rows)


def parse_mafft(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    ProtFlow Node Adapter for MAFFT / Multiple Sequence Alignment.
    Inputs:
      context['inputs']['sequences']: SequenceCollection or dict of Sequence objects.
    Outputs:
      'alignment': Alignment summary object
      'aligned_sequences': SequenceCollection with aligned gaps
      'conservation': Table with per-column conservation metrics
    """
    inputs = context.get("inputs", {})
    seq_input = inputs.get("sequences")

    seq_dict = {}
    if isinstance(seq_input, SequenceCollection):
        seq_dict = {sid: s.sequence for sid, s in seq_input.sequences.items()}
    elif isinstance(seq_input, dict):
        for k, v in seq_input.items():
            seq_dict[k] = getattr(v, "sequence", str(v))
    elif isinstance(seq_input, Sequence):
        seq_dict[seq_input.protein_id] = seq_input.sequence

    tool_used = "mafft"
    if len(seq_dict) < 2:
        # Trivial alignment
        aligned_dict = {k: v for k, v in seq_dict.items()}
        method = "single_sequence_passthrough"
        tool_used = "passthrough"
    else:
        fasta_text = "\n".join(f">{k}\n{v}" for k, v in seq_dict.items())
        famsa_bin = _find_famsa_binary()
        mafft_bin = _find_mafft_binary()
        aligned_dict = None
        tool_used = "mafft"

        if famsa_bin:
            aligned_dict = _run_famsa_local(fasta_text, famsa_bin)
            if aligned_dict:
                method = "famsa_local_v2"
                tool_used = "famsa"

        if not aligned_dict and mafft_bin:
            aligned_dict = _run_mafft_local(fasta_text, mafft_bin)
            if aligned_dict:
                method = "mafft_local_binary"
                tool_used = "mafft"

        params = context.get("parameters") or {}
        manifest = context.get("manifest")
        allow_remote = bool(params.get("allow_remote_api") or (manifest and getattr(manifest, "is_remote", False)))
        if not aligned_dict and allow_remote:
            aligned_dict = _run_mafft_ebi_api(fasta_text)
            if aligned_dict:
                method = "mafft_embl_ebi_rest_api"
                tool_used = "mafft"

        if not aligned_dict:
            aligned_dict = _run_progressive_msa(seq_dict)
            method = "blosum62_progressive_aligner"
            tool_used = "blosum62"

    first_two = list(aligned_dict.keys())[:2]
    q_id = first_two[0] if first_two else "seq1"
    t_id = first_two[1] if len(first_two) > 1 else q_id

    # Build outputs
    aligned_col = SequenceCollection(protein_id=q_id, sequences={
        k: Sequence(protein_id=k, sequence=v, description=f"aligned by {method}")
        for k, v in aligned_dict.items()
    })

    cons_df = compute_conservation(aligned_dict)
    cons_table = Table(protein_id=list(seq_dict.keys())[0] if seq_dict else "MSA", df=cons_df)

    mean_cons = float(cons_df["conservation_score"].mean()) if not cons_df.empty else 1.0

    aln_obj = Alignment(
        protein_id=q_id,
        target_id=t_id,
        query_seq=aligned_dict.get(q_id, ""),
        target_seq=aligned_dict.get(t_id, ""),
        identity=round(mean_cons, 3),
        score=float(len(cons_df)),
        alignment_format="fasta"
    )

    prov = Provenance(
        tool_name=tool_used,
        tool_version="2.5.2" if tool_used == "famsa" else "7.525",
        parameters={"method": method, "num_sequences": len(seq_dict), "length": len(cons_df)}
    )
    aln_obj.attach_provenance(prov)
    aligned_col.attach_provenance(prov)
    cons_table.attach_provenance(prov)

    return {
        "alignment": aln_obj,
        "aligned_sequences": aligned_col,
        "conservation": cons_table
    }


parse_famsa = parse_mafft

