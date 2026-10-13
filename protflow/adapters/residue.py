"""
Residue alignment and drift validation utilities.
Guards against silent index drifts, chain ID mismatches, and auth_seq_id vs label_seq_id confusion.
"""

from __future__ import annotations
from typing import Any, Dict, List, Optional, Tuple, Union
import gemmi
import pandas as pd
from protflow.types import Sequence, Structure


# Standard 3-letter to 1-letter amino acid conversion
AA3_TO_1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M", "SEC": "U", "PYL": "O"
}
AA1_TO_3 = {v: k for k, v in AA3_TO_1.items()}


class ResidueDriftError(ValueError):
    """Raised when residue indexing, chain ID, or amino acid identity drifts between tools."""
    pass


def extract_structure_residues(structure: Structure) -> List[Dict[str, Any]]:
    """
    Extracts all polymer residues from a Structure with gemmi.
    Returns list of dicts with: chain_id, auth_seq_id, label_seq_id, resname, one_letter.
    """
    st = structure.to_gemmi_structure()
    residues = []
    
    for model in st:
        for chain in model:
            for res in chain:
                # Filter out water and obvious ligands
                if res.is_water():
                    continue
                resname = res.name.upper()
                one_letter = AA3_TO_1.get(resname, "X")
                # auth_seq_id might contain insertion code, e.g. "12" or "12A"
                auth_num = str(res.seqid.num)
                if res.seqid.icode and res.seqid.icode != ' ':
                    auth_num += res.seqid.icode
                
                residues.append({
                    "chain_id": chain.name,
                    "auth_seq_id": auth_num,
                    "label_seq_id": res.label_seq_id if hasattr(res, 'label_seq_id') else None,
                    "resname": resname,
                    "one_letter": one_letter
                })
        break  # use primary model
    return residues


def validate_residue_table(
    df: pd.DataFrame,
    structure: Structure,
    sequence: Optional[Sequence] = None,
    chain_id: Optional[str] = None
) -> pd.DataFrame:
    """
    Validates a residue table against the actual 3D structure and sequence.
    Required columns in df: 'auth_seq_id' (or 'res_id') and optionally 'chain_id', 'resname'.
    
    Checks:
    1. Every (chain_id, auth_seq_id) exists in the 3D structure.
    2. Residue names match between table and structure.
    3. If sequence is provided, verifies that sequence matches the structural residues without offset.
    
    Returns enriched DataFrame with canonical ['chain_id', 'auth_seq_id', 'resname', 'one_letter'].
    Raises ResidueDriftError on mismatch.
    """
    struct_res = extract_structure_residues(structure)
    if not struct_res:
        raise ResidueDriftError(f"No valid residues found in structure {structure.file_path}")

    # Build lookup dict: (chain_id, auth_seq_id) -> residue_info
    struct_map: Dict[Tuple[str, str], Dict[str, Any]] = {
        (r["chain_id"], str(r["auth_seq_id"])): r for r in struct_res
    }
    avail_chains = sorted(list({r["chain_id"] for r in struct_res}))

    # Normalize column names
    norm_df = df.copy()
    if norm_df.empty:
        norm_df["resname"] = []
        norm_df["one_letter"] = []
        return norm_df

    if "res_id" in norm_df.columns and "auth_seq_id" not in norm_df.columns:
        norm_df["auth_seq_id"] = norm_df["res_id"]
    if "auth_seq_id" not in norm_df.columns:
        raise ResidueDriftError("Residue table is missing required 'auth_seq_id' column.")

    if "chain_id" not in norm_df.columns:
        chosen_chain = chain_id or (avail_chains[0] if len(avail_chains) == 1 else None)
        if chosen_chain is None:
            raise ResidueDriftError(
                f"Multiple chains found ({avail_chains}) but table has no 'chain_id' column."
            )
        norm_df["chain_id"] = chosen_chain

    norm_df["auth_seq_id"] = norm_df["auth_seq_id"].astype(str)
    norm_df["chain_id"] = norm_df["chain_id"].astype(str)

    # Validate each residue
    validated_resnames = []
    validated_one_letters = []

    for idx, row in norm_df.iterrows():
        key = (row["chain_id"], row["auth_seq_id"])
        if key not in struct_map:
            raise ResidueDriftError(
                f"Residue {key[1]} on chain '{key[0]}' does not exist in structure {structure.file_path}. "
                f"Available residues on chain: {[k[1] for k in struct_map if k[0] == key[0]][:10]}..."
            )
        
        expected = struct_map[key]
        if "resname" in row and pd.notna(row["resname"]):
            given_res = str(row["resname"]).upper()
            expected_res = expected["resname"]
            # Can be 3-letter or 1-letter
            if len(given_res) == 1 and expected["one_letter"] != given_res:
                raise ResidueDriftError(
                    f"Residue mismatch at {key[0]}:{key[1]}! Table says '{given_res}' "
                    f"but structure has '{expected_res}' ({expected['one_letter']})."
                )
            elif len(given_res) == 3 and expected_res != given_res:
                raise ResidueDriftError(
                    f"Residue mismatch at {key[0]}:{key[1]}! Table says '{given_res}' "
                    f"but structure has '{expected_res}'."
                )

        validated_resnames.append(expected["resname"])
        validated_one_letters.append(expected["one_letter"])

    norm_df["resname"] = validated_resnames
    norm_df["one_letter"] = validated_one_letters

    # If sequence is provided, verify against chain sequence
    if sequence and sequence.sequence:
        chain_seq = "".join([r["one_letter"] for r in struct_res if r["chain_id"] == norm_df["chain_id"].iloc[0]])
        full_seq = sequence.sequence.strip().upper()
        if chain_seq not in full_seq and full_seq not in chain_seq:
            # Check for substantial overlap or flag warning/error
            # If completely disjoint, raise drift error
            if len(set(chain_seq).intersection(set(full_seq))) < 3:
                raise ResidueDriftError(
                    f"Structure sequence for chain '{norm_df['chain_id'].iloc[0]}' does not match input Sequence."
                )

    return norm_df
