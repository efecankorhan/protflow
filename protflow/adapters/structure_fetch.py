"""AlphaFold DB / structure fetcher."""

from __future__ import annotations
import os
from pathlib import Path
from typing import Any, Dict
import pandas as pd
import requests
import gemmi

from protflow.types import Structure, ResidueScore, Provenance
from protflow.adapters.residue import AA1_TO_3


class StructureNotFoundError(RuntimeError):
    """Raised when 3D structure cannot be found in AFDB and no de novo folding engine is configured."""
    pass


def parse_alphafold_fetcher(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Fetches mmCIF and extracts residue-level pLDDT from AlphaFold DB API response.
    Supports de novo fallback when protein is missing from public database.
    """
    protein_obj = context["inputs"]["protein"]
    pid = protein_obj.protein_id
    out_dir = Path(context["output_dir"])
    params = context.get("parameters", {})

    custom_input = (
        params.get("custom_structure")
        or params.get("pdb_file")
        or params.get("structure_path")
        or getattr(protein_obj, "file_path", None)
    )

    cif_path = out_dir / f"{pid}.cif"
    cif_url = None
    mean_plddt = 90.0
    model_source = "alphafold_db"
    st_path = cif_path
    st_format = "cif"

    if custom_input and Path(custom_input).exists():
        st_path = Path(custom_input)
        st_format = "pdb" if st_path.suffix.lower() == ".pdb" else "cif"
        model_source = "custom_import"
    else:
        if isinstance(raw_data, list) and raw_data:
            entry = raw_data[0]
            cif_url = entry.get("cifUrl")
            mean_plddt = float(entry.get("globalMetricValue", 90.0))
        elif isinstance(raw_data, dict):
            if not raw_data.get("not_found"):
                cif_url = raw_data.get("cifUrl")
                mean_plddt = float(raw_data.get("globalMetricValue", 90.0))

        if cif_url and not cif_path.exists():
            resp = requests.get(cif_url, timeout=15)
            if resp.status_code == 200:
                cif_path.write_bytes(resp.content)

        # If structure is not found in public AFDB
        is_missing_from_afdb = (isinstance(raw_data, dict) and raw_data.get("not_found")) or (not cif_url)
        if is_missing_from_afdb or not cif_path.exists():
            allow_fallback = params.get("allow_synthetic_fallback", os.environ.get("PROTFLOW_ALLOW_SYNTHETIC") == "1")
            if not allow_fallback:
                raise StructureNotFoundError(
                    f"Protein '{pid}' not found in AlphaFold DB (404/400). "
                    f"Configure a de novo folding node (ColabFold/ESMFold), provide a custom PDB file, "
                    f"or pass 'allow_synthetic_fallback: true' to permit geometric mock structures."
                )

            # Explicit fallback for local simulation / test
            model_source = "synthetic_dev_mock"
            st = gemmi.Structure()
            model = gemmi.Model("1")
            chain = gemmi.Chain("A")
            seq = protein_obj.sequence or "MAG"
            for idx, aa in enumerate(seq, start=1):
                resname = AA1_TO_3.get(aa, "ALA")
                res = gemmi.Residue()
                res.name = resname
                res.seqid = gemmi.SeqId(str(idx))
                atom = gemmi.Atom()
                atom.name = "CA"
                atom.element = gemmi.Element("C")
                atom.pos = gemmi.Position(float(idx * 1.5), 0.0, 0.0)
                atom.b_iso = mean_plddt
                res.add_atom(atom)
                chain.add_residue(res)
            model.add_chain(chain)
            st.add_model(model)
            st.make_mmcif_document().write_file(str(cif_path))
            st_path = cif_path

    # Read structure to extract exact residue pLDDT scores
    st = gemmi.read_structure(str(st_path))

    score_rows = []
    for model in st:
        for chain in model:
            for res in chain:
                if res.is_water(): continue
                b_val = res[0].b_iso if len(res) > 0 else mean_plddt
                score_rows.append({
                    "protein_id": pid,
                    "chain_id": chain.name,
                    "auth_seq_id": str(res.seqid.num),
                    "resname": res.name,
                    "score_type": "plddt",
                    "score": float(b_val)
                })
        break

    score_df = pd.DataFrame(score_rows)
    struct_obj = Structure(protein_id=pid, file_path=str(st_path), format=st_format)
    residue_scores_obj = ResidueScore(protein_id=pid, df=score_df)
    if model_source in ("synthetic_dev_mock", "custom_import"):
        prov = Provenance(
            tool_name="alphafold_db",
            parameters={
                "synthetic": (model_source == "synthetic_dev_mock"),
                "model_source": model_source,
                "plddt": mean_plddt if model_source == "synthetic_dev_mock" else (score_df['score'].mean() if not score_df.empty else 85.0)
            }
        )
        struct_obj.attach_provenance(prov)
        residue_scores_obj.attach_provenance(prov)

    return {"structure": struct_obj, "residue_scores": residue_scores_obj}


