"""
Core bio-data types and provenance metadata for ProtFlow.
7 primary types: Sequence, Structure, ResidueScore, Pocket, Alignment, Annotation, Table.
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import json
import pandas as pd


@dataclass
class Provenance:
    """Provenance tracking for any generated data object."""
    tool_name: str
    tool_version: str = "unknown"
    parameters: Dict[str, Any] = field(default_factory=dict)
    confidence: Optional[Union[float, Dict[str, float]]] = None
    execution_time: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    is_remote: bool = False
    endpoint: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Provenance:
        return cls(**data)


@dataclass
class BaseBioType:
    """Base class for all bio entities, keyed by protein_id."""
    protein_id: str
    provenance: Optional[Provenance] = None

    def attach_provenance(self, prov: Provenance) -> BaseBioType:
        self.provenance = prov
        return self


@dataclass
class Sequence(BaseBioType):
    """Sequence entity (FASTA representation)."""
    sequence: str = ""
    description: str = ""

    def save_fasta(self, path: Union[str, Path]) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", encoding="utf-8") as f:
            header = f">{self.protein_id} {self.description}".strip()
            f.write(f"{header}\n{self.sequence}\n")
        return p

    @classmethod
    def from_fasta(cls, path: Union[str, Path], protein_id: Optional[str] = None) -> Sequence:
        p = Path(path)
        lines = p.read_text(encoding="utf-8").strip().splitlines()
        header = lines[0][1:].strip() if lines and lines[0].startswith(">") else ""
        seq_parts = []
        for line in (lines[1:] if lines and lines[0].startswith(">") else lines):
            if line.startswith(">"):
                break
            seq_parts.append(line.strip())
        seq = "".join(seq_parts)
        pid = protein_id or (header.split()[0] if header else p.stem)
        desc = header[len(pid):].strip() if header.startswith(pid) else header
        return cls(protein_id=pid, sequence=seq, description=desc)


@dataclass
class Structure(BaseBioType):
    """3D Structure entity (mmCIF or PDB format)."""
    file_path: str = ""
    format: str = "cif"  # cif or pdb

    def to_gemmi_structure(self):
        import gemmi
        if not self.file_path:
            raise ValueError(f"Structure '{self.protein_id}' has empty file_path (likely dropped by filter).")
        p = Path(self.file_path)
        if not p.exists():
            raise FileNotFoundError(f"Structure file not found: {self.file_path}")
        if self.format.lower() in ("cif", "mmcif") or p.suffix.lower() == ".cif":
            doc = gemmi.cif.read_file(str(p))
            return gemmi.make_structure_from_block(doc.sole_block())
        return gemmi.read_structure(str(p))

    def get_chains(self) -> List[str]:
        st = self.to_gemmi_structure()
        return [ch.name for model in st for ch in model]


@dataclass
class ResidueScore(BaseBioType):
    """
    Residue-level scores (pLDDT, conservation, SASA, etc.).
    Columns: [protein_id, chain_id, auth_seq_id, resname, score_type, score]
    Stored as Parquet.
    """
    df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(
        columns=["protein_id", "chain_id", "auth_seq_id", "resname", "score_type", "score"]
    ))

    def save_parquet(self, path: Union[str, Path]) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self.df.to_parquet(p, index=False)
        return p

    @classmethod
    def from_parquet(cls, path: Union[str, Path], protein_id: Optional[str] = None) -> ResidueScore:
        p = Path(path)
        df = pd.read_parquet(p)
        pid = protein_id or (df["protein_id"].iloc[0] if not df.empty else p.stem)
        return cls(protein_id=pid, df=df)


@dataclass
class Pocket(BaseBioType):
    """
    Pocket predictions (e.g., P2Rank, fpocket).
    Columns: [protein_id, pocket_id, score, center_x, center_y, center_z, volume, residues]
    Stored as Parquet.
    """
    df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(
        columns=["protein_id", "pocket_id", "score", "center_x", "center_y", "center_z", "volume", "residues"]
    ))

    def save_parquet(self, path: Union[str, Path]) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self.df.to_parquet(p, index=False)
        return p

    @classmethod
    def from_parquet(cls, path: Union[str, Path], protein_id: Optional[str] = None) -> Pocket:
        p = Path(path)
        df = pd.read_parquet(p)
        pid = protein_id or (df["protein_id"].iloc[0] if not df.empty else p.stem)
        return cls(protein_id=pid, df=df)


@dataclass
class Alignment(BaseBioType):
    """
    Sequence/Structure alignment record.
    Stored as Parquet or tabular data.
    """
    target_id: str = ""
    query_seq: str = ""
    target_seq: str = ""
    identity: float = 0.0
    score: float = 0.0
    evalue: float = 0.0
    alignment_format: str = "tabular"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "protein_id": self.protein_id,
            "target_id": self.target_id,
            "query_seq": self.query_seq,
            "target_seq": self.target_seq,
            "identity": self.identity,
            "score": self.score,
            "evalue": self.evalue,
            "alignment_format": self.alignment_format
        }


@dataclass
class Annotation(BaseBioType):
    """
    Feature Annotation entity (GFF3 format or tabular).
    Columns: [seqid, source, type, start, end, score, strand, phase, attributes]
    """
    records: List[Dict[str, Any]] = field(default_factory=list)

    def save_gff3(self, path: Union[str, Path]) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", encoding="utf-8") as f:
            f.write("##gff-version 3\n")
            for rec in self.records:
                attrs_str = ";".join([f"{k}={v}" for k, v in rec.get("attributes", {}).items()])
                line = (
                    f"{rec.get('seqid', self.protein_id)}\t"
                    f"{rec.get('source', '.')}\t"
                    f"{rec.get('type', 'region')}\t"
                    f"{rec.get('start', 1)}\t"
                    f"{rec.get('end', 1)}\t"
                    f"{rec.get('score', '.')}\t"
                    f"{rec.get('strand', '+')}\t"
                    f"{rec.get('phase', '.')}\t"
                    f"{attrs_str}\n"
                )
                f.write(line)
        return p

    @classmethod
    def from_gff3(cls, path: Union[str, Path], protein_id: Optional[str] = None) -> Annotation:
        p = Path(path)
        records = []
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 9:
                attrs = {}
                for item in parts[8].split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs[k] = v
                records.append({
                    "seqid": parts[0],
                    "source": parts[1],
                    "type": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": parts[5],
                    "strand": parts[6],
                    "phase": parts[7],
                    "attributes": attrs
                })
        pid = protein_id or (records[0]["seqid"] if records else p.stem)
        return cls(protein_id=pid, records=records)


@dataclass
class Table(BaseBioType):
    """General tabular data representation stored as Parquet."""
    df: pd.DataFrame = field(default_factory=pd.DataFrame)

    def save_parquet(self, path: Union[str, Path]) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        self.df.to_parquet(p, index=False)
        return p

    @classmethod
    def from_parquet(cls, path: Union[str, Path], protein_id: Optional[str] = None) -> Table:
        p = Path(path)
        df = pd.read_parquet(p)
        pid = protein_id or (df["protein_id"].iloc[0] if "protein_id" in df.columns and not df.empty else p.stem)
        return cls(protein_id=pid, df=df)


@dataclass
class SequenceCollection(BaseBioType):
    """Collection of Sequence entities (e.g. multi-FASTA / proteome)."""
    sequences: Dict[str, Sequence] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.sequences)

    def __iter__(self):
        return iter(self.sequences.values())

    def __getitem__(self, item: str) -> Sequence:
        return self.sequences[item]

    def add(self, seq: Sequence) -> SequenceCollection:
        self.sequences[seq.protein_id] = seq
        return self

    @classmethod
    def from_fasta(cls, path: Union[str, Path], collection_id: Optional[str] = None) -> SequenceCollection:
        p = Path(path)
        seqs = {}
        cur_id = None
        cur_desc = ""
        cur_lines = []
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = Sequence(protein_id=cur_id, sequence="".join(cur_lines), description=cur_desc)
                parts = line[1:].strip().split(maxsplit=1)
                cur_id = parts[0] if parts else f"seq_{len(seqs)+1}"
                cur_desc = parts[1] if len(parts) > 1 else ""
                cur_lines = []
            elif line:
                cur_lines.append(line)
        if cur_id:
            seqs[cur_id] = Sequence(protein_id=cur_id, sequence="".join(cur_lines), description=cur_desc)
        return cls(protein_id=collection_id or p.stem, sequences=seqs)


@dataclass
class StructureCollection(BaseBioType):
    """Collection of 3D Structure entities."""
    structures: Dict[str, Structure] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.structures)

    def __iter__(self):
        return iter(self.structures.values())

    def __getitem__(self, item: str) -> Structure:
        return self.structures[item]

    def add(self, st: Structure) -> StructureCollection:
        self.structures[st.protein_id] = st
        return self


# Registry mapping string type names to classes
TYPE_REGISTRY = {
    "Sequence": Sequence,
    "Sequence[]": SequenceCollection,
    "SequenceCollection": SequenceCollection,
    "Structure": Structure,
    "Structure[]": StructureCollection,
    "StructureCollection": StructureCollection,
    "ResidueScore": ResidueScore,
    "Pocket": Pocket,
    "Alignment": Alignment,
    "Annotation": Annotation,
    "Table": Table,
}

