"""Genome input, NCBI fetch, codon table and gene prediction."""

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import json
import re
import numpy as np
import pandas as pd
import requests
import gemmi

from protflow.types import (
    Sequence,
    Structure,
    ResidueScore,
    Pocket,
    Alignment,
    Annotation,
    Table,
    Provenance
)
from protflow.adapters.residue import AA3_TO_1, AA1_TO_3


# --- 1. Prodigal / ORF Gene Prediction ---

CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGA':'*',
    'TGG':'W', 'TGT':'C', 'TGC':'C'
}

COMPLEMENT_MAP = str.maketrans("ACGTURYKMSWBDHVNacgturykmswbdhvn", "TGCAAYRMKSWVHDBNtgcaayrmkswvhdbn")

def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT_MAP)[::-1]

def translate_dna(dna: str, is_cds: bool = False) -> str:
    dna = dna.upper()
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        if i == 0 and is_cds and codon in ('ATG', 'GTG', 'TTG'):
            aa.append('M')
        else:
            aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def fetch_ncbi_genome(query: str, max_bp: Optional[int] = None) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """
    Fetches a genome nucleotide sequence from NCBI E-utilities by organism name
    (e.g., 'Bacillus subtilis', 'SARS coronavirus Tor2') or accession (e.g. 'GCF_000864885.1', 'NC_000913.3').
    """
    import urllib.request
    import urllib.parse
    import ssl

    query = (query or "").strip()
    if not query:
        return None, "Arama sorgusu boş olamaz."

    ctx = ssl._create_unverified_context()
    headers = {"User-Agent": "ProtFlow/1.0 (bio-pipeline-client)"}

    try:
        # 1. Check if assembly accession (GCF_ or GCA_)
        if query.upper().startswith("GCF_") or query.upper().startswith("GCA_"):
            clean_acc = query.split()[0]
            s_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={urllib.parse.quote(clean_acc)}&retmode=json"
            req = urllib.request.Request(s_url, headers=headers)
            with urllib.request.urlopen(req, context=ctx, timeout=15) as r:
                sdata = json.loads(r.read().decode("utf-8"))
                aid_list = sdata.get("esearchresult", {}).get("idlist", [])
            if not aid_list:
                return None, f"NCBI Assembly veritabanında '{clean_acc}' bulunamadı."
            aid = aid_list[0]
            l_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=assembly&db=nuccore&id={aid}&retmode=json"
            with urllib.request.urlopen(urllib.request.Request(l_url, headers=headers), context=ctx, timeout=15) as r:
                ldata = json.loads(r.read().decode("utf-8"))
                linksets = ldata.get("linksets", [{}])[0].get("linksetdbs", [])
                nuc_ids = linksets[0].get("links", []) if linksets else []
            if not nuc_ids:
                return None, f"Assembly '{clean_acc}' için nükleotit sekansı bulunamadı."
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={nuc_ids[0]}&rettype=fasta&retmode=text"

        # 2. Check if direct nucleotide accession (NC_, NZ_, NM_, NR_)
        elif any(query.upper().startswith(p) for p in ("NC_", "NZ_", "NM_", "NR_")):
            clean_acc = query.split()[0]
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={urllib.parse.quote(clean_acc)}&rettype=fasta&retmode=text"

        # 3. Organism name search
        else:
            term = f"{query}[Organism] AND biomol_genomic[prop] AND complete[prop]"
            s_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={urllib.parse.quote(term)}&retmode=json&retmax=1"
            with urllib.request.urlopen(urllib.request.Request(s_url, headers=headers), context=ctx, timeout=15) as r:
                ids = json.loads(r.read().decode("utf-8")).get("esearchresult", {}).get("idlist", [])
            if not ids:
                s_url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={urllib.parse.quote(query)}+complete+genome&retmode=json&retmax=1"
                with urllib.request.urlopen(urllib.request.Request(s_url2, headers=headers), context=ctx, timeout=15) as r:
                    ids = json.loads(r.read().decode("utf-8")).get("esearchresult", {}).get("idlist", [])
            if not ids:
                return None, f"NCBI üzerinde '{query}' için uygun genom kaydı bulunamadı."
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ids[0]}&rettype=fasta&retmode=text"

        if max_bp:
            fetch_url += f"&seq_start=1&seq_stop={max_bp}"

        with urllib.request.urlopen(urllib.request.Request(fetch_url, headers=headers), context=ctx, timeout=45) as r:
            raw_fasta = r.read().decode("utf-8")

        lines = raw_fasta.strip().split("\n")
        header = lines[0] if lines else ""
        raw_seq = "".join(lines[1:]).replace(" ", "").replace("\r", "").upper()
        clean_seq = "".join(c for c in raw_seq if c in "ACGTUN")

        if not clean_seq:
            return None, "NCBI sekans verisi boş döndü."

        gc_count = clean_seq.count("G") + clean_seq.count("C")
        gc_pct = round((gc_count / len(clean_seq)) * 100, 1)

        org_title = query
        if header.startswith(">"):
            parts = header[1:].split(" ", 1)
            if len(parts) > 1:
                org_title = parts[1].split(",")[0].strip()

        return {
            "header": header,
            "organism": org_title,
            "sequence": clean_seq,
            "length": len(clean_seq),
            "gc_content": gc_pct
        }, None
    except Exception as e:
        return None, f"NCBI bağlantı hatası: {str(e)}"


def parse_prior_input(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Adapter for prior_input root node.
    Loads host organism profile or user-provided genome DNA, producing
    'genome' (Sequence) and 'codon_table' (Table).
    """
    params = context.get("parameters", {})
    org_id = params.get("organism", "ecoli_k12")
    org_file = Path("organisms") / f"{org_id}.json"

    org_data = {}
    if org_file.exists():
        try:
            with open(org_file, "r", encoding="utf-8") as f:
                org_data = json.load(f)
        except Exception:
            pass

    dna = params.get("genome_dna") or org_data.get("genome_dna") or "ATGGCTGGGCTGGTGGCCATCCTGGAGAAATTCGGCGACGAGCTGCGTTGGGCGCCGAAAACCTGCCACTATGGCGCGGTGATCCTGTCCTAA"
    
    # Derive protein ID cleanly from params or organism name, avoiding hardcoded E. coli fallback
    pid = params.get("protein_id") or org_data.get("id")
    if not pid:
        if org_id and org_id != "ecoli_k12":
            import re
            clean_id = re.sub(r'[^A-Za-z0-9_]', '_', org_id).strip('_').upper()
            pid = clean_id or "GENOME_ORF"
        else:
            pid = "ECOLI_MG1655_FRAG1"

    desc = org_data.get("name") or params.get("organism") or "Host genome input"

    genome_seq = Sequence(protein_id=pid, sequence=dna, description=desc)

    codon_usage = org_data.get("codon_usage", {})
    rows = []
    for aa, c_map in codon_usage.items():
        for codon, freq in c_map.items():
            rows.append({"amino_acid": aa, "codon": codon, "frequency": freq})
    codon_df = pd.DataFrame(rows) if rows else pd.DataFrame([{"amino_acid": "M", "codon": "ATG", "frequency": 1.0}])
    codon_table = Table(protein_id=pid, df=codon_df)

    return {"genome": genome_seq, "codon_table": codon_table}


def parse_or_predict_genes(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Gene finding with pyrodigal (real Prodigal): both strands, real scores.
    Genomes >= 100 kb are self-trained; shorter ones use the metagenomic models.
    Every candidate is kept; nothing is ranked or picked here.
    Outputs: proteins (all candidates), genes (one readable row per gene, incl. CDS nucleotides), annotations (GFF-like).
    """
    import pyrodigal
    from protflow.types import SequenceCollection

    genome_seq_obj = context["inputs"]["genome"]
    dna = genome_seq_obj.sequence.upper()
    pid_prefix = genome_seq_obj.protein_id or "GENE"

    # shortcut: no table 4 / translation-table detection, upgrade for Mycoplasma-like hosts
    finder = pyrodigal.GeneFinder(meta=len(dna) < 100_000)
    if len(dna) >= 100_000:
        finder.train(dna.encode())
    genes = finder.find_genes(dna.encode())
    if not genes:
        raise ValueError("pyrodigal bu dizide gen bulamadi (dizi cok kisa veya kodlayici degil).")

    rows = []
    for idx, g in enumerate(genes, 1):
        prot = g.translate().rstrip("*")
        rows.append({
            "gene_id": f"{pid_prefix}_{idx:04d}",
            "start": g.begin,
            "end": g.end,
            "strand": "+" if g.strand == 1 else "-",
            "length_aa": len(prot),
            "partial": bool(g.partial_begin or g.partial_end),
            "start_codon": g.start_type,
            "rbs_motif": g.rbs_motif,
            "gc_content": round(g.gc_cont, 3),
            "confidence": round(g.confidence(), 1),
            "score": round(g.score, 2),
            "protein": prot,
            "cds_nt": g.sequence(),
        })
    df = pd.DataFrame(rows)

    proteins = SequenceCollection(
        protein_id=genome_seq_obj.protein_id,
        sequences={r["gene_id"]: Sequence(protein_id=r["gene_id"], sequence=r["protein"],
                                          description=f"{r['start']}..{r['end']} ({r['strand']})"
                                                      + (" partial" if r["partial"] else ""))
                   for r in rows},
    )
    genes_table = Table(protein_id=genome_seq_obj.protein_id, df=df)
    genes_table.provenance = Provenance(
        tool_name="pyrodigal",
        parameters={"total_genes": len(rows), "reverse_strand": int((df["strand"] == "-").sum()),
                    "partial": int(df["partial"].sum())},
    )

    annot_obj = Annotation(protein_id=genome_seq_obj.protein_id, records=[{
        "seqid": genome_seq_obj.protein_id, "source": "pyrodigal", "type": "CDS",
        "start": r["start"], "end": r["end"], "score": str(r["score"]), "strand": r["strand"],
        "phase": "0", "attributes": {"ID": r["gene_id"], "Name": r["gene_id"], "Length": str(r["length_aa"])},
    } for r in rows])

    return {"proteins": proteins, "genes": genes_table, "annotations": annot_obj}


def parse_protein_select(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """Picks ONE protein from a candidate set by id (parameter `protein_id`); a set of one needs no id."""
    coll = context["inputs"]["proteins"]
    if not coll.sequences:
        raise ValueError("Aday protein kümesi boş.")
    pid = (context.get("parameters", {}).get("protein_id") or "").strip()
    if not pid:
        if len(coll.sequences) == 1:
            return {"protein": next(iter(coll.sequences.values()))}
        raise ValueError(f"{len(coll.sequences)} aday protein var: 'protein_id' parametresiyle birini seçin.")
    if pid not in coll.sequences:
        raise ValueError(f"'{pid}' aday kümesinde yok.")
    return {"protein": coll.sequences[pid]}


