#!/usr/bin/env python3
"""
PAH Gene Miner — Literature-based candidate gene discovery pipeline
for familial pulmonary arterial hypertension.

Searches PubMed (Entrez) and Semantic Scholar for recent PAH gene discovery
papers, extracts candidate gene symbols, scores them, and outputs structured
tables for downstream VCF filtering.

Context: WGS trio (GRCh38), affected mother + daughter, unaffected father,
parental consanguinity (relatedness=0.47). Standard PAH panel NEGATIVE.

Usage:
    python pah_gene_miner.py

Requires: _env. file with NCBI_API_KEY, S2_API_Key in the same directory,
          OR set environment variables directly.
"""

import csv
import hashlib
import json
import os
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
CACHE_DIR = SCRIPT_DIR / "cache"
OUTPUT_DIR = SCRIPT_DIR / "output"
ENV_FILE = SCRIPT_DIR / "_env."

CACHE_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

CACHE_TTL_HOURS = 24
MAX_PAPERS_PER_QUERY = 200
NCBI_DELAY = 0.4      # seconds between NCBI calls
S2_DELAY = 1.0         # seconds between Semantic Scholar calls
MAX_RETRIES = 3
BACKOFF_BASE = 2       # exponential backoff multiplier

# Standard PAH panel genes — these get penalised in scoring
STANDARD_PANEL = frozenset([
    "BMPR2", "ACVRL1", "ENG", "SMAD9", "CAV1", "KCNK3", "TBX4",
    "ATP13A3", "GDF2", "SOX17", "BMPR1B", "ABCC8", "FOXF1",
])

# False-positive abbreviations to exclude from gene extraction
FALSE_POSITIVES = frozenset([
    # Common English / article structure
    "DNA", "RNA", "PCR", "PAH", "WGS", "WES", "BMI", "SNP", "CNV",
    "NGS", "USA", "WHO", "FDA", "AND", "NOT", "THE", "FOR", "ARE",
    "BUT", "WAS", "HAS", "HAD", "CAN", "MAY", "LET", "SET", "RUN",
    "GOT", "PUT", "OLD", "NEW", "END", "WAY", "DAY", "USE", "HER",
    "HIS", "ITS", "OUR", "OUT", "OWN", "SAY", "SHE", "ALL", "HIM",
    "HOW", "MAN", "OUR", "TWO", "AGE", "SIX", "TEN", "YES", "FEW",
    "GET", "BIG", "IQR", "GBP", "EUR", "USD",
    "RESULTS", "METHODS", "CONCLUSIONS", "CONCLUSION", "BACKGROUND",
    "OBJECTIVE", "OBJECTIVES", "PURPOSE", "DESIGN", "FINDINGS",
    "INTRODUCTION", "DISCUSSION", "ABSTRACT", "STUDY", "DATA",
    "PATIENTS", "SUBJECTS", "CASES", "CONTROLS", "GROUP", "GROUPS",
    "ANALYSIS", "TABLE", "FIGURE", "SUPPLEMENTARY", "TOTAL",
    # Medical organisations / guidelines
    "NHS", "BMJ", "ATS", "ERS", "ACC", "AHA", "ESC", "IRB", "RCT",
    "ACMG", "OMIM", "HGMD", "CADD",
    # Medical abbreviations
    "OR", "CI", "HR", "SD", "IQ", "TV", "CT", "MRI", "ECG", "ICU",
    "ED", "IV", "BP", "ER", "MI", "PE", "VQ", "RV", "LV", "RA",
    "LA", "PH", "PVR", "SVR", "CO", "SV", "PAWP", "PCWP", "RAP",
    "PAP", "MPAP", "LVEF", "RVEF", "RHC", "WHO", "NYHA", "PAH",
    "PVD", "CHD", "CTEPH", "PPH", "IPAH", "HPAH", "APAH", "PVOD",
    "PCH", "ERA", "PDE", "NO", "ET", "PGI", "DL", "VS", "NS",
    "SEM", "IQR", "ROC", "AUC", "PPV", "NPV", "FDR", "LOD",
    "COPD", "ARDS", "CHF", "CKD", "ESRD", "IBD", "SLE", "RA",
    "ADHD", "ASD", "PTSD", "OCD", "CAD", "HFpEF", "HFrEF",
    "PPHN", "BPD", "RDS", "NEC", "PDA", "ASD", "VSD",
    # Bioinformatics / methods
    "QTL", "GWAS", "SNV", "MAF", "VUS", "UTR", "CDS",
    "GOF", "LOF", "HET", "HOM", "REF", "VAF",
    "IGV", "BAM", "VCF", "BED", "GTF", "GFF",
    "RPKM", "FPKM", "TPM", "CPM", "DEG", "KEGG",
    "WGCNA", "LASSO", "CRISPR", "TALEN",
    "SKAT", "BURDEN", "META", "REVEL", "SIFT",
    "PRESSO", "STROBE", "PRISMA", "MOOSE",
    # Lab / biology terms that aren't gene names
    "PASMC", "PAEC", "HUVEC", "HMEC", "VSMC", "SMC",
    "BMP", "TGF", "TGFB", "VEGF", "PDGF", "FGF",  # pathway names, not gene symbols
    "HIF", "NF", "PI3K", "MAPK", "MTOR", "AKT",     # pathway/protein family names
    "ROS", "ATP", "GTP", "NADPH", "CAMP", "CGMP",
    "EMT", "ECM", "PASMC",
    "CD4", "CD8", "NK", "DC",  # immune cell markers (2-letter combos)
    "ELISA", "FACS", "FISH", "ISH", "IHC", "HE",
    # Clinical measurement
    "BMD", "GFR", "ALT", "AST", "BNP", "CRP", "ESR", "LDH", "LP",
    # Disease / clinical terms that aren't gene names
    "COVID", "SARS2", "FEV1", "ILD", "IPF", "SSC", "GERD",
    "NASH", "NAFLD", "HCC", "RCC", "NSCLC", "SCLC",
    # Tools / databases / methods (not genes)
    "STRING", "DAVID", "GSEA", "ANNOVAR", "GATK", "STAR",
    "HISAT", "SALMON", "PICARD", "PLINK", "GCTA", "METAL",
    "MAGMA", "PASCAL", "VEGAS", "BOLT", "SAIGE", "LDSC",
    "COLOC", "TWAS", "FUSION", "FOCUS", "FINEMAP",
    # Pathway family names (not individual gene symbols)
    "NOTCH", "HEDGEHOG", "HIPPO", "DELTA",
    # Other
    "IL",  # interleukin prefix alone is too vague
    "TNF",  # cytokine abbreviation used generically
    "APOB", "LDLR",  # lipid genes often mentioned in passing
])

# Known HGNC gene prefixes — symbols starting with these are more likely real
HGNC_PREFIXES = [
    "BMPR", "ACVR", "SMAD", "TBX", "SOX", "KCN", "ATP", "GDF",
    "CAV", "ENG", "ALK", "BMP", "TGF", "FGF", "VEG", "PDG",
    "EGF", "IGF", "WNT", "SHH", "NOTCH", "JAG", "DLL",
    "FOXF", "FOX", "GATA", "NKX", "HIF", "EPH", "EPHB",
    "COL", "FBN", "FBL", "ELN", "LOX", "ADA", "FLNA",
    "PTGIS", "NOS", "EDN", "EDNR", "KCNA", "TRPC", "TRPV",
    "ABCC", "SLC", "CYP", "AGER", "AQP", "CBLN", "CCL",
    "CXCL", "IL", "TNF", "TGFB", "RUNX", "ETS", "KLF",
    "PPAR", "RHO", "RAC", "CDC", "RAS", "RAF", "MEK",
    "ERK", "AKT", "PIK", "MTOR", "JAK", "STAT",
    "HDAC", "SIRT", "KDM", "KAT", "DNMT", "TET",
    "MYH", "MYL", "ACTA", "ACTB", "TAGLN", "CNN",
    "PECAM", "CDH", "CLDN", "OCLN", "TJP",
    "PIEZO", "TRPM", "ANO", "CFTR",
    "RASA", "NF", "TSC", "PTEN", "TP53", "RB",
    "PTCH", "SMO", "GLI", "CTNNB", "APC",
    "TOP", "AQP", "APEL", "APLNR", "EDNRA", "EDNRB",
    "THBS", "ADAMTS", "MMP", "TIMP", "SERPIN",
    "SEMA", "NRP", "PLXN",
]

# Regex for valid HGNC-like symbols: 2-10 uppercase letters + optional digits
GENE_PATTERN = re.compile(r'\b([A-Z][A-Z0-9]{1,9})\b')

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def load_env():
    """Load API keys from _env. file into os.environ."""
    if ENV_FILE.exists():
        with open(ENV_FILE, "r") as f:
            for line in f:
                line = line.strip()
                if "=" in line and not line.startswith("#"):
                    key, val = line.split("=", 1)
                    os.environ.setdefault(key.strip(), val.strip())

    ncbi_key = os.environ.get("NCBI_API_KEY", "")
    s2_key = os.environ.get("S2_API_Key", "")
    if not ncbi_key:
        print("[WARN] NCBI_API_KEY not set — rate limit will be 3 req/s")
    if not s2_key:
        print("[WARN] S2_API_Key not set — Semantic Scholar may throttle")
    return ncbi_key, s2_key


def progress(msg, current=None, total=None):
    """Simple progress indicator."""
    if current is not None and total is not None:
        pct = int(100 * current / total) if total else 0
        bar = "#" * (pct // 5) + "-" * (20 - pct // 5)
        print(f"  [{bar}] {pct:3d}% | {msg}", flush=True)
    else:
        print(f"  >> {msg}", flush=True)


def cache_key(prefix, query):
    """Generate a deterministic cache filename."""
    h = hashlib.md5(query.encode()).hexdigest()[:12]
    return CACHE_DIR / f"{prefix}_{h}.json"


def load_cache(path):
    """Load cached data if file exists and is within TTL."""
    if path.exists():
        mtime = datetime.fromtimestamp(path.stat().st_mtime)
        if datetime.now() - mtime < timedelta(hours=CACHE_TTL_HOURS):
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    return None


def save_cache(path, data):
    """Save data to cache file."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)


def request_with_retry(url, params=None, headers=None, delay=0.4):
    """HTTP GET with retry and exponential backoff."""
    for attempt in range(MAX_RETRIES):
        try:
            time.sleep(delay)
            resp = requests.get(url, params=params, headers=headers, timeout=30)
            if resp.status_code == 429:
                wait = BACKOFF_BASE ** (attempt + 1)
                print(f"    [429] Rate limited, waiting {wait}s...")
                time.sleep(wait)
                continue
            resp.raise_for_status()
            return resp
        except requests.RequestException as e:
            wait = BACKOFF_BASE ** (attempt + 1)
            print(f"    [ERR] {e} — retry {attempt+1}/{MAX_RETRIES} in {wait}s")
            time.sleep(wait)
    return None


# ---------------------------------------------------------------------------
# STEP 1 — Literature Search
# ---------------------------------------------------------------------------

PUBMED_QUERIES = [
    "pulmonary arterial hypertension novel gene variant 2022:2026[dp]",
    "familial PAH whole genome sequencing candidate gene",
    "pulmonary hypertension rare variant exome sequencing",
    "PAH genetic architecture novel locus",
    "idiopathic pulmonary arterial hypertension germline variant",
    "pulmonary vascular disease gene discovery",
]

S2_QUERIES = [
    "novel gene pulmonary arterial hypertension WGS WES",
    "familial pulmonary hypertension genetic cause new gene",
    "PAH unsolved cases genomic sequencing",
]


def search_pubmed_ids(query, api_key, max_results=MAX_PAPERS_PER_QUERY):
    """Search PubMed and return list of PMIDs."""
    cpath = cache_key("pm_ids", query)
    cached = load_cache(cpath)
    if cached:
        progress(f"PubMed IDs (cached): {len(cached)} for query")
        return cached

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    all_ids = []
    retstart = 0
    batch = 100

    while retstart < max_results:
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": min(batch, max_results - retstart),
            "retstart": retstart,
            "retmode": "json",
            "sort": "relevance",
        }
        if api_key:
            params["api_key"] = api_key

        resp = request_with_retry(base, params=params, delay=NCBI_DELAY)
        if not resp:
            break

        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            break
        all_ids.extend(ids)
        retstart += batch

        total_available = int(data.get("esearchresult", {}).get("count", 0))
        if retstart >= total_available:
            break

    all_ids = list(dict.fromkeys(all_ids))  # deduplicate preserving order
    save_cache(cpath, all_ids)
    return all_ids


def fetch_pubmed_details(pmids, api_key):
    """Fetch title, abstract, year, journal for a batch of PMIDs."""
    if not pmids:
        return []

    cpath = cache_key("pm_det", ",".join(sorted(pmids)))
    cached = load_cache(cpath)
    if cached:
        return cached

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    papers = []

    # Fetch in batches of 100
    for i in range(0, len(pmids), 100):
        batch = pmids[i:i+100]
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "xml",
            "rettype": "abstract",
        }
        if api_key:
            params["api_key"] = api_key

        resp = request_with_retry(base, params=params, delay=NCBI_DELAY)
        if not resp:
            continue

        papers.extend(_parse_pubmed_xml(resp.text))

    save_cache(cpath, papers)
    return papers


def _parse_pubmed_xml(xml_text):
    """Lightweight XML parsing for PubMed efetch results (no lxml needed)."""
    papers = []

    # Split into individual articles
    articles = re.split(r'<PubmedArticle>', xml_text)
    for article in articles[1:]:  # skip first empty split
        pmid = _xml_tag(article, "PMID")
        title = _xml_tag(article, "ArticleTitle")
        abstract_parts = re.findall(
            r'<AbstractText[^>]*>(.*?)</AbstractText>', article, re.DOTALL
        )
        abstract = " ".join(abstract_parts)
        # Clean XML tags from abstract
        abstract = re.sub(r'<[^>]+>', '', abstract)
        abstract = abstract.strip()

        journal = _xml_tag(article, "Title")  # journal title
        if not journal:
            journal = _xml_tag(article, "ISOAbbreviation")

        # Year extraction
        year_match = re.search(
            r'<PubDate>.*?<Year>(\d{4})</Year>', article, re.DOTALL
        )
        year = int(year_match.group(1)) if year_match else None
        if not year:
            # Try MedlineDate
            md = re.search(r'<MedlineDate>(\d{4})', article)
            year = int(md.group(1)) if md else None

        if pmid and abstract and year and year >= 2019:
            papers.append({
                "pmid": pmid,
                "title": _clean(title) if title else "",
                "abstract": _clean(abstract),
                "year": year,
                "journal": _clean(journal) if journal else "",
                "citations": 0,  # will be updated later
                "source": "pubmed",
            })
    return papers


def _xml_tag(text, tag):
    """Extract first occurrence of a simple XML tag."""
    m = re.search(rf'<{tag}[^>]*>(.*?)</{tag}>', text, re.DOTALL)
    return m.group(1).strip() if m else None


def _clean(text):
    """Remove XML tags, collapse whitespace, fix encoding."""
    text = re.sub(r'<[^>]+>', '', text)
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def search_semantic_scholar(query, api_key, max_results=MAX_PAPERS_PER_QUERY):
    """Search Semantic Scholar and return paper dicts."""
    cpath = cache_key("s2", query)
    cached = load_cache(cpath)
    if cached:
        progress(f"S2 (cached): {len(cached)} for query")
        return cached

    base = "https://api.semanticscholar.org/graph/v1/paper/search"
    papers = []
    offset = 0
    limit = 100

    headers = {}
    if api_key:
        headers["x-api-key"] = api_key

    while offset < max_results:
        params = {
            "query": query,
            "offset": offset,
            "limit": min(limit, max_results - offset),
            "fields": "paperId,externalIds,title,abstract,year,venue,citationCount",
            "year": "2019-2026",
        }
        resp = request_with_retry(base, params=params, headers=headers, delay=S2_DELAY)
        if not resp:
            break

        data = resp.json()
        results = data.get("data", [])
        if not results:
            break

        for r in results:
            abstract = r.get("abstract") or ""
            year = r.get("year")
            pmid = (r.get("externalIds") or {}).get("PubMed", "")
            if abstract and year and year >= 2019:
                papers.append({
                    "pmid": str(pmid) if pmid else f"S2:{r.get('paperId', '')[:12]}",
                    "title": r.get("title", ""),
                    "abstract": abstract,
                    "year": year,
                    "journal": r.get("venue", ""),
                    "citations": r.get("citationCount", 0) or 0,
                    "source": "semantic_scholar",
                })
        offset += limit
        total = data.get("total", 0)
        if offset >= total:
            break

    save_cache(cpath, papers)
    return papers


def fetch_citation_counts(pmids, api_key):
    """Fetch citation counts from iCite (NCBI) for PubMed papers."""
    if not pmids:
        return {}

    cpath = cache_key("cite", ",".join(sorted(pmids[:50])))
    cached = load_cache(cpath)
    if cached:
        return cached

    counts = {}
    # iCite API accepts up to 200 PMIDs at once
    for i in range(0, len(pmids), 200):
        batch = pmids[i:i+200]
        url = "https://icite.od.nih.gov/api/pubs"
        params = {"pmids": ",".join(batch)}
        resp = request_with_retry(url, params=params, delay=0.5)
        if resp:
            data = resp.json()
            for pub in data.get("data", []):
                pid = str(pub.get("pmid", ""))
                counts[pid] = pub.get("citation_count", 0) or 0

    save_cache(cpath, counts)
    return counts


# ---------------------------------------------------------------------------
# STEP 2 — Gene Extraction
# ---------------------------------------------------------------------------

def looks_like_gene(symbol):
    """Heuristic: does this look like a real HGNC gene symbol?"""
    if symbol in FALSE_POSITIVES:
        return False
    if len(symbol) < 2 or len(symbol) > 10:
        return False
    # Must have at least one letter
    if not any(c.isalpha() for c in symbol):
        return False
    # Check known prefixes
    for prefix in HGNC_PREFIXES:
        if symbol.startswith(prefix):
            return True
    # Patterns like: letters + digits (e.g. ABCA3, BMP10)
    if re.match(r'^[A-Z]{2,6}\d{1,3}[A-Z]?$', symbol):
        return True
    # Pure letters 3-7 chars — possibly a gene but needs caution
    if re.match(r'^[A-Z]{3,7}$', symbol):
        # Very short pure-letter symbols are usually abbreviations, not genes
        if len(symbol) <= 3:
            return False
        # 4-letter pure-letter symbols: only if they match known prefixes (handled above)
        if len(symbol) == 4:
            return False  # too many false positives (MICE, RISK, RARE, etc.)
        return True
    return False


def extract_genes_from_abstract(abstract):
    """Extract HGNC-like gene symbols from an abstract."""
    candidates = GENE_PATTERN.findall(abstract)
    genes = set()
    for c in candidates:
        if looks_like_gene(c):
            genes.add(c)
    return genes


def extract_key_sentence(abstract, gene, max_chars=200):
    """Extract sentence containing the gene name, truncated."""
    sentences = re.split(r'(?<=[.!?])\s+', abstract)
    for sent in sentences:
        if re.search(r'\b' + re.escape(gene) + r'\b', sent):
            return sent[:max_chars]
    return ""


# ---------------------------------------------------------------------------
# STEP 3 — Prioritization scoring
# ---------------------------------------------------------------------------

def score_gene(gene, mentions):
    """
    Score a candidate gene based on its mention contexts.

    mentions: list of dicts with keys {pmid, year, citations, abstract, title}
    """
    score = 0

    for m in mentions:
        text = (m["abstract"] + " " + m["title"]).lower()

        # +3: causal / pathogenic / disease-causing
        if any(kw in text for kw in ["causal", "pathogenic", "disease-causing",
                                      "disease causing", "causative"]):
            score += 3

        # +2: familial / inherited / germline
        if any(kw in text for kw in ["familial", "inherited", "germline",
                                      "heritable", "hereditary"]):
            score += 2

        # +2: specifically "pulmonary arterial hypertension"
        if "pulmonary arterial hypertension" in text:
            score += 2

        # +1: per unique PMID
        score += 1

        # +1: recent paper (2023-2026)
        if m["year"] >= 2023:
            score += 1

        # +2: highly cited
        if m["citations"] > 20:
            score += 2

    # -3: standard panel penalty
    if gene in STANDARD_PANEL:
        score -= 3

    return score


# ---------------------------------------------------------------------------
# STEP 4 / 5 — Output
# ---------------------------------------------------------------------------

def write_candidate_genes(gene_data, output_path):
    """Write pah_candidate_genes.tsv."""
    # Sort by score descending
    ranked = sorted(gene_data.items(), key=lambda x: x[1]["score"], reverse=True)

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"# PAH Candidate Gene Mining Results\n")
        f.write(f"# Generated: {ts}\n")
        f.write(f"# Total unique genes: {len(ranked)}\n")
        f.write(f"# Standard panel genes penalised: {', '.join(sorted(STANDARD_PANEL))}\n")
        f.write(f"# Context: Familial PAH, WGS trio, GRCh38, consanguinity=0.47\n")
        f.write(f"#\n")

        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "RANK", "GENE_SYMBOL", "SCORE", "N_PAPERS", "YEARS",
            "CITATION_TOTAL", "PMIDS", "TITLES_BRIEF", "KEY_SENTENCES"
        ])

        for rank, (gene, info) in enumerate(ranked, 1):
            years = f"{info['min_year']}-{info['max_year']}"
            pmids_str = ",".join(info["pmids"][:15])  # cap at 15
            titles = " | ".join(t[:80] for t in info["titles"][:3])
            sents = " ||| ".join(info["key_sentences"][:3])
            writer.writerow([
                rank, gene, info["score"], info["n_papers"],
                years, info["citation_total"], pmids_str, titles, sents
            ])

    print(f"\n  Wrote {output_path} ({len(ranked)} genes)")
    return ranked


def write_papers_database(papers, output_path):
    """Write papers_database.tsv."""
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"# PAH Literature Database\n")
        f.write(f"# Generated: {ts}\n")
        f.write(f"# Total papers: {len(papers)}\n")
        f.write(f"#\n")

        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "PMID", "YEAR", "JOURNAL", "CITATIONS", "TITLE",
            "GENES_EXTRACTED", "ABSTRACT_SNIPPET"
        ])

        for p in sorted(papers, key=lambda x: x.get("year", 0), reverse=True):
            genes = ",".join(sorted(p.get("genes", set())))
            snippet = p["abstract"][:300].replace("\t", " ").replace("\n", " ")
            writer.writerow([
                p["pmid"], p["year"], p["journal"], p["citations"],
                p["title"], genes, snippet
            ])

    print(f"  Wrote {output_path} ({len(papers)} papers)")


def print_vcf_commands(ranked_genes, top_n=50):
    """Print bcftools commands for top N novel candidate genes (excluding panel)."""
    top = [g for g, info in ranked_genes if g not in STANDARD_PANEL][:top_n]

    print("\n" + "=" * 70)
    print("VCF CHECK COMMANDS — Top candidate genes")
    print("=" * 70)
    print("\n# Individual gene checks:")
    print("# (Replace trio_raw_joint.vcf.gz with your actual VCF path)\n")

    for gene in top[:50]:
        print(f'bcftools view -H trio_raw_joint.vcf.gz | grep -w "{gene}"')

    print("\n# Combined grep for all top genes at once:")
    pattern = "\\|".join(top[:50])
    print(f'bcftools view -H trio_raw_joint.vcf.gz | grep -wE "{"|".join(top[:50])}"')

    print("\n# Using bcftools with regions (requires gene BED file):")
    print("# bcftools view -R top_genes.bed trio_raw_joint.vcf.gz")

    print("\n# SnpSift filter by gene name (if annotated with SnpEff/VEP):")
    gene_list = "', '".join(top[:50])
    print(f"# SnpSift filter \"ANN[*].GENE in ['{gene_list}']\" trio_raw_joint.vcf.gz")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    start_time = time.time()
    print("=" * 70)
    print("PAH GENE MINER — Literature-Based Candidate Discovery")
    print(f"Run started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Load API keys
    ncbi_key, s2_key = load_env()

    # ------------------------------------------------------------------
    # STEP 1: Literature search
    # ------------------------------------------------------------------
    print("\n[STEP 1] Searching literature databases...")

    all_papers = {}  # pmid -> paper dict (deduplication)

    # --- PubMed ---
    print("\n  --- PubMed ---")
    all_pubmed_ids = []
    for qi, query in enumerate(PUBMED_QUERIES):
        progress(f"Query {qi+1}/{len(PUBMED_QUERIES)}: {query[:60]}...")
        ids = search_pubmed_ids(query, ncbi_key)
        all_pubmed_ids.extend(ids)
        progress(f"  Found {len(ids)} PMIDs", qi+1, len(PUBMED_QUERIES))

    # Deduplicate PMIDs
    unique_pmids = list(dict.fromkeys(all_pubmed_ids))
    print(f"\n  Total unique PubMed IDs: {len(unique_pmids)}")

    # Fetch details in batches
    print("  Fetching paper details...")
    pm_papers = fetch_pubmed_details(unique_pmids, ncbi_key)
    print(f"  Retrieved {len(pm_papers)} papers with abstracts (year >= 2019)")

    for p in pm_papers:
        all_papers[p["pmid"]] = p

    # Fetch citation counts
    print("  Fetching citation counts (iCite)...")
    pm_pmids = [p["pmid"] for p in pm_papers if not p["pmid"].startswith("S2:")]
    cite_counts = fetch_citation_counts(pm_pmids, ncbi_key)
    for p in pm_papers:
        if p["pmid"] in cite_counts:
            p["citations"] = cite_counts[p["pmid"]]
            all_papers[p["pmid"]]["citations"] = cite_counts[p["pmid"]]

    # --- Semantic Scholar ---
    print("\n  --- Semantic Scholar ---")
    for qi, query in enumerate(S2_QUERIES):
        progress(f"Query {qi+1}/{len(S2_QUERIES)}: {query[:60]}...")
        s2_papers = search_semantic_scholar(query, s2_key)
        progress(f"  Found {len(s2_papers)} papers", qi+1, len(S2_QUERIES))

        for p in s2_papers:
            # Avoid duplicates: if we have the PMID already, merge citation count
            if p["pmid"] in all_papers:
                existing = all_papers[p["pmid"]]
                existing["citations"] = max(existing["citations"], p["citations"])
            else:
                all_papers[p["pmid"]] = p

    total_papers = len(all_papers)
    print(f"\n  TOTAL UNIQUE PAPERS: {total_papers}")

    # ------------------------------------------------------------------
    # STEP 2: Gene extraction
    # ------------------------------------------------------------------
    print("\n[STEP 2] Extracting gene symbols from abstracts...")

    gene_mentions = defaultdict(list)  # gene -> list of mention contexts

    papers_list = list(all_papers.values())
    for pi, paper in enumerate(papers_list):
        genes = extract_genes_from_abstract(paper["abstract"])
        # Also check title
        genes |= extract_genes_from_abstract(paper["title"])
        paper["genes"] = genes

        for g in genes:
            gene_mentions[g].append({
                "pmid": paper["pmid"],
                "year": paper["year"],
                "citations": paper["citations"],
                "abstract": paper["abstract"],
                "title": paper["title"],
            })

        if (pi + 1) % 100 == 0:
            progress(f"Processed {pi+1}/{total_papers}", pi+1, total_papers)

    print(f"  Extracted {len(gene_mentions)} unique gene-like symbols")

    # ------------------------------------------------------------------
    # STEP 3: Scoring
    # ------------------------------------------------------------------
    print("\n[STEP 3] Scoring candidate genes...")

    gene_data = {}
    for gene, mentions in gene_mentions.items():
        sc = score_gene(gene, mentions)
        pmids = list(dict.fromkeys(m["pmid"] for m in mentions))
        years = [m["year"] for m in mentions]
        total_cites = sum(m["citations"] for m in mentions)
        titles = list(dict.fromkeys(m["title"] for m in mentions))

        # Key sentences
        key_sents = []
        seen_sents = set()
        for m in mentions:
            sent = extract_key_sentence(m["abstract"], gene)
            if sent and sent not in seen_sents:
                key_sents.append(sent)
                seen_sents.add(sent)
            if len(key_sents) >= 3:
                break

        gene_data[gene] = {
            "score": sc,
            "n_papers": len(pmids),
            "min_year": min(years),
            "max_year": max(years),
            "citation_total": total_cites,
            "pmids": pmids,
            "titles": titles,
            "key_sentences": key_sents,
        }

    # Filter: only genes with score > 0
    gene_data = {g: d for g, d in gene_data.items() if d["score"] > 0}
    print(f"  Genes with positive score: {len(gene_data)}")

    # ------------------------------------------------------------------
    # STEP 4: Output
    # ------------------------------------------------------------------
    print("\n[STEP 4] Writing output files...")

    genes_path = OUTPUT_DIR / "pah_candidate_genes.tsv"
    papers_path = OUTPUT_DIR / "papers_database.tsv"

    ranked = write_candidate_genes(gene_data, genes_path)
    write_papers_database(papers_list, papers_path)

    # ------------------------------------------------------------------
    # STEP 5: VCF commands
    # ------------------------------------------------------------------
    print_vcf_commands(ranked, top_n=50)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - start_time
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print(f"  Time elapsed: {elapsed:.1f}s")
    print(f"  Papers analysed: {total_papers}")
    print(f"  Candidate genes (score > 0): {len(gene_data)}")
    print(f"  Output: {genes_path}")
    print(f"  Output: {papers_path}")
    print("=" * 70)

    # Top 10 preview
    print("\nTOP 10 CANDIDATE GENES:")
    print(f"{'RANK':<5} {'GENE':<12} {'SCORE':<6} {'PAPERS':<8} {'YEARS':<12} {'CITES':<8}")
    print("-" * 55)
    for rank, (gene, info) in enumerate(ranked[:10], 1):
        yrs = f"{info['min_year']}-{info['max_year']}"
        panel = " *PANEL*" if gene in STANDARD_PANEL else ""
        print(f"{rank:<5} {gene:<12} {info['score']:<6} {info['n_papers']:<8} {yrs:<12} {info['citation_total']:<8}{panel}")


if __name__ == "__main__":
    main()
