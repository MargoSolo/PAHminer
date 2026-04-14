#!/usr/bin/env python3
"""
PAH Gene Miner v2.0 — Evidence-based causal gene discovery pipeline
for familial pulmonary arterial hypertension.

FUNDAMENTAL CHANGE from v1:
  v1 scored genes by FREQUENCY of mention in PAH abstracts.
  v2 scores genes by STRENGTH OF CAUSAL EVIDENCE in each abstract.

  "SMAD4 is part of the BMP pathway" (pathway context) scores ZERO.
  "FBLN2 variants were identified in 3 PAH families" scores HIGH.

Context: WGS trio (GRCh38), affected mother + daughter, unaffected father,
parental consanguinity (relatedness=0.47). Standard PAH panel NEGATIVE.

Usage:
    python pah_gene_miner_v2.py
"""

import csv
import hashlib
import json
import os
import re
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
NCBI_DELAY = 0.4
S2_DELAY = 1.0
MAX_RETRIES = 3
BACKOFF_BASE = 2

# Standard PAH panels — penalised heavily
STANDARD_PANEL = frozenset([
    "BMPR2", "ACVRL1", "ENG", "SMAD9", "CAV1", "KCNK3", "TBX4",
    "ATP13A3", "GDF2", "SOX17", "BMPR1B", "ABCC8", "FOXF1",
    "KDR", "PTGIS",
])

# Known aliases mapping to panel genes
PANEL_ALIASES = {
    "ALK1": "ACVRL1", "BMP9": "GDF2", "ENDOGLIN": "ENG",
    "PROSTACYCLIN": "PTGIS",
}

# False-positive abbreviations — NOT gene symbols
FALSE_POSITIVES = frozenset([
    # Article structure
    "RESULTS", "METHODS", "CONCLUSIONS", "CONCLUSION", "BACKGROUND",
    "OBJECTIVE", "OBJECTIVES", "PURPOSE", "DESIGN", "FINDINGS",
    "INTRODUCTION", "DISCUSSION", "ABSTRACT", "STUDY", "DATA",
    "PATIENTS", "SUBJECTS", "CASES", "CONTROLS", "GROUP", "GROUPS",
    "ANALYSIS", "TABLE", "FIGURE", "SUPPLEMENTARY", "TOTAL",
    # Common English
    "DNA", "RNA", "PCR", "PAH", "WGS", "WES", "BMI", "SNP", "CNV",
    "NGS", "USA", "WHO", "FDA", "AND", "NOT", "THE", "FOR", "ARE",
    "BUT", "WAS", "HAS", "HAD", "CAN", "MAY", "LET", "SET", "RUN",
    "GOT", "PUT", "OLD", "NEW", "END", "WAY", "DAY", "USE", "HER",
    "HIS", "ITS", "OUR", "OUT", "OWN", "SAY", "SHE", "ALL", "HIM",
    "HOW", "MAN", "TWO", "AGE", "SIX", "TEN", "YES", "FEW", "GET",
    "BIG",
    # Medical organisations / scoring systems
    "NHS", "BMJ", "ATS", "ERS", "ACC", "AHA", "ESC", "IRB", "RCT",
    "ACMG", "OMIM", "HGMD", "CADD",
    # Medical / clinical abbreviations
    "OR", "CI", "HR", "SD", "IQ", "TV", "CT", "MRI", "ECG", "ICU",
    "ED", "IV", "BP", "ER", "MI", "PE", "VQ", "RV", "LV", "RA",
    "LA", "PH", "PVR", "SVR", "CO", "SV", "PAWP", "PCWP", "RAP",
    "PAP", "MPAP", "LVEF", "RVEF", "RHC", "NYHA",
    "PVD", "CHD", "CTEPH", "PPH", "IPAH", "HPAH", "APAH", "PVOD",
    "PCH", "ERA", "PDE", "NO", "ET", "PGI", "DL", "VS", "NS",
    "SEM", "IQR", "ROC", "AUC", "PPV", "NPV", "FDR", "LOD",
    "COPD", "ARDS", "CHF", "CKD", "ESRD", "IBD", "SLE",
    "ADHD", "ASD", "PTSD", "OCD", "CAD", "PPHN", "BPD", "RDS",
    "NEC", "PDA", "VSD", "IPF", "SSC", "GERD", "NASH", "NAFLD",
    "HCC", "RCC", "NSCLC", "SCLC", "ILD", "COVID", "SARS",
    # Bioinformatics / method terms
    "QTL", "GWAS", "SNV", "MAF", "VUS", "UTR", "CDS",
    "GOF", "LOF", "HET", "HOM", "REF", "ALT", "VAF",
    "IGV", "BAM", "VCF", "BED", "GTF", "GFF",
    "RPKM", "FPKM", "TPM", "CPM", "DEG", "KEGG",
    "WGCNA", "LASSO", "CRISPR", "TALEN", "SKAT", "BURDEN",
    "META", "REVEL", "SIFT", "PRESSO", "STROBE", "PRISMA", "MOOSE",
    "STRING", "DAVID", "GSEA", "ANNOVAR", "GATK", "STAR",
    "HISAT", "SALMON", "PICARD", "PLINK", "GCTA", "METAL",
    "MAGMA", "PASCAL", "VEGAS", "BOLT", "SAIGE", "LDSC",
    "COLOC", "TWAS", "FUSION", "FOCUS", "FINEMAP",
    # Lab / biology abbreviations (not gene symbols)
    "PASMC", "PAEC", "HUVEC", "HMEC", "VSMC", "SMC",
    "ROS", "ATP", "GTP", "NADPH", "CAMP", "CGMP",
    "EMT", "ECM", "ELISA", "FACS", "FISH", "ISH", "IHC",
    "BMD", "GFR", "AST", "BNP", "CRP", "ESR", "LDH",
    "ACE", "ACE2",
    # Pathway family names (not individual gene symbols)
    "BMP", "TGF", "TGFB", "VEGF", "PDGF", "FGF", "HIF",
    "NF", "PI3K", "MAPK", "MTOR", "AKT", "RAS", "JAK", "STAT",
    "WNT", "NOTCH", "HEDGEHOG", "HIPPO", "DELTA",
    "TNF", "IL",
    # Elements / small molecules
    "NO", "CO", "CA", "NA", "MG", "FE", "ZN", "CU", "PH",
    # Immune cell labels
    "CD4", "CD8", "NK", "DC",
    # Drugs / compounds (not genes)
    "SU5416", "SUGEN", "MCT", "DMSO", "PDTC",
    # Clinical measurements / indices
    "RVSP", "PVRI", "TAPSE", "PAAT", "RVEDP", "FEV1",
    # Lab techniques
    "MLPA", "EMSA", "CHIP", "TUNEL",
    # Organisations / consortia / studies
    "NIHR", "CSTAR", "BRIDGE",
    # Bioinformatics tools
    "CIBERSORT", "DESEQ2", "DESEQ", "LIMMA", "EDGER",
    # Abbreviations / disease subtypes (not gene symbols)
    "FPAH", "MPAH", "CPAH", "DPAH",
    "ACDMPV", "GLUT1", "VEGFR2", "SARS2",
    "LADD", "VATER", "CHARGE",
    # Viral vectors (not genes)
    "AAV1", "AAV2", "AAV5", "AAV8", "AAV9",
    # Generic English words that pass length filters
    "ASSOCIATED", "VARIANTS", "VARIANT", "GENETIC", "CLINICAL",
    "IDENTIFY", "BETWEEN", "THROUGH", "WITHOUT",
    "HOWEVER", "SEVERAL", "DURING", "BEFORE",
    "AFTER", "ABOUT", "ABOVE", "BELOW",
    "WITHIN", "AMONG", "AGAINST", "AROUND",
    "TOWARD", "BEYOND", "ACROSS", "ALONG",
    "DISEASE", "SEVERE", "SIGNIFICANT", "ANALYSIS",
    "SAMPLE", "RESULT", "REPORT", "SHOWED",
    "FUNCTION", "REDUCED", "INDUCED", "TREATED",
    "CONTROL", "NORMAL", "INCREASED", "DECREASED",
    "PROTEIN", "EXPRESSION", "PATHWAY", "RECEPTOR",
    "VASCULAR", "ARTERIAL", "CARDIAC", "VENOUS",
    "TISSUE", "PLASMA", "SERUM", "BLOOD",
    "TREATMENT", "THERAPY", "RESPONSE", "OUTCOME",
    "POPULATION", "COHORT",
    # Gene family abbreviations used generically
    "BMPR", "TGFBR", "EPHB", "VEGFR", "PPAR", "ERK1",
    # HIF1 and HIF2 are complexes, actual genes are HIF1A, EPAS1
    "HIF1", "HIF2", "PHD2",  # actual gene symbol is EGLN1
])

# ---------------------------------------------------------------------------
# Gene symbol validation
# ---------------------------------------------------------------------------

GENE_PATTERN = re.compile(r'\b([A-Z][A-Z0-9]{1,9})\b')

# Suffixes common in real gene symbols
GENE_SUFFIXES = re.compile(r'^[A-Z]{2,6}\d{1,3}[A-Z]?$')
# Pattern: prefix + number (BMPR2, KCNK3, ATP13A3, NOTCH3, etc.)
GENE_LETTER_DIGIT = re.compile(r'^[A-Z]{2,7}\d{1,4}[A-Z]?\d?$')


def is_valid_gene_symbol(sym):
    """Strict HGNC symbol validation. Rejects pathway names and abbreviations."""
    if sym in FALSE_POSITIVES:
        return False
    if sym in STANDARD_PANEL:
        return True  # valid, but will be penalised in scoring
    if sym in PANEL_ALIASES:
        return True
    if len(sym) < 2 or len(sym) > 10:
        return False

    # Letters+digits pattern (most real gene symbols): BMPR2, KCNA5, TRPC6
    if GENE_LETTER_DIGIT.match(sym):
        return True

    # Pure letters: only accept if 4+ chars and not common words
    if re.match(r'^[A-Z]+$', sym):
        if len(sym) < 4:
            return False
        # Reject common English/medical words
        common = {
            "MICE", "RISK", "RARE", "TYPE", "CELL", "GENE", "MALE",
            "DOSE", "DRUG", "MILD", "LUNG", "BONE", "SKIN", "IRON",
            "LEAD", "MEAN", "HIGH", "LOSS", "ROLE", "LEFT", "LIKE",
            "MANY", "MORE", "MOST", "MUCH", "NEED", "ONLY", "PART",
            "PLAY", "POOR", "RATE", "SAME", "SUCH", "TAKE", "THAN",
            "THAT", "THEM", "THEN", "THIS", "THUS", "TIME", "TRUE",
            "USED", "VERY", "WELL", "WHEN", "WITH", "WORK", "YEAR",
            "ALSO", "BOTH", "EACH", "FROM", "HAVE", "HERE", "INTO",
            "JUST", "KNOW", "LONG", "MADE", "NEXT", "ONCE", "OVER",
            "SEEN", "SHOW", "SOME", "WILL", "BEEN", "CASE", "DONE",
            "DOWN", "EVEN", "EVER", "FORM", "FULL", "GAVE", "GOOD",
            "HALF", "HAND", "HEAD", "HELP", "HOME", "IDEA", "KEEP",
            "KIND", "LAST", "LINE", "LIVE", "LOOK", "MAKE", "MARK",
            "NAME", "NOTE", "OPEN", "PASS", "REST", "RISE", "SAID",
            "SIDE", "SIZE", "STEP", "STOP", "SURE", "TEST", "TURN",
            "UPON", "VIEW", "WANT", "WIDE", "WORD", "AREA", "BACK",
            "BASE", "BEST", "BODY", "BORN", "CAME", "COME", "COPY",
            "DEEP", "DOES", "DREW", "DROP", "EASY", "FACE", "FACT",
            "FALL", "FEEL", "FIND", "FINE", "FIRE", "FIVE", "FLAT",
            "FOUR", "FREE", "GREW", "GROW", "HELD", "HOLD", "HOLE",
            "HOUR", "HUGE", "HUMAN", "MAJOR", "NOVEL", "PANEL",
            "CASES", "GENES", "BASED", "EARLY", "FIRST", "FOUND",
            "LEVEL", "LOWER", "GROUP", "THREE", "THOSE", "UNDER",
            "USING", "WHICH", "WHILE", "OTHER", "RIGHT", "SMALL",
            "STUDY", "THEIR", "THESE", "TOTAL", "WHERE", "WHOLE",
            "YEARS", "YOUNG", "AFTER", "AMONG", "CAUSE", "CLASS",
            "CLEAR", "COULD", "DEATH", "EIGHT", "EVERY", "GIVEN",
            "GREAT", "HEART", "KNOWN", "LARGE", "LATER", "MIGHT",
            "NEVER", "OFTEN", "OLDER", "POINT", "PRIOR", "RANGE",
            "SCORE", "SHALL", "SHORT", "SINCE", "STILL", "THERE",
            "TODAY", "VALUE", "WOMEN",
            # Biological/medical multi-letter terms
            "SMAD", "GATA", "FOXF", "RUNX",  # prefixes, not gene symbols alone
            # More English / medical words
            "THESE", "THOSE", "THEIR", "THERE", "WHERE", "WHICH",
            "WHILE", "OTHER", "MIGHT", "COULD", "WOULD", "SHOULD",
            "BEING", "STILL", "SINCE", "UNTIL", "GIVEN", "GREAT",
            "FIRST", "EARLY", "LATER", "OFTEN", "POINT", "RANGE",
            "YOUNG", "OLDER", "LOWER", "UPPER", "MAJOR", "MINOR",
            "NORTH", "SOUTH", "EAST", "WEST", "WORLD", "STATE",
            "LEVEL", "POWER", "UNDER", "ABOVE", "BELOW", "CAUSE",
            "DEATH", "APPLY", "BEGIN", "BUILD", "CARRY", "CHECK",
            "CLASS", "CLEAR", "CLOSE", "COUNT", "COVER", "CROSS",
            "DRIVE", "ENTER", "EQUAL", "EVENT", "EXACT", "EXIST",
            "EXTRA", "FIELD", "FINAL", "FLOOR", "FORCE", "FRONT",
            "GREEN", "HEAVY", "INNER", "JOINT", "JUDGE", "KNOWN",
            "LARGE", "LAYER", "LEARN", "LIGHT", "LIMIT", "MATCH",
            "MODEL", "MONTH", "MOUTH", "NEVER", "NIGHT", "OCCUR",
            "ORDER", "OUTER", "PAPER", "PANEL", "PHASE", "PLACE",
            "PLAIN", "PLANT", "PLATE", "PROVE", "QUICK", "QUITE",
            "RATIO", "REACH", "RAPID", "READY", "RIGHT", "ROUND",
            "ROUGH", "SCALE", "SCENE", "SCORE", "SERVE", "SHAPE",
            "SHARE", "SHIFT", "SIGHT", "SINCE", "SIXTH", "SLEEP",
            "SMALL", "SMOKE", "SOLID", "SOLVE", "SOUND", "SPACE",
            "SPEAK", "SPEED", "SPEND", "STAGE", "STAND", "START",
            "STEEP", "STOCK", "STORE", "STORM", "STORY", "STRIP",
            "STUFF", "STYLE", "SUGAR", "SWEET", "TEETH", "THICK",
            "THIRD", "TIGHT", "TITLE", "TOTAL", "TOUCH", "TRACK",
            "TRADE", "TRAIN", "TREAT", "TRIAL", "TRUST", "TRUTH",
            "TWICE", "TYPES", "UNITY", "USUAL", "VALID", "VALUE",
            "VIRUS", "VITAL", "VOICE", "WASTE", "WATCH", "WATER",
            "WHEEL", "WHITE", "WOMEN", "WORSE", "WORST", "WORTH",
            "WRITE", "WRONG", "YIELD",
            "ACUTE", "ADULT", "BASIC", "BRIEF", "CHIEF", "DAILY",
            "DENSE", "EARLY", "FALSE", "FATAL", "FIXED", "FOCAL",
            "GROSS", "HUMAN", "IDEAL", "LOCAL", "NASAL", "NOVEL",
            "OBESE", "ONSET", "ORGAN", "PRONE", "RENAL", "SAFER",
            "TOXIC", "TUMOR", "CELLS", "GENES", "LUNGS", "LIVER",
            "HEART", "BRAIN", "ASSAY", "BASED", "FOUND", "USING",
            "STUDY", "GROUP", "AFTER", "AMONG",
        }
        if sym in common:
            return False
        return True  # 4+ letter uppercase, not in common words — possible gene

    return False


def resolve_alias(sym):
    """Map known aliases to canonical symbol, flag panel genes."""
    if sym in PANEL_ALIASES:
        return PANEL_ALIASES[sym], True  # alias, is_panel
    if sym in STANDARD_PANEL:
        return sym, True
    return sym, False


# ---------------------------------------------------------------------------
# Evidence classification NLP
# ---------------------------------------------------------------------------

# Regex patterns for causal evidence types

# HUMAN_CAUSAL: variants/mutations found in patients/families
HUMAN_CAUSAL_PATTERNS = [
    # "identified GENE as a novel gene for PAH"
    r'identified\s+{gene}\s+as\s+a?\s*(?:novel|new|candidate|causal)',
    # "GENE variants/mutations were found/identified in N PAH patients/families"
    r'{gene}\s+(?:variant|mutation|pathogenic|deleterious|truncating|missense|frameshift|splice|loss.of.function)s?\s+(?:were|was|are|is)?\s*(?:found|identified|detected|observed|reported)',
    r'(?:variant|mutation|pathogenic|deleterious)s?\s+(?:in|of)\s+{gene}\s+(?:were|was|are|is)?\s*(?:found|identified|detected|reported)',
    # "N families/patients/probands with GENE variants/mutations"
    r'\d+\s+(?:famil|patient|proband|case|individual|carrier|subject).*?{gene}\s+(?:variant|mutation)',
    r'{gene}\s+(?:variant|mutation).*?\d+\s+(?:famil|patient|proband|case)',
    # "familial PAH caused by GENE"
    r'(?:familial|heritable|hereditary|inherited).*?(?:PAH|pulmonary\s+arterial\s+hypertension).*?(?:caused?\s+by|due\s+to|attributed\s+to).*?{gene}',
    r'{gene}.*?(?:caused?|causing|responsible\s+for).*?(?:PAH|pulmonary\s+arterial\s+hypertension)',
    # "loss-of-function GENE" / "haploinsufficiency of GENE"
    r'(?:loss.of.function|haploinsufficienc|biallelic|homozygous|compound\s+heterozygous|de\s+novo)\s+(?:variant|mutation|of|in)?\s*{gene}',
    # "whole exome/genome sequencing identified/revealed GENE"
    r'(?:whole\s+)?(?:exome|genome)\s+sequencing\s+(?:identified|revealed|uncovered|detected).*?{gene}',
    r'{gene}.*?(?:identified|revealed|uncovered|detected)\s+(?:by|through|via|using)\s+(?:whole\s+)?(?:exome|genome)\s+sequencing',
    # "GENE as a new/novel risk/causal gene for PAH"
    r'{gene}\s+as\s+a\s+(?:new|novel)\s+(?:risk|causal|susceptibility|candidate)\s+gene',
    # "pathogenic variant in GENE"
    r'pathogenic\s+variant\s+in\s+{gene}',
    # "GENE was the [only/most] gene harboring..."
    r'{gene}\s+(?:was|is)\s+(?:the\s+)?(?:only|most\s+likely|candidate)\s+gene',
    # "rare variant in GENE ... PAH"
    r'rare\s+(?:variant|mutation)s?\s+in\s+{gene}.*?(?:PAH|pulmonary)',
    # "segregated with PAH/disease in the family"
    r'{gene}.*?segregat.*?(?:PAH|disease|phenotype|affected)',
]

# ANIMAL_MODEL: mouse/rat models with PH phenotype
ANIMAL_CAUSAL_PATTERNS = [
    r'{gene}\s+(?:knockout|knock.out|KO|null|deficien|haploinsufficien|conditional\s+deletion)',
    r'(?:knockout|knock.out|KO|null|deficien|deletion|ablation)\s+(?:of\s+)?{gene}',
    r'{gene}.*?(?:mice|mouse|rat|zebrafish).*?(?:develop|exhibit|show|display).*?(?:pulmonary\s+hypertension|elevated\s+RVSP|right\s+ventricular)',
    r'(?:mice|mouse|rat).*?(?:lacking|without|deficient\s+in)\s+{gene}.*?(?:pulmonary|vascular)',
    r'{gene}.*?(?:transgenic|overexpressing|knockin|knock.in).*?(?:pulmonary|vascular|PAH)',
]

# HUMAN_ASSOCIATION: GWAS / statistical association
ASSOCIATION_PATTERNS = [
    r'{gene}.*?(?:GWAS|genome.wide\s+association|association\s+study|association\s+signal)',
    r'{gene}.*?(?:associated?\s+with|risk\s+(?:locus|variant|allele|factor)).*?(?:PAH|pulmonary)',
    r'(?:locus|SNP|variant)\s+(?:near|at|in)\s+{gene}.*?(?:associated|linked|correlated)',
    r'{gene}.*?(?:odds\s+ratio|relative\s+risk|hazard\s+ratio|p\s*[<=]\s*\d)',
]

# PATHWAY: context-only (no patient/animal data) — used to SUPPRESS score
PATHWAY_PATTERNS = [
    r'{gene}\s+(?:signaling|pathway|cascade|axis)',
    r'(?:signaling|pathway|cascade|axis)\s+(?:of|through|via|involving)\s+{gene}',
    r'{gene}\s+(?:is|are|was|were)\s+(?:a\s+)?(?:member|component|part|downstream|upstream|mediator|effector|regulator)',
    r'(?:BMP|TGF|Wnt|Notch|VEGF|PDGF|Hedgehog).{0,30}{gene}.{0,30}(?:signal|pathway)',
    r'(?:we\s+used|as\s+a\s+control|internal\s+control|loading\s+control|reference\s+gene|housekeeping).*?{gene}',
    r'unlike\s+{gene}',
    r'(?:similar\s+to|compared\s+(?:to|with))\s+{gene}',
]


def _match_patterns(patterns, gene, text):
    """Check if any pattern matches in text. Returns first matching sentence."""
    text_lower = text.lower()
    gene_esc = re.escape(gene.lower())
    # Split into sentences
    sentences = re.split(r'(?<=[.!?])\s+', text)

    for pat_template in patterns:
        pat = pat_template.replace(r'{gene}', gene_esc)
        try:
            if re.search(pat, text_lower):
                # Find the sentence containing this match
                for sent in sentences:
                    if re.search(pat, sent.lower()):
                        return sent.strip()[:250]
                # Fallback: return first sentence with gene name
                for sent in sentences:
                    if gene.lower() in sent.lower():
                        return sent.strip()[:250]
                return True
        except re.error:
            continue
    return None


def classify_evidence(gene, abstract, title):
    """
    Classify the type of evidence a paper provides for a gene.

    Returns (evidence_type, causal_score, best_sentence).
    """
    full_text = abstract + " " + title
    text_lower = full_text.lower()

    # Check if gene is actually mentioned in a meaningful way
    if not re.search(r'\b' + re.escape(gene) + r'\b', full_text):
        return None, 0, ""

    # Detect if paper is a review
    is_review = bool(re.search(
        r'\breview\b|\bmeta.analysis\b|\bsystematic\s+review\b|\boverview\b',
        (title + " " + abstract[:200]).lower()
    ))

    # Try HUMAN_CAUSAL first (strongest evidence)
    match = _match_patterns(HUMAN_CAUSAL_PATTERNS, gene, full_text)
    if match:
        sent = match if isinstance(match, str) else _extract_gene_sentence(full_text, gene)
        return "HUMAN_CAUSAL", _compute_causal_score(gene, abstract, title, "HUMAN_CAUSAL", is_review), sent

    # Try ANIMAL_MODEL
    match = _match_patterns(ANIMAL_CAUSAL_PATTERNS, gene, full_text)
    if match:
        sent = match if isinstance(match, str) else _extract_gene_sentence(full_text, gene)
        return "ANIMAL_MODEL", _compute_causal_score(gene, abstract, title, "ANIMAL_MODEL", is_review), sent

    # Try ASSOCIATION
    match = _match_patterns(ASSOCIATION_PATTERNS, gene, full_text)
    if match:
        sent = match if isinstance(match, str) else _extract_gene_sentence(full_text, gene)
        return "HUMAN_ASSOCIATION", _compute_causal_score(gene, abstract, title, "HUMAN_ASSOCIATION", is_review), sent

    # Check if only pathway context
    match = _match_patterns(PATHWAY_PATTERNS, gene, full_text)
    if match:
        sent = match if isinstance(match, str) else _extract_gene_sentence(full_text, gene)
        return "PATHWAY", 0, sent

    # Gene mentioned but no clear evidence pattern — check for soft causal hints
    if _has_soft_causal_hints(gene, text_lower):
        sent = _extract_gene_sentence(full_text, gene)
        return "HUMAN_CAUSAL", _compute_causal_score(gene, abstract, title, "HUMAN_CAUSAL", is_review) - 2, sent

    # Default: pathway context (mentioned but no evidence)
    sent = _extract_gene_sentence(full_text, gene)
    if is_review:
        return "REVIEW", 0, sent
    return "PATHWAY", 0, sent


def _has_soft_causal_hints(gene, text_lower):
    """Detect softer causal language that doesn't match strict patterns."""
    gene_lower = gene.lower()
    # Gene near words like "causal", "pathogenic", "disease-causing"
    # within 80 chars
    for keyword in ["causal", "pathogenic", "disease.causing", "causative",
                     "disease gene", "susceptibility gene"]:
        idx_gene = text_lower.find(gene_lower)
        idx_kw = text_lower.find(keyword)
        if idx_gene >= 0 and idx_kw >= 0 and abs(idx_gene - idx_kw) < 100:
            return True
    return False


def _compute_causal_score(gene, abstract, title, evidence_type, is_review):
    """Compute per-paper causal score for a gene."""
    score = 0
    text_lower = (abstract + " " + title).lower()

    # +5: explicit causal language
    causal_words = ["causal", "pathogenic", "disease-causing", "disease causing",
                    "causative", "loss-of-function", "loss of function"]
    if any(w in text_lower for w in causal_words):
        score += 5

    # +3: reports human patients/families with variants
    patient_words = ["patient", "proband", "family", "families", "carrier",
                     "index case", "affected individual", "cohort"]
    variant_words = ["variant", "mutation", "deletion", "insertion",
                     "frameshift", "missense", "nonsense", "splice",
                     "truncating", "rare variant"]
    has_patients = any(w in text_lower for w in patient_words)
    has_variants = any(w in text_lower for w in variant_words)
    if has_patients and has_variants:
        score += 3

    # +3: animal model with PH phenotype
    animal_words = ["knockout", "knock-out", "null mice", "ko mice",
                    "mouse model", "transgenic", "zebrafish", "rat model",
                    "conditional deletion"]
    ph_words = ["pulmonary hypertension", "elevated rvsp",
                "right ventricular", "pulmonary vascular"]
    has_animal = any(w in text_lower for w in animal_words)
    has_ph = any(w in text_lower for w in ph_words)
    if has_animal and has_ph and evidence_type == "ANIMAL_MODEL":
        score += 3

    # +2: original research (not review)
    if not is_review:
        score += 2

    # +2: gene not in any known panel
    _, is_panel = resolve_alias(gene)
    if not is_panel:
        score += 2

    # +1: specifically "pulmonary arterial hypertension" (not just PH)
    if "pulmonary arterial hypertension" in text_lower:
        score += 1

    # review discount
    if is_review:
        score = max(0, score // 2)

    return score


def _extract_gene_sentence(text, gene):
    """Extract the best sentence containing the gene name."""
    sentences = re.split(r'(?<=[.!?])\s+', text)
    gene_sentences = [s for s in sentences if re.search(r'\b' + re.escape(gene) + r'\b', s)]

    if not gene_sentences:
        return ""

    # Prefer sentences with causal/evidence words
    evidence_words = ["identified", "causal", "pathogenic", "novel", "variant",
                      "mutation", "family", "patient", "proband", "knockout",
                      "loss-of-function", "discovered", "sequencing"]
    best = gene_sentences[0]
    best_score = 0
    for sent in gene_sentences:
        s_lower = sent.lower()
        sc = sum(1 for w in evidence_words if w in s_lower)
        if sc > best_score:
            best = sent
            best_score = sc

    return best.strip()[:250]


# ---------------------------------------------------------------------------
# HTTP / Caching utilities (same as v1)
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
        print("[WARN] NCBI_API_KEY not set")
    if not s2_key:
        print("[WARN] S2_API_Key not set")
    return ncbi_key, s2_key


def progress(msg, current=None, total=None):
    if current is not None and total is not None:
        pct = int(100 * current / total) if total else 0
        bar = "#" * (pct // 5) + "-" * (20 - pct // 5)
        print(f"  [{bar}] {pct:3d}% | {msg}", flush=True)
    else:
        print(f"  >> {msg}", flush=True)


def cache_key(prefix, query):
    h = hashlib.md5(query.encode()).hexdigest()[:12]
    return CACHE_DIR / f"{prefix}_{h}.json"


def load_cache(path):
    if path.exists():
        mtime = datetime.fromtimestamp(path.stat().st_mtime)
        if datetime.now() - mtime < timedelta(hours=CACHE_TTL_HOURS):
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    return None


def save_cache(path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)


def request_with_retry(url, params=None, headers=None, delay=0.4):
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
# STEP 1 — Literature search
# ---------------------------------------------------------------------------

PUBMED_QUERIES = [
    # From v1
    "pulmonary arterial hypertension novel gene variant 2022:2026[dp]",
    "familial PAH whole genome sequencing candidate gene",
    "pulmonary hypertension rare variant exome sequencing",
    "PAH genetic architecture novel locus",
    "idiopathic pulmonary arterial hypertension germline variant",
    "pulmonary vascular disease gene discovery",
    # New in v2 — targeted for causal evidence
    "pulmonary arterial hypertension[MeSH] AND genetic[tiab] AND novel[tiab] AND (2022:2026[dp])",
    "familial pulmonary hypertension exome sequencing gene",
    "PAH causal variant new gene family",
    "pulmonary hypertension knockout mouse gene",
    "Mendelian pulmonary hypertension whole genome",
]

S2_QUERIES = [
    # From v1
    "novel gene pulmonary arterial hypertension WGS WES",
    "familial pulmonary hypertension genetic cause new gene",
    "PAH unsolved cases genomic sequencing",
    # New in v2
    "novel causal gene familial pulmonary arterial hypertension sequencing",
    "new gene discovery pulmonary hypertension mouse model",
]


def search_pubmed_ids(query, api_key, max_results=MAX_PAPERS_PER_QUERY):
    cpath = cache_key("pm_ids", query)
    cached = load_cache(cpath)
    if cached:
        progress(f"PubMed IDs (cached): {len(cached)}")
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

    all_ids = list(dict.fromkeys(all_ids))
    save_cache(cpath, all_ids)
    return all_ids


def fetch_pubmed_details(pmids, api_key):
    if not pmids:
        return []

    cpath = cache_key("pm_det_v2", ",".join(sorted(pmids[:80])))
    cached = load_cache(cpath)
    if cached:
        return cached

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    papers = []

    for i in range(0, len(pmids), 100):
        batch_ids = pmids[i:i+100]
        params = {
            "db": "pubmed",
            "id": ",".join(batch_ids),
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
    papers = []
    articles = re.split(r'<PubmedArticle>', xml_text)
    for article in articles[1:]:
        pmid = _xml_tag(article, "PMID")
        title = _xml_tag(article, "ArticleTitle")
        abstract_parts = re.findall(
            r'<AbstractText[^>]*>(.*?)</AbstractText>', article, re.DOTALL
        )
        abstract = " ".join(abstract_parts)
        abstract = re.sub(r'<[^>]+>', '', abstract).strip()

        journal = _xml_tag(article, "Title")
        if not journal:
            journal = _xml_tag(article, "ISOAbbreviation")

        year_match = re.search(r'<PubDate>.*?<Year>(\d{4})</Year>', article, re.DOTALL)
        year = int(year_match.group(1)) if year_match else None
        if not year:
            md = re.search(r'<MedlineDate>(\d{4})', article)
            year = int(md.group(1)) if md else None

        if pmid and abstract and year and year >= 2019:
            papers.append({
                "pmid": pmid,
                "title": _clean(title) if title else "",
                "abstract": _clean(abstract),
                "year": year,
                "journal": _clean(journal) if journal else "",
                "citations": 0,
                "source": "pubmed",
            })
    return papers


def _xml_tag(text, tag):
    m = re.search(rf'<{tag}[^>]*>(.*?)</{tag}>', text, re.DOTALL)
    return m.group(1).strip() if m else None


def _clean(text):
    text = re.sub(r'<[^>]+>', '', text)
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def search_semantic_scholar(query, api_key, max_results=MAX_PAPERS_PER_QUERY):
    cpath = cache_key("s2", query)
    cached = load_cache(cpath)
    if cached:
        progress(f"S2 (cached): {len(cached)}")
        return cached

    base = "https://api.semanticscholar.org/graph/v1/paper/search"
    papers = []
    offset = 0
    limit = 100
    headers = {"x-api-key": api_key} if api_key else {}

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
    if not pmids:
        return {}
    cpath = cache_key("cite", ",".join(sorted(pmids[:50])))
    cached = load_cache(cpath)
    if cached:
        return cached

    counts = {}
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
# STEP 2+3 — Gene extraction + evidence classification + scoring
# ---------------------------------------------------------------------------

def extract_genes_from_text(text):
    """Extract valid HGNC-like gene symbols from text."""
    candidates = GENE_PATTERN.findall(text)
    genes = set()
    for c in candidates:
        if is_valid_gene_symbol(c):
            genes.add(c)
    return genes


def process_papers(all_papers):
    """
    For each paper, extract genes and classify evidence.
    Returns:
      gene_evidence: {gene: [list of evidence dicts]}
      paper_genes: {pmid: set of genes}
    """
    gene_evidence = defaultdict(list)
    paper_genes = {}

    total = len(all_papers)
    for pi, paper in enumerate(all_papers):
        genes = extract_genes_from_text(paper["abstract"])
        genes |= extract_genes_from_text(paper["title"])
        paper_genes[paper["pmid"]] = genes

        for gene in genes:
            canonical, is_panel = resolve_alias(gene)

            ev_type, causal_score, sentence = classify_evidence(
                gene, paper["abstract"], paper["title"]
            )
            if ev_type is None:
                continue

            # Citation bonus
            if paper["citations"] > 50:
                causal_score += 1
            # Recency bonus
            if paper["year"] >= 2023:
                causal_score += 1
            # Panel penalty
            if is_panel:
                causal_score -= 5

            gene_evidence[canonical].append({
                "pmid": paper["pmid"],
                "year": paper["year"],
                "citations": paper["citations"],
                "evidence_type": ev_type,
                "causal_score": causal_score,
                "sentence": sentence,
                "title": paper["title"],
                "is_panel": is_panel,
                "original_symbol": gene,
            })

        if (pi + 1) % 100 == 0:
            progress(f"Classified {pi+1}/{total}", pi+1, total)

    return gene_evidence, paper_genes


def aggregate_gene_scores(gene_evidence):
    """
    Aggregate per-paper evidence into gene-level scores.
    Apply final inclusion/exclusion filters.
    """
    gene_data = {}

    for gene, evidences in gene_evidence.items():
        _, is_panel = resolve_alias(gene)
        if is_panel:
            continue  # exclude panel genes entirely

        # Check inclusion criteria: at least one HUMAN_CAUSAL or ANIMAL_MODEL
        ev_types = {e["evidence_type"] for e in evidences}
        has_causal = bool(ev_types & {"HUMAN_CAUSAL", "ANIMAL_MODEL"})
        if not has_causal:
            continue  # pathway-only or association-only — skip

        # Aggregate
        total_score = sum(e["causal_score"] for e in evidences)
        if total_score <= 0:
            continue

        pmids = list(dict.fromkeys(e["pmid"] for e in evidences))
        years = [e["year"] for e in evidences]

        n_human = sum(1 for e in evidences if e["evidence_type"] == "HUMAN_CAUSAL")
        n_animal = sum(1 for e in evidences if e["evidence_type"] == "ANIMAL_MODEL")
        n_assoc = sum(1 for e in evidences if e["evidence_type"] == "HUMAN_ASSOCIATION")

        # Best evidence sentence: from the highest-scoring evidence entry
        best_ev = max(evidences, key=lambda e: e["causal_score"])
        best_sentence = best_ev["sentence"]

        # All evidence sentences
        all_sentences = []
        seen = set()
        for e in sorted(evidences, key=lambda x: x["causal_score"], reverse=True):
            if e["sentence"] and e["sentence"] not in seen:
                all_sentences.append(e["sentence"])
                seen.add(e["sentence"])

        gene_data[gene] = {
            "score": total_score,
            "n_human": n_human,
            "n_animal": n_animal,
            "n_assoc": n_assoc,
            "n_papers": len(pmids),
            "min_year": min(years),
            "max_year": max(years),
            "citation_total": sum(e["citations"] for e in evidences),
            "pmids": pmids,
            "best_sentence": best_sentence,
            "all_sentences": all_sentences[:5],
            "evidences": evidences,
        }

    return gene_data


# ---------------------------------------------------------------------------
# STEP 4 — Output
# ---------------------------------------------------------------------------

def write_candidate_genes(gene_data, output_path):
    ranked = sorted(gene_data.items(), key=lambda x: x[1]["score"], reverse=True)

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"# PAH Novel Gene Discovery v2.0 — Evidence-Based Scoring\n")
        f.write(f"# Generated: {ts}\n")
        f.write(f"# Inclusion: genes with HUMAN_CAUSAL or ANIMAL_MODEL evidence\n")
        f.write(f"# Exclusion: standard panel genes, pathway-only mentions\n")
        f.write(f"# Total novel genes with causal evidence: {len(ranked)}\n")
        f.write(f"#\n")

        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "RANK", "GENE", "CAUSAL_SCORE", "N_HUMAN_PAPERS",
            "N_ANIMAL_PAPERS", "N_ASSOCIATION_PAPERS",
            "BEST_EVIDENCE_SENTENCE", "ALL_PMIDS", "YEARS"
        ])

        for rank, (gene, info) in enumerate(ranked, 1):
            years = f"{info['min_year']}-{info['max_year']}"
            pmids_str = ",".join(info["pmids"][:20])
            writer.writerow([
                rank, gene, info["score"], info["n_human"],
                info["n_animal"], info["n_assoc"],
                info["best_sentence"], pmids_str, years
            ])

    print(f"\n  Wrote {output_path} ({len(ranked)} genes)")
    return ranked


def write_evidence_database(gene_evidence, output_path):
    rows = []
    for gene, evidences in gene_evidence.items():
        _, is_panel = resolve_alias(gene)
        for e in evidences:
            rows.append({
                "pmid": e["pmid"],
                "year": e["year"],
                "citations": e["citations"],
                "gene": gene,
                "evidence_type": e["evidence_type"],
                "causal_score": e["causal_score"],
                "sentence": e["sentence"],
                "title": e["title"],
            })

    rows.sort(key=lambda x: (-x["causal_score"], x["gene"]))

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"# PAH Evidence Database v2.0\n")
        f.write(f"# Generated: {ts}\n")
        f.write(f"# Total gene-paper evidence pairs: {len(rows)}\n")
        f.write(f"#\n")

        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "PMID", "YEAR", "CITATIONS", "GENE", "EVIDENCE_TYPE",
            "CAUSAL_SCORE", "EVIDENCE_SENTENCE", "TITLE"
        ])
        for r in rows:
            writer.writerow([
                r["pmid"], r["year"], r["citations"], r["gene"],
                r["evidence_type"], r["causal_score"],
                r["sentence"], r["title"],
            ])

    print(f"  Wrote {output_path} ({len(rows)} evidence pairs)")


def write_summary_report(gene_data, total_papers, output_path):
    ranked = sorted(gene_data.items(), key=lambda x: x[1]["score"], reverse=True)

    # Categorise genes
    known_pah_biology = {
        "SMAD1", "SMAD4", "SMAD5", "SMAD6", "BMP10",
        "NOTCH1", "NOTCH3", "KCNA5", "KLF2",
        "FLNA", "TGFB1", "TGFBR1", "TGFBR2",
        "EPHB4", "PTGIS", "KDR",
    }
    novel = []
    known_interesting = []
    for gene, info in ranked:
        if gene in known_pah_biology:
            known_interesting.append((gene, info))
        else:
            novel.append((gene, info))

    with open(output_path, "w", encoding="utf-8") as f:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write("=" * 72 + "\n")
        f.write("PAH GENE MINER v2.0 — SUMMARY REPORT\n")
        f.write(f"Generated: {ts}\n")
        f.write("=" * 72 + "\n\n")

        f.write(f"Total papers searched: {total_papers}\n")
        f.write(f"Total genes with causal evidence: {len(ranked)}\n")
        f.write(f"  - Truly novel (not in known PAH biology): {len(novel)}\n")
        f.write(f"  - Known PAH biology (interesting overlap): {len(known_interesting)}\n\n")

        f.write("-" * 72 + "\n")
        f.write("TOP 20 CANDIDATE GENES (by causal evidence score)\n")
        f.write("-" * 72 + "\n\n")

        for i, (gene, info) in enumerate(ranked[:20], 1):
            tag = " [KNOWN_PAH_BIOLOGY]" if gene in known_pah_biology else " [NOVEL]"
            f.write(f"{i:2d}. {gene}{tag}\n")
            f.write(f"    Score: {info['score']}  |  "
                    f"Human: {info['n_human']}  |  "
                    f"Animal: {info['n_animal']}  |  "
                    f"Assoc: {info['n_assoc']}  |  "
                    f"Papers: {info['n_papers']}  |  "
                    f"Years: {info['min_year']}-{info['max_year']}\n")
            f.write(f"    Best evidence: {info['best_sentence'][:200]}\n")
            f.write(f"    PMIDs: {', '.join(info['pmids'][:5])}\n\n")

        if novel:
            f.write("-" * 72 + "\n")
            f.write("TRULY NOVEL GENES (not previously in PAH context)\n")
            f.write("-" * 72 + "\n\n")
            for gene, info in novel[:30]:
                f.write(f"  {gene:12s}  score={info['score']:3d}  "
                        f"human={info['n_human']} animal={info['n_animal']}  "
                        f"years={info['min_year']}-{info['max_year']}\n")
                f.write(f"    {info['best_sentence'][:150]}\n\n")

        if known_interesting:
            f.write("-" * 72 + "\n")
            f.write("KNOWN PAH BIOLOGY GENES (with causal evidence)\n")
            f.write("-" * 72 + "\n\n")
            for gene, info in known_interesting:
                f.write(f"  {gene:12s}  score={info['score']:3d}  "
                        f"human={info['n_human']} animal={info['n_animal']}\n")
                f.write(f"    {info['best_sentence'][:150]}\n\n")

    print(f"  Wrote {output_path}")


def print_vcf_commands(ranked_genes):
    """Print annotated bcftools commands for top novel candidates."""
    # Filter to non-panel genes only
    top = []
    for gene, info in ranked_genes:
        _, is_panel = resolve_alias(gene)
        if not is_panel:
            top.append((gene, info))
        if len(top) >= 50:
            break

    if not top:
        print("\n  [!] No novel candidate genes found.")
        return

    print("\n" + "=" * 72)
    print("VCF CHECK COMMANDS — Novel PAH Candidate Genes")
    print("=" * 72)
    print()
    print("# Replace trio_raw_joint.vcf.gz with your actual VCF path")
    print("# These commands assume VEP/SnpEff annotations with gene names")
    print()

    for gene, info in top:
        ev_types = set()
        best_pmid = ""
        for e in info["evidences"]:
            ev_types.add(e["evidence_type"])
            if not best_pmid:
                best_pmid = e["pmid"]
        ev_str = "+".join(sorted(ev_types - {"PATHWAY", "REVIEW"}))
        print(f'# {gene} — Score:{info["score"]} — Evidence:{ev_str} — PMID:{best_pmid}')
        print(f'bcftools view -H trio_raw_joint.vcf.gz | grep -w "{gene}"')
        print()

    gene_names = [g for g, _ in top]
    pattern = "|".join(gene_names)
    print("# === Combined grep for ALL novel candidates at once ===")
    print(f'bcftools view -H trio_raw_joint.vcf.gz | grep -wE "{pattern}"')

    print()
    print("# === SnpSift filter (if annotated with SnpEff) ===")
    gene_list = "', '".join(gene_names)
    print(f"# SnpSift filter \"ANN[*].GENE in ['{gene_list}']\" trio_raw_joint.vcf.gz")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    start_time = time.time()
    print("=" * 72)
    print("PAH GENE MINER v2.0 — Evidence-Based Causal Gene Discovery")
    print(f"Run started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 72)

    ncbi_key, s2_key = load_env()

    # ------------------------------------------------------------------
    # STEP 1: Literature search
    # ------------------------------------------------------------------
    print("\n[STEP 1] Searching literature databases...")

    all_papers = {}

    # --- PubMed ---
    print("\n  --- PubMed ---")
    all_pubmed_ids = []
    for qi, query in enumerate(PUBMED_QUERIES):
        progress(f"Query {qi+1}/{len(PUBMED_QUERIES)}: {query[:55]}...")
        ids = search_pubmed_ids(query, ncbi_key)
        all_pubmed_ids.extend(ids)
        progress(f"  Got {len(ids)} PMIDs", qi+1, len(PUBMED_QUERIES))

    unique_pmids = list(dict.fromkeys(all_pubmed_ids))
    print(f"\n  Unique PubMed IDs: {len(unique_pmids)}")

    print("  Fetching paper details...")
    pm_papers = fetch_pubmed_details(unique_pmids, ncbi_key)
    print(f"  Retrieved {len(pm_papers)} papers (year >= 2019, with abstract)")

    for p in pm_papers:
        all_papers[p["pmid"]] = p

    print("  Fetching citation counts...")
    pm_pmids = [p["pmid"] for p in pm_papers if not p["pmid"].startswith("S2:")]
    cite_counts = fetch_citation_counts(pm_pmids, ncbi_key)
    for p in pm_papers:
        if p["pmid"] in cite_counts:
            p["citations"] = cite_counts[p["pmid"]]
            all_papers[p["pmid"]]["citations"] = cite_counts[p["pmid"]]

    # --- Semantic Scholar ---
    print("\n  --- Semantic Scholar ---")
    for qi, query in enumerate(S2_QUERIES):
        progress(f"Query {qi+1}/{len(S2_QUERIES)}: {query[:55]}...")
        s2_papers = search_semantic_scholar(query, s2_key)
        progress(f"  Got {len(s2_papers)} papers", qi+1, len(S2_QUERIES))

        for p in s2_papers:
            if p["pmid"] in all_papers:
                all_papers[p["pmid"]]["citations"] = max(
                    all_papers[p["pmid"]]["citations"], p["citations"]
                )
            else:
                all_papers[p["pmid"]] = p

    total_papers = len(all_papers)
    print(f"\n  TOTAL UNIQUE PAPERS: {total_papers}")

    # ------------------------------------------------------------------
    # STEP 2+3: Gene extraction, evidence classification, scoring
    # ------------------------------------------------------------------
    print("\n[STEP 2] Extracting genes and classifying evidence...")

    papers_list = list(all_papers.values())
    gene_evidence, paper_genes = process_papers(papers_list)

    total_symbols = len(gene_evidence)
    print(f"  Raw gene symbols found: {total_symbols}")

    print("\n[STEP 3] Aggregating scores and filtering...")
    gene_data = aggregate_gene_scores(gene_evidence)
    print(f"  Novel genes with causal/animal evidence: {len(gene_data)}")

    # ------------------------------------------------------------------
    # STEP 4: Output
    # ------------------------------------------------------------------
    print("\n[STEP 4] Writing output files...")

    genes_path = OUTPUT_DIR / "pah_novel_genes_v2.tsv"
    evidence_path = OUTPUT_DIR / "evidence_database_v2.tsv"
    summary_path = OUTPUT_DIR / "summary_report.txt"

    ranked = write_candidate_genes(gene_data, genes_path)
    write_evidence_database(gene_evidence, evidence_path)
    write_summary_report(gene_data, total_papers, summary_path)

    # ------------------------------------------------------------------
    # STEP 5: VCF commands
    # ------------------------------------------------------------------
    print_vcf_commands(ranked)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - start_time
    print("=" * 72)
    print("PIPELINE COMPLETE")
    print(f"  Time: {elapsed:.1f}s")
    print(f"  Papers: {total_papers}")
    print(f"  Novel genes with causal evidence: {len(gene_data)}")
    print(f"  Files: {genes_path}")
    print(f"         {evidence_path}")
    print(f"         {summary_path}")
    print("=" * 72)

    # Top 15 preview
    print("\nTOP 15 NOVEL PAH CANDIDATE GENES (evidence-based):")
    print(f"{'#':<4} {'GENE':<12} {'SCORE':<6} {'HUM':<5} {'ANI':<5} {'ASSOC':<6} {'YEARS':<12} {'BEST EVIDENCE'}")
    print("-" * 100)
    for rank, (gene, info) in enumerate(ranked[:15], 1):
        yrs = f"{info['min_year']}-{info['max_year']}"
        sent = info["best_sentence"][:55] + "..." if len(info["best_sentence"]) > 55 else info["best_sentence"]
        print(f"{rank:<4} {gene:<12} {info['score']:<6} {info['n_human']:<5} "
              f"{info['n_animal']:<5} {info['n_assoc']:<6} {yrs:<12} {sent}")


if __name__ == "__main__":
    main()
