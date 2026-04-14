# PAHminer

Evidence-based gene discovery pipeline for **Pulmonary Arterial Hypertension (PAH)**.

Mines PubMed and Semantic Scholar for causal genetic evidence beyond standard diagnostic gene panels, ranking candidates by strength of causal support rather than simple mention frequency.

## Motivation

Standard PAH gene panels (BMPR2, ACVRL1, ENG, SMAD9, TBX4, etc.) explain only ~25% of heritable cases. PAHminer systematically searches the literature for novel candidate genes with causal evidence — human genetic studies, animal models, and association data — to support variant interpretation in panel-negative families.

## How it works

```
PubMed / Semantic Scholar
        |
   structured queries (PAH + genetics/WGS/WES/variants)
        |
   abstract-level NLP scoring
        |
   evidence classification:
        - human genetic (variants in patients)
        - animal model (knockout/transgenic)
        - association (GWAS, burden tests)
        |
   composite score per gene  -->  ranked candidate list
```

**Key design decision (v2):** abstracts are scored by *strength of causal evidence*, not frequency of mention. A gene described in the context of "BMP signaling pathway" scores zero; a gene with "rare variants identified in PAH families" scores high.

## Output

| File | Description |
|------|-------------|
| `output/pah_novel_genes_v2.tsv` | Ranked novel candidate genes with evidence scores |
| `output/evidence_database_v2.tsv` | Full evidence database: gene, paper, evidence type, snippet |
| `output/papers_database.tsv` | All papers analyzed with metadata |
| `output/summary_report.txt` | Human-readable top-20 report |

## Quick start

### Prerequisites

- Python 3.8+
- `requests` library

```bash
pip install requests
```

### API keys

Create a `_env.` file in the project root:

```
NCBI_API_KEY=your_ncbi_key
S2_API_Key=your_semantic_scholar_key
GOOGLE_API_Key=your_google_key
```

NCBI API key: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

### Run

```bash
python pah_gene_miner_v2.py
```

Results are written to `output/`. Cached API responses are stored in `cache/` (24h TTL).

## Versions

- **v1** (`pah_gene_miner.py`) — frequency-based scoring; ranks genes by how often they appear in PAH abstracts.
- **v2** (`pah_gene_miner_v2.py`) — evidence-based scoring; classifies each abstract by type of causal evidence (human genetic, animal model, association) and scores accordingly. Penalizes known panel genes to surface novel candidates.

## Example results (v2)

Top novel candidates from 1,023 papers analyzed:

| Gene | Score | Human | Animal | Evidence |
|------|-------|-------|--------|----------|
| AQP1 | 72 | 4 | 2 | Rare variants in PAH cohorts |
| PDGFD | 58 | 2 | 0 | Rare deleterious variants, FDR < 0.1 |
| LDLR | 49 | 3 | 0 | Pathogenic variant in PAH patient |
| FBLN2 | 46 | 1 | 0 | Rare deleterious variants, FDR < 0.1 |
| TET2 | 45 | 4 | 0 | Somatic mutations, clonal hematopoiesis |
| FABP4 | 43 | 0 | 4 | Genetic deletion attenuates PH in mice |

## Clinical context

Designed for a WGS trio analysis (GRCh38): affected mother + daughter, unaffected father, parental consanguinity (relatedness coefficient = 0.47). Standard PAH panel negative. The pipeline generates a prioritized list of candidate genes for variant filtering in the WGS data.

## License

For research use.
