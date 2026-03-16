# From Algorithm Theory to Molecular Biology: A Computational Study of Sequence Alignment Algorithms

**Course:** CSCI 7432 – Algorithms and Data Structures  
**Institution:** Georgia Southern University | Department of Computer Science  
**Semester:** Spring 2026  
**Instructor:** Dr. Yao Xu  
**Author:** Yasmin Rocio Orduz Landazabal | yo00553@georgiasouthern.edu  
**Project Type:** Application  

---

## Overview

This project implements and analyzes two foundational sequence alignment algorithms — **Needleman-Wunsch** (global alignment) and **Smith-Waterman** (local alignment) — and examines how dynamic programming supports critical problems in molecular biology. Both algorithms are implemented from scratch in Python without relying on external bioinformatics libraries.

The project bridges algorithm theory and real biological applications by applying both methods to protein comparison, mutation detection, and short motif matching using real sequences from UniProt and NCBI.

---

## Algorithms

### Needleman-Wunsch (Global Alignment)
Introduced by Needleman & Wunsch (1970). Computes the optimal **global** alignment — forces a full end-to-end comparison of both sequences.
- **Time complexity:** O(mn)
- **Space complexity:** O(mn)
- **Best for:** Homologous proteins, sequences of similar length

### Smith-Waterman (Local Alignment)
Introduced by Smith & Waterman (1981). Finds the optimal **local** alignment — the highest-scoring region of similarity anywhere in both sequences.
- **Time complexity:** O(mn)
- **Space complexity:** O(mn)
- **Key difference:** Zero floor on all matrix cells — score never goes negative
- **Best for:** Motif detection, mutation detection, finding conserved regions

---

## Project Structure
```
Sequence-Alignment-Project/
├── src/
│   ├── __init__.py
│   ├── needleman_wunsch.py      # Global alignment implementation
│   ├── smith_waterman.py        # Local alignment implementation
│   └── scoring.py               # match=+1, mismatch=-1, gap=-2
├── tests/
│   ├── __init__.py
│   └── test_algorithms.py       # 22 formal validation tests
├── experiments/
│   ├── __init__.py
│   ├── biological_experiments.py  # 3 real biological experiments
│   ├── runtime_experiment.py      # O(mn) empirical validation
│   └── output/
│       ├── runtime_plot.png       # Generated complexity plot
│       └── runtime_results.csv    # Raw timing data
├── main.py                      # Demo script
└── README.md
```

---

## Scoring Scheme

| Parameter | Value |
|-----------|-------|
| Match | +1 |
| Mismatch | -1 |
| Gap penalty | -2 |

Configurable in `src/scoring.py`.

---

## Requirements

- Python 3.9 or higher
- matplotlib
```bash
pip install matplotlib
```

---

## How to Run

### 1. Clone the repository
```bash
git clone https://github.com/Yasmin-maker1/Sequence-Alignment-Project.git
cd Sequence-Alignment-Project
```

### 2. Run the demo
```bash
python main.py
```

### 3. Run the test suite
```bash
python tests/test_algorithms.py
```
Expected: **22/22 tests passing**

### 4. Run biological experiments
```bash
python experiments/biological_experiments.py
```

### 5. Run runtime complexity experiment
```bash
python experiments/runtime_experiment.py
```
Generates `experiments/output/runtime_plot.png` and `runtime_results.csv`

---

## Test Results

All 22 validation tests pass, covering:

| Group | Tests | What is validated |
|-------|-------|-------------------|
| NW-1 | 2 | Perfect match: score=4, no gaps |
| NW-2 | 1 | One mismatch: score=2 |
| NW-3 | 2 | One gap required: score=1 |
| NW-4 | 1 | All mismatches: score=-2 |
| NW-5 | 1 | Empty sequence: score=-6 |
| NW-6 | 2 | GATTACA/GCATGCU textbook example |
| NW-7 | 2 | Single character inputs |
| SW-1 | 1 | Perfect match: score=4 |
| SW-2 | 2 | Needle in haystack: score=4, correct region |
| SW-3 | 1 | Zero floor: score never negative |
| SW-4 | 1 | No similarity: score=0 |
| SW-5 | 2 | Finds buried region, ignores flanking |
| SW-6 | 1 | Aligned strings always equal length |
| SW-7 | 1 | SW score >= NW score always holds |
| SW-8 | 1 | Single match: score=1 |
| SW-9 | 1 | Single mismatch hits zero floor: score=0 |

---

## Biological Experiments Results

### Experiment 1 — Protein Comparison (Needleman-Wunsch)
**Human vs Horse Hemoglobin Alpha** (UniProt P69905 vs P01958)

| Metric | Result |
|--------|--------|
| Alignment score | 108 |
| Alignment length | 142 positions |
| Percent identity | 88.0% |
| Gaps in human sequence | 0 |
| Gaps in horse sequence | 0 |
```
Human      MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
           ||||.||||||||||.|||.|||||||||||||||.|||||||||||||||||||||.||
Horse      MVLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKAHG
```
→ 88% identity confirms evolutionary conservation of hemoglobin α across mammalian species.

---

### Experiment 2 — Mutation Detection (Smith-Waterman)
**BRCA1 Wild-Type vs 185delAG Mutant** (NCBI NM_007294)

| Metric | Result |
|--------|--------|
| Best local score | 145 |
| Aligned length | 151 positions |
| Gaps detected | 2 — marks exact deletion site |
```
WT     ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAA
       ||||||||||||||||||||||||||||||||||||  ||||||||||||||||||||||||
Mutant ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAA--TCATTAATGCTATGCAGAAA
```
→ The `--` marks exactly where the 2-base AG deletion occurred.

---

### Experiment 3 — Motif Matching (Smith-Waterman)
**TATA Box in Promoter Sequence**

| Metric | Result |
|--------|--------|
| Motif | TATAAAA |
| Score WITH TATA box | 7 (perfect 7/7) |
| Score WITHOUT TATA box | 2 |
| Match position | col=39 in promoter |

→ Score drops from 7 → 2 when TATA box is absent, confirming SW reliably detects biological motifs.

---

## Runtime Complexity Results

Both algorithms confirmed as **O(mn)** empirically.

| Sequence Length n | NW Time (s) | SW Time (s) | n² |
|-------------------|-------------|-------------|-----|
| 50 | 0.0006 | 0.0005 | 2,500 |
| 100 | 0.0022 | 0.0021 | 10,000 |
| 200 | 0.0100 | 0.0160 | 40,000 |
| 500 | 0.0798 | 0.0680 | 250,000 |
| 1000 | 0.3095 | 0.2830 | 1,000,000 |
| 1500 | 0.8510 | 0.7826 | 2,250,000 |

When plotted against n², both curves form a near-perfect straight line — confirming O(n²) behavior empirically.

---

## References

1. Needleman, S. B., & Wunsch, C. D. (1970). *Journal of Molecular Biology, 48*(3), 443–453.
2. Smith, T. F., & Waterman, M. S. (1981). *Journal of Molecular Biology, 147*(1), 195–197.
3. Altschul, S. F., et al. (1990). *Journal of Molecular Biology, 215*(3), 403–410.
4. Chao, J., Tang, F., & Xu, L. (2022). *Biomolecules, 12*(4), 546.
5. Zielezinski, A., et al. (2019). *Genome Biology, 20*, 144.
6. Petti, S., et al. (2023). *Bioinformatics, 39*(1), btac724.
7. Alberts, B., et al. (2019). *Essential Cell Biology* (5th ed.). W. W. Norton.

---

## GenAI Use Statement

Claude (Anthropic, claude.ai, Claude Sonnet 4.6) was used on February 12, 2026, for preliminary project development including reference identification and dataset exploration.