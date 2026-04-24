# festuca-transcriptome-annotation

Trinotate-based functional annotation pipeline for a *de novo* Trinity transcriptome assembly of *Festuca rubra*, including parallel Pfam/EggNOG/SignalP processing and Poaceae-curated BLAST, with downstream GO and KEGG background generation.

---

## Repository structure

```
festuca-transcriptome-annotation/
├── pipeline/                          # Numbered scripts (run in order)
│   ├── trinotate_00_setup.bash        # Setup guide + full run order (READ THIS FIRST)
│   ├── trinotate_01_longorfs_blastp.bash
│   ├── trinotate_02_parallel_pfam_eggnog.bash
│   ├── trinotate_02_pfam_array.pbs
│   ├── trinotate_03_transdecoder.bash
│   ├── trinotate_04_eggnog_array.pbs
│   ├── trinotate_05_tmhmm.bash
│   ├── trinotate_06_parallel_signalp6.bash
│   ├── trinotate_06_signalp6_array.pbs
│   ├── trinotate_07_blastp_rice.bash
│   ├── trinotate_08_blastp_wheat.bash
│   ├── trinotate_09_load_report.bash
│   ├── trinotate_10_merge_annotations.py
│   ├── trinotate_11_go_background.py
│   └── trinotate_12_kegg_background.py
└── independent/                       # Run independently; do not block the main pipeline
    ├── trinotate_IND_blastx.bash      # BLASTx against SwissProt
    └── trinotate_IND_rfam.bash        # Rfam cmscan (~200 h; can be parallelised)
```

Scripts ending in `_array.pbs` are PBS job array scripts called by their corresponding `_parallel_*.bash` dispatcher. All `.bash`/`.pbs` scripts target a PBS/Torque HPC cluster.

---

## Usage

1. **Read `pipeline/trinotate_00_setup.bash` first.** It contains the full setup guide (conda environment, database downloads, SignalP/TMHMM manual installation) and the step-by-step run order.
2. Edit the `USER CONFIGURATION` block at the top of each script to set your paths before running anything.
3. Run steps in the numbered order given in `trinotate_00_setup.bash`.

---

## Dependencies

All tools are installed via the `trinotate_env` conda environment except SignalP 6 and TMHMM 2 (see setup guide for manual installation).

| Tool | Version used | Purpose |
|------|-------------|---------|
| Trinotate | 3.2.2 | Functional annotation framework |
| TransDecoder | 5.5.0 | ORF prediction |
| Trinity | 2.8.5 | Transcriptome assembly (assembly pre-computed; used here for utility scripts) |
| DIAMOND | 2.0.15 | Fast protein alignment (SwissProt, Rice, Wheat) |
| HMMER / hmmscan | 3.4 | Pfam domain search |
| EggNOG-mapper | 2.1.13 | Orthology and functional annotation |
| SignalP | 6.0+h (fast model) | Signal peptide prediction |
| TMHMM | 2.0c | Transmembrane topology prediction |
| Infernal / cmscan | 1.1.4 | Rfam RNA family annotation |
| seqkit | 2.13.0 | FASTA splitting for parallel jobs |
| Python | 3.13.2 | Post-processing scripts (steps 10–12) |
| pandas | 2.2.3 | Tabular data processing in Python scripts |

> **EggNOG database:** emapperdb-5.0.2

---

## Citation

This repository accompanies the following manuscript (citation to be added upon publication):

> 

If you use this code before the manuscript is published, please cite this repository directly:

> Sammarco I. (2026). *festuca-transcriptome-annotation* [Software]. GitHub. https://github.com/iris-sammarco/festuca-transcriptome-annotation

---

## License

MIT License — see [LICENSE](LICENSE) for details.
