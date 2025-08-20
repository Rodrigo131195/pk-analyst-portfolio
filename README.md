# pk-analyst-portfolio

Turning noisy clinical PK data into decision-ready outputs — **R-first**, with a dash of Python/Jupyter.  
Focus: **NONMEM dataset generation**, **compartmental modeling assets**, **clustering (Gower & ML)**, and clean **PK visuals**.  
> All examples use **public/synthetic data**. No confidential files are included.

---

## 🗂️ Repo structure
```
.
├─ README.md
├─ LICENSE
├─ .gitignore
├─ data/
│  ├─ public/                # tiny synthetic/public demo files (no real data)
│  └─ synthetic/             # optional extra toy data
├─ derived/                  # auto-generated outputs (figures, CSVs)
├─ models/
│  ├─ run02-001.mod     
│  ├─ run02-012.mod     
│  ├─ run02-015.mod     
│  ├─ run02-016.mod     
│  ├─ run02-017.mod     
│  └─ run02-018.mod
├─ notebooks/
│  └─ 31_nonmem_dataset_builder.ipynb
├─ scripts/
  ├─ 10_visual_inspection.R             # ADA/drug time profiles & coverage
  ├─ 20_clustering_gower.R              # Gower + PAM + silhouette/MDS
  ├─ 25_clustering_feature_selection.R  # extended features + k selection
  ├─ 30_nonmem_gof_xpose.R              # GOF (DV~PRED/IPRED) + ETA (CSV/xpose)
  └─ 40_gru_autoencoder_clustering.R    # GRU embeddings → k-means clusters

```


---

## 🔁 Reproducibility
- **R**: managed via **renv** (`renv.lock`).  
- **Python**: `environment.yml` for notebook deps.  
- **Torch**: script defaults to **CPU**; set `use_gpu = TRUE` if CUDA is available.

---

## 🔒 Confidentiality
This is a **portfolio**. All inputs are synthetic/public; any resemblance to real studies is coincidental.  
No sponsor names, IDs, dates, or raw clinical tables are included.

---


## 📜 License & citation
- **License:** MIT (see `LICENSE`).  
- **Cite:** If this repo helps you, cite “Rodrigo Muñoz, *pk-analyst-portfolio* (2025)”.

---
