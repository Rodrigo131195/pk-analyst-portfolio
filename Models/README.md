# NONMEM models

This folder contains **control files (.mod)** only. No outputs (`.lst`, `.ext`, `.phi`, etc.) are committed.

Below is a short guide to the models you’ll find here and what each demonstrates.  
All naming follows the **run02-XXX.mod** convention from my internship; contents are sanitized.

---

## Files & intent

| File | Intent (what it shows) | Typical additions vs base |
|------|------------------------|---------------------------|
| `run02-001.mod` | **Base structural PK** model (oral dosing). Foundation for later runs. | Standard FOCE-I, inter-individual variability (ETA) on key PK params, log-normal residual (or similar). |
| `run02-012.mod` | **Nonlinear (Michaelis–Menten) elimination**. | Adds **VMAX, KM** to elimination; switches CL to concentration-dependent clearance via MM. |
| `run02-015.mod` | **Data-driven covariates via Gower clusters**. | Introduces categorical cluster effect (e.g., `ADA_G` or `CLUST`) on parameters (e.g., CL, V). |
| `run02-017.mod` | **Immunogenicity effect via Hill function** (on exposure/clearance). | Adds **Hill/Emax** relationship (e.g., EMAX, EC50, HILL) to capture ADA impact on PK. |
| `run02-016.mod` | **Delay/Transit compartments** version of 015. | Adds **transit/lag** structure (e.g., MTT/KTR/NN or ALAG) on absorption or ADA effect pathway. |
| `run02-018.mod` | **Everything combined**. | Integrates **Gower clusters + Hill ADA effect + delay/transit** (and MM if applicable) into the base. |

