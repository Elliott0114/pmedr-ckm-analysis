# Survey inputs

Obtain the August 2021–August 2023 files from the official
[NHANES data page](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2021-2023).
Place these eighteen transport files in `data/nhanes/`:

```text
ALB_CR_L.xpt  BIOPRO_L.xpt  BMX_L.xpt    BPQ_L.xpt    BPXO_L.xpt
DEMO_L.xpt    DIQ_L.xpt     FSQ_L.xpt    GHB_L.xpt    GLU_L.xpt
HDL_L.xpt     HIQ_L.xpt     HUQ_L.xpt    KIQ_U_L.xpt  MCQ_L.xpt
SMQ_L.xpt     TCHOL_L.xpt   TRIGLY_L.xpt
```

Obtain the SAS transport archive from the official
[2023 BRFSS data page](https://www.cdc.gov/brfss/annual_data/annual_2023.html).
Extract `LLCP2023.XPT` into `data/brfss/`.

The scripts use the complete national survey file, not optional BRFSS module
files. The state analysis includes the forty-eight participating states and the
District of Columbia with observed common cardiovascular-disease status.

`../data_sources.csv` records the file sizes and SHA-256 hashes of the source
copies used for verification. Compare the files if a later download produces a
different sample or estimate. Missing values and unknown survey responses must
remain distinct from negative responses.
