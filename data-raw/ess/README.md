# ESS source data (not tracked)

This folder is the drop point for European Social Survey (ESS) integrated data
files used to regenerate the youth-unemployment case study. **The raw ESS
microdata is license-restricted and is git-ignored on purpose** — never commit
it. Only `data-raw/fetch_ess_youth.R` (the regeneration script) and the derived,
non-identifying outputs (VPC/PCV tables, plots) are shipped.

## How to obtain the file

1. Register (free) and log in at <https://ess.sikt.no>.
2. Accept the data-use conditions.
3. Download the **integrated file** for the round(s) you want, in **CSV** or
   **Stata (.dta)** format, and place the file(s) here, e.g.:
   - `ESS11.dta`  (Round 11, 2023)
   - `ESS10.dta`, `ESS9.dta`  (optional, for pooling to thicken youth cells)

## Variables the pipeline uses

| Concept             | ESS variable(s)              |
|---------------------|------------------------------|
| Country             | `cntry`                      |
| Gender              | `gndr`                       |
| Age                 | `agea`                       |
| Main activity       | `mnactic` (+ `uempla`/`uempli`) |
| Education (highest) | `eisced`                     |
| Migration background| `brncntr`, `facntr`, `mocntr`|

Source: European Social Survey, ESS Data Portal (Sikt). Licensed CC BY-NC-SA 4.0.
