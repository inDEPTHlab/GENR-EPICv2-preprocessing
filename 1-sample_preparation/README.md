# GenR EPICv2 sample preparation

This section contains the scripts used for **sample selection** and **plate randomization** prior to EPICv2 DNA methylation profiling in the Generation R Study (funded by the HappyMums project). 

---

## Section structure

| File | Description |
|------|-------------|
| `1.1-genr_epicv2_sample_selection_report.Rmd` | Reproducible R Markdown report: eligibility filtering, sibling removal, and supplementary sample identification |
| `1.2-bridges_and_extras_selection.R` | Selection of bridge samples (for array cross-validation at birth) and supplementary *"immuno"* samples |
| `1.3-randomization.R` | Multivariate randomization of the plate layout delivered to the lab |
| `docs/index.html` | Rendered HTML output of the sample selection report |

---

## Sample selection

The goal was to identify **`3,776` samples** for EPICv2 methylation profiling: 

* `3,696` "primary" samples spanning 4 (postnatal) longitudinal time points (age 6, 10, 14, and 18 years) 
* + `80` bridges at birth.

#### Eligibility criteria

Children were included if they met **all** of the following:

- Existing DNA methylation data at birth (Illumina EPICv1 or 450K array)
- DNA (EDTA blood) available at **at least one** postnatal follow-up time point (6, 10, 14, or 18 years)
- DNA concentration ≥ 11.1 ng/µl

Children who already had 450K data at 6 and/or 10 years were excluded from re-selection at those specific time points to avoid redundancy.

#### Sibling exclusion

To avoid non-independence in downstream analyses, only one child per family was retained. When two siblings both passed eligibility, the child with the **greater number of available postnatal time points** was preferred. In the case of a tie, preference was given to the child with EPICv1 (rather than 450K) birth data; if still tied, one sibling was selected at random.

#### Supplementary samples

After DNA isolation and sex-check QC, **3,474** of the originally selected samples passed. The remaining **222** slots were filled using samples pre-isolated by the Immunology department (at 6 and 10 years) that satisfied the inclusion criteria.

#### Bridge samples

**80 samples** (from the 3,776 total) were reserved as *bridge* samples. These are birth time points that were re-profiled with EPICv2 to enable direct comparison with existing EPICv1 and 450K data. Bridge samples were selected so that all EPICv1 / 450K plates were represented **at least 2 individuals per plate**, prioritizing participants with the highest number of available follow-up time points.

---

## Plate randomization
The selected samples were distributed across **40 96-well plates** (39 full plates + 1 partially filled plate) using multivariate constrained randomization to minimize systematic confounding between key biological variables and array batch effects.

#### Method

Randomization was performed by adapting functionality from the [`Omixer`](https://bioconductor.org/packages/Omixer/) Bioconductor package. 

Briefly, we generated 1,000 candidate plate layouts and selected the one with the lowest correlation between key variables and plate/chip/position. The following variables were balanced across plates:

| Variable | Rationale |
|----------|-----------|
| `Period` (time point) | Prevent time point clustering within plates |
| `Followup` (number of time points per individual) | Balance longitudinal depth |
| `Biobank` (sample source) | Avoid batch-source confounding |
| `sex` | Distribute sex evenly |
| `age` at blood draw | Prevent age gradients across plates |
| `mom_education` | Balance socioeconomic proxies |
| `mom_smoking` | Balance prenatal exposure proxies |
| `birth_gestational_age` | Balance perinatal characteristics |
| `birth_weight` | Balance perinatal characteristics |

Samples from the same individual (i.e., multiple time points) were kept on the **same plate**, ensuring within-person longitudinal comparisons are not confounded by inter-plate variation.

## Notes

- All scripts assume a working directory containing a `/files` subdirectory with the required `.sav` input files (data management files prepared by the Generation R data management team; not included in this repository for data privacy reasons).
- Random seeds are set explicitly (`3108` for sibling selection and randomization steps) to ensure reproducibility.
- The rendered sample selection report is available at [`docs/index.html`](docs/index.html).
