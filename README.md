# Thesis Analysis Code and ACE-Bounding Algorithm

This repository contains the R code used for the statistical analyses in the master’s thesis by **Khaled Aboul Azm**. It also includes a custom column-generation algorithm for bounding the **Average Causal Effect (ACE)** using instrumental-variable assumptions.

## Repository contents

The analysis workflow consists of the following scripts:

1. **`01_Datamanagement.R`**  
   Performs the initial data management and generates the cleaned datasets.

2. **`02_Datacheck.R`**  
   Performs data-quality and consistency checks.

3. **`03_SNP_MetaA.R`**  
   Performs SNP selection and the associated filtering procedures.

4. **`04_vitGeno_1000G.R`**  
   Processes and filters the vitamin D genetic data using the 1000 Genomes reference data.

5. **`05_Instrument_strength.R`**  
   Evaluates the strength of the proposed genetic instruments, including calculation of the relevant F-statistics.

6. **Main thesis analysis script**  
   Performs the primary statistical analyses and generates the corresponding tables and figures.

The repository also contains the implementation of the custom ACE-bounding procedure described below.

## Analysis workflow

The scripts should generally be run in numerical order:

```text
01_Datamanagement.R
        ↓
02_Datacheck.R
        ↓
03_SNP_MetaA.R
        ↓
04_vitGeno_1000G.R
        ↓
05_Instrument_strength.R
        ↓
Main thesis analysis script
```

Scripts `01_Datamanagement.R` and `02_Datacheck.R` were executed on the **BISON server**. The main analysis script therefore loads their resulting datasets directly from the project data directory.

Scripts `03_SNP_MetaA.R`, `04_vitGeno_1000G.R`, and `05_Instrument_strength.R` should normally be run before the main analysis script. These prerequisite scripts are executed separately rather than sourced automatically. This prevents the unintended repetition of data-management steps or computationally intensive genetic-data processing.

When running the analysis on the **GenEpi server**, the main analysis script is intended to be run from:

```bash
~/ii2025MRVitD/analysis/Notebooks
```

## Data availability

The data used in the thesis are not included in this repository. Access may be restricted because of privacy, governance, licensing, or institutional requirements.

The scripts may therefore not run without:

- access to the relevant project data;
- the expected directory structure;
- the required R packages and external software;
- access to the relevant institutional computing environment; and
- any necessary reference data, including 1000 Genomes data.

No individual-level or otherwise restricted research data should be committed to this repository.

## Bounding the Average Causal Effect


1. **Initialization**  
   Initialize the RMP with artificial variables assigned a large objective penalty. Solve the RMP to obtain the initial dual variables.

2. **Pricing**  
   Enumerate the possible outcome response functions \(Y^X\). For each \(Y^X\), construct \(X^{G*}\) by independently selecting, for each instrument level, the exposure with the largest relevant dual contribution.

3. **Column selection**  
   Calculate the reduced cost of each candidate response type and select the candidate with the smallest reduced cost.

4. **Column addition**  
   If the minimum reduced cost is less than \(-\epsilon\), add the corresponding column to the RMP.

5. **Iteration**  
   Resolve the RMP to obtain updated dual variables and repeat the pricing procedure.

6. **Termination**  
   Stop when the minimum reduced cost is greater than or equal to \(-\epsilon\).

The procedure is applied to the relevant optimization objectives to obtain bounds on the ACE under the assumptions and constraints encoded in the model.

## Outputs

The analysis code generates the tables, figures, model results, and intermediate objects used in the thesis.

The main analysis script consolidates code originally developed across multiple analysis notebooks. It is intended to provide a simplified way to review and reproduce the analysis workflow. Consequently, not every figure is automatically saved to disk, and some intermediate objects may be retained or written for use by other scripts.

## Software requirements

The analyses were conducted in **R**. Required package libraries are loaded within the individual scripts.

Because some analyses rely on institutional servers, genetic reference data, and restricted project data, the exact computational environment may need to be recreated before running the complete workflow. 

## Author

**Khaled Aboul Azm**

This repository was created to accompany the code and methodological work presented in the author’s master’s thesis.

## Citation

If you use or adapt the code or ACE-bounding algorithm, please cite the accompanying master’s thesis and this repository. I am planning on make the bounding code its own separate R package.

## License

Unless a license file states otherwise, the code is provided for academic review and reproducibility. Permission should be obtained from the author before redistribution or reuse.
