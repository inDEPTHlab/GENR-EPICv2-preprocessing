# Longitudinal DNAm Profiling Workflow

This repositiory contains scripts for longitudinal DNA methylation (DNAm) profiling, developed for and used in The Generation R Study.  

![](images/DNAm_workflow_github.png?raw=true "Longitudinal DNAm profiling workflow")

---

## Repository structure  

1. **[Sample preparation](1-sample_preparation/README.md)**  
	*Select participant samples and allocate processing plates.*    
	1.1 [Sample selection](1-sample_preparation/1.1-genr_epicv2_sample_selection_report.Rmd)  
	1.2 [Bridge selection](1-sample_preparation/1.2-bridges_and_extras_selection.R)  
	1.3 [Randomization](1-sample_preparation/1.3-randomization.R)  
2. **[Quality Control](2-quality_control/README.md)**  
	*Reduce technical bias while preserving biological variation.*  
	2.1 [QC and functional normalization preparation](2-quality_control/2.1-qc_and_normalization_preparation.R)  
	2.2 [Functional normalization](2-quality_control/2.2-normalization.R)  
	2.4 [Outlier handling](2-quality_control/2.4-outlier_handling.R)  
3. **[Analytical considerations](3-analytical_considerations/README.md)**  
	*Address remaining technical and biological variation in downstream analyses.*    
	3.1 [Cell type proportion calculation](3-analytical_considerations/3.1-cell_type_proportions.R)  
