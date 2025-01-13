# Variable_Selection_FDR_noisy
Variable selection with FDR control for noisy data – with application to
screening metabolites that are associated with breast and colorectal cancer. For the paper, please see https://arxiv.org/abs/2310.06696.

:paperclip: Folders 

***Key_function***: myest is used for half-min imputation, myest_mi is used for multiple imputation.

***Simulation***:

•	data_gen.r is used for generating the data  
•	Missing_only. R is used for the design that data are generated from the Gaussian distribution and only contain missing values (Table 1, S1)  
•	Measurement_error_only.R is used for the design that data are generated from the Gaussian distribution and only contain measurement error (Table 2)  
•	Missing_measurement.R is used for the design that data are generated from the Gaussian distribution and contain both measurement error and missing values (Table 3, S2, S3)  
•	Empirical_data_distribution.R is used for the design that data are generated from the empirical data distribution (Table 4, S6)  
•	Missing_measurement_multiple_datasets.R is used for the design of multiple outcomes (Table S5)  
•	Oracle. R is used for the results of Oracle results (Table S4).  

***Real_Data_analysis***: 
Codes of data analysis for four platforms: GC-MS, LC-MS, Lipidyzer, NMR.

•	imputation.R: impute the datasets using multiple imputation method: min-imputation, multiple imputation  
•	method_MI.R: the algorithm for multiple imputation datasets  
•	method_min.R: the algorithm for min imputation datasets  
•	single.R: Summary the results for single outcome (either breast cancer or colorectal cancer)  
•	both.R: Summary the results for both outcomes (common features for  breast cancer and colorectal cancers)  
