This is the code to reproduce all the plots in the paper Integrative Analysis of Differentially Expressed Genes in Time-Course Multi-Omics Data with MINT-DE.
Please contact Hao Xue (hx222(at)cornell.edu) for any questions. The core function is the function select_genes in eda_function.R. It is defined as
select_genes = function(param,method,p1,p2,p3=NULL,Genes=NULL)

This function selects candidate genes given three p-values (p-values for differential expression in two modalities and p-value for correlation test across modalities).

The inputs are 

                param: parameters correspond to each method, cutoff for p-value cutoff in commonDE/Venn diagram method, top_k for number of top k selected genes in Fisher's method or Edgington's method.
                
                method: the method to choose, "V" for commonDE/Venn diagram method, "F" for Fisher's method, "E" for Edgington's method;
                
                p1: p-values for significance of differential expression in the first modality;
                
                p2: p-values for significance of differential expression in the second modality;
                
                p3: p-values for significance of correlation test.
                
flyDevAnalysis.html is the script which performs data-preprocessing and analysis to reproduce figures using fly development data in the paper. 

pseudoPlatform.html is the counterpart for synthetic data. 

Fly_development_proteomics_normalized_counts.txt is the proteomics data used in the paper.

eda_function.R contains all the auxillary functions. 
