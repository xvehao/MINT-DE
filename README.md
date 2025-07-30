This is the code to reproduce all the plots in the paper Integrative Analysis of Differentially Expressed Genes in Time-Course Multi-Omics Data with MINT-DE.
Please contact Hao Xue (hx222(at)cornell.edu) for any questions. 
The R-package can be installed directly from github by 

1) devtools::install_github("https://github.com/xvehao/MINT-DE/")

or locally by

1) downloading mintde_0.1.0.tar.gz
2) In terminal: cd the_directory_where_mintde_is_downloaded 
3) R CMD install mintde

                
flyDevAnalysis.html is the script which performs data-preprocessing and analysis to reproduce figures using fly development data in the paper. 

pseudoPlatform.html is the counterpart for synthetic data. 

Fly_development_proteomics_normalized_counts.txt is the proteomics data used in the paper.

eda_function.R contains all the auxillary functions. 

Background: Time-course multi-omics experiments have been highly informative for obtaining a comprehensive understanding of the dynamic relationships between molecules in a biological process, especially if the different profiles are obtained from the same samples. A fundamental step in analyzing time-course multi-omics data involves selecting a short list of genes or gene regions ("sites") that warrant further study. Two important criteria for site selection are the magnitude of change and the temporal dynamic consistency. However, existing methods only consider one of these criteria, while neglecting the other.

Results: In our study, we propose a framework called MINT-DE (Multi-omics INtegration of Time-course for Differential Expression analysis) to address this limitation. MINT-DE is capable of selecting sites based on summarized measures of both aforementioned aspects. We calculate evidence measures assessing the extent of differential expression for each assay and for the dynamical similarity across assays. Then based on the summary of the evidence assessment measures, sites are ranked. To evaluate the performance of MINT-DE, we apply it to analyze a time-course multi-omics dataset of Drosophila development. We compare the selection obtained from MINT-DE with those obtained from other existing methods. The analysis reveals that MINT-DE is able to identify differentially expressed time-course pairs with the highest correlations. Their corresponding genes are significantly enriched for known biological functions, as measured by gene-gene interaction networks and the Gene Ontology enrichment.

Conclusions: These findings suggest the effectiveness of MINT-DE in selecting sites that are both differentially expressed within at least one assay and temporally related across assays. This highlights the potential of MINT-DE to identify biologically important sites for downstream analysis and provide a complementarity of sites that are neglected by existing methods.

