flyDevAnalysis
================
Hao Xue  
Cornell University
2023-12-26

# Abstract

We propose a framework for site selection for time-course multi-omics
experiments based on summarized measures of both magnitude of change and
temporal dynamic consistency. Specifically, for each modality, we
calculate a $p$-value measuring the extend of differential expression
and a $p$-value measuring the dynamical similarity across modalities.
Then based on the summary of $p$-values, sites are ranked.

To evaluate our approach, we applied it to analyze a developmental
time-course multi-omics dataset (Becker et al. 2018) and compared the
selection with other existing methods. The visualization of selected
time-series pairs from different assays shows high temporal concordance.
These findings demonstrate the effectiveness of our method in selecting
sites that are both differentially expressed within a modality and
temporally related across modalities, suggesting the potential of our
method to identify biologically important sites for downstream analysis
and provide a complementarity of sites that is ignored by existing
methods.

This demo reproduces Figures 3b-d in the paper.

# Load and Preprocess Data

``` r
# rm(list=ls())
library(plotly)
library(mixOmics) # import the mixOmics library
library(ggplot2)
library(dplyr)
library(gridExtra)
library(limma)
library(edgeR)
source('eda_function.R')
```

The dataset used in this demo could be downloaded from
[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121160).
First, convert mRNA level and protein level into LFC. Use *voom*(Law et
al. 2014) for mRNA. Genes with low mRNA read count are filtered as
suggested by *voom*.

``` r
files = list.files(path="~/time_course/data/GSE121160_RAW",full.names = T)
x <- readDGE(files, columns=c(1,2), 
             labels=str_replace(str_replace(files,'/home/hx222/time_course/data/GSE121160_RAW/GSM34271.._',''),
                                '.tsv.gz','')) 
proteomics <- read.delim("~/time_course/data/flyDev/Fly_development_proteomics_normalized_counts.txt")

hours = c(00,01,02,03,04,05,06,08,10,12,14,16,18,20)
rep_idx = list()
for (i in 1:4){
  rep_idx[[i]] = seq(i,ncol(x$counts),by=4)
}
hour_idx = list()
for (j in 1:length(hours)){
  hour_idx[[j]] = ((j-1)*4+1):(4*j)
}
group <- as.factor(str_replace(colnames(x$counts),'h_.','')) 
x$samples$group <- group 
design = model.matrix((~0+factor(rep(hours,each=4))))
colnames(design) = paste0('group',hours)
contrast.matrix = makeContrasts(group20-group0, group18-group0, group16-group0,
                                group14-group0, group12-group0, group10-group0,
                                group8-group0, group6-group0, group5-group0,
                                group4-group0, group3-group0, group2-group0,
                                group1-group0, levels=design)
keep <- filterByExpr(x$counts, design)
x_filtered = x
x_filtered$counts = x$counts[keep,]
rna_voom <- voom(x_filtered, design)

rna_lfc = rna_voom$E
for (h in 1:length(hours)){
  rna_lfc[,hour_idx[[h]]] = rna_voom$E[,hour_idx[[h]]]-rna_voom$E[,hour_idx[[1]]]
}

pro_lfc = proteomics[,-c(1,2)]
for (h in 1:length(hours)){
  pro_lfc[,hour_idx[[h]]] = proteomics[,hour_idx[[h]]+2]-
    proteomics[,hour_idx[[1]]+2]
}
pro_lfc$gene_names = proteomics$gene_names
pro_lfc$gene_ids = proteomics$gene_ids
```

# MINT-DE with Edgington’s Method

Given two omics, we want to find a shortlist of genes that are
differentially expressed in at least one assay, meanwhile sharing
similar temporal pattern across two assays. Assume one has log fold
change (LFC) against the baseline, $\{x_{tg}\}$, for the first assay
(i.e., gene expression level), $g=1,...,G$, $t=0,1,...,T$, where $G$ is
the number of sites sequenced, and $T$ is the total number of time
points. From the other assay (i.e., proteomics), one has $\{y_{tg}\}$,
$g=1,...,G$, and $t=0,1,...,T$. We assume that
$x_{tg}=\mu^X_{tg}+\epsilon^X_{tg}$ and
$y_{tg}=\mu^Y_{tg}+\epsilon^Y_{tg}$, where
$\epsilon^X_{tg},\,\epsilon^X_{tg}\sim N(0,\sigma_g^2)$.

We propose a method which simultaneously considers the significance of
change and temporal concordance across assays. We

1)  compute the p-value based on the hypothesis testing whether the site
    is differentially expressed at each time point. The p-value for gene
    $g$ at time $t$ is denoted as $p_{gt}^X$. Then one takes the minimum
    of p-values across all time points, i.e.,
    $p_{g}^X=\min(p_{1g}^{X},p_{2g}^X,...,p_{Tg}^X)$.

2)  repeat (i) with $y_{gt}^{(i)}$ and get $p_{g}^Y$.

``` r
### prepare data ###
rnaFitLimma = lmFit(rna_voom,design)
rnaFitLimmaContrast = contrasts.fit(rnaFitLimma,contrast.matrix)
rnaFitLimmaContrast = eBayes(rnaFitLimmaContrast)
topRNAgenes = topTable(rnaFitLimmaContrast)
rna_lfc_padj = cbind(rna_lfc,rnaFitLimmaContrast$p.value[rownames(rna_lfc),])
colnames(rna_lfc_padj)[(ncol(rna_lfc_padj)-length(hours)+2):ncol(rna_lfc_padj)] = paste0('rna.pvalue',hours[-1],'hr')

proFitLimma = lmFit(proteomics[,-c(1,2)],design)
proFitLimmaContrast = contrasts.fit(proFitLimma,contrast.matrix)
proFitLimmaContrast = eBayes(proFitLimmaContrast)
pro = cbind(proteomics,proFitLimmaContrast$p.value)
colnames(pro)[(ncol(pro)-length(hours)+2):ncol(pro)] = paste0('pro.pvalue',hours[-1],'hr')
pro_complete = pro[complete.cases(pro),]
pro_complete = pro_complete[pro_complete$gene_ids%in%rownames(rna_lfc_padj),]
common = cbind(pro_complete,rna_lfc_padj[pro_complete$gene_ids,])

common['rna.minp'] = common %>%
  dplyr::select(contains('rna.pvalue')) %>% as.matrix() %>% rowMins(na.rm = T)
common['pro.minp'] = common %>%
  dplyr::select(contains('pro.pvalue')) %>% as.matrix() %>% rowMins(na.rm = T)
```

3)  compute the p-values of the correlation test (Spearman or Pearson
    correlation test) between $x_{.g}$ and $y_{.g}$, denoted as
    $p_{g}^{XY}$ for gene $g$.

``` r
common['cor'] = as.numeric(as.character(rowcor(common %>%
                                                 dplyr::select(contains('imputed.log2.LFQ.intensity')),
                                               common %>%
                                                 dplyr::select('00h_1':'20h_4'))))
common['cor.test.p'] = rowcortest(common %>%
                                    dplyr::select(contains('imputed.log2.LFQ.intensity')),
                                  common %>%
                                    dplyr::select('00h_1':'20h_4'))
```

4)  use Edgington’s method (Edgington 1972) to summarize the $p$-values.
    Then sites with the smallest summarized $p$-values can be picked up.
    We pick top $10$, $20$, $50$, $100$, $200$, $500$, $1000$, and
    $2000$ sums in this demo.

``` r
topks = c(10,20,50,100,200,500,1000,2000)
edgington_genes = select_genes(topks,'E',common$rna.minp,common$pro.minp,
                               common$cor.test.p,Genes=common$gene_ids)
```

In Fisher’s method, the sums of log of $p$-values are used instead of
the direct sums. The Fisher’s method can be overwhelmed by small
$p$-values, while the Edgington’s method is sensitive to large
$p$-values.

``` r
fisher_genes = select_genes(topks,'F',common$rna.minp,common$pro.minp,
                            common$cor.test.p,Genes=common$gene_ids)
```

# Other methods

## TimeOmics

*timeOmics*(Bodein et al. 2022) uses sparse partial least square (sPLS)
(Lê Cao et al. 2009) to find $H$ latent factors for given matrices
$X,Y\in\mathbb{R}^{T\times G}$ separately such that the correlation
coefficients between latent factors are maximized with sparse penalty.
The omptimization problem is: where $X_h$ and $Y_h$ are deflated data at
the $h$th step and $P_{\lambda}(u)=\sum_{g=1}^G p_{\lambda}(|u_g|)$ is a
penalty function. The total loading genes could be interpreted as the
selected genes.

According to *timeOmics*, we first filtered out genes with low
coefficient of variation (CV) and keep only the top $2500$ genes with
the largest CV. Then we apply spls with canonical mode and set the
parameter pair, the number of component and the number of loadings of
each component to be $(1,10)$, $(2,10)$, $(2,25)$, $(2,50)$, $(4,50)$,
$(5,100)$, $(5,200)$, $(5,400)$. We do not use the cross-validation
feature in the package to choose the optimal number of components since
we want to fix the number of selections.

``` r
params_comb = array(c(1,10,2,10,2,25,2,50,4,50,5,100,5,200,5,400),dim=c(2,8))
common_geneids = intersect(rownames(rna_lfc),pro_lfc$gene_ids)
rna_spls_filted = spls_filter(rna_lfc[common_geneids,unlist(rep_idx)],
                              pro_lfc[pro_lfc$gene_ids%in%common_geneids,unlist(rep_idx)])
spls_rna = rna_spls_filted$spls_X
spls_pro = rna_spls_filted$spls_Y
spls_design <- data.frame(sample = rep(1:4,each=14))
spls_res = apply_spls(spls_rna, spls_pro, params_comb,spls_design,fpath='spls_res.Rdata')
load('eda/spls_res.Rdata')
pls_list1 = spls_res$pls_list1
pls_list2 = spls_res$pls_list2
```

*timeOmics* without filtering low count genes.

``` r
params_comb = array(c(1,10,2,10,2,25,2,50,4,50,5,100,5,200,5,400),dim=c(2,8))
spls_res_woFilter = apply_spls(t(rna_lfc[common_geneids,unlist(rep_idx)]), 
                      t(pro_lfc[pro_lfc$gene_ids%in%common_geneids,unlist(rep_idx)]), 
                      params_comb,spls_design,fpath='spls_res_woFilter.Rdata')
load('eda/spls_res_woFilter.Rdata')
pls_list1_woFilter = spls_res_woFilter$pls_list1
pls_list2_woFilter = spls_res_woFilter$pls_list2
```

## commonDE Method

A straightforward approach, referred as the Venn diagram method, is to
simply select sites with magnitude of change of time-courses exceeding
specified thresholds in both assays.

In Venn diagram method, we use different $p$-value cutoff to select
genes whose $p$-vlues for differential expression are below this cutoff
in two assays simultaneously.

``` r
cutoffs = c(-26,-24.7,-22.6,-21.3,-19.6,-16.4,-13.3,-8.4)
venn_genes = select_genes(cutoffs,'V',common$rna.minp,common$pro.minp,Genes=common$gene_ids)
```

# Results Analsyis

The correlation between selected pairs of time-courses (mRNA level LFC
and protein level LFC) are computed. This reproduces Figure 3b. The
correlations, shown in Figure b, between the selected pair of
time-courses by Venn diagram method are significantly lower than those
selected by Edgington’s method.

``` r
cor_df = rbind(data.frame('cor'=comp_cor(edgington_genes[[5]]$idx),
                          'method' = 'MINT-DE'),
               data.frame('cor'=comp_cor(venn_genes[[5]]$idx),
                          'method' = 'commonDE'),
               data.frame('cor'=comp_cor(1:nrow(common)),
                          'method' = 'genome'),
               data.frame('cor'=comp_cor(unlist(pls_list1[[5]]),"gene_id"),
                          'method' = 'timeOmics1'),
               data.frame('cor'=comp_cor(unlist(pls_list2[[5]]),"protein_index"),
                          'method' = 'timeOmics2'))
cor_df$type = 'cor'
ggplot(data=cor_df)+
  geom_boxplot(aes(x=method,y=abs(cor)))+
  ylab('|Correlation|')+
  xlab('')+
  theme_bw()+
  theme(text = element_text(size = 13))+
  scale_x_discrete(limits=c('genome','MINT-DE','commonDE','timeOmics1','timeOmics2'))
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggsave("top200cor.png", device = "png", width=6,height=4)
```

``` r
cor_df = rbind(data.frame('cor'=comp_cor(edgington_genes[[5]]$idx),
                          'method' = 'MINT-DE'),
               data.frame('cor'=comp_cor(venn_genes[[5]]$idx),
                          'method' = 'commonDE'),
               data.frame('cor'=comp_cor(1:nrow(common)),
                          'method' = 'genome'),
               data.frame('cor'=comp_cor(unlist(pls_list1[[5]]),"gene_id"),
                          'method' = 'timeOmics1'),
               data.frame('cor'=comp_cor(unlist(pls_list2[[5]]),"protein_index"),
                          'method' = 'timeOmics2'),
               data.frame('cor'=comp_cor(fisher_genes[[5]]$idx),
                          'method' = 'Fisher'))
```

    ## Warning in myids[i] <- which(common$gene_ids == input[i]): number of items to
    ## replace is not a multiple of replacement length

    ## Warning in myids[i] <- which(common$gene_ids == input[i]): number of items to
    ## replace is not a multiple of replacement length

    ## Warning in myids[i] <- which(common$gene_ids == pro_lfc[input[i], "gene_ids"]):
    ## number of items to replace is not a multiple of replacement length

    ## Warning in myids[i] <- which(common$gene_ids == pro_lfc[input[i], "gene_ids"]):
    ## number of items to replace is not a multiple of replacement length

    ## Warning in myids[i] <- which(common$gene_ids == pro_lfc[input[i], "gene_ids"]):
    ## number of items to replace is not a multiple of replacement length

    ## Warning in myids[i] <- which(common$gene_ids == pro_lfc[input[i], "gene_ids"]):
    ## number of items to replace is not a multiple of replacement length

    ## Warning in myids[i] <- which(common$gene_ids == pro_lfc[input[i], "gene_ids"]):
    ## number of items to replace is not a multiple of replacement length

    ## Warning in myids[i] <- which(common$gene_ids == pro_lfc[input[i], "gene_ids"]):
    ## number of items to replace is not a multiple of replacement length

``` r
cor_df$type = 'cor'
ggplot(data=cor_df)+
  geom_boxplot(aes(x=method,y=abs(cor)))+
  ylab('|Correlation|')+
  xlab('')+
  theme_bw()+
  theme(text = element_text(size = 13))+
  scale_x_discrete(limits=c('genome','MINT-DE','commonDE','timeOmics1','timeOmics2','Fisher'))
```

    ## Warning: Removed 36 rows containing non-finite values (`stat_boxplot()`).

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
cor_df2 = rbind(data.frame('cor'=comp_cor(edgington_genes[[4]]$idx),
                          'method' = 'MINT-DE'),
               data.frame('cor'=comp_cor(venn_genes[[4]]$idx),
                          'method' = 'commonDE'),
               data.frame('cor'=comp_cor(1:nrow(common)),
                          'method' = 'genome'),
               data.frame('cor'=comp_cor(unlist(pls_list1[[4]]),"gene_id"),
                          'method' = 'timeOmics1'),
               data.frame('cor'=comp_cor(unlist(pls_list2[[4]]),"protein_index"),
                          'method' = 'timeOmics2'),
               data.frame('cor'=comp_cor(fisher_genes[[4]]$idx),
                          'method' = "MINT-DE'"))
cor_df2$type = 'cor'
ggplot(data=cor_df2)+
  geom_boxplot(aes(x=method,y=abs(cor)))+
  ylab('|Correlation|')+
  xlab('')+
  theme_bw()+
  theme(text = element_text(size = 12))+
  scale_x_discrete(limits=c('genome','MINT-DE',"MINT-DE'",'commonDE','timeOmics1','timeOmics2'))
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ggsave("top100cor.png", device = "png", width=6,height=4)
```

The range of mRNA LFC time-courses are computed. This reproduces Figure
3c. The ranges of change of selections by both Edington’s method and the
Venn diagram method are overall larger than those by , while, on the
contrary, the LFC distribution of ’s selection is similar to a randomly
chosen genes from the genome.

``` r
range_venn = comp_range(rna_lfc[venn_genes[[5]]$gene_id,])
range_edgington = comp_range(rna_lfc[edgington_genes[[5]]$gene_id,])
range_spls1 = comp_range(rna_lfc[unlist(pls_list1[[5]]),])
range_spls2 = comp_range(as.matrix(pro_lfc[unlist(pls_list2[[5]]),unlist(rep_idx)]))
range_all = comp_range(as.matrix(common%>%dplyr::select('00h_1':'20h_4')))
range_df = rbind(data_frame(range=range_spls1,
                            method='timeOmics1'),
                 data_frame(range=range_spls2,
                            method='timeOmics2'),
                 data_frame(range=range_edgington,
                            method='MINT-DE'),
                 data_frame(range=range_venn,
                            method='commonDE'),
                 data_frame(range=range_all,
                            method='genome'))
range_plot = ggplot(data=range_df)+
  geom_boxplot(aes(x=method,y=range))+
  xlab('')+
  ylab('Range')+
  theme_bw()+
  theme(text = element_text(size = 13))+
  scale_x_discrete(limits=c('genome','MINT-DE','commonDE','timeOmics1','timeOmics2'))
range_plot
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggsave("top200range.png", device = "png", width=6,height=4)
```

``` r
range_fisher = comp_range(rna_lfc[fisher_genes[[5]]$gene_id,])
range_all = comp_range(as.matrix(common%>%dplyr::select('00h_1':'20h_4')))
range_df = rbind(data_frame(range=range_spls1,
                            method='timeOmics1'),
                 data_frame(range=range_spls2,
                            method='timeOmics2'),
                 data_frame(range=range_edgington,
                            method='MINT-DE'),
                 data_frame(range=range_venn,
                            method='commonDE'),
                 data_frame(range=range_all,
                            method='genome'),
                 data_frame(range=range_fisher,
                            method='Fisher'))
range_plot = ggplot(data=range_df)+
  geom_boxplot(aes(x=method,y=range))+
  xlab('')+
  ylab('Range')+
  theme_bw()+
  theme(text = element_text(size = 13))+
  scale_x_discrete(limits=c('genome','MINT-DE','commonDE','timeOmics1','timeOmics2','Fisher'))
range_plot
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
range_fisher = comp_range(rna_lfc[fisher_genes[[4]]$gene_id,])
range_venn = comp_range(rna_lfc[venn_genes[[4]]$gene_id,])
range_edgington = comp_range(rna_lfc[edgington_genes[[4]]$gene_id,])
range_spls1 = comp_range(rna_lfc[unlist(pls_list1[[4]]),])
range_spls2 = comp_range(as.matrix(pro_lfc[unlist(pls_list2[[4]]),unlist(rep_idx)]))
range_all = comp_range(as.matrix(common%>%dplyr::select('00h_1':'20h_4')))
range_df = rbind(data_frame(range=range_spls1,
                            method='timeOmics1'),
                 data_frame(range=range_spls2,
                            method='timeOmics2'),
                 data_frame(range=range_edgington,
                            method='MINT-DE'),
                 data_frame(range=range_fisher,
                            method="MINT-DE'"),
                 data_frame(range=range_venn,
                            method='commonDE'),
                 data_frame(range=range_all,
                            method='genome'))
range_plot = ggplot(data=range_df)+
  geom_boxplot(aes(x=method,y=range))+
  xlab('')+
  ylab('Range')+
  theme_bw()+
  theme(text = element_text(size = 12))+
  scale_x_discrete(limits=c('genome','MINT-DE',"MINT-DE'",'commonDE','timeOmics1','timeOmics2'))
range_plot
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggsave("top100range.png", device = "png", width=6,height=4)
```

Uncomment to save the results.

``` r
write.csv(venn_genes[[5]]$gene_id,
          'eda/VennGeneList200.csv',
          row.names=F)
write.csv(edgington_genes[[5]]$gene_id,
          'eda/EdgingtonGeneList200.csv',
          row.names=F)
write.csv(unlist(pls_list1[[5]]),
          'eda/plsGeneList200_filtered.csv',
          row.names=F)
write.csv(rownames(rna_lfc),
          'eda/backgroundGeneList.csv',
          row.names=F)
write.csv(unlist(pls_list1_woFilter[[5]]),
          'eda/plsGeneList200_woFilter.csv',
          row.names=F)
```

To further illustrate the difference between the selection based
Edgington’s method and that based on *timeOmics*, we fit a spline with 5
degrees of freedom with the mRNA level LFC of the selected genes. We
plot the $p$-values of the fitted model vs. average LFC. This reproduce
Figure 3d. The selection based on *timeOmics* are concentrated around
the y-axis, since it can not specifically choose those genes with large
LFC. On the contrary, our method tends to choose genes more scattered
away from the y-axis, which implies those genes have larger LFC.

``` r
library(splines)
X <- ns(rep(hours,each=4),df=5)
spDesign <- model.matrix(~X)
fit <- lmFit(rna_lfc, spDesign)
fit <- eBayes(fit)
iii=5
selected_sites = rbind(data.frame(coef=topTable(fit[unlist(pls_list1[[iii]]),], coef=2:6,
                                                length(unlist(pls_list1[[iii]])))$AveExpr,
                                  p.value=topTable(fit[unlist(pls_list1[[iii]]),], coef=2:6,
                                                   length(unlist(pls_list1[[iii]])))$P.Value,
                                  method='timeOmics'),
                       data.frame(coef=topTable(fit[edgington_genes[[iii]]$gene_id,], coef=2:6,
                                                length(edgington_genes[[iii]]$gene_id))$AveExpr,
                                  p.value=topTable(fit[edgington_genes[[iii]]$gene_id,], coef=2:6,
                                                   length(edgington_genes[[iii]]$gene_id))$P.Value,
                                  method='MINT-DE'))
# all_sites = data.frame(coef=topTable(fit, coef=2:6,
#                                      number=nrow(fit))$AveExpr,
#                        p.value=topTable(fit, coef=2:6,
#                                         number=nrow(fit))$P.Value,
#                        method='all')
ggplot()+
  geom_point(data=selected_sites,
             aes(x=coef,
                 y=-log10(p.value),
                 colour = method),
                 alpha=0.5)+
  xlab('Average Expression')+
  ylab('-log(P-value)')+
  theme_bw()+
  theme(text = element_text(size = 20))
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
ggsave("VolcanoPlot.png", device = "png", width=12,height=8)
```

``` r
iii = 5
selected_sites = rbind(data.frame(coef=topTable(fit[venn_genes[[iii]]$gene_id,], coef=2:6,
                                                length(venn_genes[[iii]]$gene_id))$AveExpr,
                                  p.value=topTable(fit[venn_genes[[iii]]$gene_id,], coef=2:6,
                                                   length(venn_genes[[iii]]$gene_id))$P.Value,
                                  method='commonDE'),
                       data.frame(coef=topTable(fit[edgington_genes[[iii]]$gene_id,], coef=2:6,
                                                length(edgington_genes[[iii]]$gene_id))$AveExpr,
                                  p.value=topTable(fit[edgington_genes[[iii]]$gene_id,], coef=2:6,
                                                   length(edgington_genes[[iii]]$gene_id))$P.Value,
                                  method='MINT-DE'))
# all_sites = data.frame(coef=topTable(fit, coef=2:6,
#                                      number=nrow(fit))$AveExpr,
#                        p.value=topTable(fit, coef=2:6,
#                                         number=nrow(fit))$P.Value,
#                        method='all')
ggplot()+
  geom_point(data=selected_sites,
             aes(x=coef,
                 y=-log10(p.value),
                 colour = method),
                 alpha=0.5)+
  xlab('Average Expression')+
  ylab('-log(P-value)')+
  theme_bw()+
  theme(text = element_text(size = 20))
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggsave("VolcanoPlot_supp.png", device = "png", width=12,height=8)
```



    ```r
    plot(rna_lfc['FBgn0264975',])

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='Nrg',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0010894',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='sinu',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0015777',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='nrv2',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0033032',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='kune',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0010382',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='CycE',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0011206',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='bol',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0016070',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='smg',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0026620',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='tacc',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

### spindle organization

``` r
plot(rna_lfc['FBgn0004379',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='Klp67A',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

``` r
plot(rna_lfc['FBgn0011606',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='Klp3A',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

### transposition

``` r
plot(rna_lfc['FBgn0000146',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='aub',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

### visual behavior

``` r
plot(rna_lfc['FBgn0004646',])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
plot(1:56,pro_lfc[pro_lfc$gene_names=='ogre',-c(57,58)])
```

![](flyDevAnalysis_files/figure-gfm/unnamed-chunk-30-2.png)<!-- -->

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-becker2018quantifying" class="csl-entry">

Becker, Kolja, Alina Bluhm, Nuria Casas-Vila, Nadja Dinges, Mario
Dejung, Sergi Sayols, Clemens Kreutz, Jean-Yves Roignant, Falk Butter,
and Stefan Legewie. 2018. “Quantifying Post-Transcriptional Regulation
in the Development of Drosophila Melanogaster.” *Nature Communications*
9 (1): 4970.

</div>

<div id="ref-bodein2022timeomics" class="csl-entry">

Bodein, Antoine, Marie-Pier Scott-Boyer, Olivier Perin, Kim-Anh Lê Cao,
and Arnaud Droit. 2022. “timeOmics: An r Package for Longitudinal
Multi-Omics Data Integration.” *Bioinformatics* 38 (2): 577–79.

</div>

<div id="ref-edgington1972additive" class="csl-entry">

Edgington, Eugene S. 1972. “An Additive Method for Combining Probability
Values from Independent Experiments.” *The Journal of Psychology* 80
(2): 351–63.

</div>

<div id="ref-law2014voom" class="csl-entry">

Law, Charity W, Yunshun Chen, Wei Shi, and Gordon K Smyth. 2014. “Voom:
Precision Weights Unlock Linear Model Analysis Tools for RNA-Seq Read
Counts.” *Genome Biology* 15 (2): 1–17.

</div>

<div id="ref-le2009sparse" class="csl-entry">

Lê Cao, Kim-Anh, Pascal GP Martin, Christèle Robert-Granié, and Philippe
Besse. 2009. “Sparse Canonical Methods for Biological Data Integration:
Application to a Cross-Platform Study.” *BMC Bioinformatics* 10: 1–17.

</div>

</div>
