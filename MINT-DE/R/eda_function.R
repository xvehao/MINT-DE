library(DESeq2)
library(genomation)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(dplyr)
library(stringr)
library(ggplot2)
library(readxl)
library(ggtext)

apply_deseq2 = function(dat, hours=c(0,3,6,12,24,48,72)){
  exp_level_col = str_detect(colnames(dat),pattern='hr_')
  n_rep = sum(exp_level_col)/length(hours)
  dat_rmNA = dat[complete.cases(dat),]
  coldat = data.frame(condition=factor(rep(c("untreated","treated"),each=n_rep)))
  coldat$condition <- relevel(coldat$condition, ref = "untreated")
  lfc = data.frame(matrix(0,nrow=nrow(dat_rmNA),ncol=length(hours)))
  pvalue = padj = data.frame(matrix(1,nrow=nrow(dat_rmNA),ncol=length(hours)))
  colnames(lfc) = colnames(pvalue) = colnames(padj)=c(paste0(hours,'hr'))
  count_normalized = dat_rmNA
  ref_col = colnames(dat_rmNA)[str_detect(colnames(dat_rmNA),pattern='0hr')]
  for (h in hours[-1]){
    trt_col = colnames(dat_rmNA)[str_detect(colnames(dat_rmNA),pattern=paste0(h,'hr'))]
    dds <- DESeqDataSetFromMatrix(countData = dat_rmNA[,c(ref_col,trt_col)],
                                  colData = coldat,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res = results(dds)
    lfc[paste0(h,'hr')]=res$log2FoldChange
    pvalue[paste0(h,'hr')]=res$pvalue
    padj[paste0(h,'hr')]=res$padj
    count_normalized[,trt_col] = counts(dds)[,trt_col]
  }
  count_normalized[,ref_col] = counts(dds)[,ref_col]
  if('geneName' %in% colnames(dat_rmNA)){
    lfc['GENENAME'] = pvalue['geneName'] = padj['geneName'] = dat_rmNA$geneName
  }else{
    loc_idx = colnames(dat_rmNA)[!exp_level_col]
    lfc[,loc_idx] = pvalue[,loc_idx] = padj[,loc_idx] = dat_rmNA[,loc_idx]
  }
  return(list(lfc=lfc,pvalue=pvalue,padj=padj,count_normalized=count_normalized))
}

annotate_peak = function(peak_gr){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  peakAnno <- annotatePeak(peak_gr, TxDb=txdb, verbose=FALSE)
  peakAnno_df = data.frame(peakAnno)
  entrez <- peakAnno_df$geneId
  
  # Return the gene symbol for the set of Entrez IDs
  annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                           keys = entrez,
                                           columns = c("GENENAME"),
                                           keytype = "ENTREZID")
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  peakAnno_DF = peakAnno_df %>% 
    left_join(annotations_edb, by=c("geneId"="ENTREZID"))
  return(peakAnno_DF)
}

read_clu = function(){
  mmc3 <- read_excel("../references/atac-seq/mmc3.xlsx",
                     sheet = "2.Differential peaks-expression")
  chip = mmc3[-1,1:4]
  colnames(chip) = mmc3[1,1:4]
  atac = mmc3[-1,6:9]
  colnames(atac) = mmc3[1,6:9]
  rna = mmc3[-1,11:12]
  colnames(rna) = mmc3[1,11:12]
  return(list(chip=chip,rna=rna,atac=atac))
}

plot_lfc_KGG = function(genes,lfcdat=rnaDE$lfc,ylab="RNA lfc"){
  datFC = lfcdat[lfcdat$GENENAME %in% genes,]
  # fold change plot:
  datFC = gather (datFC, hours, FC, '0hr':'72hr', factor_key=TRUE)
  datFC$hours = gsub('hr', '', datFC$hours)
  if ('Start' %in% colnames(datFC)){
    ggplot(data=datFC, aes(x=as.numeric(as.character(hours)),
                           y=as.numeric(FC), 
                           group=Start,
                           color=GENENAME
    ))+
      geom_line(size=.5, alpha=0.4) + geom_point(size=1) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                           as.numeric(as.character(datFC$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                               as.numeric(as.character(datFC$hours)),
                                             name ="hours"))+
      ylab(ylab)+
      ggtitle(paste0(c(KEGG[j,c('ID','Description','GeneRatio','p.adjust','Count')]),collapse = ' '))+
      theme(plot.title = element_textbox_simple(
        size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
  }else{
    ggplot(data=datFC, aes(x=as.numeric(as.character(hours)),
                           y=as.numeric(FC), 
                           color=GENENAME
    ))+
      geom_line(size=.5) + geom_point(size=1) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                           as.numeric(as.character(datFC$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                               as.numeric(as.character(datFC$hours)),
                                             name ="hours"))+
      ylab(ylab)+
      ggtitle(paste0(c(KEGG[j,c('ID','Description','GeneRatio','p.adjust','Count')]),collapse = ' '))+
      theme(plot.title = element_textbox_simple(
        size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
  }
}

plot_lfc = function(genes,lfcdat=rnaDE$lfc,ylab="RNA lfc"){
  datFC = lfcdat[lfcdat$GENENAME %in% genes,]
  # fold change plot:
  datFC = gather (datFC, hours, FC, '0hr':'72hr', factor_key=TRUE)
  datFC$hours = gsub('hr', '', datFC$hours)
  if ('Start' %in% colnames(datFC)){
    ggplot(data=datFC, aes(x=as.numeric(as.character(hours)),
                           y=as.numeric(FC), 
                           group=Start,
                           color=GENENAME
    ))+
      geom_line(size=.5, alpha=0.4) + geom_point(size=1) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                           as.numeric(as.character(datFC$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                               as.numeric(as.character(datFC$hours)),
                                             name ="hours"))+
      ylab(ylab)+
      ggtitle(paste0(c(KEGG[j,c('ID','Description','GeneRatio','p.adjust','Count')]),collapse = ' '))+
      theme(plot.title = element_textbox_simple(
        size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
  }else{
    ggplot(data=datFC, aes(x=as.numeric(as.character(hours)),
                           y=as.numeric(FC), 
                           color=GENENAME
    ))+
      geom_line(size=.5) + geom_point(size=1) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                           as.numeric(as.character(datFC$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                               as.numeric(as.character(datFC$hours)),
                                             name ="hours"))+
      ylab(ylab)+
      theme(plot.title = element_textbox_simple(
        size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
  }
}

rowcor = function(A,B){
  sapply(seq.int(dim(A)[1]), function(i) {
    if (!any(c(is.na(as.numeric(A[i,])), is.na(as.numeric(B[i,]))))){
      cor(as.numeric(A[i,]), as.numeric(B[i,]))
    }
  })
}

rowcortest = function(A,B){
  sapply(seq.int(dim(A)[1]), function(i) cor.test(as.numeric(A[i,]), as.numeric(B[i,]),method='pearson')$p.value)
}

apply_deseq2_flyDev = function(dat, hour_idx, hours=c(0,3,6,12,24,48,72)){
  n_rep = length(hour_idx[[1]])
  dat_rmNA = dat[complete.cases(dat),]
  coldat = data.frame(condition=factor(rep(c("untreated","treated"),each=n_rep)))
  coldat$condition <- relevel(coldat$condition, ref = "untreated")
  lfc = data.frame(matrix(0,nrow=nrow(dat_rmNA),ncol=length(hours)))
  pvalue = padj = data.frame(matrix(1,nrow=nrow(dat_rmNA),ncol=length(hours)))
  colnames(lfc) = colnames(pvalue) = colnames(padj)=c(paste0(hours,'hr'))
  count_normalized = dat_rmNA
  ref_col = colnames(dat_rmNA)[hour_idx[[1]]]
  for (h in 2:length(hour_idx)){
    trt_col = colnames(dat_rmNA)[hour_idx[[h]]]
    dds <- DESeqDataSetFromMatrix(countData = dat_rmNA[,c(ref_col,trt_col)],
                                  colData = coldat,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res = results(dds)
    lfc[paste0(hours[h],'hr')]=res$log2FoldChange
    pvalue[paste0(hours[h],'hr')]=res$pvalue
    padj[paste0(hours[h],'hr')]=res$padj
    count_normalized[,trt_col] = counts(dds)[,trt_col]
  }
  count_normalized[,ref_col] = counts(dds)[,ref_col]
  return(list(lfc=lfc,pvalue=pvalue,padj=padj,count_normalized=count_normalized))
}

plot_pls = function(rna_lfc, comp_loadings, destination){
  pdf(file=destination, onefile=T)
  par(mfrow = c(2,2))
  for (i in 1:10){
    comp_loading = comp_loadings[[i]]
    matplot(t(rna_lfc[comp_loading, seq(1,ncol(rna_lfc),by=4)]),
            type='l',main=paste0('sam1_mrna_latent_variable',i),ylab='lfc')
    matplot(t(rna_lfc[comp_loading, seq(2,ncol(rna_lfc),by=4)]),
            type='l',main=paste0('sam2_mrna_latent_variable',i),ylab='lfc')
    matplot(t(rna_lfc[comp_loading, seq(3,ncol(rna_lfc),by=4)]),
            type='l',main=paste0('sam3_mrna_latent_variable',i),ylab='lfc')
    matplot(t(rna_lfc[comp_loading, seq(4,ncol(rna_lfc),by=4)]),
            type='l',main=paste0('sam4_mrna_latent_variable',i),ylab='lfc')
    # legend("topright", comp_loading, col=seq_len(length(comp_loading)),
    #        cex=0.8, inset=c(-.28,0), lty=seq_len(length(comp_loading)),
    #        lwd=2, xpd=TRUE)
  }
  dev.off()
}

plot_flyDev_rep = function(destination = './eda/flyCount_rna_rna.pdf',
                           mygenes=setdiff(venn_genes,Edgington_genes),
                           assay1=rna1_lfc_padj, assay2=rna2_lfc_padj){
  pdf(file=destination, onefile=T)
  par(mfrow = c(1,1))
  for (gene in mygenes){
    datFC1 = assay1[gene,]
    
    # fold change plot:
    datFC1 = gather (datFC1, timepoint, FC, lfc0hr:lfc20hr, factor_key=TRUE)
    datFC1$timepoint = 1:length(hours)
    datFC1$hours = hours
    datFC1$stage = c(rep("early", 7), rep("late", 7))
    gFC1 = ggplot(data=datFC1, aes(x=as.numeric(as.character(hours)),
                                   y=as.numeric(FC)
    ))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC1$hours)), labels =
                           as.numeric(as.character(datFC1$timepoint)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC1$hours)), labels =
                                               as.numeric(as.character(datFC1$hours)),
                                             name ="hours"))+
      ylab("RNA log fold change")+#ggtitle(firsthalf_genedescription)+
      ggtitle(gene)+
      theme(legend.position="none",
            plot.title = element_textbox_simple(
              size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
    
    datFC2 = assay2[gene,]
    
    # fold change plot:
    datFC2 = gather (datFC2, timepoint, FC, lfc0hr:lfc20hr, factor_key=TRUE)
    datFC2$timepoint = 1:length(hours)
    datFC2$hours = hours
    datFC2$stage = c(rep("early", 7), rep("late", 7))
    gFC2 = ggplot(data=datFC2, aes(x=as.numeric(as.character(hours)),
                                   y=as.numeric(FC)
    ))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC2$hours)), labels =
                           as.numeric(as.character(datFC2$timepoint)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC2$hours)), labels =
                                               as.numeric(as.character(datFC2$hours)),
                                             name ="hours"))+
      ylab("RNA log fold change")+#ggtitle(firsthalf_genedescription)+
      ggtitle(gene)+
      theme(legend.position="none",
            plot.title = element_textbox_simple(
              size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
    
    # counts plot:
    datCounts1 = rna_count[rownames(rna_count)==gene,c(rep_idx[[1]],rep_idx[[2]])]
    #colnames(datCounts) = c("gene_ID", "gene_name", paste("X", colnames(datCounts)[-c(1,2)], sep=""))
    datCounts1 = gather (datCounts1, timepoint, Counts, X00h_1:X20h_2, factor_key=TRUE)
    datCounts1$timepoint = rep(c(1:length(hours)), each=2)
    datCounts1$hours = rep(hours, 2)
    datCounts1$replicate = rep(c("A","B"),each=nrow(datCounts1)/2)
    
    gCounts1 = ggplot(data=datCounts1, aes(x=as.numeric(as.character(hours)),
                                           y=as.numeric(Counts), shape=replicate))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC1$hours)), labels =
                           as.numeric(as.character(datFC1$timepoint)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC1$hours)), labels =
                                               as.numeric(as.character(datFC1$hours)),
                                             name ="hours"))+
      ylab("RNA normalized counts")+
      ggtitle(paste0("normalized counts for both replicates ", gene))
    
    datCounts2 = rna_count[rownames(rna_count)==gene,c(rep_idx[[3]],rep_idx[[4]])]
    #colnames(datCounts) = c("gene_ID", "gene_name", paste("X", colnames(datCounts)[-c(1,2)], sep=""))
    datCounts2 = gather (datCounts2, timepoint, Counts, X00h_3:X20h_4, factor_key=TRUE)
    datCounts2$timepoint = rep(c(1:length(hours)), each=2)
    datCounts2$hours = rep(hours, 2)
    datCounts2$replicate = rep(c("A","B"),each=nrow(datCounts2)/2)
    gCounts2 = ggplot(data=datCounts2, aes(x=as.numeric(as.character(hours)),
                                           y=as.numeric(Counts), shape=replicate))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC2$hours)), labels =
                           as.numeric(as.character(datFC2$timepoint)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC2$hours)), labels =
                                               as.numeric(as.character(datFC2$hours)),
                                             name ="hours"))+
      ylab("RNA normalized counts")+
      ggtitle(paste0("normalized counts for both replicates ", gene))
    
    grid.arrange(gFC1, gFC2, gCounts1, gCounts2, nrow=2)
  }
  dev.off()
}


plot_flyDev_rna_pro = function(destination = './eda/flyCount_rna_pro.pdf',
                           myids=setdiff(venn_genes,Edgington_genes)){
  if (length(destination)>0){
    pdf(file=destination, onefile=T)
    par(mfrow = c(1,1))
  }
  for (myid in myids){
    # counts plot:
    gene = common[myid,'gene_ids']
    datCounts1 = rna_lfc[gene,]
    datCounts1 = data.frame(lfc=rna_lfc[gene,],
                            hours=rep(hours, each=4),
                            replicate=rep(c("A","B","C","D"),ncol(rna_lfc)/4))
    gCounts1 = ggplot(data=datCounts1, 
                      aes(x=as.numeric(as.character(hours)),
                          y=as.numeric(lfc), shape=replicate),
                      alpha=0.5)+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=10),
            plot.title=element_text(size=20),
            axis.text=element_text(size=20),
            legend.text=element_text(size=10),
            legend.title=element_text(size=20),
            axis.title.x=element_text(size=20, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=20, margin=margin(0,5,0,0))) +
      scale_x_continuous("Hours", breaks=as.numeric((datCounts1$hours)), labels =
                           as.numeric(as.character(datCounts1$hours)))+
      ylab("Gene Expression LFC")+
      ggtitle(paste0(gene,": ",common[myid,'gene_names']))
    
    mypro_lfc = mypro = common[myid,] %>% dplyr::select(imputed.log2.LFQ.intensity.00h_1:imputed.log2.LFQ.intensity.20h_4)
    for (h in 1:length(hours)){
      mypro_lfc[,hour_idx[[h]]] = mypro[,hour_idx[[h]]]-
        mypro[,hour_idx[[1]]]
    }
    datCounts2 = data.frame(lfc=as.numeric(mypro_lfc),
                            hours=rep(hours, each=4),
                            replicate=rep(c("A","B","C","D"),ncol(rna_lfc)/4))
    avgcor = mean(sapply(unique(datCounts2$replicate),
                         function(x) cor(datCounts2[datCounts2$replicate==x,'lfc'],
                                         datCounts1[datCounts1$replicate==x,'lfc'])))
    gCounts2 = ggplot(data=datCounts2, 
                      aes(x=hours,
                          y=lfc, 
                          shape=replicate),
                      alpha=0.5)+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=10),
            plot.title=element_text(size=20),
            axis.text=element_text(size=20),
            legend.text=element_text(size=10),
            legend.title=element_text(size=20),
            axis.title.x=element_text(size=20, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=20, margin=margin(0,5,0,0))) +
      scale_x_continuous("Hours", breaks=as.numeric((datCounts2$hours)), labels =
                           as.numeric(as.character(datCounts2$hours)))+
      ylab("Proteomics LFC")+
      ggtitle(paste0('Average Correlation: ',round(avgcor,2)))
    
    grid.arrange(gCounts1, gCounts2, nrow=2)
  }
  if (length(destination)>0){
    dev.off()
  }
}


plot_flyDev_clustered_rna = function(destination = './eda/flyCount_rna_pro.pdf',
                               mygenes=setdiff(venn_genes,Edgington_genes)){
  
  if (length(destination)>0){
    pdf(file=destination, onefile=T)
    par(mfrow = c(1,1))
  }
  plist = list()
  for (i in 1:length(mygenes)){
    # counts plot:
    gene = mygenes[i]
    datCounts1 = rna_lfc[gene,]
    datCounts1 = data.frame(lfc=rna_lfc[gene,],
                            hours=rep(hours, each=4),
                            replicate=rep(c("A","B","C","D"),ncol(rna_lfc)/4))
    gCounts1 = ggplot(data=datCounts1, aes(x=as.numeric(as.character(hours)),
                                           y=as.numeric(lfc), shape=replicate,
                                           alpha=0.5))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=10),
            axis.text=element_text(size=10),
            legend.text=element_text(size=10),
            legend.title=element_text(size=10),
            axis.title.x=element_text(size=10, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=10, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datCounts1$hours)), labels =
                           as.numeric(as.character(datCounts1$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datCounts1$hours)), labels =
                                               as.numeric(as.character(datCounts1$hours)),
                                             name ="hours"))+
      ylab("RNA lfc")+
      ggtitle(paste0("lfc for four replicates ", gene))
    plist[[i]]=gCounts1
    if (i%%4==0){
      grid.arrange(grobs=plist[(i-3):i], nrow=2,ncol=2)
    }
    if (length(mygenes)<4){
      grid.arrange(grobs=plist[1:length(mygenes)], nrow=2, ncol=2)
    }
  }
  if (length(destination)>0){
    dev.off()
  }
}


plot_flyDev_clustered_pro = function(destination = './eda/flyCount_rna_pro.pdf',
                                     mygenes=setdiff(venn_genes,Edgington_genes)){
  
  if (length(destination)>0){
    pdf(file=destination, onefile=T)
    par(mfrow = c(1,1))
  }
  plist = list()
  for (i in 1:length(mygenes)){
    # counts plot:
    gene = mygenes[i]
    datCounts1 = pro_lfc[pro_lfc$gene_names==gene,]
    #colnames(datCounts) = c("gene_ID", "gene_name", paste("X", colnames(datCounts)[-c(1,2)], sep=""))
    datCounts1 = tidyr::gather (datCounts1, timepoint, Counts, 
                         imputed.log2.LFQ.intensity.00h_1:imputed.log2.LFQ.intensity.20h_4, factor_key=TRUE)
    datCounts1$hours = rep(hours, each=4)
    datCounts1$replicate = rep(c("A","B","C","D"),nrow(datCounts1)/4)
    
    gCounts1 = ggplot(data=datCounts1, aes(x=as.numeric(as.character(hours)),
                                           y=as.numeric(Counts), shape=replicate,
                                           alpha=0.5))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datCounts1$hours)), labels =
                           as.numeric(as.character(datCounts1$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datCounts1$hours)), labels =
                                               as.numeric(as.character(datCounts1$hours)),
                                             name ="hours"))+
      ylab("RNA lfc")+
      ggtitle(paste0("lfc for four replicates ", gene))
    plist[[i]]=gCounts1
    if (i%%4==0){
      grid.arrange(grobs=plist[(i-3):i], nrow=2,ncol=2)
    }
  }
  if (length(destination)>0){
    dev.off()
  }
}

plot_flyDev_pseudo = function(destination = './eda/flyCount_pseudo.pdf',
                              myids, platform1, platform2){
  if (length(destination)>0){
    pdf(file=destination, onefile=T)
    par(mfrow = c(1,1))
  }
  for (myid in myids){
    # counts plot:
    datCounts1 = data.frame(lfc=as.numeric(platform1[myid,]),
                            hours=rep(hours, each=2),
                            replicate=rep(c("A","B"),ncol(platform1)/2))
    
    gCounts1 = ggplot(data=datCounts1, aes(x=as.numeric(as.character(hours)),
                                           y=as.numeric(lfc), shape=replicate,
                                           alpha=0.5))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datCounts1$hours)), labels =
                           as.numeric(as.character(datCounts1$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datCounts1$hours)), labels =
                                               as.numeric(as.character(datCounts1$hours)),
                                             name ="hours"))+
      ylab("RNA lfc")+
      ylim(-10,10)
    
    datCounts2 = data.frame(lfc=as.numeric(platform2[myid,]),
                            hours=rep(hours, each=2),
                            replicate=rep(c("A","B"),ncol(platform1)/2))
    gCounts2 = ggplot(data=datCounts2, aes(x=hours,
                                           y=lfc, 
                                           shape=replicate,
                                           alpha=0.5))+
      geom_line(size=1) + geom_point(size=3) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datCounts2$hours)), labels =
                           as.numeric(as.character(datCounts2$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datCounts2$hours)), labels =
                                               as.numeric(as.character(datCounts2$hours)),
                                             name ="hours"))+
      ylab("Proteomics log2LFQ")+
      ylim(-10,10)
    
    grid.arrange(gCounts1, gCounts2, nrow=2)
  }
  if (length(destination)>0){
    dev.off()
  }
}


comp_coefvar = function(input){
  return(sd(input)/mean(input))
}


apply_spls = function(X, Y, params_comb,spls_design,fpath=NULL){
  pls_list1 = pls_list2 = list()
  for (i in 1:ncol(params_comb)){
    ncomp = params_comb[1,i]
    tokeep = params_comb[2,i]
    print(paste0('running pls with ncomp: ',ncomp,' to keep: ',tokeep))
    spls.result <- spls(X,
                        Y,
                        keepX = rep(tokeep,ncomp),
                        keepY = rep(tokeep,ncomp),
                        ncomp = ncomp,
                        multilevel = spls_design,
                        mode = "canonical")
    comp_loadings1 = list()
    comp_loadings2 = list()
    if (ncomp==1){
      comp_loading1 = names(spls.result$loadings$X[abs(spls.result$loadings$X[,1])>0,])
      comp_loadings1[[1]] = comp_loading1
      comp_loading2 = names(spls.result$loadings$Y[abs(spls.result$loadings$Y[,1])>0,])
      comp_loadings2[[1]] = comp_loading2
    }else{
      for (j in 1:ncomp){
        comp_loading1 = rownames(spls.result$loadings$X[abs(spls.result$loadings$X[,j])>0,])
        comp_loadings1[[j]] = comp_loading1
        comp_loading2 = rownames(spls.result$loadings$Y[abs(spls.result$loadings$Y[,j])>0,])
        comp_loadings2[[j]] = comp_loading2
      }
    }
    pls_list1[[i]] = comp_loadings1
    pls_list2[[i]] = comp_loadings2
  }
  spls_res = list(pls_list1=pls_list1,
       pls_list2=pls_list2)
  if (length(fpath)>0){
    save(list=c('spls_res'),
         file=paste0('eda/',fpath))
  }
  return(spls_res)
}

spls_filter = function(lfc1, lfc2, tokeep=2500){
  cv1 = apply(lfc1, 1, function(x){ 
    abs(sd(x, na.rm = TRUE)/mean(x, na.rm= TRUE))
  })
  cv2 = apply(lfc2, 1, function(x){ 
    abs(sd(x, na.rm = TRUE)/mean(x, na.rm= TRUE))
  })
  
  spls_X = t(lfc1)[,rank(-cv1)<=tokeep]
  spls_Y = t(lfc2)[,rank(-cv2)<=tokeep]
  return(list(spls_X=spls_X,spls_Y=spls_Y))
}

select_genes = function(param,method,p1,p2,p3=NULL,Genes=NULL){
  selected_genes = list()
  for (k in 1:length(param)){
    if (method=='F'){
      selected_genes[[k]] = Fisher_method(param[k],p1,p2,p3,Genes=Genes)
    }else if (method=='E'){
      selected_genes[[k]] = Edgington_method(param[k],p1,p2,p3,Genes=Genes)
    }else if (method=='V'){
      selected_genes[[k]] = Venn_method(param[k],p1,p2,Genes=Genes)
    }
    print(length(selected_genes[[k]]$gene_id))
  }
  return(selected_genes)
}

summarize_res = function(selected_genes,dataset){
  myvar = mycoefvar = list()
  for (k in 1:length(selected_genes)){
    myvar[[k]] = mycoefvar[[k]] = array(NA,c(length(selected_genes[[k]]$gene_id),length(hours)-1),
                                        dimnames=list(selected_genes[[k]]$gene_id,hours[-1]))
    for (j in 2:length(hours)){
      rowdat = dataset[selected_genes[[k]]$gene_id,] %>% 
        dplyr::select(contains(paste0(hours[j],'h_'))) %>%
        as.matrix() 
      myvar[[k]][selected_genes[[k]]$gene_id,hours[j]] = rowdat %>% rowVars()
      rowsd = rowdat %>% rowSds()
      rowmean = rowdat %>% rowMeans()
      mycoefvar[[k]][selected_genes[[k]]$gene_id,hours[j]] = rowsd/rowmean
    }
  }
  sum_res = data.frame(var=unlist(lapply(1:length(myvar), 
                                         function(i) rowMeans(myvar[[i]]))),
                       cv=unlist(lapply(1:length(mycoefvar),
                                        function(i) rowMeans(mycoefvar[[i]]))),
                       cutoff=paste0('top',rep(topks,sapply(1:length(myvar),
                                                            function(i) nrow(myvar[[i]])))))
  return(sum_res)
}

summarize_spls_res = function(selected_genes,dataset){
  myvar = mycoefvar = list()
  for (k in 1:length(selected_genes)){
    myvar[[k]] = mycoefvar[[k]] = array(NA,c(length(unlist(selected_genes[[k]])),length(hours)-1),
                                        dimnames=list(unlist(selected_genes[[k]]),hours[-1]))
    for (j in 2:length(hours)){
      rowdat = dataset[unlist(selected_genes[[k]]),] %>% 
        dplyr::select(contains(paste0(hours[j],'h_'))) %>%
        as.matrix() 
      myvar[[k]][unlist(selected_genes[[k]]),hours[j]] = rowdat %>% rowVars()
      rowsd = rowdat %>% rowSds()
      rowmean = rowdat %>% rowMeans()
      mycoefvar[[k]][unlist(selected_genes[[k]]),hours[j]] = rowsd/rowmean
    }
  }
  sum_res = data.frame(var=unlist(lapply(1:length(myvar), 
                                         function(i) rowMeans(myvar[[i]]))),
                       cv=unlist(lapply(1:length(mycoefvar),
                                        function(i) rowMeans(mycoefvar[[i]]))),
                       cutoff=paste0('top',rep(topks,sapply(1:length(myvar),
                                                            function(i) nrow(myvar[[i]])))))
  return(sum_res)
}


Venn_method = function(cutoff,p1,p2,Genes){
  assay1.pcut = 10^(cutoff)#10^(-13)
  assay2.pcut = 10^(cutoff)#10^(-13)
  Venn = p1 <assay1.pcut & p2 < assay2.pcut
  selected_genes = Genes[Venn]
  return(list(gene_id=selected_genes,
              idx=which(Venn)))
}

Edgington_method = function(topk,p1,p2,p3,Genes){
  Edgington = rank(p1 + p2 + p3) < topk
  selected_genes = Genes[Edgington]
  return(list(gene_id=selected_genes,
              idx=which(Edgington)))
}
Fisher_method = function(topk,p1,p2,p3,Genes){
  Fisher = rank(log(p1) + log(p2) + log(p3)) < topk
  selected_genes = Genes[Fisher]
  return(list(gene_id=selected_genes,
              idx=which(Fisher)))
}


comp_cor = function(input, inputType="index"){
  input_size = length(input)
  if (inputType=="index"){
    myids = input
  }else if (inputType=="gene_id"){
    myids = rep(NA,input_size)
    for (i in 1:input_size){
      myids[i] = which(common$gene_ids==input[i])
    }
  }else if (inputType=="protein_index"){
    myids = rep(NA,input_size)
    for (i in 1:input_size){
      myids[i] = which(common$gene_ids==pro_lfc[input[i],"gene_ids"])
    }
  }
  res_cor = rep(NA,input_size)
  names(res_cor) = common[myids,'gene_ids']
  for (myid in myids){
    gene = common[myid,'gene_ids']
    datCounts1 = data.frame(lfc=rna_lfc[gene,],
                            hours=rep(hours, each=4),
                            replicate=rep(c("A","B","C","D"),ncol(rna_lfc)/4))
    mypro_lfc = mypro = common[myid,] %>% dplyr::select(imputed.log2.LFQ.intensity.00h_1:imputed.log2.LFQ.intensity.20h_4)
    for (h in 1:length(hours)){
      mypro_lfc[,hour_idx[[h]]] = mypro[,hour_idx[[h]]]-
        mypro[,hour_idx[[1]]]
    }
    datCounts2 = data.frame(lfc=as.numeric(mypro_lfc),
                            hours=rep(hours, each=4),
                            replicate=rep(c("A","B","C","D"),ncol(rna_lfc)/4))
    avgcor = mean(sapply(unique(datCounts2$replicate),
                         function(x) cor(datCounts2[datCounts2$replicate==x,'lfc'],
                                         datCounts1[datCounts1$replicate==x,'lfc'])))
    res_cor[gene] = avgcor
  }
  return(res_cor)
}

comp_range = function(input){
  return(rowMaxs(input)-rowMins(input))
}
