# group samples by gene
# for each gene, it will be a data.frame, including grouped copy number, expression values
# make sure you have sample_info, cn.ovl, exprs.ovl, zscore.ovl
group_samples_for_loss <- function(gene){
  data.plot <- sample_info
  cn <- cn.ovl[gene,]
  expr <- exprs.ovl[gene,]
  zscore <- zscore.ovl[gene,]
  is_loss <- sapply(cn,function(x) x == -1 | x == -2)
  cutoff <- median(expr)
  data.plot$expression <- expr
  data.plot$zscore <- zscore
  data.plot$cnv <- as.character(cn)
  data.plot$is.loss <- ifelse(is_loss==TRUE,"L","N")
  data.plot$group.cnv <- paste(data.plot$group,data.plot$is.loss,sep="")
  data.plot$group.cnv <- factor(data.plot$group.cnv,levels=c('BL','BN','CL','CN'))
  data.plot$exprLevel <- ifelse(expr>cutoff,"H",ifelse(expr<=cutoff,"L","O"))
  data.plot$group.exprLevel <- paste(data.plot$group,data.plot$exprLevel,sep="")
  data.plot$loss.exprLevel <- paste(data.plot$is.loss,data.plot$exprLevel,sep="") 
  data.plot$group.cnv.exprLevel <- paste(data.plot$group,data.plot$is.loss,data.plot$exprLevel,sep="")
  return(data.plot)
}

# group samples by gene
# for each gene, it will be a data.frame, including grouped copy number, expression values
# make sure you have sample_info, cn.ovl, exprs.ovl, zscore.ovl
group_samples_for_gain <- function(gene){
  data.plot <- sample_info
  cn <- cn.ovl[gene,]
  expr <- exprs.ovl[gene,]
  zscore <- zscore.ovl[gene,]
  is_gain <- sapply(cn,function(x) x == 1 | x == 2)
  cutoff <- median(expr)
  data.plot$expression <- expr
  data.plot$zscore <- zscore
  data.plot$cnv <- as.character(cn)
  data.plot$is.gain <- ifelse(is_gain==TRUE,"G","N")
  data.plot$group.cnv <- paste(data.plot$group,data.plot$is.gain,sep="")
  data.plot$group.cnv <- factor(data.plot$group.cnv,levels=c('BG','BN','CG','CN'))
  data.plot$exprLevel <- ifelse(expr>cutoff,"H",ifelse(expr<=cutoff,"L","O"))
  data.plot$group.exprLevel <- paste(data.plot$group,data.plot$exprLevel,sep="")
  data.plot$loss.exprLevel <- paste(data.plot$is.gain,data.plot$exprLevel,sep="") 
  data.plot$group.cnv.exprLevel <- paste(data.plot$group,data.plot$is.gain,data.plot$exprLevel,sep="")
  return(data.plot)
}
                    
# count samples by copy number and gene expression groups
count_frequency_of_cn_ex_groups_for_loss <- function(genes.plot,data.plot){

  counts <- lapply(genes.plot,function(x){
    calculate_pvalue <- function(df,groups){
      y = df[df$group.cnv.exprLevel %in% groups,]
      if( length(unique(as.vector(y$group.cnv.exprLevel))) < 2 )
          return(NA)
      freq <- table(y$group.cnv.exprLevel) < 3
      if( sum(freq) > 0 )
          return(NA)
      wilcox.test(expression ~ group.cnv.exprLevel,data=y,exact=FALSE)$p.value
    }
    calculate_foldchange <- function(df,groups){
      g1 <- as.vector(df[df$group.cnv.exprLevel==groups[1],]$expression)
      g2 <- as.vector(df[df$group.cnv.exprLevel==groups[2],]$expression)
      fc <- round(mean(g1)/mean(g2),2)
    }
    calculate_median_foldchange <- function(df,groups){
      g1 <- as.vector(df[df$group.cnv.exprLevel==groups[1],]$expression)
      g2 <- as.vector(df[df$group.cnv.exprLevel==groups[2],]$expression)
      fc <- round(median(g1)/median(g2),2)
    }
    #print(x)
    df=data.plot[[x]]; 
    y = c(nrow(df),
      sum(df$group.cnv.exprLevel=='BLH'),
      sum(df$group.cnv.exprLevel=='BLL'),
      sum(df$group.cnv.exprLevel=='BNH'), 
      sum(df$group.cnv.exprLevel=='BNL'),   
      sum(df$group.cnv.exprLevel=='CLH'),
      sum(df$group.cnv.exprLevel=='CLL'),
      sum(df$group.cnv.exprLevel=='CNH'),
      sum(df$group.cnv.exprLevel=='CNL'),
      calculate_pvalue(df,c('BLL','CNL')),
      calculate_foldchange(df,c('BLL','CNL')),
      calculate_median_foldchange(df,c('BLL','CNL'))
    )
    } ) 

  counts <- as.data.frame(matrix(unlist(counts),nrow=length(genes.plot),ncol=12,byrow=TRUE))
  rownames(counts) <- genes.plot
  colnames(counts) <- c('Total','BLH','BLL','BNH','BNL','CLH','CLL','CNH','CNL',
                      'pval.BLL.vs.CNL','fc.BLL.vs.CNL','fc2.BLL.vs.CNL')

  counts <- counts %>% mutate(
    BL = BLH + BLL,
    BLL.perc = BLL/BL,
    CLL.perc   = CLL/(CLH + CLL), 
    BNL.perc = BNL/(BNH + BNL),    
    CL.perc   = (CLL+CNL)/(CLH + CLL + CNL + CNH),
    CNL.perc   = CNL/(CNL + CNH),
    ratio1 = BLL/BL,
    ratio2 = ratio1/CNL.perc,
    Gene.symbol = genes.plot       
    ) %>% arrange(desc(BLL.perc,ratio1,ratio2))
  print(c("ratio1 = BLL/BL",
        "ratio2 = ratio1/(CNL/CN)"))
  print(nrow(counts))
  return(counts)

}


# count samples by copy number and gene expression groups
count_frequency_of_cn_ex_groups_for_gain <- function(genes.plot,data.plot){

  counts <- lapply(genes.plot,function(x){
    calculate_pvalue <- function(df,groups){
      y = df[df$group.cnv.exprLevel %in% groups,]
      if( length(unique(as.vector(y$group.cnv.exprLevel))) < 2 )
          return(NA)
      freq <- table(y$group.cnv.exprLevel) < 3
      if( sum(freq) > 0 )
          return(NA)
      wilcox.test(expression ~ group.cnv.exprLevel,data=y,exact=FALSE)$p.value
    }
    calculate_foldchange <- function(df,groups){
      g1 <- as.vector(df[df$group.cnv.exprLevel==groups[1],]$expression)
      g2 <- as.vector(df[df$group.cnv.exprLevel==groups[2],]$expression)
      fc <- round(mean(g1)/mean(g2),2)
    }
    calculate_median_foldchange <- function(df,groups){
      g1 <- as.vector(df[df$group.cnv.exprLevel==groups[1],]$expression)
      g2 <- as.vector(df[df$group.cnv.exprLevel==groups[2],]$expression)
      fc <- round(median(g1)/median(g2),2)
    }
    #print(x)
    df=data.plot[[x]]; 
    y = c(nrow(df),
      sum(df$group.cnv.exprLevel=='BGH'),
      sum(df$group.cnv.exprLevel=='BGL'),
      sum(df$group.cnv.exprLevel=='BNH'), 
      sum(df$group.cnv.exprLevel=='BNL'),   
      sum(df$group.cnv.exprLevel=='CGH'),
      sum(df$group.cnv.exprLevel=='CGL'),
      sum(df$group.cnv.exprLevel=='CNH'),
      sum(df$group.cnv.exprLevel=='CNL'),
      calculate_pvalue(df,c('BGH','CNH')),
      calculate_foldchange(df,c('BGH','CNH')),
      calculate_median_foldchange(df,c('BGH','CNH'))
    )
    } ) 

  counts <- as.data.frame(matrix(unlist(counts),nrow=length(genes.plot),ncol=12,byrow=TRUE))
  rownames(counts) <- genes.plot
  colnames(counts) <- c('Total','BGH','BGL','BNH','BNL','CGH','CGL','CNH','CNL',
                      'pval.BGH.vs.CNH','fc.BGH.vs.CNH','fc2.BGH.vs.CNH')

  counts <- counts %>% mutate(
    BG = BGH + BGL,
    BGH.perc = BGH/BG,
    CGH.perc   = CGH/(CGH + CGL), 
    BNH.perc = BNH/(BNH + BNL),    
    CNH.perc   = CNH/(CNH + CNL),
    ratio1 = BGH.perc,
    ratio2 = ratio1/CNH.perc,
    Gene.symbol = genes.plot       
    ) %>% arrange(desc(ratio1,ratio2))
  print(c("ratio1 = BGH/BG",
        "ratio2 = (BGH/BG)/(CNH/CN)"))
  print(nrow(counts))
  return(counts)

}

# plot gene exprs grouped by group.cnv 
plot_gene_exprs <- function(genes.plot,data.genes,gene2band,outdir,output_prefix){

  genes.plot <- genes.plot[genes.plot %in% names(data.genes)]
  print(length(genes.plot))
  data.plot <- lapply(genes.plot,function(x) data.genes[[x]])
  names(data.plot) <- genes.plot

  # plot expression level by group.cnv
  library(ggpubr)
  plots <- lapply(genes.plot,function(x){ 
    give.n <- function(x) return(c(y = max(x)*1.1, label = length(x)))
    p <- ggboxplot(data.plot[[x]],x="group.cnv",y="expression",title=paste(x,gene2band[x]),xlab=FALSE,legend="right") + rotate_x_text(angle = 90) + geom_jitter(aes(colour=exprLevel),width = 0.2,alpha=1,size=1) + 
                stat_summary(fun.data = give.n, geom = "text", fun = max) + labs(color="Exp.\nlevel") + ylab("Expression")
    })
  pdf(paste0(outdir,"/",output_prefix,".pdf"),10,7)
  print(ggarrange(plotlist = plots,nrow=2,ncol=3))
  dev.off()

  # zscore
  plots <- lapply(genes.plot,function(x){ 
    give.n <- function(x) return(c(y = max(x)*1.1, label = length(x)))
    p <- ggboxplot(data.plot[[x]],x="group.cnv",y="zscore",title=paste(x,gene2band[x]),xlab=FALSE,legend="right") + rotate_x_text(angle = 90) + geom_jitter(aes(colour=exprLevel),width = 0.2,alpha=1,size=1) + 
                stat_summary(fun.data = give.n, geom = "text", fun = max) + labs(color="Exp.\nlevel") + ylab("Z-score")
    })
  pdf(paste0(outdir,"/",output_prefix,".zscore.pdf"),10,7)
  print(ggarrange(plotlist = plots,nrow=2,ncol=3))
  dev.off()

}       

# plot genes on chr
make_chr_plot <- function(chrs,counts,colname,ylab="ratio",colorgroup="loss",colors=c('loss'='blue','not loss'='black')){
  library(ggpubr)
  library(scales)
  plots <- lapply(chrs,function(x){ 
           counts.plot <- counts[counts$Gene.chr==x,]
           ggscatter(counts.plot,x="coordinates",y=colname,size=0.5,color=colorgroup) + ggtitle(x) + xlab("Chromosome coordinates") + ylab(ylab) + labs(color="Cytobands") + scale_color_manual(values=colors) + scale_y_continuous(breaks= pretty_breaks())
          } )
  pdf(paste0(outdir,"/chrplots_",colname,".pdf"),10,10)
  print(ggarrange(plotlist = plots,ncol=2,nrow=4))
  dev.off()
}

# plot genes on cytobands
make_band_plot <- function(bands,counts,colname,ylab="ratio",colorgroup="loss",colors=c('loss'='blue','not loss'='black')){
  library(ggpubr)
  plots <- lapply(bands,function(x){ 
           counts.plot <- counts[counts$chr_cb==x,]
           ggscatter(counts.plot,x="coordinates",y=colname,size=0.5,color=colorgroup) + ggtitle(x) + xlab("Chromosome coordinates") + ylab(ylab) + labs(color="CNA") + scale_color_manual(values=colors) + scale_y_continuous(breaks= pretty_breaks())
          } )
  pdf(paste0(outdir,"/bandplots_",colname,".pdf"),10,10)
  print(ggarrange(plotlist = plots,ncol=2,nrow=4))
  dev.off()
}

add_cosmic_annotation <- function(table,by.x="gene",by.y="Gene.Symbol"){
  cosmic <- read.csv("data/COSMIC/Cancer_Gene_Census/cancer_gene_census.csv")
  table <- merge(table,cosmic[,c('Gene.Symbol','Name','Role.in.Cancer')],by.x=by.x,by.y=by.y,all.x=T,sort=F)
}
                      