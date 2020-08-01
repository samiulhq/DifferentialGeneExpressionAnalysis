getDEgenes <- function(countTable, group, min.cpm=2, n.min.cpm =3,FDR=0.05, batch_factor = c(1,1,1,1,1,1),lfc=0,batchEffect=TRUE,sample='unknown'){
##### This function returns a list of DE genes with associated qvalues(FDR) from multiple samples. batch effects are also considered. 
  #INPUT ARGS
  #counTable gene counts for all the samples
  #group treatment group for each sample
  #min.cmp minimum count per million a gene should have to be considered
  #n.min.cpm  number of samples that should have at least min.cpm gene expression
  #FDR false discovery rate (q value based on BH correction)
  #batch_factor = if samples were prepared at different times they may have batch efects
  #lfc  log fold change threshold for DE genes
  #batchEffect if TRUE batch_factor is used 
  #Sample : Name of the experiment 
  
   require("edgeR")
  require("baySeq")
#  require(Rsubread)
  require(RColorBrewer)
  require(gplots)
  require(gdata)
  
  ## Begin DE analysis
  len = read.table("../Data/tair10.genelength_featurecount.txt",sep="\t",fill = TRUE)
  len <- len[order(len$V1),] 
  y <- DGEList(counts = countTable, group=group, genes=data.frame(Length=len))
  #str(y)
  #dim(y)
  
  
  y$samples$lib.size <- colSums(y$counts)
  
  min.cpm <- as.double(cpm(10,mean(y$samples$lib.size))) # min.cpm is set in such way that raw counts is atleast approximately 10 accross samples. 
  # paramenter for filtering low expressed genes
  #min.cpm <- 2
  #n.min.cpm <- 3
  keep <- rowSums(cpm(y)>min.cpm) >= n.min.cpm
  table(keep)
  y <- y[keep,]
  dim(y) 
  #Counts <- Counts[keep,]

  # Check distribution of samples and library sizes
  # sizes
  y$samples$lib.size
  # plot sizes
  #barplot(y$samples$lib.size,names=colnames(y),las=2)
  
  
  # Get log2 counts per million
  logcounts <- cpm(y,log=TRUE)
  # Check distributions of samples using boxplots
 # boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
  # Let's add a blue horizontal line that corresponds to the median logCPM
  #abline(h=median(logcounts),col="blue")
 # title("Boxplots of logCPMs (unnormalized)")
  
  # TMM normalization
  y <- calcNormFactors(y, method="TMM")
  y$samples
  fname=paste("BCVnew_",names(cm)[1],".pdf",sep='')
 # pdf(fname)
  plotMDS(y, method="bcv", col=as.numeric(y$samples$group))
  plotMDS(y, col=as.numeric(y$samples$group))
  # print(distmat)
  #dev.off()
  
#  keep <- keepcol
#  y  <- y[,keep]
  y$samples$group <- relevel(y$samples$group, ref="plusFe")
  
  group_factor <- y$samples$group
  group_factor <-relevel(group_factor, ref="plusFe")
  batch_factor=factor(batch_factor)
  if(batchEffect){
  design <- model.matrix(~group_factor+batch_factor,data=y$samples)
  }
  else
  {
    design <- model.matrix(~group_factor,data=y$samples)
  }
  
  print(design)
 
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design,robust = TRUE)

  plotQLDisp(fit)

  if(batchEffect){
  #qlf <- glmQLFTest(fit,contrast = c(1,-1,rep(0,ncol(design)-2)))
  qlf <- glmQLFTest(fit,coef=2)
  print(summary(dup<-decideTestsDGE(qlf,p.value = FDR,lfc=lfc)))
  }
  else{
    qlf <- glmQLFTest(fit,coef=2)
    print(summary(dup<-decideTestsDGE(qlf,p.value = FDR,lfc=lfc)))
    
  }
  
  
  # 
  # keg <- kegga(qlf, species="At")
  # print(topKEGG(keg, sort="up"))
  # 
  ###

  all_genes <- topTags(qlf, n=Inf,sort.by = 'none')
  values <- c("Up =1","No change=0","Down=-1")

  z=rep("N",times=33602)
  
  ind1=which(dup[,1]==1) #up regulated genes
  ind2=which(dup[,1]==-1) # down regulated genes
  
  z[ind1]="U"
  z[ind2]="D"
  
  table(dup)
  # if(length(ind1)>length(ind2)){
  #   hl.col = c("blue","red")
  # }
  # else
  # {
  #   hl.col = c("red","blue")
  # }
  hl.col = c("blue","red")
  
  if(batchEffect){
  plotMD(qlf,status=z, values<-c("U","D"), hl.col = hl.col,main = paste("Differential expression (-Fe vs. +Fe) with batch effect adjusted. Sample :", sample), xlab = "Average log expression", ylab = "log-fold Change")
  }
  else{
    plotMD(qlf,status=z, values<-c("U","D"), hl.col = hl.col,main = paste("Differential expression (-Fe vs. +Fe). Sample :", sample), xlab = "Average log expression", ylab = "log-fold Change")
    
  }
  
  plotSmear(qlf, de.tags=rownames(qlf)[dup!=0])
  
  
  abline(h=c(-1, 1), col="green")
  #geneNames <- all_genes$table[,1]
  #qlfDEGenes <- geneNames[sort(c(ind1,ind2))]
  ind=which(dup!=0)
  colnames(all_genes$table)[1]<-"GeneId"
  qlfDEGenes <- all_genes$table[ind,c("GeneId","FDR")]
  df <- data.frame(qlfDEGenes)
  #qvalues for the genes
  names(df) = c("GeneId",'qval')
  return (df)
  
}