
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T) 
}

if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="OTU/ASV table in reads count [default %default]"),
    make_option(c("-d", "--design"), type="character", default="metadata.tsv",
                help="Design file or metadata [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt",
                help="Taxonomy file [default %default]"),
    make_option(c("-T", "--topNtax"), type="character", default="tax/tax_phylum.top",
                help="Top 10 phylum [default %default]"),    
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-c", "--compare"), type="character", default="KO-OE",
                help="Groups comparison [default %default]"),
    make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
                help="Threshold of P-value [default %default]"),
    make_option(c("-f", "--fdr"), type="numeric", default=0.1,
                help="Threshold of FDR [default %default]"),
    make_option(c("-o", "--output"), type="character", default="compare/",
                help="Output directory; name according to compare and  output.txt and .pdf [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=8,
                help="Figure width [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=5,
                help="Figure heidth [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  suppressWarnings(dir.create(opts$output))
  output_dir = opts$output
  

  opts$output=paste(opts$output,opts$compare, sep = "")
  

  print("Parameters are as follows. Please check it!")
  print(paste("The input data matrix file is ", opts$input,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The taxonomy file is ", opts$taxonomy,  sep = ""))
  print(paste("Top 10 phylum file is ", opts$topNtax,  sep = ""))
  print(paste("Group name is ", opts$group,  sep = ""))
  print(paste("Group compare is ", opts$compare,  sep = ""))
  print(paste("Threshold of P-value is ", opts$pvalue,  sep = ""))
  print(paste("Threshold of FDR is ", opts$fdr,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file is ", opts$output, "*.pdf/txt", sep = ""))
  print("",quote = F)
}

package_list = c("ggplot2","pheatmap","dplyr","devtools")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

package_list = c("limma","edgeR")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


dat = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "") 

design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char = "") 

design$group=design[,opts$group]

run_edgeR = function(dat,design,compare){

  group_list=strsplit(opts$compare,'-')[[1]]
  idx = design$group %in% group_list
  sub_design=design[idx,]
  sub_dat=as.matrix(dat[,rownames(sub_design)])  

  d = DGEList(counts=sub_dat,group=factor(sub_design$group))
  d = calcNormFactors(d)
  d$samples 
  
  design.mat = model.matrix(~ 0 + factor(sub_design$group))
  rownames(design.mat)=colnames(sub_dat)
  colnames(design.mat)=levels(factor(sub_design$group))
  DAO = estimateDisp(d,design.mat)
  fit = glmFit(DAO,design.mat)
  BvsA <- makeContrasts(contrasts = opts$compare, levels=design.mat)
  lrt = glmLRT(fit,contrast=BvsA)
  
  nrDAO=as.data.frame(topTags(lrt, n=nrow(dat)))
  nrDAO=as.data.frame(nrDAO)
  head(nrDAO)

  nrDAO$logFC=round(nrDAO$logFC,3)
  nrDAO$logCPM=round(nrDAO$logCPM,3)
  nrDAO$level = ifelse(nrDAO$logFC>0 & nrDAO$PValue<opts$pvalue & nrDAO$FDR < opts$fdr, "Enriched",ifelse(nrDAO$logFC<0 & nrDAO$PValue<opts$pvalue & nrDAO$FDR < opts$fdr, "Depleted","NotSig"))
  nrDAO$level=factor(nrDAO$level,levels = c("Enriched","Depleted","NotSig"))

  if (file.exists(opts$taxonomy)){
  tax = read.table(opts$taxonomy, header=T, row.names= 1, sep="\t", comment.char = "") 
  tax = tax[rownames(nrDAO),]
  nrDAO=cbind(nrDAO,tax)
  
  x=nrDAO
  x$otu=rownames(x)
  x$neglogp=-log10(x$PValue)

  x = arrange(x, Phylum, Class, Order, Family, Genus,otu)
  rownames(x) = x$otu
  
  x$tax=gsub("p__","",x$Phylum)
  if (file.exists(opts$topN)){
    topN = read.table(opts$topNtax)
    topN = as.vector(topN$V1)
  }else{
    topN=sort(c("Actinobacteria","Bacteroidetes",
            "Firmicutes","Proteobacteria"))
  }
  print(paste("The top",length(topN)[1],"phylum as legends.", sep=" "))
  print(topN)
  print("",quote = F)
  x$tax=as.vector(x$tax)
  # levels(x$tax)=c(unique(x$tax),"Low Abundance")
  if (length(unique(x$tax)) > length(topN)){
    x[!(x$tax %in% topN),]$tax = "Low Abundance"
  }
  x$otu=factor(x$otu,levels = x$otu)
  FDR = min(x$neglogp[x$level=="Enriched"])
  x$Level=x$level
  x$Phylum=x$tax
  p = ggplot(x, aes(x=otu, y=neglogp, color=Phylum, size=logCPM, shape=Level)) +
    geom_point(alpha=.7) + 
    geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
    scale_shape_manual(values=c(17, 25, 20))+
    scale_size(breaks=c(5, 10, 15)) +
    labs(x="Genus", y="-log10(P)", title=paste(group_list[1], "vs", group_list[2], sep=" ")) +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="right")
  p
  ggsave(paste(opts$output, "_manhattan.pdf", sep=""), p, 
         width = opts$width*2, height = opts$height)
  }
  
  norm = t(t(sub_dat)/colSums(sub_dat,na=T))*100

  colSums(norm)

  A_list = subset(sub_design, group %in% group_list[1])
  A_norm = norm[, rownames(A_list)]
  A_mean = as.data.frame(rowMeans(A_norm))
  colnames(A_mean)=c("MeanA")

  B_list = subset(sub_design, group %in% group_list[2])
  B_norm = norm[, rownames(B_list)]
  B_mean = as.data.frame(rowMeans(B_norm))
  colnames(B_mean)=c("MeanB")
  
  Mean = round(cbind(A_mean, B_mean, A_norm, B_norm),3)
  Mean = Mean[rownames(nrDAO),]   
  output=cbind(nrDAO[,-3],Mean)

  write.table("OTUID\t", file=paste(opts$output,"_all.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  write.table(output,file=paste0(opts$output,"_all.txt",sep=""),append = T,quote = F,sep = '\t',row.names = T)

  NoE= dim(output[output$level=="Enriched",])[1]
  NoD= dim(output[output$level=="Depleted",])[1]

  p = ggplot(output, aes(x=logFC, y=logCPM, color=level)) + 
    geom_point() + xlim(-4, 4) + theme_classic()+
    scale_colour_manual(values=c("red","green","grey")) + 
    labs(x="log2(fold change)", y="log2(count per million)", 
         title=paste(group_list[1], "vs", group_list[2], sep=" "))+ 
    annotate("text",x=-3,y=15,label=paste(NoD,sep=""))+ 
    annotate("text",x=3,y=15,label=paste(NoE,sep=""))
  p
  ggsave(paste(opts$output, "_volcano.pdf", sep=""), p, 
         width = opts$width, height = opts$height)

  output=output[output$PValue < opts$pvalue,]
  output=output[output$FDR < opts$fdr,]

  write.table("OTUID\t", file=paste(opts$output,"_sig.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  write.table(output,file=paste0(opts$output,"_sig.txt",sep=""),append = T,quote = F,sep = '\t',row.names = T)

  g=group_list
    write.table("OTUID\t", file=paste0(output_dir, g[1], "vs", g[2],"_enriched.txt"),append = F, quote = F, eol = "", row.names = F, col.names = F)
  write.table(output[output$level == "Enriched",],file=paste0(output_dir, g[1], "vs", g[2],"_enriched.txt"),append = T,quote = F,sep = '\t',row.names = T)

      write.table("OTUID\t", file=paste0(output_dir, g[2], "vs", g[1],"_enriched.txt"),append = F, quote = F, eol = "", row.names = F, col.names = F)
  write.table(output[output$level == "Depleted",],file=paste0(output_dir, g[2], "vs", g[1],"_enriched.txt"),append = T,quote = F,sep = '\t',row.names = T)
    
  if (file.exists(opts$taxonomy)){

  anno_row=data.frame(Level = output$level, 
                            Taxonomy=output$Phylum, 
                            row.names = rownames(output))

  anno_col=data.frame(Group = sub_design$group,
                      row.names = rownames(sub_design))
  
  pheatmap(norm[rownames(output),],
           scale = "row",
           cutree_rows=2,cutree_cols = 2,
           annotation_col = anno_col, annotation_row = anno_row,
           filename = paste(opts$output, "_heatmap.pdf", sep=""),
           width=opts$width, height=opts$height, 
           annotation_names_row= T,annotation_names_col=F,
           show_rownames=F,show_colnames=T,
           main = paste("Differential abundance OTUs of",group_list[1], "vs", group_list[2],sep=" "),
           fontsize=7,display_numbers=F)
}
}

run_edgeR(dat,design,opts$compare)

print(paste("Output in files in ", opts$output,"*.txt/pdf. All works done!!!", sep = ""))
print("",quote = F)
