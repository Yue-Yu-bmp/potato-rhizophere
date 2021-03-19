site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"

if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="Feature table [default %default]"),
    make_option(c("-T", "--thre"), type="numeric", default=0,
                help="Threshold of abundance, such as 0.001 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="metadata.tsv",
                help="Experiment design or sample metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Column name of group [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/otutab_mean.txt",
                help="output directory and prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))

  if (opts$output==""){
    opts$output=paste(opts$input, "_", opts$group, sep = "")}

  print(paste("Feature table: ", opts$input,  sep = ""))
  print(paste("Metadata: ", opts$design,  sep = ""))
  print(paste("Group name: ", opts$group,  sep = ""))
  print(paste("Abundance threshold: ", opts$thre,  sep = ""))
  print(paste("Output filename: ", opts$output, sep = ""))
}

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list <- c("dplyr")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

otutab = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F)

design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char = "", stringsAsFactors = F) 

design$group=design[,opts$group]

norm = t(otutab)/colSums(otutab,na=T)*100

idx = colMeans(norm) > opts$thre
HA = norm[,idx]

merge=cbind(HA, design[,c("group"),drop=F])
HA_group_mean = merge %>% group_by(group) %>% summarise_all(mean)
HA_t = as.data.frame(t(HA_group_mean))
HA_t = cbind(HA_t,  c("All", colMeans(norm)))

rownames(HA_t)[1] = "OTUID"
write.table(HA_t, file=paste(opts$output, "", sep = ""), append = FALSE, sep="\t", quote=F, row.names=T, col.names=F)

