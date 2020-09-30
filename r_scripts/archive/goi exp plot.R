# Single expression plot for genes of interest

setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-4/ARSeq_figure/")

# Library

# Import data
exp <- read.csv("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-4/ARSeq_figure/normalized_data.csv", row.names = 1, header = T)
meta <- read.csv("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/round-4/metadata_aj.csv", row.names = 1, header = T)

# Function to generate the plot
goi_plot <- function(data,meta,goi,intgroup, wes_palette="GrandBudapest1"){
  require(ggplot2)
  require(reshape2)
  require(wesanderson)

  # gene of interest
  d <- data [goi,,drop=FALSE]
  d <- melt(d)
  d <- merge(d, meta, by.x= "variable", by.y = "row.names")
  d <- d[order(d[,intgroup]),]


  # Maintain order
  #Turn your 'treatment' column into a character vector
  d$variable <- as.character(d$variable)
  #Then turn it back into a factor with the levels in the correct order
  d$variable <- factor(d$variable, levels=unique(d$variable))

  # Bar plot
  ggplot(d, aes(x=variable,y=value,fill=d[,intgroup]))+ geom_bar(position="dodge", stat="identity")+theme_classic()+
    ggtitle(goi)+ theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),
                        axis.text.x = element_text(face="bold", color="#6C7A89", 
                                                   size=14),
                        axis.text.y = element_text(face="bold", color="#6C7A89", 
                                                   size=14)) +
    scale_fill_manual(values=wes_palette(n= length(unique(d[,intgroup])), name=wes_palette))+
    #scale_fill_manual(values=c("#BFBFBF", "#6C7A89"))+
    #scale_fill_brewer(palette="Greys")+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(y = "Gene Expression (log2)")
}
goi_plot (data=exp,meta=meta,goi='HGF',intgroup='zone_aj1',wes_palette='Cavalcanti1')


# Gene of intereed collapsed
goi_collapsed <- function(data,meta,goi,intgroup, groups=NULL,palette='Dark2'){
  require(ggplot2)
  require(reshape2)
  
  # gene of interest
  d <- data [goi,,drop=FALSE]
  d <- merge(t(d), meta, by = "row.names")

  # Identify the samples of interest
  columns_to_keep <- c(goi, intgroup)
  if (is.null(groups)){
    groups = unique(d[,intgroup])
  }
  d <- d[d[,intgroup] %in% unlist(groups) , colnames(d) %in% columns_to_keep,drop=FALSE]
  
  # Aggregate values
  d[, intgroup] <- as.factor(d[, intgroup])
  d <- melt(d)
  # Calculate SE
  d = summarySE(d, measurevar="value", groupvars=c(intgroup,"variable"))
  # plot
  ggplot(d, aes(x=d[,intgroup],y=value,fill=d[,intgroup]))+ geom_bar(stat="identity")+theme_classic()+
    geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2) +
    ggtitle(goi)+ theme(plot.title = element_text(hjust = 0.5))+
    scale_fill_brewer(palette=palette)+
    #scale_fill_brewer(palette="PuBuGn")+
    theme(legend.position="none", axis.text.x = element_text(angle = 75, hjust = 1,face = "bold",size = 10))

}

goi_collapsed(data,meta,goi="MX1",intgroup="zone")

goi_line <- function(data,meta,goi,intgroup,groups=NULL,palette='Dark2',order=NULL){
  require(ggplot2)
  require(reshape2)
  
  # gene of interest
  d <- data [goi,,drop=FALSE]
  d <- merge(t(d), meta, by = "row.names")
  
  # Identify the samples of interest
  columns_to_keep <- c(goi, intgroup)
  if (is.null(groups)){
    groups = unique(d[,intgroup])
  }
  d <- d[d[,intgroup] %in% unlist(groups) , colnames(d) %in% columns_to_keep,drop=FALSE]
  
  # Aggregate values
  d[, intgroup] <- as.factor(d[, intgroup])
  d <- melt(d)
  # Calculate SE
  d = summarySE(d, measurevar="value", groupvars=c(intgroup,"variable"))
  
  # reorder the data 
  if (!is.null(order)){
    target <- unique(as.character(meta[,intgroup]))
    d <- d[order(unlist(sapply(d[,intgroup], function(x) which(target == x)))),]
    row.names(d) <- seq(1: nrow(d))
    d[, intgroup] <- as.character(d[, intgroup])
    d[, intgroup] <- factor( d[, intgroup], levels =  unique(d[, intgroup]))
  }
  
  
  # line plot
  ggplot(data=d, aes(x=d[,intgroup], y= value, group=variable)) + geom_line(aes(color=variable),size=1.2) + 
    geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2) +
    geom_point()+ 
    theme(legend.position="right", axis.text.x = element_text(angle = 75, hjust = 1,face = "bold",size = 12))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("") + ylab("TMT")
  
}

#CXCR4
#CXCL12
#CCL19
#CSTA
#JUN
#DSG1, SERPINB5, TACSTD2, 

# Function to plot average log2 expression of a gene within samples based on two groups
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
average_goi_plot <- function(data, meta, goi, intgroup, groups){
  require(ggplot2)
  require(reshape2)
  # gene of interest
  d <- data [goi,,drop=FALSE]
  d <- merge(t(d), meta, by = "row.names")

  # Identify the samples of interest
  columns_to_keep <- c(goi, intgroup)
  d <- d[d[,intgroup] %in% unlist(groups) , colnames(d) %in% columns_to_keep,drop=FALSE]
  d$grouping <- ifelse(d[,intgroup] %in% groups$A, "Group_A", "Group_B")
  d <- d[ , -which(names(d) %in% c(intgroup))]
  # Aggregate values
  d$grouping <- as.factor(d$grouping)
  d <- melt(d)
  # Convert Group_B to negative
  d$value_m <- ifelse(d$grouping == "Group_B", d$value * -1, d$value)
  # Calculate SE
  d = summarySE(d, measurevar="value_m", groupvars=c("grouping","variable"))
  
  # Maintain order
  #Turn your 'treatment' column into a character vector
  d$variable <- as.character(d$variable)
  #Then turn it back into a factor with the levels in the correct order
  d$variable <- factor(d$variable, levels=unique(d$variable))
  
  #d=aggregate(.~grouping, data=d, mean)
  #rownames(d)=d[,1]
  #d=d[,-1]
  # Convert second row to negative values
  #d[2,2:ncol(d)] <- d[2,2:ncol(d)]* -1
  
  # Plot
  g1_legend = do.call(paste, c(as.list(groups[[1]]), sep = " + "))
  g2_legend = do.call(paste, c(as.list(groups[[2]]), sep = " + "))
  ggplot(d, aes(x=variable, y=value_m, fill=grouping)) +
    geom_bar(stat="identity", position="identity")  + ylab("log2 expression") + coord_flip() +
    geom_errorbar(aes(ymin=value_m-se, ymax=value_m+se),width=.2) +
    geom_hline(yintercept = 0) + theme(panel.grid = element_blank()) +
    xlab ("") + scale_y_continuous( labels=function(x)abs(x) ) +
    scale_fill_manual(values=c("#E69F00", "#56B4E9"),
                      name="ROI's", labels=c(g1_legend, g2_legend))
}

goi <- c("CCND1", "CDKN1A", "CDK4", "CDK6", "CDK2", "CDKN1B", "STAT3", "CDKN2A", "RB1", "CTNNB1", "IFNAR1")
goi <- c("TYR", "PMEL", "MLANA", "ACSL3", "QPCT", "S100B", "CDK2", "CAPN3", "TBC1D7", "TFAP2A", "EDNRB")

goi <- c('C3',"C3AR1","MSR1")

goi <- c("IGHG2", "IGHG3", "JCHAIN", "IGKC", "IGLC1", "IGLC2", "IGLC3", "CD19", "CD21")

groups <- list(A = c("B1-tumor A","B2-tumor B","C1-tumor C","C2-tumor D","D -exophytic tumor"), 
               B= c("B3-tumor A boundary","C3-tumor D boundary","B4-brisk TILs","C4-brisk TIL area for tumor D"))
# isitu vs exophytic tumor
groups <-list(A = c("B3-tumor A boundary","C3-tumor D boundary"), B= c("D -exophytic tumor"))
# isitu vs exophytic tumor
groups <-list(A = c("B3-tumor A boundary","C3-tumor D boundary"), B= c("A1-MIS","A2-MIS+infiltrate","A3-infiltrate","A4-resolved"))

groups <-list(A = c("B3-tumor A boundary","C3-tumor D boundary"), B= c("B1-tumor A","B2-tumor B","C1-tumor C","C2-tumor D"))


groups <- list(B = c("B1-tumor A","B2-tumor B","C1-tumor C","C2-tumor D","D -exophytic tumor"), 
               A= c("B4-brisk TILs","C4-brisk TIL area for tumor D"))

average_goi_plot (data=exp,meta=meta,goi=sig$genes,intgroup='zone',groups)

# Look at MIFT program AXL program and stem cell markers



# Fold change plot
foldchange_plot <- function(foldchange, goi, order=TRUE, title= "Fold Change Plot"){
  d <- foldchange[row.names(foldchange) %in% goi,]
  if (order == TRUE){d <- d[order(d$log2FoldChange),]}
  d$orderby <- row.names(d)
  d$orderby <- factor(d$orderby , levels = d$orderby)

  # Plot
  ggplot(d, aes(x=d$orderby,y=log2FoldChange,fill=padj))+
    geom_bar(position="dodge", stat="identity")+theme_classic() +
    coord_flip() + xlab ("") +
    scale_fill_continuous(limits=c(0, 0.05))+
    ggtitle(title)
}


# Cell cycle phase
cell_cycle <- function(data, meta, sig, intgroup, groups){
  require(ggplot2)
  require(reshape2)

  # sub,set groups of intererst
  meta <- meta[meta[,intgroup] %in% groups, ]
  data <- data[ , which(names(data) %in% row.names(meta))]

  # MErge
  d = merge(data,sig, by.x="row.names", by.y="gene")
  row.names(d) <- d[,1]
  d <- d[,-1]
  d$phase <- as.factor(d$phase)
  y = aggregate(.~phase, data=d, mean)
  y = t(y)
  colnames(y) <- y[1,]
  y <- data.frame(y[-1,])
  y$G1_S <- as.numeric(y$G1_S)
  y$G2_M <- as.numeric(y$G2_M)

  # Merge with type
  m <- meta[,intgroup,drop=FALSE]
  z = merge(y,m,by="row.names")

  #Plot
  ggplot(z, aes(x=G1_S, y=G2_M,color=z[,intgroup], shape=z[,intgroup])) +
    geom_point(size=3) + theme_classic()

}


# Run cell cyle analyzer
sig <- read.csv("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/cell_cycle_sig.csv", header = T,stringsAsFactors=FALSE)
groups <- c("A1-MIS","B1-tumor A","B2-tumor B","C1-tumor C","C2-tumor D","D -exophytic tumor")

cell_cycle(data, meta, sig, intgroup="region", groups= c("A","B","C","D"))


# Correkation analysis
gene_corr <- function(data, meta, goi, intgroup, groups){
  require(HiClimR)
  require(ggplot2)

  # Subset metadata
  g1_m = meta[ meta[,intgroup] %in% groups[[1]],]
  g2_m = meta[ meta[,intgroup] %in% groups[[2]],]
  # Subset data
  g1_d = data [,colnames(data) %in% rownames(g1_m)]
  g2_d = data [,colnames(data) %in% rownames(g2_m)]
  # Remove all zero
  g1_d = g1_d[apply(g1_d[,-1], 1, function(x) !all(x==0)),]
  g2_d = g2_d[apply(g2_d[,-1], 1, function(x) !all(x==0)),]
  # Run correlation
  g1_corr = data.frame(fastCor(as.matrix(t(g1_d)), optBLAS = TRUE))
  g2_corr = data.frame(fastCor(as.matrix(t(g2_d)), optBLAS = TRUE))
  # Gene of interesr
  g1_goi <- g1_corr[,goi,drop=F]
  g2_goi <- g2_corr[,goi,drop=F]
  # Merge both groups
  final = merge(g1_goi, g2_goi, by="row.names")
  row.names(final) = final[,1]
  final = final[,-1]
  g1_legend = do.call(paste, c(as.list(groups[[1]]), sep = " + "))
  g2_legend = do.call(paste, c(as.list(groups[[2]]), sep = " + "))
  colnames(final) = c(g1_legend,g2_legend)
  # Drop the goi
  final = final[!(row.names(final) %in% goi), ]
  
  # Remove NA rowa
  final = final[complete.cases(final), ]
  
  # threshold
  final = final[final[,1] >= 0.5 | final[,1] <= -0.5 | final[,2] >= 0.5 | final[,2] <= -0.5,]
  
  ggplot(final, aes(x=final[,1], y=final[,2])) +
    geom_point(size=0.5) 
  
  # Return
  return(final)
  
}


goi = "KDM5B"
groups <-list(A = c("B3-tumor A boundary","C3-tumor D boundary"), B= c("A1-MIS","A2-MIS+infiltrate","A3-infiltrate","A4-resolved"))

final = gene_corr(data, meta, goi, intgroup="zone", groups)

# write file out
write.table(final, file="S100B_correaltion.csv", sep=',')

sig <- read.csv("/Users/aj/Dropbox (Partners HealthCare)/Data/melanoma_rarecyte/cell_cycle_genes.csv",header = F, stringsAsFactors=FALSE)
sig <- read.csv("immune_related.csv",header = T, stringsAsFactors=FALSE)

# Find correlation for a gene
library(HiClimR)
exp <- read.csv(file= "ARSeq/normalized_data_stable.csv", row.names = 1, header = T,check.names = F)
exp <- log2 (data)
g_corr = data.frame(fastCor(as.matrix(t(exp)), optBLAS = TRUE))

goi = "S"
g_goi <- g_corr[,goi,drop=F]
g_goi <- g_goi[order(g_goi[,1],decreasing = T),,drop=F]
head(g_goi, 20)
tail(g_goi, 100)
write.table(g_goi, file = "ITIH4_correlated_genes.csv", sep = ',')

goi_collapsed(data=h_data,meta=h_meta,goi=goi,intgroup="zone_1")
