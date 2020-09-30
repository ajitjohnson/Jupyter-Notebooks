# Genes of interest bar plot COVID19
library(ggpubr)

setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/NHP")

# Load data
data <- read.csv(file= "exp_cleaned.csv", header = T, check.names=F, row.names = 1)
meta <- read.csv(file= "meta_cleaned.csv", row.names = 1, header = T,check.names=F)

# Human
setwd("/Users/aj/Dropbox (Partners HealthCare)/Data/covid19/human/proteomics/")
data <- round(read.csv(file= "exp_cleaned.csv", header = T, check.names=F, row.names = 1),0)
meta <- read.csv(file= "meta_cleaned.csv", row.names = 1, header = T,check.names=F)


d4 <- read.csv(file= "archive/exp_cleaned_2.csv", header = T, check.names=F, row.names = 1)
d4_meta <- read.csv(file= "archive/meta_cleaned.csv", header = T, check.names=F, row.names = 1)

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
  #datac <- rename(datac, c("mean" = measurevar))
  #colnames(datac[,"mean"]) <- measurevar
  names(datac)[names(datac) == 'mean'] <- measurevar
  
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
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),
          axis.text.x = element_text(face="bold", color="#000000", 
                                     size=14),
          axis.text.y = element_text(face="bold", color="#000000", 
                                     size=14))+
    geom_hline(yintercept = 0) + theme(panel.grid = element_blank()) +
    xlab ("") + scale_y_continuous( labels=function(x)abs(x) ) +
    scale_fill_manual(values=c("#E69F00", "#56B4E9"),
                      name="ROI's", labels=c(g1_legend, g2_legend))
}
goi_plot <- function(data,meta,goi,intgroup, wes_palette="GrandBudapest1"){
  require(ggplot2)
  require(reshape2)
  require(wesanderson)
  require(dplyr)
  
  # gene of interest
  d <- data [goi,,drop=FALSE]
  d <- melt(d)
  dd <- merge(d, meta, by.x= "variable", by.y = "row.names")
  d <- dd %>% slice(match(d$variable, variable))
  #d <- d[order(d[,intgroup]),]
  
  
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
    #scale_fill_manual(values=wes_palette(n= length(unique(d[,intgroup])), name=wes_palette))+
    #scale_fill_manual(values=c("#BFBFBF", "#6C7A89"))+
    #scale_fill_brewer(palette="Greys")+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(y = "Counts")
}
goi_line <- function(data,meta,goi,intgroup,groups=NULL,palette='Dark2',order=NULL,yaxis_text="Log2 Expression"){
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  
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
  ggplot(data=d, aes(x=d[,intgroup], y= value, group=variable)) + geom_line(aes(color=variable),size=2.4) + 
    geom_errorbar(aes(ymin=value-se, ymax=value+se),width=0.2) + scale_colour_discrete("")+
    geom_point()+ 
    theme(legend.position="right")+ 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = "bold",size = 24),
          axis.text.y = element_text(face = "bold",size = 24),
          legend.text = element_text(face = "bold",size = 18),
          axis.title=element_text(size=18,face="bold")) +
    xlab("") + ylab(yaxis_text)

  
}
goi_collapsed <- function(data,meta,goi,intgroup, groups=NULL,palette='Dark2',order=NULL){
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
  
  # reorder the data 
  if (!is.null(order)){
    target <- unique(as.character(meta[,intgroup]))
    d <- d[order(unlist(sapply(d[,intgroup], function(x) which(target == x)))),]
    row.names(d) <- seq(1: nrow(d))
    d[, intgroup] <- as.character(d[, intgroup])
    d[, intgroup] <- factor( d[, intgroup], levels =  unique(d[, intgroup]))
  }
  
  # Maintain order
  #Turn your 'treatment' column into a character vector
  d$variable <- as.character(d$variable)
  #Then turn it back into a factor with the levels in the correct order
  d$variable <- factor(d$variable, levels=unique(d$variable))
  
  # Calculate SE
  d = summarySE(d, measurevar="value", groupvars=c(intgroup,"variable"))
  # plot
  ggplot(d, aes(x=d[,intgroup],y=value,fill=d[,intgroup]))+ geom_bar(stat="identity")+theme_classic()+
    geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2) +
    ggtitle(goi)+ theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank(),
          axis.text.x = element_text(face="bold", color="#000000", 
                                     size=14),
          axis.text.y = element_text(face="bold", color="#000000", 
                                     size=14))+
    scale_fill_brewer(palette=palette)+
    #scale_fill_brewer(palette="PuBuGn")+
    theme(legend.position="none", axis.text.x = element_text(angle = 75, hjust = 1,face = "bold",size = 10))
  
}
goi_stacked <- function(data,meta,goi,c_group='c1', p_val_pos=0,p_val_size=15,intgroup,gene_size=10,groups=NULL,order=NULL,ttest.ref.group=NULL){
  
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  require(patchwork)
  
  c1 <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
  c1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    # ggplot function
  ggplot_func <- function(data,meta,goi,intgroup,palette,groups,ttest.ref.group,gene_size,p_val_pos,p_val_size){
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
    
    # reorder the data 
    if (!is.null(order)){
      target <- unique(as.character(meta[,intgroup]))
      d <- d[order(unlist(sapply(d[,intgroup], function(x) which(target == x)))),]
      row.names(d) <- seq(1: nrow(d))
      d[, intgroup] <- as.character(d[, intgroup])
      d[, intgroup] <- factor( d[, intgroup], levels =  unique(d[, intgroup]))
    }
    
    # ggplot
    q <- ggplot(d, aes(x=d[, intgroup], y=value, fill=d[, intgroup])) + 
      #geom_boxplot() + 
      geom_bar(stat = "summary") +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
      scale_y_continuous(name=goi,limits = c(0, max(d$value)+quantile(d$value, 0.05)[[1]]), breaks= c(1, round(max(d$value)),1)) +
      theme(legend.position = "none", 
            plot.title = element_blank(),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.title.y = element_text(size = gene_size, angle = 0, vjust = 0.5, hjust = 1, face='bold'), 
            axis.text.y = element_text(size = rel(1)), 
            plot.margin= grid::unit(c(0, 0, 0, 0), "cm"),
            #panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")
      )  + scale_fill_manual(values = c1) + #scale_fill_grey(start = 0.5, end = 0.9) +
      stat_compare_means(label = "p.signif", 
                         method = "t.test", 
                         ref.group = ttest.ref.group, hide.ns = T,
                         color = '#FD5403',label.y.npc=p_val_pos, size = p_val_size
                         #position = position_nudge(y=-quantile(d$value, 1))
                         ) #position = position_nudge(y=-quantile(d$value, 0.01))
    # retuen the plot
    return(q)
    
  }
  
  goi = as.list(goi)
  
  plot_list<- purrr::map(goi, function(x) ggplot_func(data=data,meta=meta,goi=x,
                         intgroup=intgroup,palette=palette,gene_size=gene_size,
                         groups=groups,ttest.ref.group=ttest.ref.group,
                         p_val_pos=p_val_pos,p_val_size=p_val_size))
  # Plot 
  patchwork::wrap_plots(plotlist = plot_list, ncol = 1) 
  
}
goi_stacked_single <- function(data,meta,goi,c_group='c1', p_val_pos=0,p_val_size=15,intgroup,gene_size=10,groups=NULL,order=NULL,ttest.ref.group=NULL){
  
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  require(patchwork)
  
  c1 <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
  c1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # ggplot function
  ggplot_func <- function(data,meta,goi,intgroup,palette,groups,ttest.ref.group,gene_size,p_val_pos,p_val_size){
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
    
    # reorder the data 
    if (!is.null(order)){
      target <- unique(as.character(meta[,intgroup]))
      d <- d[order(unlist(sapply(d[,intgroup], function(x) which(target == x)))),]
      row.names(d) <- seq(1: nrow(d))
      d[, intgroup] <- as.character(d[, intgroup])
      d[, intgroup] <- factor( d[, intgroup], levels =  unique(d[, intgroup]))
    }
    
    # ggplot
    q <- ggplot(d, aes(x=d[, intgroup], y=value, fill=d[, intgroup])) + 
      #geom_boxplot() + 
      geom_bar(stat = "summary") +
      geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
      ggtitle(goi)+ theme(plot.title = element_text(hjust = 0.5))+
      #scale_y_continuous(name=goi,limits = c(0, max(d$value)+quantile(d$value, 0.05)[[1]]), breaks= c(1, round(max(d$value)),1)) +
      theme(legend.position = "none", 
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(face="bold", color="#000000", size=14,angle =- 90, vjust = 0.5),
            axis.text.y = element_text(face="bold", color="#000000", size=14),
            #axis.text.x = element_blank(), 
            #axis.ticks.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            #axis.title.y = element_text(size = gene_size, angle = 0, vjust = 0.5, hjust = 1, face='bold'), 
            plot.margin= grid::unit(c(0, 0, 0, 0), "cm"),
            #panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")
      )  + scale_fill_manual(values = c1) + #scale_fill_grey(start = 0.5, end = 0.9) +
      stat_compare_means(label = "p.signif", 
                         method = "t.test", 
                         ref.group = ttest.ref.group, hide.ns = T,
                         color = '#FD5403',label.y.npc=p_val_pos, size = p_val_size
                         #position = position_nudge(y=-quantile(d$value, 1))
      ) #position = position_nudge(y=-quantile(d$value, 0.01))
    # retuen the plot
    return(q)
    
  }
  
  goi = as.list(goi)
  
  plot_list<- purrr::map(goi, function(x) ggplot_func(data=data,meta=meta,goi=x,
                                                      intgroup=intgroup,palette=palette,gene_size=gene_size,
                                                      groups=groups,ttest.ref.group=ttest.ref.group,
                                                      p_val_pos=p_val_pos,p_val_size=p_val_size))
  # Plot 
  patchwork::wrap_plots(plotlist = plot_list, ncol = 1) 
  
}
goi_fc <- function(goi, data, meta, intgroup, fc_base, groups=NULL, order=NULL){
  require(magrittr)
  require(dplyr)
  require(ComplexHeatmap)
  require(circlize)
  
  # gene of interest
  d <- data [row.names(data) %in% goi,,drop=FALSE]
  d <- d[goi,]
  m <- meta[,intgroup,drop=F]
  d <- merge(t(d), m, by = "row.names")
  row.names(d) <- d[,1]
  d <- d[,-1]
  
  # Identify the samples of interest
  columns_to_keep <- c(goi, intgroup)
  if (is.null(groups)){
    groups = unique(d[,intgroup])
  }
  d <- d[d[,intgroup] %in% unlist(groups) , colnames(d) %in% columns_to_keep,drop=FALSE]
  
  # reorder the data 
  if (!is.null(order)){
    target <- unique(as.character(meta[,intgroup]))
    d <- d[order(unlist(sapply(d[,intgroup], function(x) which(target == x)))),]
    #row.names(d) <- seq(1: nrow(d))
    d[, intgroup] <- as.character(d[, intgroup])
    d[, intgroup] <- factor( d[, intgroup], levels =  unique(d[, intgroup]))
  }
  
  # Calculate the mean for foldchnage calculation
  # Aggregate values
  d[, intgroup] <- as.factor(d[, intgroup])
  d_collapsed <- aggregate(.~d[, intgroup], data=d, mean)
  row.names(d_collapsed) <- d_collapsed[,1]
  d_collapsed <- d_collapsed[,-1]
  d_collapsed <- data.frame(t(d_collapsed[ , -which(names(d_collapsed) %in% intgroup)]))
  # Calculate foldchange
  fc <- data.frame(d_collapsed %>% sapply(`/`, d_collapsed[,fc_base]))
  row.names(fc) <- row.names(d_collapsed)
  # drop the fc_base column
  fc <- fc[ , -which(names(fc) %in% fc_base), drop=F]
  
  # Calculate the significant foldchanges
  
  
  pval_calculator <- function(goi){
    x = d[,c(goi,intgroup)]
    p_value <- c()
    for (i in unique(x[, intgroup])){
      # subset the groups of interest
      base <- x[x[, intgroup] == fc_base,][,1] 
      contrast <- x[x[, intgroup] == i,][,1] 
      # t-test
      p_val <- t.test(base,contrast)
      p_val <- p_val$p.value
      names(p_val) <- i
      p_value <- c(p_value,p_val)
    }
    return(p_value)
  }
  
  # genees of interest
  all_pvalues <- lapply(goi, pval_calculator)
  names(all_pvalues) <- goi
  # convert to dataframe
  all_pvalues <- data.frame(t(data.frame(all_pvalues)))
  # drop the fc_base column
  all_pvalues <- all_pvalues[ , -which(names(all_pvalues) %in% fc_base), drop=F]
  # rearrange the row names to match fc
  all_pvalues <- all_pvalues[row.names(fc),,drop=F]
  # Keep only significant values
  all_pvalues[all_pvalues > 0.05] <- NA
  all_pvalues[all_pvalues < 0.05] <- 1
  
  # multiply the fc and p value to get only the significant values
  single_mat <- fc * all_pvalues
  single_mat[single_mat > 1] <- 'Sig Up'
  single_mat[single_mat < 1] <- 'Sig Down'
  # Convert NA to 0
  single_mat <- single_mat %>% mutate_all(funs(ifelse(is.na(.), 'Non-Sig', .)))
  row.names(single_mat) <- row.names(all_pvalues)
  
  #single_mat = single_mat[!single_mat$Duvelisib %in% c('Non-Sig'),,drop=F]
  
  # Set col
  #h1_col = colorRamp2(c(-0.9, 0, 0.9), c("#0000FF","#FFFFFF", "#FF0000"))
  h1_col = structure(c('red','blue','#F8F8F8'), names = c('Sig Up','Sig Down','Non-Sig')) # black, red, green, blue
  
  # Heatmap of the final dataframe
  Heatmap(as.matrix(single_mat), cluster_columns = F, cluster_rows = T, col = h1_col,
          show_row_names = TRUE, na_col = "white",name='Sig FC',
          rect_gp = gpar(col = "black", lwd = 1))
  

  
}
sig_heatmap <- function(sig, data, meta, intgroup, fc_base=NULL, order=NULL, method='heatmap'){
  require(ComplexHeatmap)
  require(circlize)
  require(RColorBrewer)
  
  # clean up the signature and find the overlap between data
  use_sig <- sig[sig %in% row.names(data)]
  
  # Subset the data
  data_subset <- data[row.names(data) %in% use_sig, ]
  
  # Sacle the data for heatmap
  data_subset_scaled = t(scale(t(data_subset)))
  # color
  h1_col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B6AAF",  '#55B0AE', "#F8F6B8","#F5A15B","#B11E4B"))
  #Set Column annotation 
  col_ann <- data.frame(as.character(meta[,intgroup]))
  colnames(col_ann) <- c('group')
  col_color <- brewer.pal(length(unique(col_ann[,1])),"Dark2") #BrBG
  if (length(col_color) > length(unique(col_ann[,1]))){col_color = col_color[1:length(unique(col_ann[,1]))]}
  names(col_color) <- unique(col_ann[,1])
  col_color <- list(group = col_color)
  col_Ann <- HeatmapAnnotation(df=col_ann, col=col_color, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
  # Heatmap
  h2 <- goi_fc (goi=use_sig, data=data, meta=meta, intgroup=intgroup, fc_base=fc_base, order=order)
  h1 <- Heatmap(as.matrix(data_subset_scaled), col = h1_col, cluster_columns = F, cluster_rows = T, show_row_names = T, name='z score',top_annotation=col_Ann,row_names_gp = gpar(fontsize = 8))
  draw(h1+h2)
  
}


goi = c('MAPK1')
goi_stacked (ndata,meta,goi,gene_size=12,intgroup='zone',groups=NULL,order="order",ttest.ref.group='control')

goi_fc (goi, ndata, meta, intgroup='group', fc_base='NHP_c', order='order')


# normadata
ndata <- read.csv(file= "ARSeq/normalized_data.csv", header = T, check.names=F, row.names = 1)

goi <- c('DDX3X','MX2','MX1')
groups <- list(A = c("uninfected"), B= c("covid_d2"))

average_goi_plot (g_data,g_meta,goi,intgroup='region',groups)
goi_plot (g_data,g_meta,goi='MX1',intgroup='region',wes_palette='Moonrise2')

goi_collapsed(data,meta,goi='C3',intgroup='group')
goi_collapsed(g_ndata,g_meta,goi='IL10',intgroup='group')

m = goi[1:5]

#goi_plot (data=d4,meta=d4_meta,goi='CTSB',intgroup='post-infecton',wes_palette='Moonrise2')


# line plot


