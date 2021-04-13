# Collection of plotting scripts


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
    labs(y = "Log2")
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
goi_collapsed_jitter <- function(data,meta,goi,intgroup, groups=NULL){
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

  # use ggpubr
  ggbarplot(d, x = intgroup, y = "value", add = c("mean_se", "jitter"),fill = intgroup,width=1)+
    theme(plot.title = element_text(hjust = 0.5,face="bold", color="#000000", size=14),
          axis.title.x = element_blank(),
          axis.text.x = element_text(face="bold", color="#000000", size=14,angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(face="bold", color="#000000", size=14),
          axis.title.y = element_text(face="bold", color="#000000", size=14))+
    ylab("expression (log2)") + ggtitle(goi)

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

  # Convert nan to 1
  is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
  single_mat[is.nan(single_mat)] <- NA

  is.nan.data.frame <- function(y)
  do.call(cbind, lapply(y, is.infinite))
  single_mat[is.nan(single_mat)] <- NA

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
          rect_gp = gpar(col = "black", lwd = 1),width = unit(2, "cm"),
          column_names_gp = gpar(fontsize = 10))



}
sig_heatmap <- function(sig, data, meta, cluster_rows=T, intgroup, fc_base=NULL, order=NULL, groups=NULL, method='heatmap'){
  require(ComplexHeatmap)
  require(circlize)
  require(RColorBrewer)

  # clean up the signature and find the overlap between data
  use_sig <- sig[sig %in% row.names(data)]

  # Subset the data
  data_subset <- data[row.names(data) %in% use_sig, ]

  # Identify the samples of interest if only subgroups are required
  if (!is.null(groups)){
    soi = row.names(meta[meta[,intgroup] %in% groups,])
    data_subset <- data_subset[,colnames(data_subset) %in% soi]
    meta_subset <- meta[row.names(meta) %in% soi,]
  }

  # Sacle the data for heatmap
  data_subset_scaled = t(scale(t(data_subset)))
  # color
  h1_col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B6AAF",  '#55B0AE', "#F8F6B8","#F5A15B","#B11E4B"))
  #Set Column annotation
  col_ann <- data.frame(as.character(meta_subset[,intgroup]))
  colnames(col_ann) <- c('group')
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols = gg_color_hue(length(unique(col_ann[,1])))
  names(cols) <- unique(col_ann[,1])
  col_color <- list(group = cols)

  #col_color <- brewer.pal(length(unique(col_ann[,1])),"Dark2") #BrBG
  #if (length(col_color) > length(unique(col_ann[,1]))){col_color = col_color[1:length(unique(col_ann[,1]))]}
  #names(col_color) <- unique(col_ann[,1])
  #col_color <- list(group = col_color)
  col_Ann <- HeatmapAnnotation(df=col_ann, col=col_color, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
  # Heatmap
  h2 <- goi_fc (goi=use_sig, data=data, meta=meta_subset, intgroup=intgroup, fc_base=fc_base, order=order)
  h1 <- Heatmap(as.matrix(data_subset_scaled), heatmap_width = unit(0.5, "npc"), col = h1_col, show_column_names=F, cluster_columns = F, cluster_rows = cluster_rows, show_row_names = T, name='z score',top_annotation=col_Ann,row_names_gp = gpar(fontsize = 8))
  draw(h1+h2)

}
heatmap <- function(sig, data, meta, intgroup, order=NULL, method='heatmap'){
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
  #col_color <- col_color[-1]
  col_Ann <- HeatmapAnnotation(df=col_ann, col=col_color, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
  # Heatmap
  h1 <- Heatmap(as.matrix(data_subset_scaled), heatmap_width = unit(0.5, "npc"), col = h1_col, show_column_names=F, cluster_columns = F, cluster_rows = T, show_row_names = T, name='z score',top_annotation=col_Ann,row_names_gp = gpar(fontsize = 10))
  draw(h1)

}
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
# correlation between two genes
two_gene_corr <- function(data,meta,goi1,goi2,intgroup,groups=NULL){
  require(ggplot2)

  # gene of interest
  d <- data [row.names(data) %in% c(goi1,goi2),,drop=FALSE]
  d <- merge(t(d), meta, by = "row.names")

  # Identify the samples of interest
  columns_to_keep <- c(goi1,goi2, intgroup)
  if (is.null(groups)){
    groups = unique(d[,intgroup])
  }
  d <- d[d[,intgroup] %in% unlist(groups) , colnames(d) %in% columns_to_keep,drop=FALSE]

  # Draw the plot
  ggplot(d, aes(x=d[,goi1], y=d[,goi2])) + geom_point(aes(color=d[,intgroup])) +
    theme_classic() + geom_smooth(method=lm,fullrange=TRUE)+
    labs(x = goi1, y= goi2, color=intgroup)

}

# PCA plot
arseq.pca.plot = function(dds, intgroup="arseq.group", ntop=500, pc.a= 1, pc.b = 2, returnData=FALSE,wes_palette="GrandBudapest1"){
  print("Performing PCA analysis")
  require(DESeq2)
  require(ggrepel)
  require(ggplot2)
  require(wesanderson)
  # Normalize data
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  # calculate the variance for each gene
  rv <- rowVars(vsd@assays@data[[1]])

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(vsd@assays@data[[1]][select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(vsd@colData))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(vsd@colData[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    vsd@colData[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(pca$x[,pc.a], pca$x[,pc.b], group=group, intgroup.df, name=colnames(vsd))
  names(d)[1] <- paste("PC",pc.a,sep = "")
  names(d)[2] <- paste("PC",pc.b,sep = "")


  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc.a:pc.b]
    return(d)
  }

  ggplot(data=d, aes_string(x=paste("PC",pc.a,sep = ""), y=paste("PC",pc.b,sep = ""), color="group")) + geom_point(size=2) +
    xlab(paste0(paste("PC",pc.a,sep = ""),": ",round(percentVar[pc.a] * 100),"% variance")) +
    ylab(paste0(paste("PC",pc.b,sep = ""),": ",round(percentVar[pc.b] * 100),"% variance")) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(face="bold", color="#1a1a1b",
                                     size=14),
          axis.text.y = element_text(face="bold", color="#1a1a1b",
                                     size=14)) +
    #scale_color_manual(values=wes_palette(n= 5, name=wes_palette))+
    geom_text_repel(aes(label = .data$name),size = 5) +
    coord_fixed() + ggtitle("Principal component analysis (PCA) Plot")+
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")
}
# ssgsea analysis
ssgsea_plot <- function(data,signature,meta=NULL,intgroup=NULL, ssgsea_score=NULL, custom_color=NULL, cluster_columns=F, show_column_names=F, padding = unit(c(2, 2, 2, 60), "mm")){

  # load required libraries
  require(GSEABase)
  require(GSVA)
  require(circlize)
  require(ComplexHeatmap)
  require(RColorBrewer)

  if(is.null(ssgsea_score)){
    # Make a copy of the group type
    group <- signature[,1,drop=FALSE]
    colnames(group) <- c('cluster')
    # Drop the groups from signature
    signature <- signature[,-1]
    # transpose the signature
    signature <- data.frame(t(signature))

    # Convert signature into a named list
    named_list <- function(df){
      final_list <- c()
      for (i in 1: ncol(df)){
        tmp <- as.character(df[,i])
        tmp <- list(tmp[nchar(tmp) > 1])
        names(tmp) <- colnames(df[,i,drop=F])
        final_list <- c(final_list, tmp)
      }
      return(final_list)
    }
    x = named_list (signature)

    # Run the enrichment
    gbm_es <- gsva(as.matrix(data), x,method='ssgsea',ssgsea.norm=T)
  } else {
    signature <- signature [row.names(signature) %in% row.names(ssgsea_score),]
    group <- signature[,1,drop=FALSE]
    colnames(group) <- c('cluster')
    gbm_es <- ssgsea_score
  }



  # Viz using heatmap with metadata
  ssgsea_scaled = t(scale(t(gbm_es)))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#4B6AAF",  '#55B0AE', "#F8F6B8","#F5A15B","#B11E4B"))
  h1_col = colorRamp2(c(-1, 0, 1), c("#F8F6B8","#F5A15B","#B11E4B"))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#f2e9e4",  '#c9ada7', "#9a8c98","#4a4e69","#22223b"))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#e0e1dd",  '#778da9', "#415a77","#1b263b","#0d1b2a"))
  #h1_col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#f8f9fa",  '#dee2e6', "#ced4da","#6c757d","#343a40"))

  # converge grouping with ssgsea
  group <- group[row.names(group) %in% row.names(gbm_es), , drop = FALSE]

  if(!is.null(meta)){
    # Add column annotation (Sample group)
    col_ann <- data.frame(as.character(meta[,intgroup]))
    colnames(col_ann) <- c('Groups')
    # Create color pallete
    colourCount <- length(unique(col_ann[,1]))
    getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
    col_color <- getPalette(colourCount)
    if (!is.null(custom_color)){col_color <- custom_color}
    names(col_color) <- unique(col_ann[,1])
    col_color <- list(Groups = col_color)
    col_Ann <- HeatmapAnnotation(df=col_ann, col=col_color, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
  }

  # Heatmap
  if(is.null(meta)){
    h1 <- Heatmap(as.matrix(ssgsea_scaled), col = h1_col, border=T, cluster_columns = F,
                  rect_gp = gpar(col = "#22223b", lwd = 1), show_column_names = show_column_names,
                  row_names_gp = gpar(fontsize = 12), cluster_rows = T,
                  row_split = factor(group$cluster, levels = unique(group$cluster)),
                  show_row_names = TRUE, name='ssGSEA score')
    draw(h1, heatmap_legend_side = "left", annotation_legend_side = "left", padding = padding)

  } else {
    h1 <- Heatmap(as.matrix(ssgsea_scaled), col = h1_col, border=T, cluster_columns = cluster_columns,
                  rect_gp = gpar(col = "#22223b", lwd = 1), show_column_names = show_column_names,
                  row_names_gp = gpar(fontsize = 12), cluster_rows = T, top_annotation=col_Ann,
                  row_split = factor(group$cluster, levels = unique(group$cluster)),
                  show_row_names = TRUE, name='ssGSEA score')
    draw(h1, heatmap_legend_side = "left", annotation_legend_side = "left", padding = padding)

  }

}
# ssGSEA DE
ssgsea_score <- function(data,signature, ssgsea_save=FALSE){

  # load required libraries
  require(GSEABase)
  require(GSVA)
  require(circlize)
  require(ComplexHeatmap)
  require(RColorBrewer)

  #data = g_ndata
  #signature = custom_sig
  #meta = g_meta
  #intgroup = 'roi'
  #group_of_interest = 'MIS'
  # make a copy
  real_signature <- signature

  # Make a copy of the group type
  group <- signature[,1,drop=FALSE]
  colnames(group) <- c('cluster')
  # Drop the groups from signature
  signature <- signature[,-1]
  # transpose the signature
  signature <- data.frame(t(signature))

  # Convert signature into a named list
  named_list <- function(df){
    final_list <- c()
    for (i in 1: ncol(df)){
      tmp <- as.character(df[,i])
      tmp <- list(tmp[nchar(tmp) > 1])
      names(tmp) <- colnames(df[,i,drop=F])
      final_list <- c(final_list, tmp)
    }
    return(final_list)
  }
  x = named_list (signature)

  # Run the enrichment
  gbm_es <- gsva(as.matrix(data), x,method='ssgsea',ssgsea.norm=T)
  gbm_es <- data.frame(gbm_es)

  # Save aggregate
  if(ssgsea_save == TRUE){
    write.csv(gbm_es, file = 'ssgsea_score.csv')
  }

  # return modified signature
  return(gbm_es)

}
# ssGSEA DE finder
ssgsea_de <- function(ssgsea_score, meta,intgroup,group_of_interest,de_cutoff=0,ssgsea_aggragate=NULL, ssgsea_aggragate_save=FALSE){

  gbm_es <- ssgsea_score
  gbm_es_raw <- ssgsea_score
  gbm_es <- data.frame(t(scale(t(gbm_es))))

  if(is.null(ssgsea_aggragate)){
    # idnetify signatures that are DE's
    gbm <- data.frame(t(gbm_es))
    gbm <- merge(gbm, meta[,intgroup,drop=F], by = 'row.names')
    row.names(gbm) <- gbm[,1]
    gbm <- gbm[,-1]
    gbm[,intgroup] <- as.factor(gbm[,intgroup])
    names(gbm)[names(gbm) == intgroup] <- 'GROUP'
    # aggregate
    gbm=aggregate(.~GROUP, data=gbm, mean)
    row.names(gbm) <- gbm[,1]
    gbm <- gbm[,-1]
    # Save aggregate
    if(ssgsea_aggragate_save == TRUE){
      write.csv(gbm, file = 'ssgsea_aggregate_score.csv')
    }

  } else {
    gbm <- ssgsea_aggragate
  }

  # Find columns with mena > 0
  gbm_group <- gbm[group_of_interest,,drop=F]
  sigs <- gbm_group[which(gbm_group > de_cutoff)]

  # Subset ssgsea score
  gbm_final <- data.frame(gbm_es_raw[row.names(gbm_es_raw) %in% colnames(sigs),])


  # return modified signature
  return(gbm_final)

}

# SSGSEA maaping to binary clinical phenotypes
ssgsea_binary_phenotype_corr <- function(data,signature,meta,min.sz=5){

  # pre-requisite
  # pass the signature as a dataframe (example: directly from GSEA website)
  # the phenotype needs to be binary (0-not active, 1- active), NA is allowed. Samples with NA will be filtered

  # Sig:
  # signature <-  read.csv(file= "/Volumes/SSD/Dropbox (Partners HealthCare)/Data/covid19/signatures/msigdb.v7.2.symbols.csv", header = F, row.names=1, check.names = F)

  # usage
  # df <- ssgsea_binary_phenotype_corr(data=ndata,signature=signature,meta=clinical_meta,min.sz=5)

  # load required libraries
  require(GSEABase)
  require(GSVA)

  # Make a copy of the group type
  group <- signature[,1,drop=FALSE]
  colnames(group) <- c('cluster')
  # Drop the groups from signature
  signature <- signature[,-1]
  # transpose the signature
  signature <- data.frame(t(signature))

  # Convert signature into a named list
  named_list <- function(df){
    final_list <- c()
    for (i in 1: ncol(df)){
      tmp <- as.character(df[,i])
      tmp <- list(tmp[nchar(tmp) > 1])
      names(tmp) <- colnames(df[,i,drop=F])
      final_list <- c(final_list, tmp)
    }
    return(final_list)
  }
  x = named_list (signature)

  # Run the enrichment
  gbm_es <- gsva(as.matrix(data), x,method='ssgsea',ssgsea.norm=T, min.sz = min.sz)
  gbm_es <- data.frame(gbm_es)
  colnames(gbm_es) <- colnames(data)

  # Perform a t-test to check if the enrichment score is different between clinical groups
  pval = list()
  for (i in 1: ncol(meta)){
    print(paste ('Processing', colnames(meta[i])))
    # Find the samples that are postive or negative for the given phenotype
    pos = row.names(subset(meta, meta[,i,drop=F] == 1))
    neg = row.names(subset(meta, meta[,i,drop=F] == 0))
    # calculate p value
    pval[i] = list(apply(gbm_es, 1, function(x) t.test(x[pos],x[neg])$p.value))
  }
  # collapse the list to a dataframe
  df <- data.frame(matrix(unlist(pval), nrow=nrow(gbm_es), byrow=F),stringsAsFactors=FALSE)
  row.names(df) <- row.names(gbm_es)
  colnames(df) <- colnames(meta)

  # return
  return(df)

}

# Volcano plot
arseq.volcano.plot <- function(deg,pCutoff=0.05,FCcutoff=1,colCustom=keyvals,selectLab=NULL){
  print("Generating a volcano plot between the constrast groups")
  require(EnhancedVolcano)
  deg.volcano <- data.frame(deg)[complete.cases(data.frame(deg)),]
  EnhancedVolcano(deg.volcano,
                  lab = rownames(deg.volcano),
                  selectLab=selectLab,
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  colCustom = colCustom,
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,
                  gridlines.minor = FALSE,
                  gridlines.major = FALSE,
                  border = 'full',
                  labSize = 4.0,
                  pointSize = 2.0,
                  colAlpha = 1,
                  #legend=c('NS','3 Log2 FC','Adj P-value (<0.005)','Adj P-value (<0.005) & 3 Log2 FC'),
                  legendPosition = 'right',
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  colConnectors = 'grey30',
                  legendIconSize = 3.0)
}

# Accessories
deg <- read.csv(file= "ARSeq/Microthrombi/Differential expression/Y vs N.csv", header = T, check.names=F, row.names = 1)
deg.volcano <- data.frame(deg)[complete.cases(data.frame(deg)),]
# Run the function (https://github.com/ajitjohnson/arseq/blob/master/R/arseq.volcano.plot.R)
keyvals <- ifelse(deg.volcano$log2FoldChange < -1 & deg.volcano$padj < 0.05, '#AA735B',
                  ifelse(deg.volcano$log2FoldChange > 1  & deg.volcano$padj < 0.05, '#AA735B',
                         ifelse(deg.volcano$log2FoldChange < -1  & deg.volcano$padj > 0.05, '#92c5de',
                                ifelse(deg.volcano$log2FoldChange > 1  & deg.volcano$padj > 0.05, '#92c5de','#D3D3D3'))))
keyvals[is.na(keyvals)] <- '#D3D3D3'
names(keyvals)[keyvals == '#AA735B'] <- '1 Log2 FC & Padj < 0.05'
names(keyvals)[keyvals == '#92c5de'] <- '1 Log2 FC'
names(keyvals)[keyvals == '#D3D3D3'] <- 'NS'

# Plots
#arseq.volcano.plot(deg,selectLab=selectLab)
arseq.volcano.plot(deg)


# loop through a set of markers and save a PDF file
goi = read.csv("/Users/aj/Dropbox (Partners HealthCare)/Data/PCA/signatures/all_stained_markers.csv", header = F)
goi = goi$V1

somePDFPath = paste(getwd(), '/pickseq.pdf', sep="")
pdf(file=somePDFPath)
for (i in goi){
  try(plot(goi_collapsed_jitter (data=p_ndata,meta=p_meta,goi=i,intgroup='roi')))
}
dev.off()

somePDFPath = paste(getwd(), '/geomx.pdf', sep="")
pdf(file=somePDFPath)
for (i in goi){
  try(plot(goi_collapsed_jitter (data=g_ndata,meta=g_meta,goi=i,intgroup='roi')))
}
dev.off()

# correla


# subset ndata
subset_meta <- p_meta[p_meta$roi %in% c('IT','ET'),]
subset_ndata <- p_ndata[,colnames(p_ndata) %in% row.names(subset_meta)]

# Find correlation for a gene
library(HiClimR)
exp <- p_ndata#subset_ndata
g_corr = data.frame(fastCor(as.matrix(t(exp)), optBLAS = TRUE))
goi = "S100B"
g_goi <- g_corr[,goi,drop=F]
g_goi <- g_goi[order(g_goi[,1],decreasing = T),,drop=F]
head(g_goi, 20)
tail(g_goi, 20)
write.table(g_goi, file = "S100B_correlated_genes.csv", sep = ',')

g_goi['CD8A',]
g_goi['CD24',]

# Normal PCA
pca <- prcomp(data.frame(t(log2(g_data))), scale. = T)
autoplot(pca, data = g_meta, colour = 'zone', size = 3) + theme_classic()+ geom_text_repel(aes(label = colnames(g_data)),size = 3)
