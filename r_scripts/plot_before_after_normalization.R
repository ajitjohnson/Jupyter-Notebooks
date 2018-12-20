# Function to vizulaize expression profile before and after normalization
# To use this function, it requires two data frames

plotba <- function(x,y) {
	require(ggplot2)
	require(reshape2)
	require(cowplot)
	# Before normalization 
	p1 <- ggplot(data = melt(log2(x)), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))+
		  theme_classic()+
		  theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))
	# After normalization
	p2 <- ggplot(data = melt(y), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))+
		  theme_classic()+
		  theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))
	return(plot_grid(p1, p2, labels = "AUTO"))
	}