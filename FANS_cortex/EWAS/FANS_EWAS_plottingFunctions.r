

## Plotting functions specific to FANS EWAS results ##

library(ggplot2)
library(viridis)
library(dplyr)
options(scipen=999) #remove scientific notation


ESscatter_fetalVSadult <- function(df, age.group, type='Point', bins=120){

# scatter
if(type=='Point'){
p <- {if(age.group=='Fetal'){
		ggplot(df %>% arrange(Significance), aes(x=Fetal, y=Adult, colour=Significance))+
		xlim(min(df$Fetal)-0.001, max(df$Fetal)+0.001)+
		xlab("Fetal\nEffect Size")+
		ylab("Adult\nEffect Size")
	}else{ # if age.group == 'Adult' | 'All'
		ggplot(df %>% arrange(Significance), aes(x=Adult, y=Fetal, colour=Significance))+
		ylim(min(df$Fetal)-0.001, max(df$Fetal)+0.001)+
		ylab("Fetal\nEffect Size")+
		xlab("Adult\nEffect Size")
	}}+ 
	geom_point(size=2)+
	geom_vline(xintercept=0, linetype="dashed", color = "darkgrey", linewidth=1)+
	geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth=1)+
	ggtitle(ttl)+
	theme_minimal()+
	theme(legend.position = c(0.10, 0.15), legend.text = element_text(size=15), legend.title = element_text(size=15))+ #for Cell
	theme(axis.text=element_text(size=19), axis.title=element_text(size=22))+
	theme(plot.title = element_text(size=23))+
	scale_colour_manual(values = c('Not significant'="grey",  'Fetal'=viridis(20)[11], 'Adult'=viridis(20)[6], 'Both'=viridis(20)[16])) #grey, dark purple, greeny-blue, yellow
}

# heatmap
if(type=='2D'){

p <- {if(age.group=='Fetal'){
		ggplot(df, aes(x=Fetal, y=Adult))+
		xlim(min(df$Fetal)-5, max(df$Fetal)+5)+
		xlab("Fetal\nEffect Size")+
		ylab("Adult\nEffect Size")
	}else{ # if age.group == 'Adult' | 'All'
		ggplot(df, aes(x=Adult, y=Fetal))+
		ylim(min(df$Fetal)-5, max(df$Fetal)+5)+
		ylab("Fetal\nCell-type effect size (%)")+
		xlab("Adult\nCell-type effect size (%)")
	}}+ 
		geom_vline(xintercept=0, linetype="dashed", color = "darkgrey", linewidth=1)+
		geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth=1)+
		geom_bin2d(bins=bins)+
		ggtitle(ttl)+
		theme_minimal()+
		theme(legend.position = c(0.1, 0.85), legend.text = element_text(size=15), legend.title = element_text(size=15))+
		theme(axis.text=element_text(size=19), axis.title=element_text(size=22))+
		theme(plot.title = element_text(size=23))
}

show(p)
}
	


ESscatter_AgeCellSpecific <- function(df, age.group, Xlim=NULL, Ylim=NULL){

p <- ggplot(df %>% arrange(Significance), aes(x=Neuronal, y=Non.neuronal, colour=Significance))+
		{if(!is.null(Xlim)){ xlim(Xlim) }}+
		{if(!is.null(Ylim)){ ylim(Ylim) }}+
		xlab("Neuronal\nAge Effect Size")+
		ylab("Non-neuronal\nAge Effect Size")+
		geom_point(size=2)+
		geom_vline(xintercept=0, linetype="dashed", color = "darkgrey", linewidth=1)+
		geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth=1)+
		theme_minimal()+
		theme(legend.position = c(0.15, 0.85), legend.text = element_text(size=15), legend.title = element_text(size=15))+
		theme(axis.text=element_text(size=19), axis.title=element_text(size=22), plot.title=element_text(size=23))+
		scale_colour_manual(values = c('Not significant'="grey",  'Neuronal'=viridis(20)[11], 'Non-neuronal'=viridis(20)[6], 'Both'=viridis(20)[16])) #grey, dark purple, greeny-blue, yellow

show(p)
}