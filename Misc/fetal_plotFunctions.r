
library(ggplot2)

age_scatter <- function(res, pheno, betas, i, age.column){

	# extract info
	probe <- rownames(res)[i]													# the i_th probe in res
	betas.plot <- (betas[probe,])*100											# corresponding DNAm values for probe. Multiplied by 100 for %
	age <- pheno[,age.column] 														# age column of sample sheet
	data <- data.frame(Betas=betas.plot, Age=age)								# combine DNAm (y-axis) and age (x-axis) into df

	# creating the title using probe ID, gene name (if applicable) and chr
	if(res[i,'Gene']!=''){														# if gene is present...
		ttl <- paste0(probe, ' - ', res[i,'Gene'], ' - CHR ',res[i,'CHR'])			# combine probe ID, gene name and chromosome
	}else{																		# if there's no gene name...
		ttl <- paste0(probe, ' - CHR ', res[i,'CHR'])								# just probe ID and chr
	}
	ttl <- paste0(ttl, "\np = ", signif(res[i,'P.Age'], digits=3))				# append the p-value for that probe to the title as a new line

	# plot
	ggplot(data, aes(x=Age, y=Betas)) + 
		geom_point(size=3)+
		ggtitle(ttl)+
		xlab("Age (PCW)")+
		ylab("DNA methylation (%)")+
		ylim(0,100)+
		stat_smooth(method = lm, formula = y ~ x,se=T) +
		theme_minimal()+
		theme(axis.text=element_text(size=19), axis.title=element_text(size=20))+
		theme(plot.title = element_text(size=21))
}


group_boxPlot <- function(res, pheno, betas, i, p.column, group.column, group.plot.name=NULL, xlab="", ylab="DNA methylation (%)", colours){

	# create plotting df
	probe <- rownames(res)[i]
	betas.plot <- (betas[probe,])*100
	p.value <- signif(res[i,paste(p.column)], digits=3)
	df <- data.frame(Betas=betas.plot, Group=pheno[,paste0(group.column)])

	if(is.null(group.plot.name)){
		group.plot.name <- group.column
	}
	
	# create plot title
	if(res[i,'Gene']!=''){ # if probe is annotated to a gene
			ttl <- paste0(probe, ' - ', res[i,'Gene'], ' - CHR ',res[i,'CHR'])
	}else{
		ttl <- paste0(probe, ' - CHR ', res[i,'CHR'])
	}
	ttl <- paste0(ttl,"\np = ",p.value)	# append p-value to title


	ggplot(df, aes(x=Group, y=Betas, fill=Group)) + 
		geom_boxplot(outlier.shape=NA)+
		geom_jitter(height = 0, width = 0.1)+
		ggtitle(ttl)+
		xlab(xlab)+
		ylab(ylab)+
		ylim(0,100)+
		theme_minimal()+
		#theme(legend.position = "none")+
		theme(axis.text=element_text(size=19), axis.title=element_text(size=20))+
		theme(plot.title = element_text(size=21))+
		scale_fill_manual(values=colours, name=group.plot.name)
}


ageByGroup_scatter <- function(res, pheno, betas, i, age.column, p.column, group.column, group.plot.name=NULL, xlab="Age (pcw)", ylab="DNA methylation (%)", colours, Xlim=NULL, Ylim=NULL){
	probe <- rownames(res)[i]
	betas.plot <- (as.numeric(betas[probe,]))*100
	p.value <- signif(res[i,paste(p.column)], digits=3)
	age <- pheno[,age.column]
	df <- data.frame(Betas=betas.plot, Age=age, Group=pheno[,paste0(group.column)])

	if(is.null(Xlim)){
		Xlim <- c(min(df$Age), max(df$Age))
	}
	if(is.null(Ylim)){
		Ylim <- c(0,100)
	}
	if(is.null(group.plot.name)){
		group.plot.name <- group.column
	}
	
	# creating the title using probe ID, gene name (if applicable) and chr
	if(res[i,'Gene']!=''){														# if gene is present...
		ttl <- paste0(probe, ' - ', res[i,'Gene'], ' - CHR ',res[i,'CHR'])			# combine probe ID, gene name and chromosome
	}else{																		# if there's no gene name...
		ttl <- paste0(probe, ' - CHR ', res[i,'CHR'])								# just probe ID and chr
	}
	ttl <- paste0(ttl,"\np = ",p.value)

	ggplot(df, aes(x=Age, y=Betas, col=Group)) + 
		geom_point(size=3)+
		ggtitle(ttl)+
		xlab(xlab)+
		ylab(ylab)+
		scale_x_continuous(limits=c(Xlim),expand=c(0.1,0.1))+ # using scale continuous as usual lim() arg expands beyond desired age range
		scale_y_continuous(limits=c(Ylim),expand=c(0.1,0.1),breaks=c(0,50,100))+
		stat_smooth(method = lm, formula = y ~ x,se=T)+
		theme_minimal()+
		theme(axis.text=element_text(size=19), axis.title=element_text(size=20))+
		theme(plot.title = element_text(size=21))+
		scale_color_manual(values=colours, name=group.plot.name)
}