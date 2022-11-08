##################################
#    Load required libraries     #
##################################
library(phyloseq)
library(vegan) # adonis
library(writexl)
library(ggplot2)
library(iNEXT)
library(readr)
library(dplyr) #select
library(metacoder)#filter_taxa function in conflict with phyloseq
library(writexl)                          
library(tidyr) #unite select
library(plyr)
##################################
colors<-list();
colors[["ISS"]]<- "red"
colors[["ISS_G"]]<- "blue"
colors[["final LAR"]]<- "#FF7F00"
colors[["final LAR_G"]]<- "#666600"
colors[["t.2 - 8 Weeks LAR"]]<- "#FF4000"
colors[["t.1 - 4 Weeks LAR"]]<- "#FF4000"
colors[["t.0 - Pre Launch LAR"]]<- "#FF4000"
colors[["t.2 - 8 Weeks LAR_G"]]<- "#006600"
colors[["t.1 - 4 Weeks LAR_G"]]<- "#006600"
colors[["t.0 - Pre Launch LAR_G"]]<- "#006600"
colors[["baseline"]]<- "#444444"
load("physeq.Rdata")
load("meta.Rdata")

ps=physeq
########### Comparison ###########
obj<-6
key1<-"Compare"
key2<-"Compare"
key3<-"Compare"
key4<-"Compare"

comp1<-"ISS"
comp1a<-"ISS"
comp1filename<-"ISS"
comp1title<-"ISS"

comp2<-"ISS_G"
comp2a<-"ISS_G"
comp2filename<-"ISS_G"
comp2title<-"ISS_G"

comp3<-"final LAR"
comp3a<-"final LAR"
comp3filename<-"final_LAR"
comp3title<-"LAR"

comp4<-"final LAR_G"
comp4a<-"final LAR_G"
comp4filename<-"final_LAR_G"
comp4title<-"LAR_G"

##################################
#     Top 20 Genera Box Plot     #
##################################
#rare=56000
#library(phyloseq)
#library(ggplot2)
#
#load("physeq.Rdata")
#load("meta.Rdata")
#ps=physeq
##################################
mycols=c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black",
"gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
"gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1",
"steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
"darkorange4", "brown")

##################################
#        FB Ratio Plot           #
##################################
#rare=56000
#library(phyloseq)
#library(writexl)
#library(ggplot2)
#library(plyr)
##################################
#load("../20200408_ISS_LAR/physeq.Rdata")
#load("../20200408_ISS_LAR/meta.Rdata")
#ps=physeq
##################################
make_ps <- function(){
	script=paste0("subset_samples(ps,( TT == \"baseline_final\" |",key1," == comp1 | ",key2," == comp2 |",key3," == comp3 |",key4," == comp4))")
	pss=eval(parse(text=script));
	pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
	#pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
	#pssfpr<<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
	pssfpr<<-pssf
}
##################################
make_fbratio<-function(){
	pssfpr_phylum=tax_glom(pssfpr,taxrank="Phylum")
	phylum_tax=tax_table(pssfpr_phylum)[,"Phylum"]
	phylum_otu=otu_table(pssfpr_phylum)
	rownames(phylum_otu)=phylum_tax[,1]
	phylum_otu_tax_per=phylum_otu/rep(colSums(phylum_otu),each=nrow(phylum_otu))
	t_phylum_otu_tax_per=t(phylum_otu_tax_per)
	fb=NULL
	for (i in 1:length(t_phylum_otu_tax_per[,1])){fb=rbind(fb,t_phylum_otu_tax_per[i,])}
	fb=cbind("SampleID"=rownames(t_phylum_otu_tax_per),fb)
	fbdf=as.data.frame(fb)
	meta.fb=meta[meta$X %in% fbdf$SampleID,]
	names(meta.fb)[1]="SampleID"
	fbdfm<<-merge(fbdf,meta.fb,by="SampleID")
	fbratio<<-data.frame(
		'fbratio'  =log2(as.numeric(as.character(fbdfm[,which(colnames(fbdfm)=="Firmicutes")]))/as.numeric(as.character(fbdfm[,which(colnames(fbdfm)=="Bacteroidetes")]))),
		Treatment=fbdfm[,which(colnames(fbdfm)=="Compare")]
	)
	#ttt=t.test(fbratio ~ Treatment, data=fbratio)
	# www=wilcox.test(fbratio ~ Treatment, data=fbratio)
	fbratio$Treatment=factor(fbratio$Treatment) #reset factors/levels
	color1=colors[[levels(fbratio$Treatment)[1]]]
	color2=colors[[levels(fbratio$Treatment)[2]]]	
	color3=colors[[levels(fbratio$Treatment)[3]]]
	color4=colors[[levels(fbratio$Treatment)[4]]]
	color5=colors[[levels(fbratio$Treatment)[5]]]
	levels(fbratio$Treatment)[which(levels(fbratio$Treatment)==comp1a)]<-comp1title
	levels(fbratio$Treatment)[which(levels(fbratio$Treatment)==comp2a)]<-comp2title
	pdf(paste0("FB_ratios_boxplot_",comp1filename,"_vs_",comp2filename,"_vs_",comp3filename,"_vs_",comp4filename,"_vs_baseline.pdf"),width=5,height=5)
	ppp=ggplot(fbratio,aes(x=Treatment,y=fbratio,color=Treatment))
	ppp=ppp+stat_boxplot(geom='errorbar',width=0.5)+geom_boxplot()
	ppp=ppp+theme(legend.position="none")+geom_jitter(width=0.1)
	#ppp=ppp+ggtitle("Firmicutes/Bacteroidetes Ratios",subtitle=paste0("p values: ",round(ttt$p.value,4)," (t-test); ",round(www$p.value,4)," (Wilcoxon)"))
	ppp=ppp+labs(y="Log2 F/B Ratios")+theme(axis.text=element_text(size=10))
	droplevels(fbratio$Treatment)
	ppp=ppp+scale_fill_manual(values=c(color1, color2,color3,color4,color5))+scale_color_manual(values=c(color1, color2,color3,color4,color5))
	print(ppp)
	dev.off()
	svg(paste0("FB_ratios_boxplot_",comp1filename,"_vs_",comp2filename,"_vs_",comp3filename,"_vs_",comp4filename,"_vs_baseline.svg"),width=5,height=5)	
	print(ppp)
	dev.off()

	filename=paste0("FB_ratios_significance_",comp1filename,"_vs_",comp2filename,"_vs_",comp3filename,"_vs_",comp4filename,"_vs_baseline.txt")
	www=pairwise.wilcox.test(fbratio$fbratio,fbratio$Treatment,p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=F)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(fbratio$fbratio,fbratio$Treatment,p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(fbratio$fbratio,fbratio$Treatment,p.adj='holm')
	write("\n######## p values adjustment method: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(fbratio$fbratio,fbratio$Treatment,p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(fbratio$fbratio,fbratio$Treatment,p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	
	write("\n",file=filename,append=T)
	kruskal=kruskal.test(fbratio ~ Treatment, data=fbratio)
	capture.output(kruskal,file = filename,append=T)
	write("\n######## Wilcoxon Test ########",file=filename,append=T)


	# t.df=data.frame(
	# "t-test"=c(
	# "t value","fegree freedom","p value",
	# "95% confidence interval (lower)",
	# "95% confidence interval (higher)",
	# names(ttt$estimate[1]),
 	# names(ttt$estimate[2]),
 	# "stderr",
 	# "alternative hypothesis",
 	# "method",
 	# "data"
	# ),
 	# "value"=c(
 	# ttt$statistic,
 	# ttt$parameter,
 	# ttt$p.value,
 	# ttt$conf.int[1],
 	# ttt$conf.int[2],
 	# ttt$estimate[1],
 	# ttt$estimate[2],
 	# ttt$stderr,
 	# "true difference in means is not equal to 0",
 	# ttt$method,
 	# ttt$data
	# ))
	# w.df=data.frame(
	# "wilcox test"=
	# c(
	# "wilcox statistic","p value",
 	# "alternative hypothesis",
 	# "method",
 	# "data"
	# ),
 	# "value"
 	# =c(
 	# www$statistic,
 	# www$p.value,
 	# " true location shift is not equal to 0",
 	# www$method,
 	# www$data
	# ))	
	# write_xlsx(list('t.test'=t.df,'wilcox'=w.df),paste0("FB_Ratios_Statistics_",comp1filename,"_vs_",comp2filename,".xlsx"))
}
##################################
make_fbstat<-function(){
	phylum_per<-list(
	tmp01=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp1a),which(colnames(fbdfm)=="Actinobacteria" )])),
	tmp02=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp2a),which(colnames(fbdfm)=="Actinobacteria" )])),
	tmp03=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp1a),which(colnames(fbdfm)=="Bacteroidetes"  )])),
	tmp04=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp2a),which(colnames(fbdfm)=="Bacteroidetes"  )])),
	tmp05=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp1a),which(colnames(fbdfm)=="Deferribacteres")])),
	tmp06=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp2a),which(colnames(fbdfm)=="Deferribacteres")])),
	tmp07=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp1a),which(colnames(fbdfm)=="Firmicutes"     )])),
	tmp08=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp2a),which(colnames(fbdfm)=="Firmicutes"     )])),
	tmp09=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp1a),which(colnames(fbdfm)=="Proteobacteria" )])),
	tmp10=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp2a),which(colnames(fbdfm)=="Proteobacteria" )])),
	tmp11=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp1a),which(colnames(fbdfm)=="Tenericutes"    )])),
	tmp12=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp2a),which(colnames(fbdfm)=="Tenericutes"    )])),
	tmp13=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp1a),which(colnames(fbdfm)=="Verrucomicrobia")])),
	tmp14=as.numeric(as.character(fbdfm[which(fbdfm$Compare==comp2a),which(colnames(fbdfm)=="Verrucomicrobia")]))
	)
	names(phylum_per)=
	c(
	paste0("Actinobacteria" ,"_",comp1filename),
	paste0("Actinobacteria" ,"_",comp2filename),
	paste0("Bacteroidetes"  ,"_",comp1filename),
	paste0("Bacteroidetes"  ,"_",comp2filename),
	paste0("Deferribacteres","_",comp1filename),
	paste0("Deferribacteres","_",comp2filename),
	paste0("Firmicutes"     ,"_",comp1filename),
	paste0("Firmicutes"     ,"_",comp2filename),
	paste0("Proteobacteria" ,"_",comp1filename),
	paste0("Proteobacteria" ,"_",comp2filename),
	paste0("Tenericutes"    ,"_",comp1filename),
	paste0("Tenericutes"    ,"_",comp2filename),
	paste0("Verrucomicrobia","_",comp1filename),
	paste0("Verrucomicrobia","_",comp2filename)
	)
	lll=ldply(phylum_per,rbind)
	colnames(lll)[1]="Phylum/Sample"
	write_xlsx(lll,paste0("Phylum_abundance_",comp1filename,"_vs_",comp2filename,"_vs_",comp3filename,"_vs_",comp4filename,"_vs_baseline.xlsx"))
}
##################################
make_ps()
make_fbratio()
make_fbstat()
