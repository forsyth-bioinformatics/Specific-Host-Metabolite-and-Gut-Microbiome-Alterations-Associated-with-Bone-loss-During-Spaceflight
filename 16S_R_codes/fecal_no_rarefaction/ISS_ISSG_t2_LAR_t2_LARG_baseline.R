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
rare=40000
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

comp3<-"t.2 - 8 Weeks LAR"
comp3a<-"t.2 - 8 Weeks LAR"
comp3filename<-"t2_LAR"
comp3title<-"End-Point LAR"

comp4<-"t.2 - 8 Weeks LAR_G"
comp4a<-"t.2 - 8 Weeks LAR_G"
comp4filename<-"t2_LAR_G"
comp4title<-"End-Point LAR_G"

##############   ASV   ###########
fld="."
level="ASV"
ps=physeq
##################################
script=paste0("subset_samples(ps,( TT == \"baseline_final\" | ",key1," == comp1 | ",key2," == comp2 | ",key3," == comp3 | ",key4," == comp4))")
pss=eval(parse(text=script));
pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
#pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
#pssfpr<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
pssfpr=pssf
################################## alpha box plot
p_rich=plot_richness(pssfpr,x=key1,measures=c("Observed","Shannon","Simpson"), title="Alpha Diversity",color=key1) + geom_boxplot() + theme(plot.title = element_text(hjust = 0.5)) + xlab("Group")
svg(paste0(fld,"/alpha_diversity_box_plots_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".svg"), width=8, height=4)
p_rich
dev.off()
pdf(paste0(fld,"/alpha_diversity_box_plots_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".pdf"), width=8, height=4)
p_rich
dev.off()
##################################
richness <- estimate_richness(pssfpr,measures=c("Observed","Shannon","Simpson"))
write.table(richness,file=paste0("Alpha_Diversity_Indices_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt"),quote = F,sep="\t",col.names=NA)

idx="Observed"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," ",level," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


idx="Shannon"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


idx="Simpson"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


##################################




############## Species ###########
fld="."
level="Species"
ps=tax_glom(physeq,taxrank="Species")
##################################
script=paste0("subset_samples(ps,( TT == \"baseline_final\" |",key1," == comp1 | ",key2," == comp2 | ",key3," == comp3 | ",key4," == comp4))")
pss=eval(parse(text=script));
pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
pssfpr<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
################################## alpha box plot
p_rich=plot_richness(pssfpr,x=key1,measures=c("Observed","Shannon","Simpson"), title="Alpha Diversity",color=key1) + geom_boxplot() + theme(plot.title = element_text(hjust = 0.5)) + xlab("Group")
svg(paste0(fld,"/alpha_diversity_box_plots_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".svg"), width=8, height=4)
p_rich
dev.off()
pdf(paste0(fld,"/alpha_diversity_box_plots_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".pdf"), width=8, height=4)
p_rich
dev.off()
##################################
richness <- estimate_richness(pssfpr,measures=c("Observed","Shannon","Simpson"))
write.table(richness,file=paste0("Alpha_Diversity_Indices_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt"),quote = F,sep="\t",col.names=NA)

idx="Observed"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," ",level," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


idx="Shannon"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


idx="Simpson"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


##################################
##############  Genus  ###########
fld="."
level="Genus"
ps=tax_glom(physeq,taxrank="Genus")
##################################
script=paste0("subset_samples(ps,( TT == \"baseline_final\" |",key1," == comp1 | ",key2," == comp2 | ",key3," == comp3 | ",key4," == comp4))")
pss=eval(parse(text=script));
pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
pssfpr<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
################################## alpha box plot
p_rich=plot_richness(pssfpr,x=key1,measures=c("Observed","Shannon","Simpson"), title="Alpha Diversity",color=key1) + geom_boxplot() + theme(plot.title = element_text(hjust = 0.5)) + xlab("Group")
svg(paste0(fld,"/alpha_diversity_box_plots_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".svg"), width=8, height=4)
p_rich
dev.off()
pdf(paste0(fld,"/alpha_diversity_box_plots_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".pdf"), width=8, height=4)
p_rich
dev.off()
##################################
richness <- estimate_richness(pssfpr,measures=c("Observed","Shannon","Simpson"))
write.table(richness,file=paste0("Alpha_Diversity_Indices_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt"),quote = F,sep="\t",col.names=NA)

idx="Observed"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," ",level," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


idx="Shannon"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


idx="Simpson"
anova = aov(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
tukeyHSD=TukeyHSD(anova)
kruskal=kruskal.test(richness[[idx]] ~ sample_data(pssfpr)[[key1]])
wilcox=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]], p.adj = "bonf")

filename=paste0("Alpha_Group_Significance_by_",idx,"_",comp1filename,".vs.",comp2filename,".vs.",comp3filename,".vs.",comp4filename,".vs.","baseline_",level,".txt")
write(paste0(idx," Diversity Group Significance"),file=filename)
write("\n######## ANOVA Test ########",file=filename,append=T)
capture.output(anova,file = filename,append=T)
write("\n######## tukeyHSD Test ########",file=filename,append=T)
capture.output(tukeyHSD,file = filename,append=T)
write("\n######## Kruskal-Wallis Rank Sum Test ########",file=filename,append=T)
capture.output(kruskal,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BH')
	write("\n######## p values adjustment method: Benjamini & Hochberg = FDR ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='bonferroni')
	write("\n######## p values adjustment method: Bonferroni ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)		
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='holm')
	write("\n######## p values adjustment methfbratio$fbratio,fbratio$Treatment,od: Holm ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='hochberg')
	write("\n######## p values adjustment method: Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)
	www=pairwise.wilcox.test(richness[[idx]] , sample_data(pssfpr)[[key1]],p.adj='BY')
	write("\n######## p values adjustment method: Benjamini & Hochberg ########",file=filename,append=T)
	capture.output(www,file = filename,append=T)	


##################################
