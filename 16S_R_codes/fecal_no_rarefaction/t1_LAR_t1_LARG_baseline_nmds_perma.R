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

comp1<-"t.1 - 4 Weeks LAR"
comp1a<-"t.1 - 4 Weeks LAR"
comp1filename<-"t1_LAR"
comp1title<-"Post-Earth Return LAR"

comp2<-"t.1 - 4 Weeks LAR_G"
comp2a<-"t.1 - 4 Weeks LAR_G"
comp2filename<-"t1_LAR_G"
comp2title<-"Post-Earth Return LAR_G"


##############   ASV   ###########
fld="."
level="ASV"
ps=physeq
##################################
script=paste0("subset_samples(ps,( TT == \"baseline_pre\" | ",key1," == comp1 | ",key2," == comp2))")
pss=eval(parse(text=script));
pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
#pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
#pssfpr<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
pssfpr=pssf
##################################
pssfpr.ord=phyloseq::ordinate(pssfpr,"NMDS","bray")
treatment_order=c(comp1a,comp2a,"baseline")
p_pssfpr=phyloseq::plot_ordination(pssfpr,pssfpr.ord,type="samples",color="Compare") + labs(color = "Treatment") + geom_point(size=2) + stat_ellipse()  
p_pssfpr$data$Compare<-as.character(p_pssfpr$data$Compare)
p_pssfpr$data$Compare<-factor(p_pssfpr$data$Compare,levels=treatment_order)
pdf(paste0("NMDS_Bray_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".pdf"),width=6,height=6)
color1=colors[[levels(p_pssfpr$data$Compare)[1]]]
color2=colors[[levels(p_pssfpr$data$Compare)[2]]]
color3=colors[["baseline"]]
ppp=p_pssfpr + scale_color_manual(values=c(color1,color2,color3),labels=c(comp1title,comp2title,"Baseline")) + ggtitle(paste0("NMDS - ",level))
print(ppp)
dev.off()
svg(paste0("NMDS_Bray_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".svg"),width=6,height=6)
print(ppp)
dev.off()
bray_dist=phyloseq::distance(pssfpr, method="bray",type="samples")
dataframe=data.frame(sample_data(pssfpr))
pnova=adonis(bray_dist ~ Compare, data = dataframe)
write_xlsx(pnova$aov.tab,paste0("Permanova_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".xlsx"))
##################################



############## Species ###########
fld="."
level="Species"
ps=tax_glom(physeq,taxrank="Species")
##################################
script=paste0("subset_samples(ps,( TT == \"baseline_pre\" |",key1," == comp1 | ",key2," == comp2))")
pss=eval(parse(text=script));
pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
pssfpr<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)

##################################
pssfpr.ord=phyloseq::ordinate(pssfpr,"NMDS","bray")
treatment_order=c(comp1a,comp2a,"baseline")
p_pssfpr=phyloseq::plot_ordination(pssfpr,pssfpr.ord,type="samples",color="Compare") + labs(color = "Treatment") + geom_point(size=2) + stat_ellipse()  
p_pssfpr$data$Compare<-as.character(p_pssfpr$data$Compare)
p_pssfpr$data$Compare<-factor(p_pssfpr$data$Compare,levels=treatment_order)
pdf(paste0("NMDS_Bray_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".pdf"),width=6,height=6)
color1=colors[[levels(p_pssfpr$data$Compare)[1]]]
color2=colors[[levels(p_pssfpr$data$Compare)[2]]]
color3=colors[["baseline"]]
ppp=p_pssfpr + scale_color_manual(values=c(color1,color2,color3),labels=c(comp1title,comp2title,"Baseline")) + ggtitle(paste0("NMDS - ",level))
print(ppp)
dev.off()
svg(paste0("NMDS_Bray_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".svg"),width=6,height=6)
print(ppp)
dev.off()
bray_dist=phyloseq::distance(pssfpr, method="bray",type="samples")
dataframe=data.frame(sample_data(pssfpr))
pnova=adonis(bray_dist ~ Compare, data = dataframe)
write_xlsx(pnova$aov.tab,paste0("Permanova_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".xlsx"))
##################################

##############  Genus  ###########
fld="."
level="Genus"
ps=tax_glom(physeq,taxrank="Genus")
##################################
script=paste0("subset_samples(ps,( TT == \"baseline_pre\" |",key1," == comp1 | ",key2," == comp2))")
pss=eval(parse(text=script));
pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
pssfpr<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
##################################
pssfpr.ord=phyloseq::ordinate(pssfpr,"NMDS","bray")
treatment_order=c(comp1a,comp2a,"baseline")
p_pssfpr=phyloseq::plot_ordination(pssfpr,pssfpr.ord,type="samples",color="Compare") + labs(color = "Treatment") + geom_point(size=2) + stat_ellipse()  
p_pssfpr$data$Compare<-as.character(p_pssfpr$data$Compare)
p_pssfpr$data$Compare<-factor(p_pssfpr$data$Compare,levels=treatment_order)
pdf(paste0("NMDS_Bray_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".pdf"),width=6,height=6)
color1=colors[[levels(p_pssfpr$data$Compare)[1]]]
color2=colors[[levels(p_pssfpr$data$Compare)[2]]]
color3=colors[["baseline"]]
ppp=p_pssfpr + scale_color_manual(values=c(color1,color2,color3),labels=c(comp1title,comp2title,"Baseline")) + ggtitle(paste0("NMDS - ",level))
print(ppp)
dev.off()
svg(paste0("NMDS_Bray_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".svg"),width=6,height=6)
print(ppp)
dev.off()
bray_dist=phyloseq::distance(pssfpr, method="bray",type="samples")
dataframe=data.frame(sample_data(pssfpr))
pnova=adonis(bray_dist ~ Compare, data = dataframe)
write_xlsx(pnova$aov.tab,paste0("Permanova_",comp1filename,".vs.",comp2filename,".vs.baseline_",level,".xlsx"))
##################################
