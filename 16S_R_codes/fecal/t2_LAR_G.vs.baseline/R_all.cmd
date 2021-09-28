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
load("../physeq.Rdata")
load("../meta.Rdata")
ps=physeq
rare=56000
########### Comparison ###########
obj<-13
key1<-"Compare"
key2<-"TT"
comp1<-"t.2 - 8 Weeks LAR_G"
comp1a<-"t.2 - 8 Weeks LAR_G"
comp1filename<-"t2_LAR_G"
comp1title<-"End-Point LAR_G"
comp2<-"baseline_final"
comp2a<-"baseline"
comp2filename<-"baseline"
comp2title<-"Baseline"

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
level="Genus"
##################################
make_ps <- function(){
	script=paste0("subset_samples(ps,( ",key1," == comp1 | ",key2," == comp2))")
	pss=eval(parse(text=script));
	pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
	pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
	pssfpr<<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
}
##################################
make_melt <- function(cutoff=0.01){
	psbargenus=tax_glom(pssfpr,"Genus")
	psbarper=transform_sample_counts(psbargenus,function(x) x / sum (x))
	ucount=0
	ucount1=0
	cutoff1=0
	i=1
	j=3
	target=21
	while(ucount!=target){
		ps_melt<-psmelt(psbarper)
		ps_melt$Genus<-as.character(ps_melt$Genus)
		ps_melt$Genus[ps_melt$Abundance <=cutoff]<-"Other"
		ps_melt$Genus<-gsub("\\[","",ps_melt$Genus)
		ps_melt$Genus<-gsub("\\]","",ps_melt$Genus)
		ucount=length(unique(ps_melt$Genus))
		print(paste0(cutoff,"==",cutoff1,"****",ucount,"==",ucount1))
		cutoff1=cutoff
		if(ucount < target){
			if(ucount1 > target){#revese direction and go by one more decimal
				j=j+1;
			}
			script=paste0(i,"e-",j)
			incr=eval(parse(text=script));
			cutoff=cutoff-incr			
		}else if(ucount > target){
			if(ucount1 <target){#revese direction and go by one more decimal
				j=j+1;
			}
			script=paste0(i,"e-",j)
			incr=eval(parse(text=script));
			cutoff=cutoff+incr
		}
		#i=i+1
		ucount1=ucount
	}
	levels(ps_melt$Compare)[which(levels(ps_melt$Compare)==comp1a)]<-comp1title
	levels(ps_melt$Compare)[which(levels(ps_melt$Compare)==comp2a)]<-comp2title	
	ps21_melt<<-ps_melt
	topct<<-ucount-1
}
##################################
make_plot <- function(){
	psmeltplot<-ggplot(data=ps21_melt,aes(x=Sample, y=Abundance, fill=Genus))+ ylab("Relative Abundance") +facet_wrap(~Compare, scales = "free_x",nrow=2) + geom_bar(aes(), stat="identity") + theme(axis.text.x = element_blank(),axis.ticks.x=element_blank()) +  guides(fill=guide_legend(ncol=1)) + scale_fill_manual(values=mycols)
	pdf(paste0("Top_",topct,"_bar_plot_",comp1filename,"_vs_",comp2filename,"_",level,".pdf"),width=8,height=8)
	print(psmeltplot)
	dev.off()
}
##################################
make_ps()
make_melt()
make_plot()
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
	script=paste0("subset_samples(ps,( ",key1," == comp1 | ",key2," == comp2))")
	pss=eval(parse(text=script));
	pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
	pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
	pssfpr<<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
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
	ttt=t.test(fbratio ~ Treatment, data=fbratio)
	www=wilcox.test(fbratio ~ Treatment, data=fbratio)
	fbratio$Treatment=factor(fbratio$Treatment) #reset factors/levels
	color1=colors[[levels(fbratio$Treatment)[1]]]
	color2=colors[[levels(fbratio$Treatment)[2]]]	
	levels(fbratio$Treatment)[which(levels(fbratio$Treatment)==comp1a)]<-comp1title
	levels(fbratio$Treatment)[which(levels(fbratio$Treatment)==comp2a)]<-comp2title
	pdf(paste0("FB_ratios_boxplot_",comp1filename,"_vs_",comp2filename,".pdf"),width=5,height=5)
	ppp=ggplot(fbratio,aes(x=Treatment,y=fbratio,color=Treatment))
	ppp=ppp+stat_boxplot(geom='errorbar',width=0.5)+geom_boxplot()
	ppp=ppp+theme(legend.position="none")+geom_jitter(width=0.1)
	ppp=ppp+ggtitle("Firmicutes/Bacteroidetes Ratios",subtitle=paste0("p values: ",round(ttt$p.value,4)," (t-test); ",round(www$p.value,4)," (Wilcoxon)"))
	ppp=ppp+labs(y="Log2 F/B Ratios")+theme(axis.text=element_text(size=10))
	droplevels(fbratio$Treatment)
	ppp=ppp+scale_fill_manual(values=c(color1, color2))+scale_color_manual(values=c(color1, color2))
	print(ppp)
	dev.off()
	t.df=data.frame(
	"t-test"=c(
	"t value","fegree freedom","p value",
	"95% confidence interval (lower)",
	"95% confidence interval (higher)",
	names(ttt$estimate[1]),
 	names(ttt$estimate[2]),
 	"stderr",
 	"alternative hypothesis",
 	"method",
 	"data"
	),
 	"value"=c(
 	ttt$statistic,
 	ttt$parameter,
 	ttt$p.value,
 	ttt$conf.int[1],
 	ttt$conf.int[2],
 	ttt$estimate[1],
 	ttt$estimate[2],
 	ttt$stderr,
 	"true difference in means is not equal to 0",
 	ttt$method,
 	ttt$data
	))
	w.df=data.frame(
	"wilcox test"=
	c(
	"wilcox statistic","p value",
 	"alternative hypothesis",
 	"method",
 	"data"
	),
 	"value"
 	=c(
 	www$statistic,
 	www$p.value,
 	" true location shift is not equal to 0",
 	www$method,
 	www$data
	))	
	write_xlsx(list('t.test'=t.df,'wilcox'=w.df),paste0("FB_Ratios_Statistics_",comp1filename,"_vs_",comp2filename,".xlsx"))
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
	write_xlsx(lll,paste0("Phylum_abundance_",comp1filename,"_vs_",comp2filename,".xlsx"))
}
##################################
make_ps()
make_fbratio()
make_fbstat()
##############   ASV   ###########
level="ASV"
ps=physeq
##################################
#    NMDS, iNEXT, MetacodeR      #
##################################
#rare=56000
#library(phyloseq)
#library(vegan) # adonis
#library(writexl)
#library(ggplot2)
#library(iNEXT)
#library(readr)
#library(dplyr) #select
#library(metacoder)#filter_taxa function in conflict with phyloseq
#library(writexl)                          
#library(tidyr) #unite select
#library(plyr)
#colors<-list();
#colors[["ISS"]]<- "red"
#colors[["ISS_G"]]<- "blue"
#colors[["final LAR"]]<- "#00661A"
#colors[["final LAR_G"]]<- "#666600"
#colors[["t.2 - 8 Weeks LAR"]]<- "#0000B2"
#colors[["t.1 - 4 Weeks LAR"]]<- "#D93600"
#colors[["t.0 - Pre Launch LAR"]]<- "#7F00FF"
#colors[["t.2 - 8 Weeks LAR_G"]]<- "#006DD9"
#colors[["t.1 - 4 Weeks LAR_G"]]<- "#FF8000"
#colors[["t.0 - Pre Launch LAR_G"]]<- "#9673FF"
#colors[["baseline"]]<- "#444444"
##################################
#load("../physeq.Rdata")
#load("../meta.Rdata")
#load("../inextout.ASV.Rdata")
##################################
#ps=physeq
#ps=tax_glom(physeq,taxrank="Species")
#ps=tax_glom(physeq,taxrank="Genus")
#level="ASV";
#level="Species";
#level="Genus";
##################################
make_ps <- function(){
	script=paste0("subset_samples(ps,( ",key1," == comp1 | ",key2," == comp2))")
	pss=eval(parse(text=script));
	pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
	pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
	pssfpr<<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
}

make_nmds_perma <- function(){
	pssfpr.ord=phyloseq::ordinate(pssfpr,"NMDS","bray")
	treatment_order=c(comp1a,comp2a)
	p_pssfpr=phyloseq::plot_ordination(pssfpr,pssfpr.ord,type="samples",color="Compare") + labs(color = "Treatment") + geom_point(size=2) + stat_ellipse()  
	p_pssfpr$data$Compare<-as.character(p_pssfpr$data$Compare)
	p_pssfpr$data$Compare<-factor(p_pssfpr$data$Compare,levels=treatment_order)
	pdf(paste0(level,"/NMDS_Bray_",comp1filename,".vs.",comp2filename,"_",level,".pdf"),width=6,height=6)
	color1=colors[[levels(p_pssfpr$data$Compare)[1]]]
	color2=colors[[levels(p_pssfpr$data$Compare)[2]]]
	ppp=p_pssfpr + scale_color_manual(values=c(color1,color2),labels=c(comp1title,comp2title)) + ggtitle(paste0("Fecal ",comp1title," & ",comp2title," NMDS - ",level))
	print(ppp)
	dev.off()
	bray_dist=phyloseq::distance(pssfpr, method="bray",type="samples")
	dataframe=data.frame(sample_data(pssfpr))
	pnova=adonis(bray_dist ~ Compare, data = dataframe)
	write_xlsx(pnova$aov.tab,paste0(level,"/Permanova_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
}

make_inext <- function(){
	pssfprfile=paste0(level,"/pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	pssfprfile1=paste0(      "pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	write.table(otu_table(pssfpr),pssfprfile,sep="\t",quote=FALSE)
	system(paste0("../make_inext_table_ave.pl ",pssfprfile1," Compare ",level))
	
	inextfile=paste0(level,"/inext_list_",pssfprfile1);
	inext.in=read_tsv(inextfile)
	inext.out=iNEXT(as.list(inext.in),q=c(0,1,2),datatype="abundance",nboot=100)
	script=paste0("inextout.",obj,".",level," <<- inext.out")
	eval(parse(text=script));
	script=paste0("save(inextout.",obj,".",level,",file='inextout.",obj,".",level,".Rdata')")
	eval(parse(text=script));

}

make_gginext <- function(){
	script=paste0("inext.out<-inextout.",obj,".",level)
	#script=paste0("inext.out<-inextout.1.",level)	
	eval(parse(text=script))
	test<-sort(c(comp1title,comp2title))
	if(test[1]==comp1title){
		color1=colors[[comp1a]]
		color2=colors[[comp2a]]
	}else{
		color1=colors[[comp2a]]
		color2=colors[[comp1a]]	
	}
#	color1=colors[[levels(inext.out$DataInfo$site)[1]]]
#	color2=colors[[levels(inext.out$DataInfo$site)[2]]]
	names(inext.out$iNextEst)[which(names(inext.out$iNextEst)==comp1a)]<-comp1title
	names(inext.out$iNextEst)[which(names(inext.out$iNextEst)==comp2a)]<-comp2title
	ggiNEXTout=ggiNEXT(inext.out, type=1, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Sample-size based rarefaction and extrapolation sampling curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Number of Sequences") + ylab(paste0(level," Diversity")) + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type1_",level,".pdf"),width=6,height=8)
	print(ppp)
	dev.off()
	ggiNEXTout=ggiNEXT(inext.out, type=2, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Sample completeness curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Number of Sequences") + ylab("Sample Coverage") + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type2_",level,".pdf"),width=6,height=4)
	print(ppp)
	dev.off()
	ggiNEXTout=ggiNEXT(inext.out, type=3, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Coverage based rarefaction and extrapolation sampling curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Sample Coverage") + ylab(paste0(level," Diversity")) + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type3_",level,".pdf"),width=6,height=8)
	print(ppp)
	dev.off()
}

inext_stat <- function(){
	pssfprfile=paste0(level,"/pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	system(paste0("../make_inext_table.pl ",pssfprfile," ",comp1filename,".vs.",comp2filename," ",level))
	inext.all=read_tsv(paste0(level,"/inext_",comp1filename,".vs.",comp2filename,".txt"),cols(.default = col_double()),col_names=T)
	inext.all.out<-iNEXT(as.list(inext.all),q=c(0,1,2),datatype="abundance",nboot=0,size=rare)
	inext.all.out_iNextEst<-inext.all.out$iNextEst
	inext.q123=NULL
	for(i in 1:length(inext.all.out_iNextEst)) {
		inext.q123=rbind(inext.q123,
			cbind("name"=names(inext.all.out_iNextEst[i]),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==0),qD),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==1),qD),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==2),qD)
			)
		)
	}
	
	meta1=meta[meta$X %in% inext.q123$name,]
	names(meta1)[1]="SampleID"
	colnames(inext.q123)=c("SampleID","qD0","qD1","qD2")
	inext.q123m=merge(inext.q123,meta1,by="SampleID")
	ttest_all.qD0=     t.test(qD0 ~ Compare, data=inext.q123m)
	ttest_all.qD1=     t.test(qD1 ~ Compare, data=inext.q123m)
	ttest_all.qD2=     t.test(qD2 ~ Compare, data=inext.q123m)
	wtest_all.qD0=wilcox.test(qD0 ~ Compare, data=inext.q123m)
	wtest_all.qD1=wilcox.test(qD1 ~ Compare, data=inext.q123m)
	wtest_all.qD2=wilcox.test(qD2 ~ Compare, data=inext.q123m)
	
	ttest.rows=c("t statistic",
	"degrees of freedom",  
	"p.value",    
	"confidence interval 1",  
	"confidence interval 2",
	"estimated mean 1",  
	"estimated mean 1",
	"stderr",     
	"alternative hypothesis",
	"t test method")
	
	ttest_all_qD0=c(
	ttest_all.qD0$statistic,  
	ttest_all.qD0$parameter,  
	ttest_all.qD0$p.value,    
	ttest_all.qD0$conf.int,   
	ttest_all.qD0$estimate,   
	ttest_all.qD0$stderr,     
	ttest_all.qD0$alternative,
	ttest_all.qD0$method)
	
	ttest_all_qD1=c(
	ttest_all.qD1$statistic,  
	ttest_all.qD1$parameter,  
	ttest_all.qD1$p.value,    
	ttest_all.qD1$conf.int,   
	ttest_all.qD1$estimate,   
	ttest_all.qD1$stderr,     
	ttest_all.qD1$alternative,
	ttest_all.qD1$method)
			
	ttest_all_qD2=c(
	ttest_all.qD2$statistic,  
	ttest_all.qD2$parameter,  
	ttest_all.qD2$p.value,    
	ttest_all.qD2$conf.int,   
	ttest_all.qD2$estimate,   
	ttest_all.qD2$stderr,     
	ttest_all.qD2$alternative,
	ttest_all.qD2$method)
	
	ttest_all=data.frame("Hill_Number"=ttest.rows,qD0=ttest_all_qD0,qD1=ttest_all_qD1,qD2=ttest_all_qD2)
	
	write_xlsx(ttest_all, paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_T_test_",level,".xlsx"))
	write_xlsx(inext.q123,paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	
	wtest.rows=c("w statistic",
	"p.value",    
	"alternative hypothesis",
	"test method")
	
	wtest_all_qD0=c(
	wtest_all.qD0$statistic,  
	wtest_all.qD0$p.value,"true location shift is not equal to 0",   
	wtest_all.qD0$method)
	
	wtest_all_qD1=c(
	wtest_all.qD1$statistic,  
	wtest_all.qD1$p.value,"true location shift is not equal to 0",   
	wtest_all.qD1$method)
	
	wtest_all_qD2=c(
	wtest_all.qD2$statistic,  
	wtest_all.qD2$p.value,"true location shift is not equal to 0",   
	wtest_all.qD2$method)
	
	wtest_all=data.frame("Hill_Number"=wtest.rows,qD0=wtest_all_qD0,qD1=wtest_all_qD1,qD2=wtest_all_qD2)
	write_xlsx(wtest_all,paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_Wilcox_ranksum_test_",level,".xlsx"))
}

################################## metacoder

plot_body_site_diff <- function(x, site_1, site_2, output_name, seed = 1, outformat="svg", level="OTUs") {
	if(level == "Genus"){
	level1="Genera"
	}else{
	level1=level
	}
	set.seed(seed)
	x %>%
	mutate_obs("tax_prop", abundance = rowMeans(x$data$tax_prop[x$data$sample_data$sample_id])) %>%
	filter_taxa(abundance >= 0.001, reassign_obs = FALSE) %>%
	filter_taxa(taxon_names != "", reassign_obs = FALSE) %>% # Some taxonomic levels are not named
	filter_obs("diff_table", treatment_1 %in% c(site_1, site_2), treatment_2 %in% c(site_1, site_2)) %>%
	heat_tree(node_size_axis_label = paste0("Number of ",level1),
	node_size = n_obs,
	node_size_range=c(0.02,.06),
	node_label_size_range=c(0.02,.04),
	node_color_axis_label = "Log 2 ratio",
	node_color = log2_median_ratio,
	node_color_range = diverging_palette(),
	node_color_trans = "linear",
	initial_layout = "re", layout = "da",
	node_color_interval = c(-5, 5),
	edge_color_interval = c(-5, 5),
	node_label = taxon_names,
	output_file = paste0(output_name, "--", site_1, ".vs.", site_2, "_",level,".", outformat)
	)
}

make_metacoder <- function(){
	#mc=parse_phyloseq(pssfpr)
	mc=parse_phyloseq(ps)
	mc$data$otu_prop=calc_obs_props(mc,data="otu_table",cols=mc$data$sample_data$sample_id)
	mc$data$tax_prop=calc_taxon_abund(mc,data="otu_prop",cols=mc$data$sample_data$sample_id)
	comcol=mc$data$sample_data$Compare
	mc$data$diff_table <- compare_groups(mc, data='tax_prop', cols= mc$data$sample_data$sample_id, groups= comcol,combinations = list(c(comp1a,comp2a)))
	mc <- mutate_obs(mc, "diff_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"), log2_median_ratio = ifelse(wilcox_p_value < 0.05 | is.na(wilcox_p_value), log2_median_ratio, 0))
	outname<-paste0(level,"/MetacodeR")
	plot_body_site_diff(mc,comp1a,comp2a,outname,1,"pdf",level)
	taxonomy=classifications(mc)
	difftable=cbind(taxonomy,mc$data$diff_table)
	write_xlsx(difftable,paste0(level,"/MetacodeR_diff_table_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	taxprotable=cbind(taxonomy,mc$data$tax_prop)
	write_xlsx(taxprotable,paste0(level,"/MetacodeR_tax_proportion_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	taxonomy=unite(mc$data$tax_data[,3:8], Taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";")
	otuprotable=cbind(taxonomy,mc$data$otu_prop)
	write_xlsx(otuprotable,paste0(level,"/MetacodeR_otu_proportion_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
}


################################## perform by calling functions
make_metacoder()
make_ps()
make_nmds_perma()
make_inext()
make_gginext()
inext_stat()

############## Species ###########
level="Species"
ps=tax_glom(physeq,taxrank="Species")
##################################
#    NMDS, iNEXT, MetacodeR      #
##################################
#rare=56000
#library(phyloseq)
#library(vegan) # adonis
#library(writexl)
#library(ggplot2)
#library(iNEXT)
#library(readr)
#library(dplyr) #select
#library(metacoder)#filter_taxa function in conflict with phyloseq
#library(writexl)                          
#library(tidyr) #unite select
#library(plyr)
#colors<-list();
#colors[["ISS"]]<- "red"
#colors[["ISS_G"]]<- "blue"
#colors[["final LAR"]]<- "#00661A"
#colors[["final LAR_G"]]<- "#666600"
#colors[["t.2 - 8 Weeks LAR"]]<- "#0000B2"
#colors[["t.1 - 4 Weeks LAR"]]<- "#D93600"
#colors[["t.0 - Pre Launch LAR"]]<- "#7F00FF"
#colors[["t.2 - 8 Weeks LAR_G"]]<- "#006DD9"
#colors[["t.1 - 4 Weeks LAR_G"]]<- "#FF8000"
#colors[["t.0 - Pre Launch LAR_G"]]<- "#9673FF"
#colors[["baseline"]]<- "#444444"
##################################
#load("../physeq.Rdata")
#load("../meta.Rdata")
#load("../inextout.ASV.Rdata")
##################################
#ps=physeq
#ps=tax_glom(physeq,taxrank="Species")
#ps=tax_glom(physeq,taxrank="Genus")
#level="ASV";
#level="Species";
#level="Genus";
##################################
make_ps <- function(){
	script=paste0("subset_samples(ps,( ",key1," == comp1 | ",key2," == comp2))")
	pss=eval(parse(text=script));
	pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
	pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
	pssfpr<<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
}

make_nmds_perma <- function(){
	pssfpr.ord=phyloseq::ordinate(pssfpr,"NMDS","bray")
	treatment_order=c(comp1a,comp2a)
	p_pssfpr=phyloseq::plot_ordination(pssfpr,pssfpr.ord,type="samples",color="Compare") + labs(color = "Treatment") + geom_point(size=2) + stat_ellipse()  
	p_pssfpr$data$Compare<-as.character(p_pssfpr$data$Compare)
	p_pssfpr$data$Compare<-factor(p_pssfpr$data$Compare,levels=treatment_order)
	pdf(paste0(level,"/NMDS_Bray_",comp1filename,".vs.",comp2filename,"_",level,".pdf"),width=6,height=6)
	color1=colors[[levels(p_pssfpr$data$Compare)[1]]]
	color2=colors[[levels(p_pssfpr$data$Compare)[2]]]
	ppp=p_pssfpr + scale_color_manual(values=c(color1,color2),labels=c(comp1title,comp2title)) + ggtitle(paste0("Fecal ",comp1title," & ",comp2title," NMDS - ",level))
	print(ppp)
	dev.off()
	bray_dist=phyloseq::distance(pssfpr, method="bray",type="samples")
	dataframe=data.frame(sample_data(pssfpr))
	pnova=adonis(bray_dist ~ Compare, data = dataframe)
	write_xlsx(pnova$aov.tab,paste0(level,"/Permanova_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
}

make_inext <- function(){
	pssfprfile=paste0(level,"/pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	pssfprfile1=paste0(      "pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	write.table(otu_table(pssfpr),pssfprfile,sep="\t",quote=FALSE)
	system(paste0("../make_inext_table_ave.pl ",pssfprfile1," Compare ",level))
	
	inextfile=paste0(level,"/inext_list_",pssfprfile1);
	inext.in=read_tsv(inextfile)
	inext.out=iNEXT(as.list(inext.in),q=c(0,1,2),datatype="abundance",nboot=100)
	script=paste0("inextout.",obj,".",level," <<- inext.out")
	eval(parse(text=script));
	script=paste0("save(inextout.",obj,".",level,",file='inextout.",obj,".",level,".Rdata')")
	eval(parse(text=script));

}

make_gginext <- function(){
	script=paste0("inext.out<-inextout.",obj,".",level)
	#script=paste0("inext.out<-inextout.1.",level)	
	eval(parse(text=script))
	test<-sort(c(comp1title,comp2title))
	if(test[1]==comp1title){
		color1=colors[[comp1a]]
		color2=colors[[comp2a]]
	}else{
		color1=colors[[comp2a]]
		color2=colors[[comp1a]]	
	}
#	color1=colors[[levels(inext.out$DataInfo$site)[1]]]
#	color2=colors[[levels(inext.out$DataInfo$site)[2]]]
	names(inext.out$iNextEst)[which(names(inext.out$iNextEst)==comp1a)]<-comp1title
	names(inext.out$iNextEst)[which(names(inext.out$iNextEst)==comp2a)]<-comp2title
	ggiNEXTout=ggiNEXT(inext.out, type=1, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Sample-size based rarefaction and extrapolation sampling curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Number of Sequences") + ylab(paste0(level," Diversity")) + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type1_",level,".pdf"),width=6,height=8)
	print(ppp)
	dev.off()
	ggiNEXTout=ggiNEXT(inext.out, type=2, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Sample completeness curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Number of Sequences") + ylab("Sample Coverage") + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type2_",level,".pdf"),width=6,height=4)
	print(ppp)
	dev.off()
	ggiNEXTout=ggiNEXT(inext.out, type=3, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Coverage based rarefaction and extrapolation sampling curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Sample Coverage") + ylab(paste0(level," Diversity")) + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type3_",level,".pdf"),width=6,height=8)
	print(ppp)
	dev.off()
}

inext_stat <- function(){
	pssfprfile=paste0(level,"/pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	system(paste0("../make_inext_table.pl ",pssfprfile," ",comp1filename,".vs.",comp2filename," ",level))
	inext.all=read_tsv(paste0(level,"/inext_",comp1filename,".vs.",comp2filename,".txt"),cols(.default = col_double()),col_names=T)
	inext.all.out<-iNEXT(as.list(inext.all),q=c(0,1,2),datatype="abundance",nboot=0,size=rare)
	inext.all.out_iNextEst<-inext.all.out$iNextEst
	inext.q123=NULL
	for(i in 1:length(inext.all.out_iNextEst)) {
		inext.q123=rbind(inext.q123,
			cbind("name"=names(inext.all.out_iNextEst[i]),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==0),qD),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==1),qD),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==2),qD)
			)
		)
	}
	
	meta1=meta[meta$X %in% inext.q123$name,]
	names(meta1)[1]="SampleID"
	colnames(inext.q123)=c("SampleID","qD0","qD1","qD2")
	inext.q123m=merge(inext.q123,meta1,by="SampleID")
	ttest_all.qD0=     t.test(qD0 ~ Compare, data=inext.q123m)
	ttest_all.qD1=     t.test(qD1 ~ Compare, data=inext.q123m)
	ttest_all.qD2=     t.test(qD2 ~ Compare, data=inext.q123m)
	wtest_all.qD0=wilcox.test(qD0 ~ Compare, data=inext.q123m)
	wtest_all.qD1=wilcox.test(qD1 ~ Compare, data=inext.q123m)
	wtest_all.qD2=wilcox.test(qD2 ~ Compare, data=inext.q123m)
	
	ttest.rows=c("t statistic",
	"degrees of freedom",  
	"p.value",    
	"confidence interval 1",  
	"confidence interval 2",
	"estimated mean 1",  
	"estimated mean 1",
	"stderr",     
	"alternative hypothesis",
	"t test method")
	
	ttest_all_qD0=c(
	ttest_all.qD0$statistic,  
	ttest_all.qD0$parameter,  
	ttest_all.qD0$p.value,    
	ttest_all.qD0$conf.int,   
	ttest_all.qD0$estimate,   
	ttest_all.qD0$stderr,     
	ttest_all.qD0$alternative,
	ttest_all.qD0$method)
	
	ttest_all_qD1=c(
	ttest_all.qD1$statistic,  
	ttest_all.qD1$parameter,  
	ttest_all.qD1$p.value,    
	ttest_all.qD1$conf.int,   
	ttest_all.qD1$estimate,   
	ttest_all.qD1$stderr,     
	ttest_all.qD1$alternative,
	ttest_all.qD1$method)
			
	ttest_all_qD2=c(
	ttest_all.qD2$statistic,  
	ttest_all.qD2$parameter,  
	ttest_all.qD2$p.value,    
	ttest_all.qD2$conf.int,   
	ttest_all.qD2$estimate,   
	ttest_all.qD2$stderr,     
	ttest_all.qD2$alternative,
	ttest_all.qD2$method)
	
	ttest_all=data.frame("Hill_Number"=ttest.rows,qD0=ttest_all_qD0,qD1=ttest_all_qD1,qD2=ttest_all_qD2)
	
	write_xlsx(ttest_all, paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_T_test_",level,".xlsx"))
	write_xlsx(inext.q123,paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	
	wtest.rows=c("w statistic",
	"p.value",    
	"alternative hypothesis",
	"test method")
	
	wtest_all_qD0=c(
	wtest_all.qD0$statistic,  
	wtest_all.qD0$p.value,"true location shift is not equal to 0",   
	wtest_all.qD0$method)
	
	wtest_all_qD1=c(
	wtest_all.qD1$statistic,  
	wtest_all.qD1$p.value,"true location shift is not equal to 0",   
	wtest_all.qD1$method)
	
	wtest_all_qD2=c(
	wtest_all.qD2$statistic,  
	wtest_all.qD2$p.value,"true location shift is not equal to 0",   
	wtest_all.qD2$method)
	
	wtest_all=data.frame("Hill_Number"=wtest.rows,qD0=wtest_all_qD0,qD1=wtest_all_qD1,qD2=wtest_all_qD2)
	write_xlsx(wtest_all,paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_Wilcox_ranksum_test_",level,".xlsx"))
}

################################## metacoder

plot_body_site_diff <- function(x, site_1, site_2, output_name, seed = 1, outformat="svg", level="OTUs") {
	if(level == "Genus"){
	level1="Genera"
	}else{
	level1=level
	}
	set.seed(seed)
	x %>%
	mutate_obs("tax_prop", abundance = rowMeans(x$data$tax_prop[x$data$sample_data$sample_id])) %>%
	filter_taxa(abundance >= 0.001, reassign_obs = FALSE) %>%
	filter_taxa(taxon_names != "", reassign_obs = FALSE) %>% # Some taxonomic levels are not named
	filter_obs("diff_table", treatment_1 %in% c(site_1, site_2), treatment_2 %in% c(site_1, site_2)) %>%
	heat_tree(node_size_axis_label = paste0("Number of ",level1),
	node_size = n_obs,
	node_size_range=c(0.02,.06),
	node_label_size_range=c(0.02,.04),
	node_color_axis_label = "Log 2 ratio",
	node_color = log2_median_ratio,
	node_color_range = diverging_palette(),
	node_color_trans = "linear",
	initial_layout = "re", layout = "da",
	node_color_interval = c(-5, 5),
	edge_color_interval = c(-5, 5),
	node_label = taxon_names,
	output_file = paste0(output_name, "--", site_1, ".vs.", site_2, "_",level,".", outformat)
	)
}

make_metacoder <- function(){
	#mc=parse_phyloseq(pssfpr)
	mc=parse_phyloseq(ps)
	mc$data$otu_prop=calc_obs_props(mc,data="otu_table",cols=mc$data$sample_data$sample_id)
	mc$data$tax_prop=calc_taxon_abund(mc,data="otu_prop",cols=mc$data$sample_data$sample_id)
	comcol=mc$data$sample_data$Compare
	mc$data$diff_table <- compare_groups(mc, data='tax_prop', cols= mc$data$sample_data$sample_id, groups= comcol,combinations = list(c(comp1a,comp2a)))
	mc <- mutate_obs(mc, "diff_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"), log2_median_ratio = ifelse(wilcox_p_value < 0.05 | is.na(wilcox_p_value), log2_median_ratio, 0))
	outname<-paste0(level,"/MetacodeR")
	plot_body_site_diff(mc,comp1a,comp2a,outname,1,"pdf",level)
	taxonomy=classifications(mc)
	difftable=cbind(taxonomy,mc$data$diff_table)
	write_xlsx(difftable,paste0(level,"/MetacodeR_diff_table_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	taxprotable=cbind(taxonomy,mc$data$tax_prop)
	write_xlsx(taxprotable,paste0(level,"/MetacodeR_tax_proportion_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	taxonomy=unite(mc$data$tax_data[,3:8], Taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";")
	otuprotable=cbind(taxonomy,mc$data$otu_prop)
	write_xlsx(otuprotable,paste0(level,"/MetacodeR_otu_proportion_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
}


################################## perform by calling functions
make_metacoder()
make_ps()
make_nmds_perma()
make_inext()
make_gginext()
inext_stat()

##############  Genus  ###########
level="Genus"
ps=tax_glom(physeq,taxrank="Genus")
##################################
#    NMDS, iNEXT, MetacodeR      #
##################################
#rare=56000
#library(phyloseq)
#library(vegan) # adonis
#library(writexl)
#library(ggplot2)
#library(iNEXT)
#library(readr)
#library(dplyr) #select
#library(metacoder)#filter_taxa function in conflict with phyloseq
#library(writexl)                          
#library(tidyr) #unite select
#library(plyr)
#colors<-list();
#colors[["ISS"]]<- "red"
#colors[["ISS_G"]]<- "blue"
#colors[["final LAR"]]<- "#00661A"
#colors[["final LAR_G"]]<- "#666600"
#colors[["t.2 - 8 Weeks LAR"]]<- "#0000B2"
#colors[["t.1 - 4 Weeks LAR"]]<- "#D93600"
#colors[["t.0 - Pre Launch LAR"]]<- "#7F00FF"
#colors[["t.2 - 8 Weeks LAR_G"]]<- "#006DD9"
#colors[["t.1 - 4 Weeks LAR_G"]]<- "#FF8000"
#colors[["t.0 - Pre Launch LAR_G"]]<- "#9673FF"
#colors[["baseline"]]<- "#444444"
##################################
#load("../physeq.Rdata")
#load("../meta.Rdata")
#load("../inextout.ASV.Rdata")
##################################
#ps=physeq
#ps=tax_glom(physeq,taxrank="Species")
#ps=tax_glom(physeq,taxrank="Genus")
#level="ASV";
#level="Species";
#level="Genus";
##################################
make_ps <- function(){
	script=paste0("subset_samples(ps,( ",key1," == comp1 | ",key2," == comp2))")
	pss=eval(parse(text=script));
	pssf=phyloseq::filter_taxa(pss, function(x) length(which(x >= 10))>=2 | sum(x) >= 100, TRUE)
	pssfp=prune_samples(sample_sums(pssf)>=rare, pssf)
	pssfpr<<-rarefy_even_depth(pssfp,rngseed=1,sample.size=rare,replace=F)
}

make_nmds_perma <- function(){
	pssfpr.ord=phyloseq::ordinate(pssfpr,"NMDS","bray")
	treatment_order=c(comp1a,comp2a)
	p_pssfpr=phyloseq::plot_ordination(pssfpr,pssfpr.ord,type="samples",color="Compare") + labs(color = "Treatment") + geom_point(size=2) + stat_ellipse()  
	p_pssfpr$data$Compare<-as.character(p_pssfpr$data$Compare)
	p_pssfpr$data$Compare<-factor(p_pssfpr$data$Compare,levels=treatment_order)
	pdf(paste0(level,"/NMDS_Bray_",comp1filename,".vs.",comp2filename,"_",level,".pdf"),width=6,height=6)
	color1=colors[[levels(p_pssfpr$data$Compare)[1]]]
	color2=colors[[levels(p_pssfpr$data$Compare)[2]]]
	ppp=p_pssfpr + scale_color_manual(values=c(color1,color2),labels=c(comp1title,comp2title)) + ggtitle(paste0("Fecal ",comp1title," & ",comp2title," NMDS - ",level))
	print(ppp)
	dev.off()
	bray_dist=phyloseq::distance(pssfpr, method="bray",type="samples")
	dataframe=data.frame(sample_data(pssfpr))
	pnova=adonis(bray_dist ~ Compare, data = dataframe)
	write_xlsx(pnova$aov.tab,paste0(level,"/Permanova_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
}

make_inext <- function(){
	pssfprfile=paste0(level,"/pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	pssfprfile1=paste0(      "pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	write.table(otu_table(pssfpr),pssfprfile,sep="\t",quote=FALSE)
	system(paste0("../make_inext_table_ave.pl ",pssfprfile1," Compare ",level))
	
	inextfile=paste0(level,"/inext_list_",pssfprfile1);
	inext.in=read_tsv(inextfile)
	inext.out=iNEXT(as.list(inext.in),q=c(0,1,2),datatype="abundance",nboot=100)
	script=paste0("inextout.",obj,".",level," <<- inext.out")
	eval(parse(text=script));
	script=paste0("save(inextout.",obj,".",level,",file='inextout.",obj,".",level,".Rdata')")
	eval(parse(text=script));

}

make_gginext <- function(){
	script=paste0("inext.out<-inextout.",obj,".",level)
	#script=paste0("inext.out<-inextout.1.",level)	
	eval(parse(text=script))
	test<-sort(c(comp1title,comp2title))
	if(test[1]==comp1title){
		color1=colors[[comp1a]]
		color2=colors[[comp2a]]
	}else{
		color1=colors[[comp2a]]
		color2=colors[[comp1a]]	
	}
#	color1=colors[[levels(inext.out$DataInfo$site)[1]]]
#	color2=colors[[levels(inext.out$DataInfo$site)[2]]]
	names(inext.out$iNextEst)[which(names(inext.out$iNextEst)==comp1a)]<-comp1title
	names(inext.out$iNextEst)[which(names(inext.out$iNextEst)==comp2a)]<-comp2title
	ggiNEXTout=ggiNEXT(inext.out, type=1, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Sample-size based rarefaction and extrapolation sampling curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Number of Sequences") + ylab(paste0(level," Diversity")) + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type1_",level,".pdf"),width=6,height=8)
	print(ppp)
	dev.off()
	ggiNEXTout=ggiNEXT(inext.out, type=2, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Sample completeness curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Number of Sequences") + ylab("Sample Coverage") + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type2_",level,".pdf"),width=6,height=4)
	print(ppp)
	dev.off()
	ggiNEXTout=ggiNEXT(inext.out, type=3, facet.var="order") + facet_wrap(~order, scales="free", ncol=1) + theme_bw(base_size=12)
	ggiNEXTout= ggiNEXTout + ggtitle(paste0("Coverage based rarefaction and extrapolation sampling curves:\nFecal ",comp1title," vs ",comp2title," - ",level," level"))
	ggiNEXTout= ggiNEXTout + xlab("Sample Coverage") + ylab(paste0(level," Diversity")) + scale_color_manual(values=c(color1,color2))
	ppp=ggiNEXTout
	pdf(paste0(level,"/iNEXT_Hill_Numbers_",comp1filename,".vs.",comp2filename,"_type3_",level,".pdf"),width=6,height=8)
	print(ppp)
	dev.off()
}

inext_stat <- function(){
	pssfprfile=paste0(level,"/pssfpr_",comp1filename,".vs.",comp2filename,"_",level,".txt")
	system(paste0("../make_inext_table.pl ",pssfprfile," ",comp1filename,".vs.",comp2filename," ",level))
	inext.all=read_tsv(paste0(level,"/inext_",comp1filename,".vs.",comp2filename,".txt"),cols(.default = col_double()),col_names=T)
	inext.all.out<-iNEXT(as.list(inext.all),q=c(0,1,2),datatype="abundance",nboot=0,size=rare)
	inext.all.out_iNextEst<-inext.all.out$iNextEst
	inext.q123=NULL
	for(i in 1:length(inext.all.out_iNextEst)) {
		inext.q123=rbind(inext.q123,
			cbind("name"=names(inext.all.out_iNextEst[i]),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==0),qD),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==1),qD),
				select(filter(inext.all.out_iNextEst[[i]], m==rare & order==2),qD)
			)
		)
	}
	
	meta1=meta[meta$X %in% inext.q123$name,]
	names(meta1)[1]="SampleID"
	colnames(inext.q123)=c("SampleID","qD0","qD1","qD2")
	inext.q123m=merge(inext.q123,meta1,by="SampleID")
	ttest_all.qD0=     t.test(qD0 ~ Compare, data=inext.q123m)
	ttest_all.qD1=     t.test(qD1 ~ Compare, data=inext.q123m)
	ttest_all.qD2=     t.test(qD2 ~ Compare, data=inext.q123m)
	wtest_all.qD0=wilcox.test(qD0 ~ Compare, data=inext.q123m)
	wtest_all.qD1=wilcox.test(qD1 ~ Compare, data=inext.q123m)
	wtest_all.qD2=wilcox.test(qD2 ~ Compare, data=inext.q123m)
	
	ttest.rows=c("t statistic",
	"degrees of freedom",  
	"p.value",    
	"confidence interval 1",  
	"confidence interval 2",
	"estimated mean 1",  
	"estimated mean 1",
	"stderr",     
	"alternative hypothesis",
	"t test method")
	
	ttest_all_qD0=c(
	ttest_all.qD0$statistic,  
	ttest_all.qD0$parameter,  
	ttest_all.qD0$p.value,    
	ttest_all.qD0$conf.int,   
	ttest_all.qD0$estimate,   
	ttest_all.qD0$stderr,     
	ttest_all.qD0$alternative,
	ttest_all.qD0$method)
	
	ttest_all_qD1=c(
	ttest_all.qD1$statistic,  
	ttest_all.qD1$parameter,  
	ttest_all.qD1$p.value,    
	ttest_all.qD1$conf.int,   
	ttest_all.qD1$estimate,   
	ttest_all.qD1$stderr,     
	ttest_all.qD1$alternative,
	ttest_all.qD1$method)
			
	ttest_all_qD2=c(
	ttest_all.qD2$statistic,  
	ttest_all.qD2$parameter,  
	ttest_all.qD2$p.value,    
	ttest_all.qD2$conf.int,   
	ttest_all.qD2$estimate,   
	ttest_all.qD2$stderr,     
	ttest_all.qD2$alternative,
	ttest_all.qD2$method)
	
	ttest_all=data.frame("Hill_Number"=ttest.rows,qD0=ttest_all_qD0,qD1=ttest_all_qD1,qD2=ttest_all_qD2)
	
	write_xlsx(ttest_all, paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_T_test_",level,".xlsx"))
	write_xlsx(inext.q123,paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	
	wtest.rows=c("w statistic",
	"p.value",    
	"alternative hypothesis",
	"test method")
	
	wtest_all_qD0=c(
	wtest_all.qD0$statistic,  
	wtest_all.qD0$p.value,"true location shift is not equal to 0",   
	wtest_all.qD0$method)
	
	wtest_all_qD1=c(
	wtest_all.qD1$statistic,  
	wtest_all.qD1$p.value,"true location shift is not equal to 0",   
	wtest_all.qD1$method)
	
	wtest_all_qD2=c(
	wtest_all.qD2$statistic,  
	wtest_all.qD2$p.value,"true location shift is not equal to 0",   
	wtest_all.qD2$method)
	
	wtest_all=data.frame("Hill_Number"=wtest.rows,qD0=wtest_all_qD0,qD1=wtest_all_qD1,qD2=wtest_all_qD2)
	write_xlsx(wtest_all,paste0(level,"/Hill_numbers_",comp1filename,".vs.",comp2filename,"_Wilcox_ranksum_test_",level,".xlsx"))
}

################################## metacoder

plot_body_site_diff <- function(x, site_1, site_2, output_name, seed = 1, outformat="svg", level="OTUs") {
	if(level == "Genus"){
	level1="Genera"
	}else{
	level1=level
	}
	set.seed(seed)
	x %>%
	mutate_obs("tax_prop", abundance = rowMeans(x$data$tax_prop[x$data$sample_data$sample_id])) %>%
	filter_taxa(abundance >= 0.001, reassign_obs = FALSE) %>%
	filter_taxa(taxon_names != "", reassign_obs = FALSE) %>% # Some taxonomic levels are not named
	filter_obs("diff_table", treatment_1 %in% c(site_1, site_2), treatment_2 %in% c(site_1, site_2)) %>%
	heat_tree(node_size_axis_label = paste0("Number of ",level1),
	node_size = n_obs,
	node_size_range=c(0.02,.06),
	node_label_size_range=c(0.02,.04),
	node_color_axis_label = "Log 2 ratio",
	node_color = log2_median_ratio,
	node_color_range = diverging_palette(),
	node_color_trans = "linear",
	initial_layout = "re", layout = "da",
	node_color_interval = c(-5, 5),
	edge_color_interval = c(-5, 5),
	node_label = taxon_names,
	output_file = paste0(output_name, "--", site_1, ".vs.", site_2, "_",level,".", outformat)
	)
}

make_metacoder <- function(){
	#mc=parse_phyloseq(pssfpr)
	mc=parse_phyloseq(ps)
	mc$data$otu_prop=calc_obs_props(mc,data="otu_table",cols=mc$data$sample_data$sample_id)
	mc$data$tax_prop=calc_taxon_abund(mc,data="otu_prop",cols=mc$data$sample_data$sample_id)
	comcol=mc$data$sample_data$Compare
	mc$data$diff_table <- compare_groups(mc, data='tax_prop', cols= mc$data$sample_data$sample_id, groups= comcol,combinations = list(c(comp1a,comp2a)))
	mc <- mutate_obs(mc, "diff_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"), log2_median_ratio = ifelse(wilcox_p_value < 0.05 | is.na(wilcox_p_value), log2_median_ratio, 0))
	outname<-paste0(level,"/MetacodeR")
	plot_body_site_diff(mc,comp1a,comp2a,outname,1,"pdf",level)
	taxonomy=classifications(mc)
	difftable=cbind(taxonomy,mc$data$diff_table)
	write_xlsx(difftable,paste0(level,"/MetacodeR_diff_table_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	taxprotable=cbind(taxonomy,mc$data$tax_prop)
	write_xlsx(taxprotable,paste0(level,"/MetacodeR_tax_proportion_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
	taxonomy=unite(mc$data$tax_data[,3:8], Taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";")
	otuprotable=cbind(taxonomy,mc$data$otu_prop)
	write_xlsx(otuprotable,paste0(level,"/MetacodeR_otu_proportion_",comp1filename,".vs.",comp2filename,"_",level,".xlsx"))
}


################################## perform by calling functions
make_metacoder()
make_ps()
make_nmds_perma()
make_inext()
make_gginext()
inext_stat()

