setwd("C:/Users/edoardo/Desktop/Varese/TCGA Data/BRCA")
FileName<-"BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
a<-read.table(paste0(getwd(),"/data/",FileName),header = T,sep = "\t",stringsAsFactors = FALSE)
#extrapolate only the real normalized counts and transform in numbers
dfGeneExp<-as.data.frame(sapply(FUN = as.numeric,X = a[2:nrow(a),2:ncol(a)]))
sapply(FUN = class,X = dfGeneExp)
########################################################
#I will use the log2(norm_counts+1) as expression value#
########################################################
#dfGeneExplog<-log2(dfGeneExp)
samp<-names(a)
#even from sample label beware to remove the headers
SampName<-samp[2:length(samp)]
#verify that the # of column of dfgeneexp is the same of the number of element in sampname
ncol(dfGeneExp)==length(SampName)
#extarpolate the gene names (entrez|locusID)
genes<-as.character(a[,1][-1])
GeneName<-as.data.frame(t(as.data.frame(strsplit(genes,split = "|",fixed = T))))
#I have noticed in a couple of cases the entrez gene is not defined hence in these cases I keep the lucus ID
GeneName$name[GeneName[,1]=="?"]<-paste0(GeneName[GeneName[,1]=="?",1],"|",GeneName[GeneName[,1]=="?",2])
GeneName$name[GeneName[,1]!="?"]<-as.character(GeneName[GeneName[,1]!="?",1])
#finally the rigth labeling
GenesName<-GeneName$name
#verify that all the dfgeneex has the same number of rows of the length of genesname 
nrow(dfGeneExp)==length(GenesName)

tab<-t(as.data.frame(strsplit(samp,fixed = T,split = ".")[-1]))
#define a table to identify the composition of the sample, 
#remember that 11 is healthy, 01 is cancer
#(beware 10 is helthy blood sample not a good comparison respect tumor tissue. the rigth comparison should be cancer tissue with healthy tissue)
table(tab[,4])
#define the index for control and case
CaseIndex<-unname(which(tab[,4]=="01A"|tab[,4]=="01B"))
ControlIndex<-unname(which(tab[,4]=="11A"|tab[,4]=="11B"))
MetaIndex<-unname(which(tab[,4]=="06A"|tab[,4]=="06B"))
#start setting up the df
dfRN<-cbind.data.frame("patience"=tab[,3])
dfRN$patience<-as.character(dfRN$patience)
dfRN$sample[CaseIndex]<-"01"
dfRN$sample[ControlIndex]<-"11"
dfRN$sample[MetaIndex]<-"06"
#add the normalized counts
row<-which(GenesName=="RNASET2")
dfRN$norm_counts<-unlist(dfGeneExp[row,])
dfRN$log2Exp<-log2(dfRN$norm_counts+1)
#########################################
#########################################
#########################################
#identification of the matched sample
#verify whether all the sample are unique
dup<-which(table(tab[,3])>1)
tab[,3][122]
duplex<-{}
for(i in names(dup)){
  print(i)
  ind<-which(tab[,3]==i)+1
  duplex[[i]]<-samp[ind]
  
}
lengths(duplex)
names(which(lengths(duplex)==3))

#metodo laborioso per far diventare la lista un df ma con piu controllo sul carattere
m<-max(lengths(duplex))
df<-{}
for(i in 1:length(duplex)){
  print(i)
  r<-c(duplex[[i]],rep(0,m-length(duplex[[i]])))
  df<-rbind(df,r)
}
#fast metod to transform the list in a df but long to make back character and to remove na
#library(plyr)
#df<-ldply(duplex,rbind)
#df$`3`<-as.character(df$`3`)
#df[is.na(df[,4]),4]<-"0"
##########################################
##########################################
##########################################

#########################################
#########################################
#load the clinical data
FileName<-"BRCA.clin.merged.csv"
Clinical<-read.csv(paste0(getwd(),"/data/",FileName),header = T,sep = ",")
#extrapolate the patience code
patience<-toupper(as.character(unlist(Clinical[1348,])[-1]))
length(patience)
#histological type
type<-as.character(unlist(Clinical[1149,])[-1])
#there is also a histollogical type from the section I prefered to use the first label (there are only small differences)
#type2<-as.character(unlist(Clinical[3703,])[-1])
#notice that all the patience are unique
#length(unique(patience))==length(patience)
#extarpolate the stage code
stage<-as.character(unlist(Clinical[1460,])[-1])
table(stage)
#i don't know exactly but it might be that the stage x means that the stage is not defined
length(stage)

#recapitulate the composition of the df
#table(tab[,4])
#patseq<-unname(tab[,3])
#pat01<-patseq[CaseIndex]
#pat06<-patseq[MetaIndex]
#pat11<-patseq[ControlIndex]
#what patience are different from the clinical and the rnaseq?
#length(pat01)==length(unique(pat01)) there are unique patience
#match(pat01,patience)
#all the patience in pat01 are present in patience (clinical data)

#length(patience)>length(pat01), which are the missing?
length(unique(patience))==length(patience)
#there are no replicate names

#what are the patience with clinical data but not sequenced with the label 01
#patience[which(is.na(match(patience,pat01))==T,arr.ind = T)]
#are they present in any of the other patience data?
#match(patience[which(is.na(match(patience,pat01))==T,arr.ind = T)],pat06)
#match(patience[which(is.na(match(patience,pat01))==T,arr.ind = T)],pat01)
#match(patience[which(is.na(match(patience,pat01))==T,arr.ind = T)],pat11)
#they are not present in the df of the rnaseq data.
#setup a df with the staging and the patience
clindf<-data.frame(patience,stage,type,stringsAsFactors = FALSE)
#remove the na in order to avoid wrong subsetting
clindf[is.na(clindf)]<-"none"
#########################################################################################
#########################################################################################
#merge the information from the two df
#match the sequenced patience with the clinical data
sum(is.na(match(dfRN$patience[dfRN$sample=="01"],clindf$patience)))
#ok all the sequenced patience are in the clinical data
dfRN$stage[dfRN$sample=="01"]<-clindf$stage[match(dfRN$patience[dfRN$sample=="01"],clindf$patience)]
#check that the sample are correctly matched
s<-sample(dfRN$patience[dfRN$sample=="01"],10)
rig1<-match(s,dfRN$patience[dfRN$sample=="01"])
dfRN$stage[dfRN$sample=="01"][rig1]
clindf$stage[match(s,clindf$patience)]
#comparison
clindf$stage[match(s,clindf$patience)]==dfRN$stage[dfRN$sample=="01"][rig1]
###################################################
#merge the histological type
dfRN$histo_type[dfRN$sample=="01"]<-clindf$type[match(dfRN$patience[dfRN$sample=="01"],clindf$patience)]
#check again the matching
s2<-sample(dfRN$patience[dfRN$sample=="01"],10)
rig2<-match(s2,dfRN$patience[dfRN$sample=="01"])
dfRN$histo_type[dfRN$sample=="01"][rig2]
clindf$type[match(s2,clindf$patience)]
#comparison
clindf$type[match(s2,clindf$patience)]==dfRN$histo_type[dfRN$sample=="01"][rig2]

###################################################################################
#so far I cannot find a staging or histo type for meta, I just don't analyze them #
###################################################################################

################################################
#statistic for the comparison of type and stage#
################################################
library(ggplot2)
g<-ggplot(dfRN[dfRN$sample=="01"& dfRN$stage!="none"& dfRN$stage!="stage x",],aes(y=log2Exp,as.factor(stage)))+geom_boxplot()
g+xlab("")+ylab("log2Exp")+ggtitle("BRCA_RNASET2")
g2<-ggplot(dfRN[dfRN$sample=="01"& dfRN$histo_type!="none",],aes(y=log2Exp,as.factor(histo_type)))+geom_boxplot()
g2+xlab("")+ylab("log2Exp")+ggtitle("BRCA_RNASET2")+theme(axis.text.x = element_text(angle = 45, hjust = 1))

a<-aov(log2Exp~stage,dfRN[dfRN$sample=="01"& dfRN$stage!="none"& dfRN$stage!="stage x",])
s1<-summary(a)
a1<-aggregate(log2Exp~stage,dfRN[dfRN$sample=="01"& dfRN$stage!="none"& dfRN$stage!="stage x",],
              FUN = function(x) c("mean"=mean(x),"sd"=sd(x)) )
tab1<-table(dfRN$stage[dfRN$sample=="01"])


b<-aov(log2Exp~histo_type,dfRN[dfRN$sample=="01"& dfRN$histo_type!="none",])
s2<-summary(b)
a2<-aggregate(log2Exp~histo_type,dfRN[dfRN$sample=="01"& dfRN$histo_type!="none",],
              FUN = function(x) c("mean"=mean(x),"sd"=sd(x)) )
tab2<-table(dfRN$histo_type[dfRN$sample=="01"])
t2<-TukeyHSD(b)

write.csv(a1,paste0(getwd(),"/log~stage.csv"))
write.csv(a2,paste0(getwd(),"/log~type.csv"))
capture.output(s1,file="anova_log~stage.txt")
capture.output(tab1,file="tab_log~stage.txt")
capture.output(s2,file="anova_log~type.txt")
capture.output(tab2,file="tab_log~type.txt")
write.csv(t2$histo_type,paste0(getwd(),"/TukeyHSD_log~type.csv"))
################################
#statistic for normal vs cancer#
################################
dfn_vs_c<-dfRN[dfRN$sample=="01"|dfRN$sample=="11",]
ttest<-t.test(log2Exp~sample,dfn_vs_c)
g3<-ggplot(dfn_vs_c,aes(y=log2Exp,factor(sample,levels=c("11","01"),labels = c("control","cancer"))))+geom_boxplot()
g3+xlab("")+ylab("log2Exp")+ggtitle("BRCA_RNASET2")
capture.output(ttest,file="ttest_normal~cancer.txt")
tab3<-table(dfn_vs_c$sample)
capture.output(tab3,file="tab_normal~cancer.txt")
