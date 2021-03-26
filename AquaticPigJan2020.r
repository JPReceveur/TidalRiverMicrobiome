#Delaware Aquatic Pig
#From https://gist.github.com/grabear/018e86413b19b62a6bb8e72a9adba349
#Parse Silva function
parse_taxonomy_silva_128 <- function(char.vec){
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec = parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(Rank1="D_0__Unassigned", Rank2="D_1__Unassigned", Rank3="D_2__Unassigned", Rank4="D_3__Unassigned",
                  Rank5="D_4__Unassigned", Rank6="D_5__Unassigned", Rank7="D_6__Unassigned")
  }
  # Define the meaning of each prefix according to SILVA taxonomy
  Tranks = c(D_0="Kingdom", D_1="Phylum", D_2="Class", D_3="Order", D_4="Family", D_5="Genus", D_6="Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
  if( length(ti) == 0L ){
    warning(
      "No silva prefixes were found. \n",
      "Consider using parse_taxonomy_delfault() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec = char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Format character vectors for Ambiguous taxa
    if( length(ti) < 7 ){
      for (key in names(char.vec)) {
        if ( char.vec[key] == "Ambiguous_taxa" ) {
          tax_no <- (as.numeric(substr(key, 5, 5)) - 1)
          char.vec[key] = sprintf("D_%s__Ambiguous_taxa", tax_no)
        }
      }
      # Reset the trimmed indicies if Ambiguous taxa
      ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
    }
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec = gsub("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks = Tranks[substr(char.vec[ti], 1, 3)]
    # Replace, being sure to avoid prefixes notK present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
  }
  return(taxvec)
}



library(caTools)
library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(ggpubr)
library(Rmisc)
library(multcompView)
library(randomForestCI)

set.seed(10)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#000000","#CC79A7","#F0E442")
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
biom=import_biom("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Aquatic Pig\\AquaticPigWTax6.25.19.biom",parseFunction=parse_taxonomy_silva_128)
tax_table(biom)
#tax_table(biom) <- tax_table(biom)[,-c(5:10,14)]#remove dummy ranks

metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Aquatic Pig\\AquaticPigMetadataWDiversity6.25.19.tsv",header = TRUE)
#levels(metadata$DecompStage)<-c("SubmergedFresh","EarlyFloating","Floatingdecay","Advancedfloatingdecay","SunkenRemains")
metadata$DecompStage = factor(metadata$DecompStage, levels = c("SubmergedFresh","EarlyFloating","Floatingdecay","Advancedfloatingdecay","SunkenRemains")) #fixes x-axis labels
head(metadata)
tree=read_tree("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Aquatic Pig\\tree6.26.19.nwk")
SondeData<-read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Aquatic Pig\\AquaticPigSondeData.txt",header=TRUE)

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$id
physeq=merge_phyloseq(biom,sampdat,tree)
physeq


# #removing singletons and phyla not present in 10% of samples (Singletons already removed in QIIME)
# physeq_phy <- tax_glom(physeq, taxrank = 'Phylum')
# physeq_phyla <- filter_taxa(physeq_phy, function(x) sum(x > 1) > (0.10*length(x)), TRUE)
# physeq_phyla_rel <- transform_sample_counts(physeq_phyla, function(OTU) OTU/sum(OTU) )
# 






SondeData<-read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Aquatic Pig\\AquaticPigSondeData.txt",header=TRUE)
head(SondeData)
SondeData<-subset(SondeData, Time > 4)# REmove first four hours (burn-in)

SondeData$Time<-SondeData$Time/24
SalinityPlot<-ggplot(SondeData,aes(x=Time,y=Sal))+geom_line(color="black",size=0.5)+xlab("Day")+ylab("Salinity (ppt)")
SalinityPlot
dev.off()
tiff("Figures/SalinityPlot.tiff", width = 74, height = 74, units = 'mm', res = 1200)
SalinityPlot
dev.off()
#Temperature Figure
TempPlot<-ggplot(SondeData, aes(x = Time, y = Temp))+geom_line(size=0.5)+
  ylab("Temperature (°C)")+xlab("Day")+ scale_color_manual(values=cbPalette)+ theme(legend.title=element_blank())

TempPlot


dev.off()
tiff("TempFig.tiff", width = 84, height = 84, units = 'mm', res = 1200)
TempPlot
dev.off()

#####Combine Temp and Salinity plots
dev.off()
tiff("Figures/TempAndSalCombined.tiff", width = 74, height = 150, units = 'mm', res = 1200)
ggarrange(TempPlot,SalinityPlot,
          labels = c("a", "b"), nrow = 2)
dev.off()


#TAble of top overall taxa

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
PhylumAll=tax_glom(GPr, "Phylum")

PhylumLevel = filter_taxa(PhylumAll, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower tha 3%
FamilyAll=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(FamilyAll, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower tha 3%
GenusAll=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GenusAll, function(x) mean(x) > 2e-2, TRUE) #filter out any taxa lower tha 2%

df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata

df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata

#Phylum level plot
PhylumLevel3 = filter_taxa(PhylumAll, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower tha 1%

df <- psmelt(PhylumLevel3)
df$Abundance=df$Abundance*100

compare_means(Abundance~DecompStage, data=df,group.by="Phylum",method = "kruskal.test",p.adjust.method="fdr")
compare_means(Abundance ~ DecompStage, data = df, group.by = "Phylum", p.adjust.method = "fdr")

Trtdata <- ddply(df, c("Phylum","DecompStage"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
#SigList<-length(unique(Trtdata$Phylum))
#SigLetters2<-vector(length=NComparisons)
#vec<-unlist(lst)
Means=compare_means(Abundance ~ DecompStage, data = df, group.by = "Phylum", p.adjust.method = "fdr",method="wilcox.test")
#Means
# Multicomp letters for Single family
# Bacteroidetes<- subset(Means,Phylum=="Bacteroidetes")
# Hyphenated<-as.character(paste0(Bacteroidetes$group1,"-",Bacteroidetes$group2))
# difference<-Bacteroidetes$p.adj
# names(difference)<-Hyphenated
# 
# Letters<-multcompLetters(difference)
# Letters

#Multicomp letters for all families
# for (i in levels(Means$Phylum)){
#   Tax<-i
#   TaxAbundance<-subset(Means,Phylum==i )
#   Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
#   difference<-TaxAbundance$p.adj
#   names(difference)<-Hyphenated
#   Letters<-multcompLetters(difference)
#   #print(Letters)
#   SigList[i]<-Letters
#   
# }
# SigList
# vec<-unlist(SigList)
# vec<-vec[-1]

vec<-c("a","a","b","b","b",
       "ab","a","ab","ab","b",
       "ab","ab","a","ab","b",
       "ab","a","bc","c","c")
# 
 KWResults<-c("Kruskal-Wallis,\n P-adj < 0.001","KW, P-adj = 0.042","KW, P-adj = 0.035","    P-adj < 0.001")
# 
# 
dat_text <- data.frame(
  Phylum = c("Bacteroidetes", "Epsilonbacteraeota", "Firmicutes","Proteobacteria"),
  label   = c("Kruskal-Wallis,\n P-adj < 0.001", "KW, \n P-adj = 0.042", "KW, \n P-adj = 0.035","    P-adj \n< 0.001"),
  DecompStage     = c(3.6, 3.6, 3.6,4.5),
  mean     = c(75,75,75,75)
)


PhylumPlotDecompStage=ggplot(Trtdata, aes(x=DecompStage,y=mean))+geom_bar(aes(fill = DecompStage),colour="black", stat="identity")+
  facet_grid(~Phylum)+xlab("Decomp Stage")+ylab("Relative Abundance (> 3%, SEM)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank())+facet_wrap(~Phylum)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+geom_text(aes(x=DecompStage, y=mean+se+10,label=vec))+
  scale_fill_manual(values=cbPalette)+ theme(legend.position = "none")+
  scale_x_discrete(labels=c("SubmergedFresh" = "SF","EarlyFloating"= "EF","Floatingdecay"="FD","Advancedfloatingdecay"="AFD","SunkenRemains"="SR"))+
  geom_text(data=dat_text,aes(label=label,x=DecompStage,y=mean),size=3.5)#+

  
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
PhylumPlotDecompStage
#theme_set(theme_bw(base_size = 11)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


dev.off()
tiff("Figures/PhylumPlotDecompStage.tiff", width = 74, height = 74, units = 'mm', res = 600)
PhylumPlotDecompStage
dev.off()


#Family level plot
FamilyLevel3 = filter_taxa(FamilyAll, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower tha 3%

df <- psmelt(FamilyLevel3)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family","DecompStage"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata

compare_means(Abundance~DecompStage, data=df,group.by="Family",method = "kruskal.test",p.adjust.method="fdr")

NComparisons<-length(unique(metadata$DecompStage))*length(unique(Trtdata$Family))
SigList<-length(unique(Trtdata$Family))
SigLetters2<-vector(length=NComparisons)
#vec<-unlist(lst)
Means=compare_means(Abundance ~ DecompStage, data = df, group.by = "Family", p.adjust.method = "fdr")
Means




# for (i in levels(Means$Family)){
#   Tax<-i 
#   TaxAbundance<-subset(Means,Family==i )
#   Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
#   difference<-TaxAbundance$p.adj
#   names(difference)<-Hyphenated
#   Letters<-multcompLetters(difference)
#   #print(Letters)
#   SigList[i]<-Letters
#   
# }
# vec2<-unlist(SigList)
# vec2<-vec2[-1]

vec2<-c("a","ab","bc","bc","c",
       "a","b","c","bc","b",
       "ab","a","b","b","b",
       "ab","a","b","bc","c",
       "a","ab","b","ab","a",
       "a","a","a","ab","b",
       "a","a","a","a","a")

#unique(Trtdata$Family)
dat_text <- data.frame(
  Family = c("Aeromonadaceae", "Bacteroidaceae", "Burkholderiaceae","Chromobacteriaceae","Clostridiaceae 1","Moraxellaceae","Rhodocyclaceae"),
  label   = c("\n KW,P-adj \n< 0.001","P-adj =\n 0.001","P-adj = \n 0.003","P-adj \n< 0.001","P-adj =\n0.005","P-adj \n< 0.001","P-adj =\n 0.037"),
  DecompStage     = c(3.6,1.6,4.3,3.6,1.6,4.4,4),
  mean     = c(28,28,25,25,25,25,25)
)


KWResults<-c("Kruskal-Wallis,\n P-adj < 0.001","P-adj = 0.001","P-adj = 0.003","P-adj < 0.001","P-adj = 0.005","P-adj < 0.001","P-adj = 0.037")
FamilyPlotDecompStage=ggplot(Trtdata, aes(x=DecompStage,y=mean))+geom_bar(aes(fill = DecompStage),colour="black", stat="identity")+
  facet_grid(~Family)+xlab("Decomp Stage")+ylab("Relative Abundance (> 3%, SEM)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank())+facet_wrap(~Family)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+geom_text(aes(x=DecompStage, y=mean+se+5,label=vec2))+
  scale_fill_manual(values=cbPalette)+ theme(legend.position = "none")+
  scale_x_discrete(labels=c("SubmergedFresh" = "SF","EarlyFloating"= "EF","Floatingdecay"="FD","Advancedfloatingdecay"="AFD","SunkenRemains"="SR"))#+
  geom_text(data=dat_text,aes(label=label,x=DecompStage,y=mean),size=2.8)

#+annotate("text", label = KWResults, size = 2.5, x = 3.6, y = 33)
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
FamilyPlotDecompStage
theme_set(theme_bw(base_size = 7.5)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/FamilyPlotDecompStage.tiff", width = 74, height = 74, units = 'mm', res = 600)
FamilyPlotDecompStage
dev.off()

#Shannon Diversity
kruskal.test(shannon~ DecompStage, data=metadata) #fig.width = 8, fig.height = 8 for line above to change size
posthoc<-compare_means(shannon ~ DecompStage, data = metadata, p.adjust.method = "fdr")


Hyphenated<-as.character(paste0(posthoc$group1,"-",posthoc$group2))
difference<-posthoc$p.format
names(difference)<-Hyphenated

LettersShannon<-multcompLetters(difference)

#posthoc <- HSD.test(av,"DecompStage",group=TRUE)

stats<-summarySE(metadata,measurevar="shannon",groupvars=c("DecompStage",na.rm=TRUE))
ShannonDecomp<-ggplot(stats,aes(x=DecompStage,y=shannon,fill=DecompStage))+geom_bar(stat="identity",colour="black")+
  geom_errorbar(aes(ymin=shannon-se,ymax=shannon+se),color="black")+
  geom_text(aes(x=DecompStage, y=shannon+se+0.5,label=LettersShannon$Letters), position=position_dodge(width=0.9), size=3.5,color="black")+ #Write in labels from posthoc
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size=7.5))+ scale_fill_manual(values=cbPalette)+xlab("")+ylab("Shannon Diversity (SEM)")+labs(fill="Decomposition Stage")+guides(fill=FALSE)+
  scale_x_discrete(labels=c("SubmergedFresh" = "Submerged\n Fresh","EarlyFloating"= "Early\n Floating","Floatingdecay"="Floating\n Decay","Advancedfloatingdecay"="Advanced\n Floating","SunkenRemains"="Sunken\n Remains")) + 
  annotate("text", label = expression(paste("KW,",chi^2, "= 18.7, P-adj < 0.001")), size = 3.5, x = 2.5, y = 8.5)#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))
theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

ShannonDecomp

dev.off()
tiff("Figures/ShannonPlotDecompStage.tiff", width = 74, height = 74, units = 'mm', res = 600)
ShannonDecomp
dev.off()

#Faith PD
kruskal.test(faith_pd~ DecompStage, data=metadata) #fig.width = 8, fig.height = 8 for line above to change size
posthoc<-compare_means(faith_pd ~ DecompStage, data = metadata, p.adjust.method = "fdr")


Hyphenated<-as.character(paste0(posthoc$group1,"-",posthoc$group2))
difference<-posthoc$p.format
names(difference)<-Hyphenated

Letters<-multcompLetters(difference)

#posthoc <- HSD.test(av,"DecompStage",group=TRUE)

stats<-summarySE(metadata,measurevar="faith_pd",groupvars=c("DecompStage",na.rm=TRUE))
faith_pdDecomp<-ggplot(stats,aes(x=DecompStage,y=faith_pd,fill=DecompStage))+geom_bar(stat="identity",colour="black")+
  geom_errorbar(aes(ymin=faith_pd-se,ymax=faith_pd+se),color="black")+
  geom_text(aes(x=DecompStage, y=faith_pd+se+4,label=Letters$Letters), position=position_dodge(width=0.9), size=3.5,color="black")+ #Write in labels from posthoc
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size=7.5))+ scale_fill_manual(values=cbPalette)+xlab("")+ylab("Faith's Phylogenetic Diversity (SEM)")+labs(fill="Decomposition Stage")+guides(fill=FALSE)+
  scale_x_discrete(labels=c("SubmergedFresh" = "Submerged\n Fresh","EarlyFloating"= "Early\n Floating","Floatingdecay"="Floating\n Decay","Advancedfloatingdecay"="Advanced\n Floating","SunkenRemains"="Sunken\n Remains"))+
  annotate("text", label = expression(paste("KW,",chi^2, "= 19.3, P-adj < 0.001")), size = 3.5, x = 2.4, y = 55)#+ theme(legend.justification=c(0.05,0.95), legend.position=c(0.05,0.95))

dev.off()
tiff("Figures/faith_pdPlotDecompStage.tiff", width = 74, height = 74, units = 'mm', res = 600)
faith_pdDecomp
dev.off()



PhylumPlotDecompStage
FamilyPlotDecompStage
ShannonDecomp
faith_pdDecomp
theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/Figure2NoP.tiff", width = 174, height = 174, units = 'mm', res = 600)
ggarrange(PhylumPlotDecompStage, FamilyPlotDecompStage,ShannonDecomp,faith_pdDecomp,
          labels = c("a","b","c","d"),
          ncol = 2, nrow = 2)
dev.off()

###################
#Beta diversity
###################
GPdist=phyloseq::distance(physeq, "jaccard")
beta=betadisper(GPdist, sample_data(physeq)$DecompStage)
permutest(beta)
boxplot(beta)


adonis(GPdist ~ DecompStage, as(sample_data(physeq), "data.frame"))
#Submerged Fresh vs Early Floating
SFVSEF<-subset_samples(physeq,DecompStage=="SubmergedFresh"|DecompStage=="EarlyFloating")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))

#Submerged Fresh vs Floating Decay
SFVSEF<-subset_samples(physeq,DecompStage=="SubmergedFresh"|DecompStage=="Floatingdecay")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))


#Submerged Fresh vs advancedFloating Decay
SFVSEF<-subset_samples(physeq,DecompStage=="SubmergedFresh"|DecompStage=="Advancedfloatingdecay")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))


#Submerged Fresh vs Sunken Remains
SFVSEF<-subset_samples(physeq,DecompStage=="SubmergedFresh"|DecompStage=="SunkenRemains")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))

#Early Floating vs FloatingDecay
SFVSEF<-subset_samples(physeq,DecompStage=="EarlyFloating"|DecompStage=="Floatingdecay")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))

#Early Floating vs Advanced Floating Decay
SFVSEF<-subset_samples(physeq,DecompStage=="EarlyFloating"|DecompStage=="Advancedfloatingdecay")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))

#Early Floating vs Sunken Remains
SFVSEF<-subset_samples(physeq,DecompStage=="EarlyFloating"|DecompStage=="SunkenRemains")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))

# Floating Decay vs Advanced Floating Decay
SFVSEF<-subset_samples(physeq,DecompStage=="Floatingdecay"|DecompStage=="Advancedfloatingdecay")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))


# Floating Decay vs Sunken Remains
SFVSEF<-subset_samples(physeq,DecompStage=="Floatingdecay"|DecompStage=="SunkenRemains")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))


# Advanced Floating Decay vs Sunken Remains
SFVSEF<-subset_samples(physeq,DecompStage=="Advancedfloatingdecay"|DecompStage=="SunkenRemains")
GPdist=phyloseq::distance(SFVSEF, "jaccard")
adonis(GPdist ~ DecompStage, as(sample_data(SFVSEF), "data.frame"))




############
#PCoA by decomp stage
##############33
# metadata2<-metadata
# metadata2$DecompStage <- gsub('SubmergedFresh', 'Submerged Fresh', metadata2$DecompStage)
# metadata2$DecompStage <- gsub('EarlyFloating', 'Early Floating', metadata2$DecompStage)
# metadata2$DecompStage <- gsub('Floatingdecay', 'Floating Decay', metadata2$DecompStage)
# metadata2$DecompStage <- gsub('Advancedfloatingdecay', 'Advanced Floating Decay', metadata2$DecompStage)
# metadata2$DecompStage <- gsub('SunkenRemains', 'Sunken Remains', metadata2$DecompStage)
# metadata2$DecompStage = factor(metadata2$DecompStage, levels = c("Submerged Fresh","Early Floating","Floating Decay","Advanced Floating Decay","Sunken Remains")) #fixes x-axis labels

metadata2<-metadata
metadata2$DecompStage <- gsub('SubmergedFresh', 'SF', metadata2$DecompStage)
metadata2$DecompStage <- gsub('EarlyFloating', 'EF', metadata2$DecompStage)
metadata2$DecompStage <- gsub('Floatingdecay', 'FD', metadata2$DecompStage)
metadata2$DecompStage <- gsub('Advancedfloatingdecay', 'AFD', metadata2$DecompStage)
metadata2$DecompStage <- gsub('SunkenRemains', 'SR', metadata2$DecompStage)
metadata2$DecompStage = factor(metadata2$DecompStage, levels = c("SF","EF","FD","AFD","SR")) #fixes x-axis labels



sampdat2=sample_data(metadata2)
sample_names(sampdat2)=metadata2$id
physeq2=merge_phyloseq(biom,sampdat2,tree)


sample_data(physeq2)

set.seed(23)
ord=ordinate(physeq2,"PCoA", "jaccard")
ordplot=plot_ordination(physeq2, ord,"samples", color="DecompStage",shape="DecompStage")+geom_point(size=3.5)+scale_color_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot<-ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = DecompStage))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position="bottom",legend.background = element_rect(color="black"),legend.title = element_blank())+
  scale_shape_manual(values=c(20,15,16,17,18))
ordplot
theme_set(theme_bw(base_size = 13)+theme(legend.title=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#  labels=c("Submerged Fresh", "Early Floating", "Floating Decay","Advanced Floating Decay","Sunken Remains")+scale_fill_manual(values=cbPalette)++scale_colour_manual(values=cbPalette)
dev.off()

tiff("Figures/PCoADecompStage.tiff", width = 100, height = 100, units = 'mm', res = 600)
ordplot
dev.off()

##############
#Random Forest
##############

ForestData=GenusAll#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
response <- as.factor(sample_data(ForestData)$DecompStage)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:10, ]#

#imp.20

destroyX = function(es) {
  f = es
  for (row in c(1:nrow(f))){ #for each column in dataframe
    if (startsWith(row.names(f)[row], "X") == TRUE)  { #if starts with 'X' ..
      row.names(f)[row] <- substr(row.names(f)[row], 2, 100) #get rid of it
    }
    assign(deparse(substitute(es)), f, inherits = TRUE)
  }
}
destroyX(imp.20)
row.names(imp.20)


otunames <- row.names(imp.20)
r <- rownames(tax_table(ForestData)) %in% otunames
otunames
PredictorTable<-kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification
PredictorTable

GenusRandomForestSubset = subset_taxa(GenusAll, row.names(tax_table(GenusAll))%in% otunames)
GenusRandomForestSubset

df <- psmelt(GenusRandomForestSubset)
df$Abundance=df$Abundance*100

dfSorted<- df[with(df, order(Genus)), ]
head(dfSorted)

Trtdata <- ddply(dfSorted, c("Genus", "DecompStage"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
TrtdataSorted<- Trtdata[with(Trtdata, order(Genus)), ]
TrtdataSorted
cdataplot=ggplot(Trtdata, aes(x=DecompStage,y=mean))+geom_bar(aes(fill = DecompStage),colour="black", stat="identity")+ facet_grid(~Genus)+xlab("Decomposition Stage")+ylab("Relative Abundance (%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+facet_wrap(~Genus,scales = "free_y")+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+scale_fill_manual(values=cbPalette)
cdataplot

compare_means(Abundance ~ DecompStage, data = dfSorted, group.by = "Genus", p.adjust.method = "fdr",method="kruskal.test")

Means<-compare_means(Abundance ~ DecompStage, data = dfSorted, group.by = "Genus", p.adjust.method = "fdr")
Means
Means=compare_means(Abundance ~ DecompStage, data = dfSorted, group.by = "Genus", p.adjust.method = "fdr")

unique(Means$Genus)
SigList<-length(unique(Trtdata$Genus))

#SigList
for (i in unique(Means$Genus)){
  Tax<-i
  TaxAbundance<-subset(Means,Genus==i )
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  print(Letters)
  SigList[i]<-Letters
  
}

vec<-unlist(SigList)

vec<-vec[-1]
#vec
#Uncultured bacteria is from Gemmatimonadaceae
#Sva0081 group is from Desulfobacteraceae
#MM2 is from Methylophilaceae



TrtdataSorted2<-TrtdataSorted
TrtdataSorted2$Genus <- gsub('Sva0081 sediment group', 'Sva0081 group', TrtdataSorted2$Genus)
TrtdataSorted2$Genus <- gsub('uncultured', 'uncultured \nGemmatimonadaceae', TrtdataSorted2$Genus)
TrtdataSorted2$Genus <- gsub('MM2', 'MM2 \n(Methylophilaceae)', TrtdataSorted2$Genus)


#mean+se+(mean+se*.1)
GenusIndicators=ggplot(TrtdataSorted2, aes(x=DecompStage,y=mean))+geom_bar(aes(fill = DecompStage),colour="black", stat="identity")+ facet_wrap(~Genus,scale="free_y")+
  ylab("Relative Abundance (%, SEM)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+ scale_fill_manual(values=cbPalette)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.justification=c(0.05,0.95), legend.position=c(0.52,0.15))+ theme(legend.title = element_blank()) +
  scale_x_discrete(labels=c("SubmergedFresh" = "SF","EarlyFloating"= "EF","Floatingdecay"="FD","Advancedfloatingdecay"="AFD","SunkenRemains"="SR"))+
  theme(legend.position = "none",strip.text.x = element_text(size = 9))+
  geom_text(aes(x=DecompStage, y=c(9,28,23,3,2,
                                 30.5,3,2.2,1.7,2.1,
                                 5.9,2.8,.6,.65,.55,
                                 .55,3.1,.75,.25,0.25,
                                 .35,.8,4.5,.8,.75,
                                 .43,.17,.1,.45,.8,
                                 .65,.25,.15,.55,1.3,
                                 .55,0.1,0.15,0.54,1.15,
                                 .1,.065,.05,.1,.46,
                                 15.5,5,2.5,1.5,1)),label=vec,size=4)
GenusIndicators
theme_set(theme_bw(base_size = 13.4)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/GenusLevelIndicators.tiff", width = 174, height = 174, units = 'mm', res = 600)
GenusIndicators
dev.off()




##########
#Random Forest OTUlevel
##########
library(ranger)
physeq
set.seed(154)

otu <- as.data.frame(t(otu_table(physeq)))


otu$id <- rownames(otu)
meta_sa <- metadata %>% select(id, ADH,Location)
otu <- merge(meta_sa, otu, by = 'id')

otu <- otu[,-1] #Remove id from predictors
otu <-otu[,-2] #Remove sample location from predictors
names(otu) <- make.names(names(otu))





sample = sample.split(otu, SplitRatio = .75)
train = subset(otu, sample == TRUE)
test  = subset(otu, sample == FALSE)

features<-setdiff(names(otu),"Day")
mTune<-tuneRF(
  x=otu[features],
  y=otu$ADH,
  ntreeTry = 5000,
  mtryStart = 5,
  stepFactor = 1.5,
  improve = 0.01,
  trace=F
)




m1 <- randomForest(
  formula = ADH ~ .,
  data    = otu,
  ntree= 500, keep.inbag=T,importance=T)

print(m1)


TreeFile<-data.frame(seq(1,500,by=1),m1$mse)
colnames(TreeFile)<-c("Tree","MSE")
head(TreeFile)
ErrorByTree<-ggplot(TreeFile,aes(x=Tree,y=MSE))+geom_line()+xlab("Tree")+ylab("Mean Squared Error")
ErrorByTree
#head(TreeFile)
#varImpPlot(m1,type=1)

X <- test
var_hat <- randomForestInfJack(m1, otu, calibrate = TRUE)

#head(var_hat)
plot(var_hat)

df <- data.frame(y = otu, var_hat)
df <- mutate(df, se = sqrt(var.hat))
#head(df)
(sqrt(m1$mse))[500] #RMSE at final tree


p1 <- ggplot(df, aes(x = y.ADH, y = y.hat))
p1<-p1 + geom_errorbar(aes(ymin=y.hat-se, ymax=y.hat+se),color="grey", width=.1) +
  geom_point() +
  geom_abline(intercept=0, slope=1, linetype=2)+
  xlab("True ADH") +
  ylab("Predicted ADH (ASV, \u00b1 95% CIs)")+annotate("text", label = "% Var explained = 80.8\n RMSE = 1943.2", size = 4, x = 4500, y = 11000)#+geom_jitter(width = 0.5)
p1

physeq

importance<-importance(m1, type=1)
head(importance)
imp <- data.frame(predictors = rownames(importance), importance)
imp

imp.sort <- arrange(imp, desc(X.IncMSE))
imp.sort
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
imp.20
TopPredictors<-ggplot(imp.20,aes(x=predictors,y=X.IncMSE))+geom_bar(stat="identity")+ylab("% Increase in MSE")+xlab("Top Predictors")+coord_flip()


#############33
#Random Forest Genus Level
##############
#Genus Level
set.seed(147)
Genus <- as.data.frame(t(otu_table(GenusAll)))

Genus$id <- rownames(Genus)
meta_sa <- metadata %>% select(id, ADH,Location)
Genus <- merge(meta_sa, Genus, by = 'id')

Genus <- Genus[,-1] #Remove id from predictors
Genus <-Genus[,-2] #Remove sample location from predictors
names(Genus) <- make.names(names(Genus))



features<-setdiff(names(Genus),"ADH")
mTune<-tuneRF(
  x=Genus[features],
  y=Genus$ADH,
  ntreeTry = 50,
  mtryStart = 2,
  stepFactor = 1.5,
  improve = 1,
  trace=F
)
tune

mGenus <- randomForest(
  formula = ADH ~ .,
  data    = Genus,
  ntree= 500, keep.inbag=T,importance=T)
#head(Genus)
print(mGenus)



TreeFile<-data.frame(seq(1,500,by=1),sqrt(mGenus$mse))
colnames(TreeFile)<-c("Tree","MSE")
head(TreeFile)
ErrorByTreeGenus<-ggplot(TreeFile,aes(x=Tree,y=MSE))+geom_line()+xlab("Tree")+ylab("RSME (Genus, ADH)")
ErrorByTreeGenus

var_hat <- randomForestInfJack(mGenus, Genus, calibrate = TRUE)

#head(var_hat)
plot(var_hat)

dfGenus <- data.frame(y = Genus, var_hat)
dfGenus <- mutate(dfGenus, se = sqrt(var.hat))
#head(df)
mGenus
(sqrt(mGenus$mse))[500] #RMSE at final tree


pGenus <- ggplot(dfGenus, aes(x = y.ADH, y = y.hat))
pGenus<-pGenus + geom_errorbar(aes(ymin=y.hat-se, ymax=y.hat+se),position = "dodge",color="grey", width=.1) +
  geom_point() +
  geom_abline(intercept=0, slope=1, linetype=2)+
  xlab("True ADH") +
  ylab("Predicted ADH (Genus, \u00b1 95% CIs)")+annotate("text", label = "% Var explained = 77.2\n RMSE = 2114.2", size = 4, x = 4500, y = 11000)#+geom_jitter(width = 0.5)
pGenus

#############33
#Random Forest Family Level
##############
#Family Level
set.seed(147)
Family <- as.data.frame(t(otu_table(FamilyAll)))

Family$id <- rownames(Family)
meta_sa <- metadata %>% select(id, ADH,Location)
Family <- merge(meta_sa, Family, by = 'id')

Family <- Family[,-1] #Remove id from predictors
Family <-Family[,-2] #Remove sample location from predictors
names(Family) <- make.names(names(Family))


mFamily <- randomForest(
  formula = ADH ~ .,
  data    = Family,
  ntree= 500, keep.inbag=T,importance=T)

print(mFamily)
var_hat <- randomForestInfJack(mFamily, Family, calibrate = TRUE)

#head(var_hat)
plot(var_hat)

df <- data.frame(y = Family, var_hat)
df <- mutate(df, se = sqrt(var.hat))
#head(df)
mFamily
(sqrt(mFamily$mse))[500] #RMSE at final tree

pFamily <- ggplot(df, aes(x = y.ADH, y = y.hat))
pFamily<-pFamily + geom_errorbar(aes(ymin=y.hat-se, ymax=y.hat+se),position = "dodge",color="grey", width=.1) +
  geom_point() +
  geom_abline(intercept=0, slope=1, linetype=2)+
  xlab("True ADH") +
  ylab("Predicted ADH (Family, \u00b1 95% CIs)")+annotate("text", label = "% Var explained = 77.3\n RMSE = 2111.9", size = 4, x = 4500, y = 11000)#+geom_jitter(width = 0.5)
pFamily


#############33
#Random Forest Phylum Level
##############
#Phylum Level
set.seed(147)
Phylum <- as.data.frame(t(otu_table(PhylumAll)))

Phylum$id <- rownames(Phylum)
meta_sa <- metadata %>% select(id, ADH,Location)
Phylum <- merge(meta_sa, Phylum, by = 'id')

Phylum <- Phylum[,-1] #Remove id from predictors
Phylum <-Phylum[,-2] #Remove sample location from predictors
names(Phylum) <- make.names(names(Phylum))


mPhylum <- randomForest(
  formula = ADH ~ .,
  data    = Phylum,
  ntree= 500, keep.inbag=T,importance=T)

print(mPhylum)

var_hat <- randomForestInfJack(mPhylum, Phylum, calibrate = TRUE)
var_hat
#head(var_hat)
plot(var_hat)

df <- data.frame(y = Phylum, var_hat)
df <- mutate(df, se = sqrt(var.hat))
#head(df)
df
mPhylum
(sqrt(mPhylum$mse))[500] #RMSE at final tree


pPhylum <- ggplot(df, aes(x = y.ADH, y = y.hat))
pPhylum<-pPhylum + geom_errorbar(aes(ymin=y.hat-se, ymax=y.hat+se),position = "dodge",color="grey", width=.1) +
  geom_abline(intercept=0, slope=1, linetype=2)+geom_point()+
  xlab("True ADH") +
  ylab("Predicted ADH (Phylum, \u00b1 95% CIs)")+annotate("text", label = "% Var explained = 79.67\n RMSE = 1997.1", size = 4, x = 4500, y = 11000)#+geom_jitter(width = 0.5)
pPhylum

#############3
#Join together ADH plots
#############
pPhylum
pFamily
pGenus
p1
theme_set(theme_bw(base_size = 11)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dev.off()
tiff("Figures/ADHRandomForest.tiff", width = 174, height = 174, units = 'mm', res = 1200)
ggarrange(pPhylum,pFamily,pGenus,p1,
          labels = c("a", "b","c","d"), ncol = 2,nrow=2)
dev.off()

