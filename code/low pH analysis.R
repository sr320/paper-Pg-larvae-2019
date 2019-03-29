setwd('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/no pH 8.2')
#load biostats package
source('~/Documents/Multivariate Stats/biostats.R')

geo.ab<-read.csv('ABACUS_cdhit0.9_lowpH_output.csv', header=T, row.names=1)

library(dplyr)
adjnsaf<-select(geo.ab, contains('ADJNSAF'))

#keep only proteins with at least 2 unique peptides
nsaf.uniq<-cbind(adjnsaf, geo.ab$ALL_NUMPEPSUNIQ)
twopeps<-subset(nsaf.uniq, select=X2017_JULY_24_GEOLARV_G10_2_ADJNSAF:X2017_JULY_24_GEOLARV_G7_3_ADJNSAF, nsaf.uniq[,49]>1)
#6357 proteins

#keep only proteins that are not contaminants
twopeps$protein<-row.names(twopeps)
geo.prot<-subset(twopeps, grepl(paste('TRINITY', collapse="|"), twopeps$protein))
#6328 proteins

colnames(geo.prot)<-sub("X2017_JULY_24_GEOLARV_", "", colnames(geo.prot))
colnames(geo.prot)<-sub("_ADJNSAF", "", colnames(geo.prot))

pH7.1cols<-cbind(geo.prot$G6_2,geo.prot$G6_3,geo.prot$G7_2,geo.prot$G7_3,geo.prot$G12_2,geo.prot$G12_3,geo.prot$G13_2,geo.prot$G13_3,geo.prot$G18_2,geo.prot$G18_3,geo.prot$G19_2,geo.prot$G19_3,geo.prot$G22,geo.prot$G22_2,geo.prot$G_23,geo.prot$G23_2,geo.prot$G28,geo.prot$G28_2,geo.prot$G29,geo.prot$G29_2,geo.prot$G34,geo.prot$G34_2,geo.prot$G35,geo.prot$G35_2)

pH7.5cols<-cbind(geo.prot$G4_2,geo.prot$G4_3,geo.prot$G5_2,geo.prot$G5_3,geo.prot$G10_2,geo.prot$G10_3,geo.prot$G11_2,geo.prot$G11_3,geo.prot$G16_2,geo.prot$G16_3,geo.prot$G17_2,geo.prot$G17_3,geo.prot$G20,geo.prot$G20_2,geo.prot$G_21,geo.prot$G21_2,geo.prot$G26,geo.prot$G26_2,geo.prot$G27,geo.prot$G27_2,geo.prot$G32,geo.prot$G32_2,geo.prot$G33,geo.prot$G33_2)

SumNSAFpH7.1<-rowSums(pH7.1cols)
SumNSAFpH7.5<-rowSums(pH7.5cols)

geo.file<-cbind(geo.prot, SumNSAFpH7.1, SumNSAFpH7.5)

uniprot<-read.csv('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/geo_larvae_OA_background.fasta.results.csv', header=F)
pH7.5<-read.csv('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/no pH 8.2/cluster assignments pH 7.5.csv', header=T)
pH7.1<-read.csv('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/no pH 8.2/cluster assignments pH 7.1.csv', header=T)

names(pH7.5)[names(pH7.5)=='X']<-'protein'
names(pH7.5)[names(pH7.5)=='x']<-'clusterpH7.5'

names(pH7.1)[names(pH7.1)=='X']<-'protein'
names(pH7.1)[names(pH7.1)=='x']<-'clusterpH7.1'

names(uniprot)[names(uniprot)=='V1']<-'protein'
names(uniprot)[names(uniprot)=='V2']<-'uniprotID'

uniprot.uniq<-unique(uniprot)

geo.merge1<-merge(x=geo.file, y=pH7.1, by='protein', all.x=T)
geo.merge2<-merge(x=geo.merge1, y=pH7.5, by='protein', all.x=T)
geo.merge3<-merge(x=geo.merge2, y=uniprot.uniq, by='protein', all.x=T)

write.csv(geo.merge3, 'Geoduck Larvae Summary File.csv', quote=F, row.names=F)

#NMDS
library(vegan)

geo.t<-t(geo.prot[,1:48])
geo.tra<-(geo.t+1)
geo.tra<-data.trans(geo.tra, method='log', plot=F)

nmds.techreps<-metaMDS(geo.tra, distance='bray', k=2, trymax=100, autotransform=F)
#Run 20 stress 0.1581164

ordiplot(nmds.techreps, choices=c(1,2), type='text', display='sites', cex=0.5)

#average tech reps
G35<-cbind(geo.prot$G35, geo.prot$G35_2)

G35.avg<-rowMeans(G35)

LowpH2<-data.frame(G4.avg, G5.avg, G6.avg, G7.avg, G10.avg, G11.avg, G12.avg, G13.avg,G16.avg, G17.avg, G18.avg, G19.avg, G20.avg, G21.avg, G22.avg, G23.avg, G26.avg, G27.avg, G28.avg, G29.avg, G32.avg, G33.avg, G34.avg, G35.avg)
rownames(LowpH2)<-rownames(geo.prot)

#day 6 "#C7E9B4"
#day 8 "#7FCDBB"
#day 10 "#41B6C4"
#day 12 "#1D91C0"
#day 14 "#225EA8"
#day 17 "#0C2C84"
#pH 7.5 squares 22, gray outline
#pH 7.1 triangles 24, gray outline

geo.t2<-t(LowpH2)
geo.tra2<-(geo.t2+1)
geo.tra2<-data.trans(geo.tra2, method='log', plot=F)

nmds.avg<-metaMDS(geo.tra2, distance='bray', k=2, trymax=100, autotransform=F)
#Run 20 stress 0.1160174
fig.avg<-ordiplot(nmds.avg, choices=c(1,2), type='none', display='sites')

points(fig.avg, 'sites', bg=c(rep("#C7E9B4", 4), rep("#7FCDBB",4), rep("#41B6C4",4), rep("#1D91C0",4), rep("#225EA8",4), rep("#0C2C84",4)), col='grey76', pch=c(22, 22, 24, 24, 22, 22, 24, 24, 22, 22, 24, 24, 22, 22, 24, 24, 22, 22, 24, 24, 22, 22, 24, 24), cex=1.2)
legend(x=0.042, y=0.075, legend=c('Day 6', 'Day 8', 'Day 10', 'Day 12', 'Day 14', 'Day 17', 'pH 7.5', 'pH 7.1'), pch=c(19, 19, 19, 19, 19, 19, 15, 17), col=c("#C7E9B4","#7FCDBB","#41B6C4","#1D91C0", "#225EA8", "#0C2C84",rep('grey76', 2)))

LowpH.row2<-data.stand(geo.t2, method='total', margin='row', plot=F)
LowpH.d2<-vegdist(LowpH.row2, 'bray')
pH.low2<-c(rep('7.5', 2), rep('7.1', 2),rep('7.5', 2), rep('7.1', 2),rep('7.5', 2), rep('7.1', 2),rep('7.5', 2), rep('7.1', 2), rep('7.5', 2), rep('7.1', 2), rep('7.5', 2), rep('7.1', 2))
low.anosim2<-anosim(LowpH.d2, grouping=pH.low2)
summary(low.anosim2)
#ANOSIM statistic R: 0.008207 
      #Significance: 0.338 

day.low2<-c(rep('6',4), rep('8',4), rep('10',4), rep('12',4), rep('14',4), rep('17',4))
low4.anosim<-anosim(LowpH.d2, grouping=day.low2)
summary(low4.anosim)
#ANOSIM statistic R: 0.6579 
      #Significance: 0.001

LowpH.eigen2<-envfit(nmds.avg$points, geo.tra2, perm=1000)