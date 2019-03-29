#Cluster plots for pH 7.5
pH7.5<-read.csv('NSAF pH 7.5.csv', header=T, row.names=1)

Day6.7.5<-cbind(pH7.5$G4.avg, pH7.5$G5.avg)
Day8.7.5<-cbind(pH7.5$G10.avg, pH7.5$G11.avg)
Day10.7.5<-cbind(pH7.5$G16.avg, pH7.5$G17.avg)
Day12.7.5<-cbind(pH7.5$G20.avg, pH7.5$G21.avg)
Day14.7.5<-cbind(pH7.5$G26.avg, pH7.5$G27.avg)
Day17.7.5<-cbind(pH7.5$G32.avg, pH7.5$G33.avg)

D6pH7.5<-rowMeans(Day6.7.5)
D8pH7.5<-rowMeans(Day8.7.5)
D10pH7.5<-rowMeans(Day10.7.5)
D12pH7.5<-rowMeans(Day12.7.5)
D14pH7.5<-rowMeans(Day14.7.5)
D17pH7.5<-rowMeans(Day17.7.5)

pH7.5.df<-data.frame(D6pH7.5, D8pH7.5, D10pH7.5,D12pH7.5, D14pH7.5, D17pH7.5)
rownames(pH7.5.df)<-rownames(pH7.5)
names(pH7.5.df)[names(pH7.5.df)=='D17pH7.5']<-'17'

pH7.5.df$sum<-rowSums(pH7.5.df)
pH7.5.df<-subset(pH7.5.df, sum>0, select=-sum)
pH7.5.bray<-vegdist(pH7.5.df, method='bray')
pH7.5.clust<-hclust(pH7.5.bray, method='average')
plot(pH7.5.clust, labels=F)
rect.hclust(pH7.5.clust, h=0.5)

prot.clust.7.5<-read.csv('cluster assignments pH 7.5.csv', header=T, row.names=1)
clust.nsaf.7.5<-merge(x=prot.clust.7.5, y=pH7.5.df, by='row.names', all.x=T)

library(ggthemes)
library(reshape)
library(ggplot2)

melt7.5<-melt(clust.nsaf.7.5, id.vars=c('Row.names', 'x'))
ggplot(melt7.5, aes(x=variable, y=value, group=Row.names)) +geom_line(alpha=0.1) +theme_bw() +facet_wrap(~x, scales='free_y') +labs(x='Day', y='Averaged Normalized Spectral Abundance Factor') +theme(axis.text.x = element_text(angle=90,size=6)) +geom_smooth(aes(group=1),method="loess", se=FALSE, span=0.6)

#cluster plots for pH 7.1
pH7.1<-read.csv('NSAF pH 7.1.csv', header=T, row.names=1)

Day6.7.1<-cbind(pH7.1$G6.avg, pH7.1$G7.avg)
Day8.7.1<-cbind(pH7.1$G12.avg, pH7.1$G13.avg)
Day10.7.1<-cbind(pH7.1$G18.avg, pH7.1$G19.avg)
Day12.7.1<-cbind(pH7.1$G22.avg, pH7.1$G23.avg)
Day14.7.1<-cbind(pH7.1$G28.avg, pH7.1$G29.avg)
Day17.7.1<-cbind(pH7.1$G34.avg, pH7.1$G35.avg)

D6pH7.1<-rowMeans(Day6.7.1)
D8pH7.1<-rowMeans(Day8.7.1)
D10pH7.1<-rowMeans(Day10.7.1)
D12pH7.1<-rowMeans(Day12.7.1)
D14pH7.1<-rowMeans(Day14.7.1)
D17pH7.1<-rowMeans(Day17.7.1)

pH7.1.df<-data.frame(D6pH7.1, D8pH7.1, D10pH7.1,D12pH7.1, D14pH7.1, D17pH7.1)
rownames(pH7.1.df)<-rownames(pH7.1)
names(pH7.1.df)[names(pH7.1.df)=='D17pH7.1']<-'17'

pH7.1.df$sum<-rowSums(pH7.1.df)
pH7.1.df<-subset(pH7.1.df, sum>0, select=-sum)


prot.clust.7.1<-read.csv('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/no pH 8.2/cluster assignments pH 7.1.csv', header=T, row.names=1)
clust.nsaf.7.1<-merge(x=prot.clust.7.1, y=pH7.1.df, by='row.names', all.x=T)

melt7.1<-melt(clust.nsaf.7.1, id.vars=c('Row.names', 'x'))
ggplot(melt7.1, aes(x=variable, y=value, group=Row.names)) +geom_line(alpha=0.1) +theme_bw() +facet_wrap(~x, scales='free_y') +labs(x='Day', y='Averaged Normalized Spectral Abundance Factor') +theme(axis.text.x = element_text(angle=90,size=6)) +geom_smooth(aes(group=1),method="loess", se=FALSE, span=0.6)