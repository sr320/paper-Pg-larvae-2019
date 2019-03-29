setwd('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/no pH 8.2')

geo.size<-read.csv('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/size.csv', header=T)

size.lowpH<-subset(geo.size, Treatment<8.2)

size.early<-subset(size.lowpH, Age.Days<12)
size.late<-subset(size.lowpH, Age.Days>10)

size.con1.1<-subset(size.early, Conical==1)
size.con3.1<-subset(size.early, Conical==3)
size.con5.1<-subset(size.early, Conical==5)
size.con6.1<-subset(size.early, Conical==6)
size.con7.1<-subset(size.early, Conical==7)
size.con8.1<-subset(size.early, Conical==8)
size.con1.2<-subset(size.late, Conical==1)
size.con3.2<-subset(size.late, Conical==3)
size.con5.2<-subset(size.late, Conical==5)
size.con6.2<-subset(size.late, Conical==6)
size.con7.2<-subset(size.late, Conical==7)
size.con8.2<-subset(size.late, Conical==8)

plot(x=size.con5.1$Age.Days, y=size.con5.1$Average.Size, ylim=c(100,300), ylab='Average Size (µm)', xlab='Age (Days)', type='l', col='orange')
lines(x=size.con6.1$Age.Days, y=size.con6.1$Average.Size, col='orange')
lines(x=size.con7.1$Age.Days, y=size.con7.1$Average.Size, col='dodgerblue')
lines(x=size.con8.1$Age.Days, y=size.con8.1$Average.Size, col='dodgerblue')
legend(x=6, y=300, legend=c('pH 7.5', 'pH 7.1'), lty=c(1, 1), col=c('dodgerblue', 'orange'))

plot(x=size.con5.2$Age.Days, y=size.con5.2$Average.Size, ylim=c(100,300), ylab='', xlab='Age (Days)', type='l', col='orange')
lines(x=size.con6.2$Age.Days, y=size.con6.2$Average.Size, col='orange')
lines(x=size.con7.2$Age.Days, y=size.con7.2$Average.Size, col='dodgerblue')
lines(x=size.con8.2$Age.Days, y=size.con8.2$Average.Size, col='dodgerblue')


size.con5<-subset(size.lowpH, Conical==5)
size.con6<-subset(size.lowpH, Conical==6)
size.con7<-subset(size.lowpH, Conical==7)
size.con8<-subset(size.lowpH, Conical==8)

plot(x=size.con5$Age.Days, y=size.con5$Average.Size, ylim=c(100,300), ylab='Average Size (µm)', xlab='Age (Days)', type='l', col='orange')
lines(x=size.con6$Age.Days, y=size.con6$Average.Size, col='orange')
lines(x=size.con7$Age.Days, y=size.con7$Average.Size, col='dodgerblue')
lines(x=size.con8$Age.Days, y=size.con8$Average.Size, col='dodgerblue')
legend(x=6, y=300, legend=c('pH 7.5', 'pH 7.1'), lty=c(1, 1), col=c('dodgerblue', 'orange'))

#Chi square
#numbers in matrix are average numbers of larvae captured on each screen size
geoday10=matrix(c(0.222E6, 0.278E6, 0.555E6, 0.356E6, 1.15E6, 1.1E6, 0.434E6, 0.489E6), nrow=2)
colnames(geoday10)=c(80, 90, 100, 118)
rownames(geoday10)=c(7.5, 7.1)

#divide each entry by sum of entries in the column (i.e. sum for each screen size)
prop.table(geoday10, 2)

barplot(prop.table(geoday10, 2), beside=T, legend.text=T, ylim=c(0,1), ylab='Proportions')

numb.av<-read.csv('~/Documents/genome_sciences_postdoc/geoduck larvae/Emma analysis/no pH 8.2/larval numbers averaged.csv', header=T)
av.d10<-subset(numb.av, Day==10, select=-Day)
d10.mat<-as.matrix(av.d10[,2:3])
rownames(d10.mat)<-av.d10$Min.Size
source('~/Documents/Multivariate Stats/biostats.R', chdir = TRUE)
d10.std<-data.stand(d10.mat, method='total')
d10.std.mat<-as.matrix(d10.std)
barplot(d10.std.mat, legend=rownames(d10.std.mat), xlim=c(0,ncol(d10.std.mat)+1), args.legend=list(x=ncol(d10.std.mat)+1, y=max(colSums(d10.std.mat)), bty='n'))

av.d17<-subset(numb.av, Day==17, select=-Day)
d17.std<-data.stand(av.d17, method='total')
d17.mat<-as.matrix(d17.std[,2:3])
rownames(d17.mat)<-av.d17$Min.Size
barplot(d17.mat, legend=rownames(d17.mat), xlim=c(0,ncol(d17.mat)+2), args.legend=list(x=ncol(d17.mat)+2, y=max(colSums(d17.mat)), bty='n'))

#sum rows and columns
x=addmargins(geoday10)

 80     90     100    118     Sum
7.5 222000 555000 1150000 434000 2361000
7.1 278000 356000 1100000 489000 2223000
Sum 500000 911000 2250000 923000 4584000

p.pH7.5<-2.361E6/4.58E6
p.pH7.1<-2.223E6/4.58E6

#expected distributions
screen80<-5E5*c(p.pH7.5, p.pH7.1)
screen90<-9.11E5*c(p.pH7.5, p.pH7.1)
screen100<-2.25E6*c(p.pH7.5, p.pH7.1)
screen118<-9.23E5*c(p.pH7.5, p.pH7.1)

Expected=cbind(screen80, screen90, screen100, screen118)
rownames(Expected)=c(7.5, 7.1)

chi.sq.d10<-sum((geoday10-Expected)^2/Expected)
#0.09898174
#degrees of freedom = (2-1)(4-1) = 3
df=3
1-pchisq(chi.sq.d10,df)

res.day10<-chisq.test(geoday10)
Pearsons Chi-squared test

data:  geoday10
X-squared = 50021, df = 3, p-value < 2.2e-16

#day 17
geoday17=matrix(c(0, 0.07E6, 0, 0, 0, 0.1545E6, 0, 0, 0, 0.3915E6,0, 0,  0.022E6, 0.2445E6, 0.111E6, 0.367E6, 0.0077E6, 0, 0.054E6, 0.162E6, 0, 0, 0.033E6, 0.033E6), nrow=4)
colnames(geoday17)=c(140, 150, 160, 180, 200, 225)
rownames(geoday17)=c(7.5, 7.1, '7.5to8.2', '7.1to8.2')

prop.table(geoday17, 2)
barplot(prop.table(geoday17,2), beside=T, legend.text=T, ylim=c(0,1), ylab='Proportions')

res.day17<-chisq.test(geoday17)
Pearsons Chi-squared test

data:  geoday17
X-squared = 1044200, df = 15, p-value < 2.2e-16

#day 19
geoday19<-matrix(c(0,0.145E6, 0, 0,0, 0.467E6, 0, 0, 0, 0.756E6,0.089E6, 0, 0, 0.277E6,  0.0462E6, 0.146E6,  0.033E6, 0.08E6,  0.0267E6, 0.113E6, 0.056E6, 0, 0.1E6, 0.244E6, 0, 0, 0.09E6, 0.075E6), nrow=4)

res.day19<-chisq.test(geoday19)
Pearsons Chi-squared test

data:  geoday19
X-squared = 2012700, df = 18, p-value < 2.2e-16

colnames(geoday19)=c(150, 160, 180, 200, 225, 236, 250)
rownames(geoday19)=c(7.5, 7.1, '7.5to8.2', '7.1to8.2')
prop.table(geoday19, 2)