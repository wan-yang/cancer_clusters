# PLOT eFigure 8

dir_home = '~/Documents/MS/Cancer_clusters/Epidemiology/model_code/' # CHANGE TO YOUR OWN LOCAL DIRECTORY

dir_data=paste0(dir_home,'data/')
dir_code=paste0(dir_home,'scripts/')
dir_proj=paste0(dir_home,'proj_cancer_clusters/')

if(!file.exists(dir_proj)) dir.create(dir_proj)

library(UScancer); library(stringr); library(xlsx); 
library(ggplot2); library(beepr);
library(data.table); library(magrittr)
library(segmented); library(RColorBrewer)
library(cluster); library(dendextend); 

source(paste0(dir_code,'Fn_cancer_extract_subtree.R'))
source(paste0(dir_code,'Fn_formatting.R'))

years=1973:2015; 
ages= 25:39;
sexes=c('Men','Women')
sex=1
da.Men=read.csv(paste0(dir_data,'table_by.sites_stdInciRate_a',min(ages),'-',max(ages),'_',sexes[sex],'.csv'))
sex=2
da.Women=read.csv(paste0(dir_data,'table_by.sites_stdInciRate_a',min(ages),'-',max(ages),'_',sexes[sex],'.csv'))

siteTable=read.csv(paste0(dir_data,'siteTable.csv'),stringsAsFactors = F) %>% data.table()
cols=c('pink','cornflowerblue','brown','black','darkgrey','orange','red','darkgreen','magenta')
names(cols)=unique(siteTable$system)

key.cans=c('lung.and.bronchus',"larynx","urinary.bladder")
inci.w=da.Women[,key.cans]
inci.m=da.Men[,key.cans]
inci.w=inci.w[,names(sort(inci.w[nrow(inci.w),]))]
inci.m=inci.m[,names(sort(inci.m[nrow(inci.m),]))]
can.cor=NULL;
for(a in names(inci.w)){
  can.cor=append(can.cor,cols[siteTable[site==a]$system])
}
ymax=max(inci.w,inci.m); ymin=min(inci.w,inci.m)

pdf(paste0(dir_proj,'figS_cp.can.trends.cluster2inMen_','_a',min(ages),'-',max(ages),'.pdf'),width = 6,height=2.5)
par(mfcol=c(1,2),mar=c(0,2.5,0,4),oma=c(2.5,.5,1.5,.5),mgp=c(1.3,.3,0),tck=-.02,cex=.8,cex.lab=.9,cex.axis=.9)
matplot(x=years,y=inci.m, ylim = c(ymin,ymax), ylab='',xlab='',col=can.cor,type='l',lty=1,log='y',xaxt='n',lwd=1,xpd=T)
matpoints(years, inci.m,col=can.cor,pch=1:ncol(inci.m),cex=.5,lwd=1,log='y',xpd=T)
ys=as.numeric(unlist(inci.m[nrow(inci.m),]))
y.diff=diff(ys)
y.diff[which(y.diff<ymax/10)]=ymax/10
ys=ys[1]+c(0,cumsum(y.diff))
points(rep(max(years)+2.5,ncol(inci.m)),ys,col=can.cor,pch=1:ncol(inci.m),cex=.7,log='y',lwd=1,xaxt='n',xpd=T)
text(rep(max(years)+.5 + 1.3,ncol(inci.m)),ys,fn_format.site(gsub("\\.",' ',colnames(inci.m))),pos=4,col=can.cor,xpd=T)
axis(1)
mtext('Year',outer=F,side=1,line=1.2,cex=.8,adj=.5)
mtext('Log Incidence Rate (per 100,000)',outer = F,side=2,line=1.3,cex=.8,adj=.5,xpd=T)
mtext('(A) Men',outer = F,side=3,line = .1,cex = .8,adj = .02)
# for women
matplot(years,inci.w,ylim = c(.02,5),ylab='',xlab='',col=can.cor,type='l',lty=1,log='y',lwd=1,xaxt='n',xpd=T)
matpoints(years,inci.w,col=can.cor,pch=1:ncol(inci.w),cex=.5,log='y',lwd=1,xaxt='n',xpd=T)
ys=as.numeric(unlist(inci.w[nrow(inci.w),]))
y.diff=diff(ys)
y.diff[which(y.diff<ymax/10)]=ymax/10
ys=ys[1]+c(0,cumsum(y.diff))
points(rep(max(years)+2.5,ncol(inci.w)),ys,col=can.cor,pch=1:ncol(inci.w),cex=.7,log='y',lwd=1,xaxt='n',xpd=T)
text(rep(max(years)+.5 + 1.3,ncol(inci.w)),ys,fn_format.site(gsub("\\.",' ',colnames(inci.w))),pos=4,col=can.cor,xpd=T)
axis(1)
mtext('Year',outer=F,side=1,line=1.2,cex=.8,adj=.5)
mtext('Log Incidence Rate (per 100,000)',outer = F,side=2,line=1.2,cex=.8,adj=.5,xpd=NA)
mtext('(B) Women',outer = F,side=3,line = .1,cex = .8,adj = .02)
dev.off()
