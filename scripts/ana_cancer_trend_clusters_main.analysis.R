# to analyze cancer trends and clusters/systematic relation among cancer sites
# Main analysis
# 5/28/19 

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
ages= 25:39 # main study

sexes=c('Men','Women')

sys.names=c('Breast','Digestive','Digestive','Female Genital',
            'Blood','Male Genital','Respiratory','Urinary','All Other Sites')
names(sys.names)=tolower(c('BREAST','COLRECT','DIGOTHR','FEMGEN','LYMYLEUK','MALEGEN','RESPIR','URINARY','OTHER'))

siteTable=read.csv(paste0(dir_data,'siteTable.csv'),stringsAsFactors = F) %>% data.table()


sexes=c('Men','Women')

# Compute the corr
for(sex in 1:2){
  TDA=read.csv(paste0(dir_data,'table_by.sites_stdInciRate_a',min(ages),'-',max(ages),'_',sexes[sex],'.csv'))
  tda=TDA[,-1]
  # exclude cance with 0 incidence
  idx=which(colSums(tda)==0)
  if(length(idx)>0) tda=tda[,-idx]
  corr=matrix(0,ncol(tda),ncol(tda))
  for(i in 1:(ncol(tda)-1)){
    for(j in (i+1):ncol(tda)){
      corr[i,j]=corr[j,i]=cor(tda[,i],tda[,j]) # 1-cor(tda[,i],tda[,j]) # distance
    }
  }
  colnames(corr)=rownames(corr)=colnames(tda)
  write.csv(corr,paste0(dir_proj,'table_by.sites_corr_a',min(ages),'-',max(ages),'_',sexes[sex],'.csv'),row.names = T)
  
}

sex=1  # 1 (for men) and 2 (for women)
da.inci=read.csv(paste0(dir_data,'table_by.sites_stdInciRate_a',min(ages),'-',max(ages),'_',sexes[sex],'.csv'))
corr=read.csv(paste0(dir_proj,'table_by.sites_corr_a',min(ages),'-',max(ages),'_',sexes[sex],'.csv'),row.names = 1)
# exclude miscellaneous
corr=corr[,!grepl('miscellaneous',colnames(corr))]
corr=corr[!grepl('miscellaneous',rownames(corr)),]
distance=1-corr
d.corr=as.dist(distance)
method='average' # 'complete' # 'single';#  'median'; # ; # 'centroid'; # 'average'
hc.corr=hclust(d.corr,method=method)

cut.off=.4
hclust.trees=f.get.subtrees(hc.corr,h=cut.off)
tree=hclust.trees;
tr.size=function(tree){length(tree$height)}
tr.height=function(tree){mean(tree$height)} # rank trees based on similarity (level of confidence)
sizes=sapply(tree,tr.size)
bigtrees=which(sizes>1)
bigtrees=bigtrees[order(sapply(tree,tr.height)[bigtrees])]
wholetree=f.get.subtrees(hc.corr,h=3)

# Plot Figures 1 (for men) and 2 (for women)
# with diff symbols
pdf(paste0(dir_proj,'fig.diff.symbo_cluster_',method,'_cut',cut.off,'_',sexes[sex],'_a',min(ages),'-',max(ages),'.pdf'),width = 6,height=length(bigtrees)*ifelse(length(bigtrees)>2,1.5,ifelse(length(bigtrees)==2,1.8,3)))
par(mfcol=c(length(bigtrees),2),lheight=.7,xpd=T,mar=c(0,0,0,4.5),oma=c(2.5,2.5,1.5,.5),mgp=c(1.3,.3,0),tck=-.02,cex=.8,cex.lab=.9,cex.axis=.9)
cc=palette(brewer.pal(n=11,name='Spectral')); # colfunc=colorRampPalette(cc);cols=colfunc(length(unique(siteTable$system)))
cols=c('pink',palette(brewer.pal(n=8,name='Dark2'))); 
cols=c('pink','cornflowerblue','brown','black','darkgrey','orange','red','darkgreen','magenta')
names(cols)=unique(siteTable$system)
cols1=c('cyan2','cyan3'); cols2=c('grey30','grey50');  cnt=cnt2=1
for(tr in 1:length(bigtrees)){ # plot incidence
  tree1=tree[[bigtrees[tr]]]
  inci=da.inci[,tree1$labels]
  inci=inci[,names(sort(inci[nrow(inci),]))]
  tree1.cor=NULL;
  for(a in names(inci)){
    tree1.cor=append(tree1.cor,cols[siteTable[site==a]$system])
  }
  ymax=max(inci)
  matplot(years,inci,ylab='',col=tree1.cor,type='l',lty=1,pch=1:ncol(inci),log='y',lwd=1,xaxt='n',xpd=T)
  matpoints(years,inci,col=tree1.cor,pch=1:ncol(inci),cex=.5,log='y',lwd=1,xaxt='n',xpd=T)
  ys=as.numeric(unlist(inci[nrow(inci),]))
  y.diff=diff(ys)
  y.diff[which(y.diff<ymax/10)]=ymax/10
  ys=ys[1]+c(0,cumsum(y.diff))
  if(sex==2 & cut.off==.4) ys=c(.25,.5,1,1.5,2.8,4,6,8.5,25)
  if(sex==2 & cut.off==.35) ys=c(.5,1.3,2.8,4,6,8.5,25)
  if(sex==1 &  tr==1 &  cut.off==.4) ys=c(.4,2.5,3.5,5,6.7)
  if(sex==1 &  tr==3 &  cut.off==.4) ys=c(.3,.8,12)
  points(rep(max(years)+2.5,ncol(inci)),ys,col=tree1.cor,pch=1:ncol(inci),cex=.7,log='y',lwd=1,xaxt='n',xpd=T)
  text(rep(max(years)+.5 + 1.3,ncol(inci)),ys,fn_format.site(gsub("\\.",' ',colnames(inci))),pos=4,col=tree1.cor,xpd=T)
}
axis(1)
mtext('Year',outer=F,side=1,line=1.2,cex=.8,adj=.5)
mtext('Log Incidence Rate (per 100,000)',outer = T,side=2,line=1.3,cex=.8,adj=.5,xpd=T)
mtext('(A) Incidence rate over time',side=3,outer=T,line=.1,adj=0,cex=.8); 
for(tr in 1:length(bigtrees)){
  tree1=tree[[bigtrees[tr]]]
  inci=da.inci[,tree1$labels[tree1$order]]
  
  ori.order=tree1$order
  # tree1.cor=cols[siteTable[site %in% tree1$labels[tree1$order]]$system] # messed up b/c sorting
  tree1.cor=NULL;
  for(a in tree1$labels[tree1$order]){
    tree1.cor=append(tree1.cor,cols[siteTable[site==a]$system])
  }
  ori.labels=gsub("\\.",' ',tree1$labels)
  ori.labels=fn_format.site(ori.labels)
  tree1$labels=ori.labels # gsub("\\.",' ',tree1$labels) # format the cancer names
  tree1=as.dendrogram(tree1); # somehow the order is messed up
  labels_colors(tree1)=tree1.cor # reg.colors[tree1.cor]; # set the color
  labels(tree1)=ori.labels[ori.order];  # reset the labels in order
  
  plot(tree1,horiz = T,edgePar = list(lwd=.5,col=cols2[tr%%2+1]),xlim=c(cut.off+.05,0),axes=F,xpd=T)
  # legend('left',adj=c(.8,.5),legend = paste0('C',tr),lty=NULL,cex=1,bty='n')
}
axis(1,at=seq(cut.off,0,by=-.05),labels = 1-seq(cut.off,0,by=-.05),cex.axis=.9)
mtext('Correlation',outer=F,side=1,line=1.2,cex=.8,adj=.5)
mtext('(B) Hierarchical clustering',side=3,outer=T,line=.1,adj=.75,cex=.8); 
dev.off()

# Plot eFigures 1 (for men) and 2 (for women)
pdf(paste0(dir_proj,'fig_cluster.all_',method,'_cut',cut.off,'_',sexes[sex],'_a',min(ages),'-',max(ages),'.pdf'),width = 6,height=7.5)
par(mfcol=c(1,1),lheight=.7,xpd=T,mar=c(0,0,0,17),oma=c(2.5,.5,.3,.3),mgp=c(1.5,.3,0),tck=-.02,cex=.8,cex.lab=.9,cex.axis=.9)
cc=palette(brewer.pal(n=11,name='Spectral')); # colfunc=colorRampPalette(cc);cols=colfunc(length(unique(siteTable$system)))
cols=c('pink',palette(brewer.pal(n=8,name='Dark2'))); 
cols=c('pink','cornflowerblue','brown','black','darkgrey','orange','red','darkgreen','magenta')
names(cols)=unique(siteTable$system)
cols1=c('cyan2','cyan3'); cols2=c('grey30','grey50');  cnt=cnt2=1
{
  tree1=wholetree[[1]]
  inci=da.inci[,tree1$labels[tree1$order]]
  
  ori.order=tree1$order
  # tree1.cor=cols[siteTable[site %in% tree1$labels[tree1$order]]$system] # messed up b/c sorting
  tree1.cor=NULL;
  for(a in tree1$labels[tree1$order]){
    tree1.cor=append(tree1.cor,cols[siteTable[site==a]$system])
  }
  ori.labels=gsub("\\.",' ',tree1$labels)
  ori.labels=fn_format.site2(ori.labels)
  tree1$labels=ori.labels # gsub("\\.",' ',tree1$labels) # format the cancer names
  tree1=as.dendrogram(tree1); # somehow the order is messed up
  labels_colors(tree1)=tree1.cor # reg.colors[tree1.cor]; # set the color
  labels(tree1)=ori.labels[ori.order];  # reset the labels in order
  
  plot(tree1,horiz = T,edgePar = list(lwd=.5,col=cols2[tr%%2+1]),
       # xlim=c(cut.off+.05,0),
       axes=F,xpd=T)
  # legend('left',adj=c(.8,.5),legend = paste0('C',tr),lty=NULL,cex=1,bty='n')
}
axis(1,at=seq(1,0,by=-.1),labels = 1-seq(1,0,by=-.1),cex.axis=.9)
abline(v=cut.off,lwd=.5,lty=5)
mtext('Correlation',outer=F,side=1,line=1.2,cex=.8,adj=.5)
mtext('Hierarchical clustering',side=3,outer=T,line=-1.5,adj=0,cex=.8); 
dev.off()

