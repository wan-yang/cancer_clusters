# Joint point anlaysis
# to check the trends in colorectal, kidney, and thyroid
# Plot eFigure 7


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
cans = c("colon.and.rectum","kidney.and.renal.pelvis","thyroid")
can.names=c('Colorectal','Kidney','Thyroid')
names(can.names)=cans

pdf(paste0(dir_proj,'FigS_joinpoint.pdf'),width = 6,height=7.5)
par(mfcol=c(length(cans),2),xaxt='n',xpd=T,mar=c(0,2.5,0,.5),oma=c(2,.5,.8,.5),mgp=c(1.2,.2,0),tck=-.02,cex=.9,cex.lab=.9,cex.axis=.85)
for (sex in c('Men','Women')){
  for(can in cans){
    tda=get(paste0('da.',sex))
    tda=tda[,c('year',can)] %>% data.table() #
    colnames(tda)=c('Year','StdInciRate')
    tda$LogRate = log(tda$StdInciRate) # log it first!
    # plot(tda$Year,tda$StdInciRate,pch=20)
    # test diff models and select based on BIC?
    psi1=median(years); 
    psi2=quantile(years,.33);
    psi3=quantile(years,.67);
    psi4=quantile(years,.25);
    psi5=quantile(years,.75);
    psi6=quantile(years,c(.33,.67)); 
    psi7=quantile(years,c(.25,.5)); 
    psi8=quantile(years,c(.5,.75)); 
    psi9=quantile(years,c(.17,.33));
    psi10=quantile(years,c(.67,.83));
    psi11=quantile(years,c(.22,.45,.67)); 
    psi12=quantile(years,c(.33,.55,.78)); 
    psi13=quantile(years,c(.25,.5,.75)); 
    BICs=numeric(13)
    out<-lm(LogRate~Year, data=tda) 
    for(mi in 1:length(BICs)){
      set.seed(10)
      this.psi=get(paste0('psi',mi))
      out.seg<- segmented(out,seg.Z=~Year,psi=this.psi)  # , psi=c(1985,2005)
      BICs[mi]=BIC(out.seg)
      # print( your_breakpoint <- round( as.vector( out.seg$psi[, "Est." ] ) ) )
    }
    BICs
    idx.best=which.min(BICs)
    this.psi=get(paste0('psi',idx.best))
    # out<-lm(LogRate~Year, data=tda) 
    out.seg<-segmented(out,seg.Z=~Year,psi=this.psi)  # , psi=c(1985,2005)
    
    ## the slopes of the three segments....
    slope(out.seg)
    summary(out.seg)
    # figuring out the breakpoint year was the purpose of this joinpoint analysis.
    ( your_breakpoint <- round( as.vector( out.seg$psi[, "Est." ] ) ) )
    # obtain the annual percent change (APC=) estimates for each time point
    APCs=as.matrix(round(slope(out.seg, APC = T)[[1]],1))
    
    plot(out.seg,conf.level=.95,shade=T,lwd=2,ylab='Log Incidence Rate (per 100,000)')
    points(tda$Year,tda$LogRate,pch='x', cex=.6)
    # points(out.seg, link=FALSE, col=2, cex=.6)
    abline(v=your_breakpoint, lty=2, col='grey50')
    x=c(min(years),your_breakpoint,max(years))
    for(ip in 1:nrow(APCs)){
      mtext(paste0(round(x[ip],0),'-',round(x[ip+1]-1,0),': APC=',APCs[ip,1],' (',APCs[ip,2],', ',APCs[ip,3],')'),
            line=-1.2 - ip *.92, side=3,adj =.05, cex=.75
      )
    }
    mtext(paste0(sex,': ',can.names[can]),side = 3,line=-1.1,cex=.8,adj=.05)
  } # cancer
  # axis(1,xpd=NA)
  xx=seq(1975,2015,by=10)
  # text(xx,y=rep(min(tda$LogRate),length(xx)),pos=1,labels=xx,cex=.75,xpd=NA)
  if (sex=='Men')  text(xx,y=rep(.608,length(xx)),pos=1,labels=xx,cex=.75,xpd=NA)
  if (sex=='Women') text(xx,y=rep(2.07,length(xx)),pos=1,labels=xx,cex=.75,xpd=NA)
  mtext('Year',side=1,line=.6,cex=.8,adj=.5)
} # sex
dev.off()

