######
#created Sep 2013, M. Beck
#creates interpolation grids for TB hirsch model
#used to create heat maps and to get normalization data

######
#get interpolation grids

######
#salinity grid
#includes beta estimates for weighted regressions for each obs
#includes back-transformation from Moyer et al. 2012

rm(list=ls())
source('M:/r_code/EPC/epc_dat.r')

sal.div<-20 #no. of divisions, range is different for each segment

segs<-unique(tb.dat$seg)

mods.out <- NULL

strt<-Sys.time()

for(seg in segs){
  
  dat.in<-tb.dat[tb.dat$seg==seg,]
  
  #salinity grid data specific to each segment
  sal.grid<-seq(min(dat.in$sal.ref),max(dat.in$sal.ref),length=sal.div)

  seg.out<-NULL
  
  for(row in 1:nrow(dat.in)){
    
    row.out<-NULL
    
    ref.in<-dat.in[row,]
    
    # log
    cat(as.character(seg), nrow(dat.in) - row,'\n')
    
    for(sal in sal.grid){
      
      ref.in$sal.ref<-sal
      ref.wts<-wt.fun(ref.in,dat.in,wt.vars=c('month.num','year','sal.ref'))
      
      # data to predict
      pred.dat <- data.frame(sal.ref=sal,dec.time=ref.in$dec.time)
      
      # crq model, estimates all quants
      mod <- quantreg::crq(
        Surv(Chla_ugl, not_cens, type = "left") ~ 
          dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time), 
        weights = ref.wts,
        data = dat.in, 
        method = "Portnoy"
        )
      
      # sometimes crq fucks up
      test <- try({coef(mod)})
      if('try-error' %in% class(test)){
        err_out <- rep(NA, 9)
        row.out<-rbind(row.out, c(err_out))
        next
      }
        
      # fitted coefficients for each model
      parms <- coef(mod, c(0.1, 0.5, 0.9))
      
      # predicted values by quantile model coefficients
      fits <- sapply(1:3, function(x){
        with(pred.dat, 
          parms[1, x] + parms[2, x] * dec.time + parms[3, x] * sal.ref + parms[4, x] * sin(2*pi*dec.time) + parms[5, x] * cos(2*pi*dec.time)
        )
      })
      names(fits) <- paste0('fit.', c('lo', 'md', 'hi'))
      
      # back transformed predicted values
      bt.fits <- exp(fits)
      names(bt.fits) <- paste0('bt.', c('lo', 'md', 'hi'))
      
      # model parameters for sal.ref
      betas <- coef(mod, c(0.1, 0.5, 0.9))['sal.ref', ]
      names(betas) <- paste0('b.', c('lo', 'md', 'hi'))
      
      #append to row out for each unique sal
      row.out<-rbind(
        row.out,
        c(fits, bt.fits, betas)
        )
      
      }
    
    wt.fits<-suppressWarnings(data.frame(
      year=ref.in$year,
      month.num=ref.in$month.num,
      date.f=ref.in$date.f,
      dec.time=ref.in$dec.time,
      seg,
      sal.grid,
      row.out
      ))
    
    seg.out<-rbind(seg.out,wt.fits)
    
    }
  
  mods.out <- rbind(mods.out, seg.out)
  
  }
  
Sys.time() - strt

sal.grd <- mods.out

save(sal.grd,file='M:/wq_models/EPC/interp_grids/salwt_grd.RData')

#####
#month and yr. grids have to be redone as of Nov. 5th for back-transformation estimates....

#####
#month grid

rm(list=ls())
source('M:/r_code/EPC/epc_dat.r')

#get interpolation grids
mo.grid<-sort(unique(tb.dat$month.num))

mods.out<-NULL

segs<-unique(tb.dat$seg)

strt<-Sys.time()

for(seg in segs){
  
  dat.in<-tb.dat[tb.dat$seg==seg,]
  
  seg.out<-NULL
  
  for(row in 1:nrow(dat.in)){
    
    wt.fit.md<-NULL
    wt.fit.hi<-NULL
    wt.fit.lo<-NULL

    ref.in<-dat.in[row,]
    
    cat(as.character(seg), nrow(dat.in) - row,'\n')
    flush.console()
    
    for(mo in mo.grid){
      
      #reset month and dec.time to eval month
      ref.in$month.num<-mo
      ref.in$dec.time<-as.numeric(ref.in$year) + mo
      ref.wts<-wt.fun(ref.in,dat.in,wt.vars=c('month.num','year','sal.ref'))
      
      #data to use for model predictions
      pred.dat<-data.frame(sal.ref=ref.in$sal.ref,dec.time=ref.in$dec.time)

      #OLS wtd model
      wt.mod.md<-lm(
        Chla_ugl~dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time),
        weights=ref.wts,
        data=dat.in
        )
      
      #OLS wtd predict
      wt.pred.md<-predict(wt.mod.md,newdata=pred.dat)
      
      if(class(try({
        
        wt.mod.hi<-rq(
          Chla_ugl~dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time),
          weights=ref.wts,
          tau=0.9,
          data=dat.in
          )
        
        wt.pred.hi<-predict(wt.mod.hi,newdata=pred.dat)
        
        }))=='try-error'){ wt.pred.hi<-NA }
      
      #quantile wtd model, 0.1, exception for error
      if(class(try({
        
        wt.mod.lo<-rq(
          Chla_ugl~dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time),
          weights=ref.wts,
          tau=0.1,
          data=dat.in
          )
        
        wt.pred.lo<-predict(wt.mod.lo,newdata=pred.dat)
        
        }))=='try-error'){ wt.pred.lo<-NA }
      
      wt.fit.md<-c(wt.fit.md,wt.pred.md)
      wt.fit.hi<-c(wt.fit.hi,wt.pred.hi)
      wt.fit.lo<-c(wt.fit.lo,wt.pred.lo)      
      
      }
    
    wt.fits<-data.frame(
      year=ref.in$year,
      month.num=dat.in[row,]$month.num,
      date.f=ref.in$date.f,
      dec.time=dat.in[row,]$dec.time,
      seg,
      mo.grid,
      wt.fit.md,
      wt.fit.hi,
      wt.fit.lo
      )
    
    seg.out<-rbind(seg.out,wt.fits)
    
    }
  
  mods.out<-rbind(mods.out,seg.out)
  
  }

cat(Sys.time()-strt,'\n')

mo.grd<-mods.out

save(mo.grd,file='M:/wq_models/EPC/interp_grids/mowt_grd.RData')

######
#annual grid

rm(list=ls())
source('M:/r_code/EPC/epc_dat.r')

#get interpolation grids
yr.grid<-sort(unique(tb.dat$year))

mods.out<-NULL

segs<-unique(tb.dat$seg)

strt<-Sys.time()

for(seg in segs){
  
  dat.in<-tb.dat[tb.dat$seg==seg,]
  
  seg.out<-NULL
  
  for(row in 1:nrow(dat.in)){
    
    wt.fit.md<-NULL
    wt.fit.hi<-NULL
    wt.fit.lo<-NULL
    
    ref.in<-dat.in[row,]
    
    cat(as.character(seg), nrow(dat.in) - row,'\n')
    flush.console()
    
    for(yr in yr.grid){
      
      #reset year and dec.time to eval year
      ref.in$year<-yr
      ref.in$dec.time<-as.numeric(ref.in$year) + ref.in$month.num
      ref.wts<-wt.fun(ref.in,dat.in,wt.vars=c('month.num','year','sal.ref'))
      
      #data to use for model predictions
      pred.dat<-data.frame(sal.ref=ref.in$sal.ref,dec.time=ref.in$dec.time)

      #OLS wtd model
      wt.mod.md<-lm(
        Chla_ugl~dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time),
        weights=ref.wts,
        data=dat.in
        )
      
      #OLS wtd predict
      wt.pred.md<-predict(wt.mod.md,newdata=pred.dat)
      
      if(class(try({
        
        wt.mod.hi<-rq(
          Chla_ugl~dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time),
          weights=ref.wts,
          tau=0.9,
          data=dat.in
          )
        
        wt.pred.hi<-predict(wt.mod.hi,newdata=pred.dat)
        
        }))=='try-error'){ wt.pred.hi<-NA }
      
      #quantile wtd model, 0.1, exception for error
      if(class(try({
        
        wt.mod.lo<-rq(
          Chla_ugl~dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time),
          weights=ref.wts,
          tau=0.1,
          data=dat.in
          )
        
        wt.pred.lo<-predict(wt.mod.lo,newdata=pred.dat)
        
        }))=='try-error'){ wt.pred.lo<-NA }
      
      wt.fit.md<-c(wt.fit.md,wt.pred.md)
      wt.fit.hi<-c(wt.fit.hi,wt.pred.hi)
      wt.fit.lo<-c(wt.fit.lo,wt.pred.lo)      
      
      }
    
    wt.fits<-data.frame(
      year=dat.in[row,]$year,
      month.num=dat.in[row,]$month.num,
      date.f=ref.in$date.f,
      dec.time=dat.in[row,]$dec.time,
      seg,
      yr.grid,
      wt.fit.md,
      wt.fit.hi,
      wt.fit.lo
      )
    
    seg.out<-rbind(seg.out,wt.fits)
    
    }
  
  mods.out<-rbind(mods.out,seg.out)
  
  }

cat(Sys.time()-strt,'\n')

yr.grd<-mods.out

save(yr.grd,file='M:/wq_models/EPC/interp_grids/yrwt_grd.RData')

######
#grid for all normalization variables: year, month, and salinity
#can be used to create monstrous 3D interp grid
yr.grid<-sort(unique(tb.dat$year))
mo.grid<-sort(unique(tb.dat$month.num))
sal.div<-50
sal.grid<-seq(min(tb.dat$sal.ref),max(tb.dat$sal.ref),length=sal.div)

int.grid<-expand.grid(yr.grid,mo.grid,sal.grid,stringsAsFactors=F)
names(int.grid)<-c('year','month.num','sal.grid')

######
#plot interpolation grids

###
#by salinity

ylabs<-expression(paste('chl  ',italic(a),' (',italic(mu),'g',l^-1,')'))

#min, max sal.ref vals to plot....
lim.vals<-aggregate(
  sal.ref~month.num+seg,
  FUN=function(x) cbind(quantile(x,0.05),quantile(x,0.95)),
  data=epc.est
  )
names(lim.vals)[2]<-'seg'

to.plo<-merge(int.grd,lim.vals,by=c('month.num','seg'),all.x=T)

#interp grid removing extreme values
p<-ggplot(to.plo,aes(x=dec.time,y=sal.grid)) + 
  geom_tile(aes(fill=exp(wt.fit.md)),width=0.1) + #adjust this to fill gaps
#   scale_fill_brewer(type='div',palette = 'BuGn') +
   scale_fill_gradient2(name=ylabs,low='blue',mid='lightgreen',high='red',midpoint=20) +
   geom_line(aes(x=dec.time,y=sal.ref[,2])) +
   geom_line(aes(x=dec.time,y=sal.ref[,1])) +
  facet_wrap(~seg,nrow=2,ncol=2) +
  theme_bw() +
  scale_x_continuous(
    breaks=seq(1974,2012,by=2),
    name='Date',
    expand = c(0,0)
    ) + 
  scale_y_continuous(name='Proportion freshwater',expand = c(0,0)) +
  theme(
    axis.text.x=element_text(angle = 90, vjust=0.5,hjust=1)
    )

pdf('C:/Users/mbeck/Desktop/sal_grd.pdf',width=11,height=6.5,family='serif')
print(p)
dev.off()

###
#by year

to.plo<-yr.grd

p<-ggplot(to.plo,aes(x=dec.time,y=as.numeric(as.character(yr.grid)))) + 
  geom_tile(aes(fill=exp(as.numeric(wt.fit.md))),width=0.2) + #adjust this to fill gaps
  scale_fill_gradient2(name=ylabs,low='blue',mid='lightgreen',high='red',midpoint=20) +
  facet_wrap(~seg,nrow=2,ncol=2) +
  theme_bw() +
  scale_x_continuous(
    breaks=seq(1974,2012,by=2),
    name='Date',
    expand = c(0,0)
    ) + 
  scale_y_continuous(
    name='Estimated annual condition',
    breaks=seq(1974,2012,by=2),
    expand = c(0,0)
    ) +
  theme(
    axis.text.x=element_text(angle = 90, vjust=0.5,hjust=1)
    ) 

pdf('C:/Users/mbeck/Desktop/yr_grd.pdf',width=11,height=6.5,family='serif')
print(p)
dev.off()

###
#by month

to.plo<-mo.grd

p<-ggplot(to.plo,aes(x=as.numeric(dec.time),y=mo.grid)) + 
  geom_tile(aes(fill=exp(as.numeric(wt.fit.md))),height=1/12,width=0.25) + #adjust this to fill gaps
  scale_fill_gradient2(name=ylabs,low='blue',mid='lightgreen',high='red',midpoint=20) +
  facet_wrap(~seg,nrow=2,ncol=2) +
  theme_bw() +
  scale_x_continuous(
    name='Date',
    breaks=seq(1974,2012,by=2),
    expand = c(0,0)
    ) + 
  scale_y_continuous(
    breaks=seq(min(to.plo$mo.grid),max(to.plo$mo.grid),length=12),
    labels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),
    name='Estimated monthly condition',
    expand = c(0,0)
    ) +
  theme(
    axis.text.x=element_text(angle = 90, vjust=0.5,hjust=1)
    )

pdf('C:/Users/mbeck/Desktop/mo_grd.pdf',width=11,height=6.5,family='serif')
print(p)
dev.off()

