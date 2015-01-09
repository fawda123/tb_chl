#created Aug 2013 by M. Beck
#Hirsch model for EPC chlorophyll data
#estimated data are obtained - from interpolation grid and actual data

######
#import data

rm(list=ls())
source('R/epc_dat.r')

######

#get estimated data from interpolation grids
#add obs column to sal.grd to indicate which is matching

obs<-numeric(nrow(sal.grd))

segs<-unique(tb.dat$seg)

fit<-NULL

strt<-Sys.time()

for(seg in segs){
  
  dat.in<-tb.dat[tb.dat$seg==seg,]
  
  seg.fit<-NULL

  for(row in 1:nrow(dat.in)){
  
    cat(paste(seg,row,'of',nrow(dat.in)),'\n')
    flush.console()
    
    #observation to normalize
    ref.in<-dat.in[row,]
    
    #get inteprolation grid data
    grd.in<-sal.grd[sal.grd$seg == seg & sal.grd$dec.time==ref.in$dec.time,]
    
    #matching closest values in grid to actual input
    grd.mtch<-grd.in[which.min(abs(grd.in$sal.grid-ref.in$sal.ref)),c('fit.md','fit.hi','fit.lo','b.md','b.hi','b.lo','bt.md','bt.hi','bt.lo')]
    
    obs[as.numeric(row.names(grd.mtch))]<-1
    
    #output for segment
    seg.fit<-rbind(seg.fit,cbind(ref.in,grd.mtch))
    
    }
  
  #output for all
  fit<-rbind(fit,seg.fit)
  
  }

Sys.time()-strt

epc.est<-fit

save(epc.est,file='data/epc_est.RData')

#save obs to sal.grd
sal.grd<-data.frame(
  sal.grd[,c('seg','year','month.num','dec.time','sal.grid')],
  obs,
  sal.grd[,c('fit.md','fit.hi','fit.lo','b.md','b.hi','b.lo','bt.md','bt.hi','bt.lo')]
  )

save(sal.grd,file='data/salwt_grd.RData')

######
#get actual model estimates, not from interpolation grid
#output appended to mods.out

rm(list = ls())

load('data/epc_tb_dat.RData')
source('R/epc_dat.r')

segs<-unique(tb.dat$seg)

mods.out<-matrix(ncol=12,nrow=nrow(tb.dat))

for(seg in segs){

  dat.in<-tb.dat[tb.dat$seg==seg,]
  
  vecs<-c('fit.lo','fit.md','fit.hi','nl.lo','nl.md','nl.hi','res.lo','res.md','res.hi','bt.lo','bt.md','bt.hi')
  sapply(vecs,function(x) assign(x,numeric(nrow(dat.in)),envir=.GlobalEnv))

  strt<-Sys.time()
  for(row in 1:nrow(dat.in)){
    
    ref.in<-dat.in[row,]

    ref.wts<-wt.fun(ref.in,dat.in,wt.vars=c('month.num','year','sal.ref'))
    
    cat(paste(row,' of ',nrow(dat.in),', ',sum(ref.wts>0),' > 0',sep=''),'\n')
    flush.console()
    
    # crq model, estimates all quants
    mod <- quantreg::crq(
      Surv(Chla_ugl, not_cens, type = "left") ~ 
        dec.time + sal.ref + sin(2*pi*dec.time) + cos(2*pi*dec.time), 
      weights = ref.wts,
      data = dat.in, 
      method = "Portnoy"
      )
    
    # null model for model evaluation, it's specific to the weighting
    nl.mod <- quantreg::crq( # null model
      Surv(Chla_ugl, not_cens, type = "left") ~ 1, 
      weights = ref.wts, 
      data = dat.in, 
      method = "Portnoy"
      )
    
    # sometimes crq fucks up
    test <- try({coef(mod)})
    if('try-error' %in% class(test)){
      fit.lo[row] <- NA
      fit.md[row] <- NA
      fit.hi[row] <- NA
      nl.lo[row] <- NA
      nl.md[row] <- NA
      nl.hi[row] <- NA
      res.lo[row] <- NA
      res.md[row] <- NA
      res.hi[row] <- NA
      bt.lo[row] <- NA
      bt.md[row] <- NA
      bt.hi[row] <- NA
      next
    }
      
    # estimate from null model, it is the coefficient which is just the intercept
    nl <- coef(nl.mod, c(0.1, 0.5, 0.9))
    names(nl) <- paste0('nl.', c('lo', 'md', 'hi'))
    
    # fitted coefficients for each model
    parms <- coef(mod, c(0.1, 0.5, 0.9))
    
    # predicted values by quantile model coefficients
    fits <- sapply(1:3, function(x){
      with(dat.in[row, ], 
        parms[1, x] + parms[2, x] * dec.time + parms[3, x] * sal.ref + parms[4, x] * sin(2*pi*dec.time) + parms[5, x] * cos(2*pi*dec.time)
      )
    })
    names(fits) <- paste0('fit.', c('lo', 'md', 'hi'))
    
    # back transformed predicted values
    bt.fits <- exp(fits)
    names(bt.fits) <- paste0('bt.', c('lo', 'md', 'hi'))
    
    # residuals
    res <- dat.in[row, 'Chla_ugl'] - fits
    names(res) <- paste0('res.', c('lo', 'md', 'hi'))
      
    fit.lo[row] <- fits['fit.lo']
    fit.md[row] <- fits['fit.md']
    fit.hi[row] <- fits['fit.hi']
    nl.lo[row] <- nl['nl.lo']
    nl.md[row] <- nl['nl.md']
    nl.hi[row] <- nl['nl.hi']
    res.lo[row] <- res['res.lo']
    res.md[row] <- res['res.md']
    res.hi[row] <- res['res.hi']
    bt.lo[row] <- exp(fits['fit.lo'])
    bt.md[row] <- exp(fits['fit.md'])
    bt.hi[row] <- exp(fits['fit.hi'])
    
    }
    
  cat(Sys.time()-strt,'\n')
  
  mods.out[tb.dat$seg==seg,]<-cbind(sapply(vecs,get))
  
  }

out<-data.frame(tb.dat,mods.out,stringsAsFactors=F)
names(out)[grep('X',names(out))]<-vecs
epc.est.act<-out

save(epc.est.act,file='data/epc_est_act.RData')
