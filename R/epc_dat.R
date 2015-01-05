######
#data import for Tampa Bay
#tb.dat is used for models, it is raw chlorophyll data aggregated by month

######
library(maptools)
library(raster)
library(rgdal)
library(scales)
library(reshape)
library(RColorBrewer)
library(ggplot2)
library(quantreg)
library(rgeos)
library(gridExtra)
library(Metrics)
library(gtable)
library(Hmisc)
library(plyr)

source('R/epc_funs.r')

######
#import data

#station wq data, named 'tb.chl'
data(epc_tb_chl)

#segment data, primary data frame for analysis, names 'tb.dat'
#note that wq data (except salinity) are log-transformed
data(epc_tb_dat)

#salinity interpolation grid, named 'sal.grd'
data(salwt_grd)

#estimated data using salinity grid, named 'epc.est'
#model predictions are log-transformed
data(epc_est)

#estimated data using obs, named 'epc.est.act'
#model predictions are log-transformed
#contains residuals (obs-pred, + is overfit), also log-transformed
data(epc_est_act)

#normalized model data by salinity, named 'sal.nrm'
data(sal_nrm)

#load seagrass data, named 'sg.dat'
data(epc_seagrass)

#load ENSO data, named enso
data(ENSO)


