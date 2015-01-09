# created 8/5/2013, updated 1/7/2015
# script combines raw data from Hillsborough water atlas
# output is two files, one has all data combined 
#and the other only includes chl and salinity data

#may have to download middle bay again w/o cockroach, I removed it but output fomrat might be diff

library(maptools)
library(rgdal)
library(scales)
library(reshape)
library(RColorBrewer)

file.names<-c('tb_lower.txt','tb_middle.txt','tb_hillsborough.txt',
  'tb_old.txt','tb_station94.txt')

#out.long<-NULL
out.shrt<-NULL

for(fl.nm in file.names){
	
	tmp.fl<-read.table(paste0('M:/wq_models/EPC/', fl.nm),header=T,sep='\t',stringsAsFactors=F)
	#out.long<-rbind(out.long,tmp)

  #subset parameters of interest
	tmp<-tmp.fl[grepl('Chla_ugl|Salinity_ppt',tmp.fl$Parameter),]

  #unique columns for month, day, and year
	times<-do.call('rbind',strsplit(tmp$SampleDate,' '))[,1]
	times<-data.frame(do.call('rbind',strsplit(times,'/')),stringsAsFactors=F)
	names(times)<-c('month','day','year')
	tmp<-data.frame(times,tmp,stringsAsFactors=T)
	
  #create month vector for names
	month.con<-data.frame(
		int=seq(1,12),
		month.name=c('January','February','March','April','May','June','July','August',
			'September','October','November','December')
		)
	tmp<-merge(tmp,month.con,by.x='month',by.y='int',all.x=T)
	
  #aggregate by month
	tmp<-aggregate(Result_Value~month+month.name+year+StationID+WBodyID+Actual_Latitude+Actual_Longitude+Parameter,data=tmp,FUN=function(x) median(x, na.rm = T))#function(x) c(mean(x),length(x)))
	
  #create unique columns for each parameter, will have to change for >2 parms
  tmp<-split(tmp,tmp$Parameter)
  parm.names<-names(tmp)
  name.sel<-unique(unlist(lapply(tmp,names)))
  name.sel<-name.sel[!name.sel %in% c('Parameter','Result_Value')]
  tmp<-merge(tmp[[1]],tmp[[2]],by=name.sel)#removes stations w/o both chl and salinity
  tmp<-tmp[,grep('Param',names(tmp),invert=T)] #remove columns with parm.names\
  names(tmp)[grep('Result',names(tmp))]<-parm.names
  
  tmp<-tmp[with(tmp,order(year,as.numeric(month))),]
	
	out.shrt<-rbind(out.shrt,tmp)
	
	}

#change station 94 WBID to lower tampa Bay
out.shrt[out.shrt$StationID=='=94','WBodyID']<-'20011'

tb.chl <- out.shrt
save(tb.chl, file = 'data/epc_tb_chl.RData')

######
#code for modifying tb chlorophyll data

rm(list = ls())

library(maptools)
library(rgdal)

#station wq data
load('data/epc_tb_chl.RData')

#station locations
tb.crds<-unique(data.frame(tb.chl[,names(tb.chl) %in% c('Actual_Longitude','Actual_Latitude','StationID')]))
tb.crds<-data.frame(tb.crds,WBID=NA)
tb.crds2<-coordinates(data.frame(tb.crds$Actual_Longitude,tb.crds$Actual_Latitude))             
tb.crds<-SpatialPointsDataFrame(tb.crds2,data=tb.crds)
#writeSpatialShape(tb.crds,'C:/Users/mbeck/Desktop/tb_sta.shp')

#wbid shapefile for estuaries, subset Tampa Bay
tb<-readShapeSpatial('M:/GIS/WBID_Run44_repr.shp')
tb<-tb[tb$GROUP_NAME %in% 'Tampa Bay',]
tb<-tb[grep('1558*|8049*',tb$WBID),]

#add column to tb.chl for wbid (smaller segs within WBodyID)
#taken from spatial join of tb.crds and tb
sel.vec<-over(tb.crds,geometry(tb))
tb.crds$WBID<-tb$WBID[sel.vec]
tb.chl<-merge(
  tb.chl,
  tb.crds@data[,names(tb.crds) %in% c('StationID','WBID')],
  by='StationID'
)

#convert bay segment ('WBodyID') to character name
tb.chl$WBodyID<-factor(
  tb.chl$WBodyID,
  levels=c('20005','20007','20010','20011'),
  labels=c('Hillsborough Bay','Middle Tampa Bay','Old Tampa Bay',
    'Lower Tampa Bay')
)

#reorder levels for month.name, will be preserved in dput
tb.chl$month.name<-factor(tb.chl$month.name,
  levels=unique(tb.chl$month.name)[order(as.numeric(unique(tb.chl$month)))]
)

#create date column as date object
#assumes day as first day of each month
month.f<-sprintf('%02d',as.numeric(tb.chl$month)) #add zero, two digits max
date.f<-paste(month.f,'01',tb.chl$year,sep='/')
date.f<-as.Date(date.f,'%m/%d/%Y')
tb.chl<-data.frame(tb.chl,date.f)

#create salinity data calibrated by ref
#ref is station 94, unless it't not found, then maximum is used
#note that you still get negative numbers if 93 is not maximum

sal.ref<-numeric(nrow(tb.chl))

uni.dates<-unique(tb.chl$date.f)
for(uni.date in 1:length(uni.dates)){
  
  print(uni.dates[uni.date])
  flush.console()
  
  tmp<-tb.chl[tb.chl$date.f == uni.dates[uni.date],]
  
  tmp.ref<-tmp[tmp$StationID=='=94','Salinity_ppt']
  
  if(length(tmp.ref)==0) tmp.ref<-max(tmp[,'Salinity_ppt'])
  
  tmp.sal<-pmax(0,1-tmp$Salinity_ppt/tmp.ref)
  
  sal.ref[tb.chl$date.f==uni.dates[uni.date]]<-tmp.sal
  
  }

tb.chl<-data.frame(tb.chl,sal.ref)

#get load data, merge with chlorophyll data
tb.load<-read.csv('M:/wq_models/EPC/TBEP_loads_85_07.csv',header=T)

load.sel<-names(tb.load) %in% c('Year','Month','TN.Load..kg.month.','Bay.Name')
tb.chl<-merge(tb.chl,tb.load[load.sel],
  by.x=c('month','year','WBodyID'),
  by.y=c('Month','Year','Bay.Name'),
  all.x=T
  )

names(tb.chl)[names(tb.chl)=='TN.Load..kg.month.']<-'TN_kg_mo'
names(tb.chl)[names(tb.chl)=='WBodyID']<-'seg'

##
#create decimal dates, invidual months and year plus months
#julian date as fraction of 365
month.num<-(strptime(tb.chl$date.f,'%Y-%m-%d')$yday+1)
lp.days<-month.num %in% c(61,92, 122,153,183,214,245,275,306,336)
month.num[lp.days]<-month.num[lp.days] - 1
month.num<-month.num/365
tb.chl$month.num<-month.num
#decimal time (year and month)
tb.chl$dec.time<-as.numeric(tb.chl$year) + tb.chl$month.num

##
#remove 2013 data because it's incomplete for whole year
#remove 1972, 1973 because no chla data
tb.chl<-tb.chl[!tb.chl$year %in% c('2013'),]

##
#detection limit for values before 2006 was 2.4
#these were set to 1.2 as suggested by Gilbert 1987
#lapply(split(tb.chl,tb.chl$year),function(x) min(x$Chla_ugl))
# tb.chl$Chla_ugl[as.numeric(tb.chl$year)<2006 & tb.chl$Chla_ugl<=2.4]<-1.2

#save output
save(tb.chl,file='data/epc_tb_chl.RData')

######
#create tb.dat from tb.chl
#tb.dat is same as tb.chl except data are aggregated by major segment
#tb.chl has data by stations

rm(list = ls())

load('data/epc_tb_chl.RData')

tb.dat<-aggregate(
  cbind(Chla_ugl,Salinity_ppt,sal.ref,TN_kg_mo)~year+month.name+date.f+month.num+dec.time+seg,
  data=tb.chl,
  median,
  na.action='na.pass'
)
tb.dat$Chla_ugl<-log(tb.dat$Chla_ugl)
tb.dat$TN_kg_mo<-log(tb.dat$TN_kg_mo)

save(tb.dat,file='data/epc_tb_dat.RData')

######
#get nutrient data for TN and TP, append to tb.dat

#may have to download middle bay again w/o cockroach, I removed it but output fomrat might be diff

rm(list=ls())

library(maptools)
library(rgdal)
library(scales)
library(reshape)
library(RColorBrewer)

file.names<-c('tb_lower.txt','tb_middle.txt','tb_hillsborough.txt',
  'tb_old.txt','tb_station94.txt')

out.shrt<-NULL

for(fl.nm in file.names){
  
	tmp.fl<-read.table(paste0('M:/wq_models/EPC/', fl.nm),header=T,sep='\t',stringsAsFactors=F)

  #subset parameters of interest
	tmp<-tmp.fl[grepl('TN_ugl|TP_ugl',tmp.fl$Parameter),]

  #unique columns for month, day, and year
	times<-do.call('rbind',strsplit(tmp$SampleDate,' '))[,1]
	times<-data.frame(do.call('rbind',strsplit(times,'/')),stringsAsFactors=F)
	names(times)<-c('month','day','year')
	tmp<-data.frame(times,tmp,stringsAsFactors=T)
	
  #create month vector for names
	month.con<-data.frame(
		int=seq(1,12),
		month.name=c('January','February','March','April','May','June','July','August',
			'September','October','November','December')
		)
	tmp<-merge(tmp,month.con,by.x='month',by.y='int',all.x=T)
	
  #aggregate by month
	tmp<-aggregate(Result_Value~month+month.name+year+StationID+WBodyID+Actual_Latitude+Actual_Longitude+Parameter,data=tmp,FUN=function(x) median(x, na.rm = T))#function(x) c(mean(x),length(x)))
	
  #create unique columns for each parameter, will have to change for >2 parms
  tmp<-split(tmp,tmp$Parameter)
  parm.names<-names(tmp)
  name.sel<-unique(unlist(lapply(tmp,names)))
  name.sel<-name.sel[!name.sel %in% c('Parameter','Result_Value')]
  tmp<-merge(tmp[[1]],tmp[[2]],by=name.sel,all=T)
  tmp<-tmp[,grep('Param',names(tmp),invert=T)] #remove columns with parm.names\
  names(tmp)[grep('Result',names(tmp))]<-parm.names
  
  tmp<-tmp[with(tmp,order(year,as.numeric(month))),]
	
	out.shrt<-rbind(out.shrt,tmp)
	
	}

#change station 94 WBID to lower tampa Bay
out.shrt[out.shrt$StationID=='=94','WBodyID']<-'20011'

NP.dat<-out.shrt

#convert WBodyID to names
NP.dat$WBodyID<-factor(
  NP.dat$WBodyID,
  levels=c('20005','20007','20010','20011'),
  labels=c('Hillsborough Bay','Middle Tampa Bay','Old Tampa Bay',
    'Lower Tampa Bay')
  )
names(NP.dat)[names(NP.dat) %in% 'WBodyID']<-'seg'

#aggregate by station
NP.dat<-aggregate(
  cbind(TN_ugl,TP_ugl)~year+month.name+seg,
  data=NP.dat,
  median,
  na.action='na.pass'
)

#log transform P and N concentrations
NP.dat$TP_ugl<-log(NP.dat$TP_ugl)
NP.dat$TN_ugl<-log(NP.dat$TN_ugl)

load('data/epc_tb_dat.RData')

tb.dat<-merge(tb.dat,NP.dat,by=c('year','seg','month.name'),all.x=T)

save(tb.dat,file='data/epc_tb_dat.RData')

######
# add censoring information

rm(list = ls())

load('data/epc_tb_dat.RData')

# data frame of censored values by year, same by segment
cens_min <- data.frame(
  year = unique(tb.dat$year), 
  lim = log(c(rep(2.4, 37), rep(3.1, 2)))
)
# combine with aggregated tampa data
tb.dat <- merge(tb.dat, cens_min, by = 'year')

# floor chl dat at limit if not
# not_cens is vector of TF incidating if values are censored (T is not censored)
tb.dat$Chla_ugl <- with(tb.dat, pmax(Chla_ugl, lim))
tb.dat$not_cens <- with(tb.dat, Chla_ugl > lim)

save(tb.dat, file = 'data/epc_tb_dat.RData')

######
# extra stuff, don't need to run for producing data

# evaluation of detection limits in observed data by stations
# don't need to do this for data processing but used to find annual detection limits that are used in processing 

setwd('M:/wq_models/EPC')

file.names <- c('tb_lower.txt', 'tb_middle.txt', 'tb_hillsborough.txt',
  'tb_old.txt', 'tb_station94.txt')

##
# find detection limits by year, segment
out.lim <- NULL

for(fl.nm in file.names){
  
  tmp.fl <- read.table(fl.nm,header=T,sep='\t',stringsAsFactors=F)

  tmp <- tmp.fl[grepl('Chla_ugl',tmp.fl$Parameter),]
  
  yrs <- strsplit(tmp$SampleDate, ' ')
  yrs <- do.call('rbind', yrs)[, 1]
  yrs <- as.Date(yrs, format = '%m/%d/%Y')
  tmp$yrs <- format(yrs, '%Y')
  nondetect <- tmp[tmp$QACode %in% 'U', ]  

  out <- unique(nondetect[, c('WaterBodyName', 'yrs', 'Parameter', 'Result_Value', 
    'QACode')])
  
  # get prop of obs at det limit by year for the segment
  prop <- plyr::ddply(tmp, 
    .variable = c('yrs'), 
    .fun = function(x) sum(x$QACode %in% 'U')/nrow(x)
    )
  out <- merge(prop, out, by = 'yrs', all.x = T)
  
  out.lim <- rbind(out.lim, out)
  
}

# data frame of proportion of obs that are censored and censoring value by segment and year
out.lim <- data.frame(out.lim, row.names = c(1:nrow(out.lim)))
names(out.lim)[names(out.lim) %in% 'V1'] <- 'prop'





