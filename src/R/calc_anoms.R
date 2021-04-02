# Description: Calculate the climatology and meteorological anomalies quickly
# with the data.table package. Considerable RAM (~60+ GB) is required. 

library(arrow)
library(tidyverse); 
library(data.table); library(lubridate);
setDTthreads(threads=parallel::detectCores()-3)

#*******************************************************************************
# Get data  ---------------------------------------------------------
#*******************************************************************************
tmp <- arrow::read_parquet("/home/sami/srifai@gmail.com/work/research/data_general/Oz_misc_data/ARD_ndvi_aclim_2021-02-28.parquet")

# data.table 
tmp <- setDT(tmp) # OR: tmp <- as.data.table(tmp)
tmp <- tmp %>% select(-x.x,-x.y,-y.x,-y.y) %>% rename(x=x_vi, y=y_vi)
tmp <- tmp[, pet:=ifelse(pet<10,10,pet)]
tmp <- tmp[, id := .GRP, by=.(x,y)]

# Subsetting to just one type of vegetation class
# tmp_vc <- tmp[date==ymd("2000-01-01"),.(x,y,vc)]
# vec_vc <- unique(tmp$vc) %>% sort
# tmp_vc <- tmp_vc[vc=="Eucalypt Open Forests"]
# tmp <- tmp[tmp_vc,on=.(x,y)] # subset to just the selected vegetation class
#*******************************************************************************
#* END SECTION
#*******************************************************************************

#*******************************************************************************
# Calculate Climatology --------------------------------------------------------
#*******************************************************************************
# calc norms 
tmp <- tmp[, `:=`(month = month(date))] # create month
tmp <- tmp[, `:=`(year = year(date))]   # create year
tmp <- tmp[,`:=`("pe" = precip/pet)]
# norms_ndvi <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
#                   .("ndvi_u" = mean(ndvi_mcd,na.rm=TRUE), 
#                     "ndvi_sd" = sd(ndvi_mcd,na.rm=TRUE)),
#                   by=.(x,y,month)] # joining on x,y,month
norms_p <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
               .("precip_u" = mean(precip,na.rm=TRUE), 
                 "precip_sd" = sd(precip,na.rm=TRUE)),
               by=.(x,y,month)]
norms_vpd <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
               .("vpd15_u" = mean(vpd15,na.rm=TRUE), 
                 "vpd15_sd" = sd(vpd15,na.rm=TRUE)),
               by=.(x,y,month)]
norms_pet <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                 .(pet_u = mean(pet,na.rm=TRUE), 
                   pet_sd = sd(pet,na.rm=TRUE)),
                 by=.(x,y,month)]
norms_pe <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                .(pe_u = mean(pe,na.rm=TRUE), 
                  pe_sd = sd(pe,na.rm=TRUE)),
                by=.(x,y,month)]
norms_tmax <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                .(tmax_u = mean(tmax,na.rm=TRUE), 
                  tmax_sd = sd(tmax,na.rm=TRUE)),
                by=.(x,y,month)]
norms_tmin <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                  .(tmin_u = mean(tmin,na.rm=TRUE), 
                    tmin_sd = sd(tmin,na.rm=TRUE)),
                  by=.(x,y,month)]

norms_map <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                 .("ap" = sum(precip)),
                 by=.(x,y,year)][,.("map"=mean(ap), 
                                    "ap_sd"=sd(ap)),by=.(x,y)]
norms_mapet <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                   .("apet" = sum(pet)),
                   by=.(x,y,year)][,.("mapet"=mean(apet), 
                                      "apet_sd"=sd(apet)),by=.(x,y)]
norms_mape <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                  .("ape" = mean(pe)),
                  by=.(x,y,year)][,.("mape"=mean(ape), 
                                     "ape_sd"=sd(ape)),by=.(x,y)]
norms_matmax <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                 .("atmax" = mean(tmax)),
                 by=.(x,y,year)][,.("matmax"=mean(atmax), 
                                    "atmax_sd"=sd(atmax)),by=.(x,y)]
norms_matmin <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                    .("atmin" = mean(tmin)),
                    by=.(x,y,year)][,.("matmin"=mean(atmin), 
                                       "atmin_sd"=sd(atmin)),by=.(x,y)]
norms_mavpd <- tmp[date>=ymd('1982-01-01')&date<=ymd("2011-12-31"), # filter to ref period
                    .("avpd" = mean(vpd15)),
                    by=.(x,y,year)][,.("mavpd15"=mean(avpd), 
                                       "avpd15_sd"=sd(avpd)),by=.(x,y)]

# join all the data frames ***
norms <- norms_p[norms_pet, on=.(x,y,month)] # join data.tables
norms <- norms[norms_vpd, on=.(x,y,month)] # join data.tables
norms <- norms[norms_pe, on=.(x,y,month)] # join data.tables
norms <- norms[norms_tmax, on=.(x,y,month)] # join data.tables
norms <- norms[norms_tmin, on=.(x,y,month)]
norms <- norms[norms_map, on=.(x,y)]
norms <- norms[norms_mapet, on=.(x,y)]
norms <- norms[norms_mape, on=.(x,y)]
norms <- norms[norms_matmax, on=.(x,y)]
norms <- norms[norms_matmin, on=.(x,y)]
norms <- norms[norms_mavpd, on=.(x,y)]

# norms <- norms[norms_ndvi, on=.(x,y,month)]
tmp <- norms[tmp, on=.(x,y,month)]
rm(norms); 
gc(verbose = T, reset = T, full = T)
#*******************************************************************************
#* END SECTION
#*******************************************************************************

#*******************************************************************************
# Calculate the anomalies ------
#*******************************************************************************
tmp <- tmp[, `:=`(#ndvi_anom = ndvi_mcd - ndvi_u, 
                  precip_anom = precip-precip_u,  # calc raw anomaly 
                  pet_anom = pet-pet_u, 
                  pe_anom = pe-pe_u, 
                  tmax_anom = tmax-tmax_u, 
                  tmin_anom = tmin-tmin_u, 
                  vpd15_anom = vpd15 - vpd15_u)]
tmp <- tmp[, `:=`(#ndvi_anom_sd = ndvi_anom/ndvi_sd,
                  precip_anom_sd = precip_anom/precip_sd,  # calc sd anomaly 
                  pet_anom_sd = pet_anom/pet_sd, 
                  pe_anom_sd = pe_anom/pe_sd, 
                  tmax_anom_sd = tmax_anom/tmax_sd, 
                  tmin_anom_sd = tmin_anom/tmin_sd, 
                  vpd15_anom_sd = vpd15_anom/vpd15_sd)]
#*******************************************************************************
#* END SECTION
#*******************************************************************************

#*******************************************************************************
# Calculate the multi-year anomalies ------
#*******************************************************************************
# calculate the rolling 12-month sums 
# tmp <- tmp[order(x,y,date)][, ndvi_12mo := frollmean(ndvi_mcd,n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, precip_12mo := frollsum(precip,n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pet_12mo := frollsum(pet,n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pe_12mo := frollmean(pe,n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, tmax_anom_12mo := frollapply(tmax_anom,FUN=max,
                                                           n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, tmin_anom_12mo := frollapply(tmin_anom,FUN=max,
                                                           n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, vpd15_12mo := frollapply(vpd15,FUN=mean,
                                                          n = 12,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, vpd15_anom_12mo := frollapply(vpd15_anom,FUN=mean,
                                                           n = 12,fill = NA,align='right'), by=.(x,y)]
# tmp <- tmp[order(x,y,date)][, ndvi_anom_12mo := frollmean(ndvi_anom,n = 12,fill = NA,align='right'), by=.(x,y)]

gc(verbose = T, reset = T, full = T)

# calculate rolling 3-month anomaly
tmp <- tmp[order(x,y,date)][, precip_anom_3mo := frollsum(precip_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pet_anom_3mo := frollsum(pet_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pe_anom_3mo := frollmean(pe_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, tmax_anom_3mo := frollmean(tmax_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, tmin_anom_3mo := frollmean(tmin_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, vpd15_anom_3mo := frollmean(vpd15_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
gc()

# calculate rolling 6-month anomaly
tmp <- tmp[order(x,y,date)][, precip_anom_6mo := frollsum(precip_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pet_anom_6mo := frollsum(pet_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pe_anom_6mo := frollmean(pe_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, tmax_anom_6mo := frollmean(tmax_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, vpd15_anom_6mo := frollmean(vpd15_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
gc()

# # rolling 2 year sums
# tmp <- tmp[order(x,y,date)][, precip_24mo := frollsum(precip,n = 24,fill = NA,align='right'), by=.(x,y)]
# tmp <- tmp[order(x,y,date)][, pet_24mo := frollsum(pet,n = 24,fill = NA,align='right'), by=.(x,y)]
# tmp <- tmp[order(x,y,date)][, pe_24mo := frollsum(pe,n = 24,fill = NA,align='right'), by=.(x,y)]
# tmp <- tmp[order(x,y,date)][, vpd15_24mo := frollmean(vpd15,n = 24,fill = NA,align='right'), by=.(x,y)]
# gc()

# rolling 3 year sums
tmp <- tmp[order(x,y,date)][, precip_36mo := frollsum(precip,n = 36,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pet_36mo := frollsum(pet,n = 36,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, pe_36mo := frollmean(pe,n = 36,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, vpd15_36mo := frollmean(vpd15,n = 36,fill = NA,align='right'), by=.(x,y)]
gc()

# # rolling 4 year sums
# tmp <- tmp[order(x,y,date)][, precip_48mo := frollsum(precip,n = 48,fill = NA,align='right'), by=.(x,y)]
# tmp <- tmp[order(x,y,date)][, pet_48mo := frollsum(pet,n = 48,fill = NA,align='right'), by=.(x,y)]
# tmp <- tmp[order(x,y,date)][, pe_48mo := frollsum(pe,n = 48,fill = NA,align='right'), by=.(x,y)]
# tmp <- tmp[order(x,y,date)][, vpd15_48mo := frollmean(vpd15,n = 48,fill = NA,align='right'), by=.(x,y)]
# gc()

# calc anoms of rolling sums
tmp <- tmp[, `:=`(precip_anom_12mo = precip_12mo-map)]
tmp <- tmp[, `:=`(pet_anom_12mo = pet_12mo-mapet)]
tmp <- tmp[, `:=`(pe_anom_12mo = pe_12mo-mape)]
tmp <- tmp[, `:=`(vpd15_anom_12mo = vpd15_12mo-mavpd15)]
gc()

# tmp <- tmp[, `:=`(precip_anom_24mo = precip_24mo-2*map)]
# tmp <- tmp[, `:=`(pet_anom_24mo = pet_24mo-2*mapet)]
# tmp <- tmp[, `:=`(pe_anom_24mo = pe_24mo-2*mape)]
# tmp <- tmp[, `:=`(vpd15_anom_24mo = vpd15_24mo-mavpd15)]
# gc()

tmp <- tmp[, `:=`(precip_anom_36mo = precip_36mo-3*map)]
tmp <- tmp[, `:=`(pet_anom_36mo = pet_36mo-3*mapet)]
tmp <- tmp[, `:=`(pe_anom_36mo = pe_36mo-3*mape)]
tmp <- tmp[, `:=`(vpd15_anom_36mo = vpd15_36mo-mavpd15)]
gc()
 
# tmp <- tmp[, `:=`(precip_anom_48mo = precip_48mo-4*map)]
# tmp <- tmp[, `:=`(pet_anom_48mo = pet_48mo-4*mapet)]
# tmp <- tmp[, `:=`(pe_anom_48mo = pe_48mo-4*mape)]
# tmp <- tmp[, `:=`(vpd15_anom_48mo = vpd15_48mo-mavpd15)]
# gc()
 
# tmp <- tmp[order(x,y,date)][, tmax_anom_24mo := frollmean(tmax_anom,n = 24,fill = NA,align='right'), by=.(x,y)]
tmp <- tmp[order(x,y,date)][, tmax_anom_36mo := frollmean(tmax_anom,n = 36,fill = NA,align='right'), by=.(x,y)]
gc(verbose = T, reset = T, full = T)
#*******************************************************************************
#* END SECTION
#*******************************************************************************


#*******************************************************************************
#* Add season and hydro year -----
#*******************************************************************************
vec_dates <- data.table(date=sort(unique(tmp$date))) %>% 
  .[,quarter:=quarter(date)] %>% 
  mutate(q = case_when(quarter==1~"DJF",
                       quarter==2~"MAM",
                       quarter==3~"JJA",
                       quarter==4~"SON")) %>% 
  mutate(season = factor(q,
                         levels=c("DJF","MAM","JJA","SON"), 
                         ordered=T)) %>% 
  select(date,season) %>% 
  mutate(hydro_year = year(date+months(1)))

tmp <- tmp[vec_dates,on=.(date)]
gc(verbose = T, reset = T, full = T)
#*******************************************************************************
#* END SECTION
#*******************************************************************************



#*******************************************************************************
# ATTACH NVIS ---------------------------------------------------------
#*******************************************************************************
library(sf); library(stars)
library(tidyverse); 
library(data.table); library(lubridate);
library(dtplyr)
setDTthreads(threads=8)
nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif") %>% 
  set_names('veg_class')
codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt",
                         fwf_widths(c(2,100)), skip = 1) %>%
  set_names(c("veg_class","veg_class_descrip")) %>%
  mutate(vc = as.factor(veg_class_descrip))
base_coords <- unique(tmp[,.(x,y)])
base_coords <- st_as_sf(base_coords,crs=st_crs(4326),coords=c("x","y"))
nvis <- st_extract(nvis, pts = base_coords)
nvis <- bind_cols(nvis %>% st_coordinates() %>% as.data.table,nvis$veg_class) %>% set_names(c("x","y","veg_class"))
nvis <- merge(nvis,codes,by='veg_class')
nvis <- nvis[veg_class<=15]
tmp <- merge(tmp,nvis,by=c("x","y"))
rm(base_coords)
#*******************************************************************************
#* END SECTION
#*******************************************************************************



arrow::write_parquet(tmp, sink="/home/sami/scratch/ARD_ndvi_aclim_anoms.parquet",
              compression = 'snappy')
gc()
