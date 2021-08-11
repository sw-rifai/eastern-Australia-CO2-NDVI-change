# library(brms)
library(mgcv); library(gratia); library(nls.multstart)
library(tidyverse)
library(data.table); setDTthreads(threads = 8)
library(lubridate); 
library(dtplyr);
library(nls.multstart)
options(mc.cores=parallel::detectCores()-3) 
set.seed(333)
# IMPORT ###################################################################

# vegetation index record
vi <- arrow::read_parquet("../data_general/MCD43/MCD43_AVHRR_NDVI_hybrid_2020-10-12.parquet" 
                          # col_select = c("x","y","date",
                          #                "ndvi_c","ndvi_mcd","ndvi_hyb", 
                          #                "evi2_hyb","evi2_mcd","sz")
) %>% 
  as.data.table()
vi <- vi %>% lazy_dt() %>% 
  mutate(ndvi_hyb_e1 = coalesce(ndvi_mcd_nm_pred, NA_real_),
         ndvi_hyb_e2 = coalesce(ndvi_mcd, NA_real_)) %>% 
  mutate(ndvi_hyb = coalesce(ndvi_hyb_e2, ndvi_hyb_e1)) %>% 
  mutate(ndvi_hyb = ifelse(between(ndvi_hyb,0,1),ndvi_hyb,NA_real_)) %>% 
  as.data.table()
norms_vi <- vi[,`:=`(month=month(date))] %>% 
  .[,.(ndvi_u = mean(ndvi_hyb,na.rm=TRUE), 
       ndvi_sd = sd(ndvi_hyb,na.rm=TRUE)),keyby=.(x,y,month)]
vi <- norms_vi[vi,on=.(x,y,month)] %>% 
  .[,`:=`(ndvi_anom = ndvi_hyb - ndvi_u)] %>% 
  .[,`:=`(ndvi_anom_sd = ndvi_anom/ndvi_sd)]

# vi <- mcf[,.(x,y,date,soil,gv,npv)][vi,on=.(x,y,date)]

# Load climate data 
dat <- arrow::read_parquet("/home/sami/scratch/ARD_ndvi_aclim_anoms.parquet",
                           col_select = c(
                             "date", "hydro_year", "id",
                             # "season",
                             "precip",  "precip_anom", 
                             "precip_anom_3mo","precip_anom_36mo",
                             "precip_anom_12mo",
                             "map",
                             "precip_12mo","precip_36mo",
                             "tmax","tmax_anom","tmax_anom_sd", "matmax",
                             "tmin","tmin_anom",
                             "vpd15","vpd15_anom","vpd15_anom_sd","mavpd15",
                             "vpd15_12mo",
                             "vpd15_u",
                             "pet","mapet","pet_anom","pet_anom_3mo","pet_u","pet_sd",
                             "pet_anom_sd", "pet_12mo","pet_36mo",
                             "pe","mape",
                             # "ndvi_u",
                             # "ndvi_anom",
                             # "ndvi_anom_12mo",
                             # "ndvi_anom_sd",
                             # "ndvi_mcd",
                             'vc','veg_class',
                             'month',
                             "x", "y", "year")) %>% 
  as.data.table() %>% 
  .[is.infinite(mape)==F]
dat <- dat[order(x,y,date)][, vpd15_12mo := frollmean(vpd15,n = 12,fill = NA,align='right',na.rm=TRUE), by=.(x,y)]
dat <- dat[,`:=`(pe = precip/pet, 
                 pe_12mo = precip_12mo/pet_12mo)]
dat <- merge(dat, 
             vi,
             by=c("x","y","date"), 
             all=TRUE,allow.cartesian=TRUE)
dat <- dat[order(x,y,date)][, ndvi_3mo := frollmean(ndvi_hyb,n = 3,fill = NA,align='center',na.rm=TRUE), by=.(x,y)]
# rm(vi); 
gc(full=TRUE)

# Attach season
dat[,`:=`(year=year(date),month=month(date))] %>%
  .[,`:=`(season = case_when(month%in%c(3:5)~'MAM',
                             month%in%c(6:8)~'JJA',
                             month%in%c(9:11)~'SON',
                             month%in%c(12,1,2)~'DJF'))]
dat[,`:=`(season = factor(season, levels = c('SON','DJF','MAM','JJA'),ordered = TRUE))]
dat[,`:=`(hydro_year=year(date+months(1)))]
dat[,`:=`(epoch=ifelse(date>=ymd("2001-01-01"),1,0))]

# FILTER TO LON >= 140 !!! **********
dat <- dat[x>=140]

coords_keep <- dat %>% lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(nobs_total = sum(is.na(ndvi_hyb)==F)) %>% 
  ungroup() %>% 
  as.data.table()
dat <- merge(dat, coords_keep, by=c("x","y"))

# attach CO2
mlo <- readr::read_table("../data_general/CO2_growth_rate/co2_mm_mlo_20200405.txt", 
                         skip = 72, col_names = F) %>% 
  set_names(
    c("year","month","ddate","co2_avg","co2_int","co2_trend","ndays")
  ) %>% 
  mutate(date = ymd(paste(year,month,1))) %>% 
  select(date,co2_int,co2_trend) %>% 
  as.data.table()
dat <- merge(mlo,dat,by="date")
dat <- dat[#is.na(ndvi_3mo)==F & # why remove the na's from ndvi_3mo?
             is.na(co2_int)==F]
center_co2 <- mean(dat$co2_int)
dat <- dat[,`:=`(cco2=co2_int-center_co2)]
gc()
dat <- dat[is.na(vc)==F]
dat <- dat[str_detect(vc,"Forests") | 
             str_detect(vc, "Eucalypt") |
             str_detect(vc, "Rainforests")]
dat <- dat[x>= 140] # FILTER TO LON >= 140
dat <- dat[ndvi_hyb>0][ndvi_anom_sd > -3.5 & ndvi_anom_sd < 3.5]
dat[,`:=`(pe_anom_12mo = pe_12mo - mape)]
dat[,`:=`(epoch = ifelse(date<ymd("2000-12-31"),'avhrr','modis'))]
dat <- dat %>% mutate(epoch = as_factor(epoch), 
                      season = factor(season, levels=c("SON","DJF","MAM","JJA"))) %>% 
  as.data.table()

#Filter out the black summer
dat <- dat[date<=ymd("2019-08-01")]
ldat <- dat %>% lazy_dt()
#*******************************************************************************

#split test & train ------------------------------------------------------------
set.seed(321)
train_dat <- dat[mape<1.5][#is.na(ndvi_3mo)==F & 
  is.na(pe_12mo)==F][sample(.N, 4e6)]
test_dat <- dat[mape<1.5][#is.na(ndvi_3mo)==F & 
  is.na(pe_12mo)==F][sample(.N, 4e6)]
gc(full = T)
#*******************************************************************************


# Fit models --------------------------------------------------------------
n4_son <- train_dat[season=='SON'] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape)+ 
                                 B3*(cco2*pe_anom_12mo/mape) + 
                                 B4*as.numeric(epoch),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0,B4=-0.1),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001,B4=0.1),
                               # supp_errors = 'Y',
                               na.action = na.omit)
n4_djf <- train_dat[season=='DJF'] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape) + 
                                 B4*as.numeric(epoch),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0,B4=-0.1),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001,B4=0.1),
                               # supp_errors = 'Y',
                               na.action = na.omit)
n4_mam <- train_dat[season=='MAM'] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape) + 
                                 B4*as.numeric(epoch),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0,B4=-0.1),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001,B4=0.1),
                               # supp_errors = 'Y',
                               na.action = na.omit)
n4_jja <- train_dat[season=='JJA'] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape) + 
                                 B4*as.numeric(epoch),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0,B4=-0.1),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001,B4=0.1),
                               # supp_errors = 'Y',
                               na.action = na.omit)


e1_n4_son <- train_dat[season=='SON'][date<=ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape)+ 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)
e1_n4_djf <- train_dat[season=='DJF'][date<=ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)
e1_n4_mam <- train_dat[season=='MAM'][date<=ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)
e1_n4_jja <- train_dat[season=='JJA'][date<=ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)

e2_n4_son <- train_dat[season=='SON'][date>ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape)+ 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)
e2_n4_djf <- train_dat[season=='DJF'][date>ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)
e2_n4_mam <- train_dat[season=='MAM'][date>ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)
e2_n4_jja <- train_dat[season=='JJA'][date>ymd("2000-12-31")] %>% 
  nls.multstart::nls_multstart(ndvi_3mo ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(pe_anom_12mo/mape) + 
                                 B2*(cco2*mape) + 
                                 B3*(cco2*pe_anom_12mo/mape),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001),
                               # supp_errors = 'Y',
                               na.action = na.omit)




# Evaluate Models ---------------------------------------------------------

gofs <- expand_grid(r2=NA_real_,rmse=NA_real_,season=c('SON','DJF','MAM','JJA'),
            epoch=c("all","modis","avhrr")) %>% 
  as.data.table()

# merged record
yardstick::rsq_trad_vec(test_dat[season=='SON']$ndvi_3mo,
                        estimate=predict(n4_son,newdata=test_dat[season=="SON"])) -> gofs[epoch=='all'&season=='SON']$r2
yardstick::rsq_trad_vec(test_dat[season=='DJF']$ndvi_3mo,
                        estimate=predict(n4_djf,newdata=test_dat[season=="DJF"])) -> gofs[epoch=='all'&season=='DJF']$r2
yardstick::rsq_trad_vec(test_dat[season=='MAM']$ndvi_3mo,
                        estimate=predict(n4_mam,newdata=test_dat[season=="MAM"])) -> gofs[epoch=='all'&season=='MAM']$r2
yardstick::rsq_trad_vec(test_dat[season=='JJA']$ndvi_3mo,
                        estimate=predict(n4_jja,newdata=test_dat[season=="JJA"])) -> gofs[epoch=='all'&season=='JJA']$r2

yardstick::rmse_vec(test_dat[season=='SON']$ndvi_3mo,
                    estimate=predict(n4_son,newdata=test_dat[season=="SON"]))  -> gofs[epoch=='all'&season=='SON']$rmse
yardstick::rmse_vec(test_dat[season=='DJF']$ndvi_3mo,
                    estimate=predict(n4_djf,newdata=test_dat[season=="DJF"])) -> gofs[epoch=='all'&season=='DJF']$rmse
yardstick::rmse_vec(test_dat[season=='MAM']$ndvi_3mo,
                    estimate=predict(n4_mam,newdata=test_dat[season=="MAM"])) -> gofs[epoch=='all'&season=='MAM']$rmse
yardstick::rmse_vec(test_dat[season=='JJA']$ndvi_3mo,
                    estimate=predict(n4_jja,newdata=test_dat[season=="JJA"])) -> gofs[epoch=='all'&season=='JJA']$rmse

# just avhrr record
yardstick::rsq_trad_vec(test_dat[season=='SON']$ndvi_3mo,
                        estimate=predict(e1_n4_son,newdata=test_dat[season=="SON"])) -> gofs[epoch=='avhrr'&season=='SON']$r2
yardstick::rsq_trad_vec(test_dat[season=='DJF']$ndvi_3mo,
                        estimate=predict(e1_n4_djf,newdata=test_dat[season=="DJF"])) -> gofs[epoch=='avhrr'&season=='DJF']$r2
yardstick::rsq_trad_vec(test_dat[season=='MAM']$ndvi_3mo,
                        estimate=predict(e1_n4_mam,newdata=test_dat[season=="MAM"])) -> gofs[epoch=='avhrr'&season=='MAM']$r2
yardstick::rsq_trad_vec(test_dat[season=='JJA']$ndvi_3mo,
                        estimate=predict(e1_n4_jja,newdata=test_dat[season=="JJA"])) -> gofs[epoch=='avhrr'&season=='JJA']$r2

yardstick::rmse_vec(test_dat[season=='SON']$ndvi_3mo,
                    estimate=predict(e1_n4_son,newdata=test_dat[season=="SON"]))  -> gofs[epoch=='avhrr'&season=='SON']$rmse
yardstick::rmse_vec(test_dat[season=='DJF']$ndvi_3mo,
                    estimate=predict(e1_n4_djf,newdata=test_dat[season=="DJF"])) -> gofs[epoch=='avhrr'&season=='DJF']$rmse
yardstick::rmse_vec(test_dat[season=='MAM']$ndvi_3mo,
                    estimate=predict(e1_n4_mam,newdata=test_dat[season=="MAM"])) -> gofs[epoch=='avhrr'&season=='MAM']$rmse
yardstick::rmse_vec(test_dat[season=='JJA']$ndvi_3mo,
                    estimate=predict(e1_n4_jja,newdata=test_dat[season=="JJA"])) -> gofs[epoch=='avhrr'&season=='JJA']$rmse

# just avhrr record
yardstick::rsq_trad_vec(test_dat[season=='SON']$ndvi_3mo,
                        estimate=predict(e2_n4_son,newdata=test_dat[season=="SON"])) -> gofs[epoch=='modis'&season=='SON']$r2
yardstick::rsq_trad_vec(test_dat[season=='DJF']$ndvi_3mo,
                        estimate=predict(e2_n4_djf,newdata=test_dat[season=="DJF"])) -> gofs[epoch=='modis'&season=='DJF']$r2
yardstick::rsq_trad_vec(test_dat[season=='MAM']$ndvi_3mo,
                        estimate=predict(e2_n4_mam,newdata=test_dat[season=="MAM"])) -> gofs[epoch=='modis'&season=='MAM']$r2
yardstick::rsq_trad_vec(test_dat[season=='JJA']$ndvi_3mo,
                        estimate=predict(e2_n4_jja,newdata=test_dat[season=="JJA"])) -> gofs[epoch=='modis'&season=='JJA']$r2

yardstick::rmse_vec(test_dat[season=='SON']$ndvi_3mo,
                    estimate=predict(e2_n4_son,newdata=test_dat[season=="SON"]))  -> gofs[epoch=='modis'&season=='SON']$rmse
yardstick::rmse_vec(test_dat[season=='DJF']$ndvi_3mo,
                    estimate=predict(e2_n4_djf,newdata=test_dat[season=="DJF"])) -> gofs[epoch=='modis'&season=='DJF']$rmse
yardstick::rmse_vec(test_dat[season=='MAM']$ndvi_3mo,
                    estimate=predict(e2_n4_mam,newdata=test_dat[season=="MAM"])) -> gofs[epoch=='modis'&season=='MAM']$rmse
yardstick::rmse_vec(test_dat[season=='JJA']$ndvi_3mo,
                    estimate=predict(e2_n4_jja,newdata=test_dat[season=="JJA"])) -> gofs[epoch=='modis'&season=='JJA']$rmse

gofs <- gofs %>% mutate(epoch = case_when(epoch=='all'~'Merged 1982-2019', 
                                  epoch=='modis'~'MODIS 2001-2019',
                                  epoch=='avhrr'~'AVHRR 1982-2000')) %>% 
  as.data.table()

# Predict Models ----------------------------------------------------------
n4_preds <- expand_grid(season=unique(train_dat$season),
                        co2 = seq(min(dat$co2_int),max(dat$co2_int),length.out=100),
                        mape = seq(0.05,1.5,length.out = 200), 
                        pct_anom = c(0), 
                        epoch = 2) %>% 
  mutate(pe_anom_12mo = 0.01*pct_anom*mape) %>%
  # mutate(pe_12mo = pe_anom_12mo+mape) %>% 
  mutate(cco2 = co2-center_co2)

n4_preds <- expand_grid(season=unique(train_dat$season),
                        co2 = seq(min(dat$co2_int),max(dat$co2_int),length.out=100),
                        mape = seq(0.05,1.5,length.out = 200), 
                        pct_anom = c(0), 
                        epoch = 2) %>% 
  mutate(pe_anom_12mo = 0.01*pct_anom*mape) %>%
  # mutate(pe_12mo = pe_anom_12mo+mape) %>% 
  mutate(cco2 = co2-center_co2)

n4_preds <- bind_rows(
  n4_preds %>% filter(season=='SON') %>% mutate(pred = predict(n4_son, newdata=.)),
  n4_preds %>% filter(season=='DJF') %>% mutate(pred = predict(n4_djf, newdata=.)),
  n4_preds %>% filter(season=='MAM') %>% mutate(pred = predict(n4_mam, newdata=.)),
  n4_preds %>% filter(season=='JJA') %>% mutate(pred = predict(n4_jja, newdata=.)))



e1_n4_preds <- expand_grid(season=unique(train_dat$season),
                           co2 = seq(min(dat$co2_int),
                                     unique(train_dat[date==ymd("2000-12-01")]$co2_int),
                                     length.out=100),
                           mape = seq(0.05,1.5,length.out = 200), 
                           pct_anom = c(0), 
                           epoch = 2) %>% 
  mutate(pe_anom_12mo = 0.01*pct_anom*mape) %>%
  # mutate(pe_12mo = pe_anom_12mo+mape) %>% 
  mutate(cco2 = co2-center_co2)
e2_n4_preds <- expand_grid(season=unique(train_dat$season),
                           co2 = seq(unique(train_dat[date==ymd("2001-01-01")]$co2_int),
                                     max(dat$co2_int),length.out=100),
                           mape = seq(0.05,1.5,length.out = 200), 
                           pct_anom = c(0), 
                           epoch = 2) %>% 
  mutate(pe_anom_12mo = 0.01*pct_anom*mape) %>%
  # mutate(pe_12mo = pe_anom_12mo+mape) %>% 
  mutate(cco2 = co2-center_co2)

n4_preds <- bind_rows(
  n4_preds %>% filter(season=='SON') %>% mutate(pred = predict(n4_son, newdata=.)),
  n4_preds %>% filter(season=='DJF') %>% mutate(pred = predict(n4_djf, newdata=.)),
  n4_preds %>% filter(season=='MAM') %>% mutate(pred = predict(n4_mam, newdata=.)),
  n4_preds %>% filter(season=='JJA') %>% mutate(pred = predict(n4_jja, newdata=.))) %>% 
  mutate(epoch = "1982-2019")
e1_n4_preds <- bind_rows(
  e1_n4_preds %>% filter(season=='SON') %>% mutate(pred = predict(e1_n4_son, newdata=.)),
  e1_n4_preds %>% filter(season=='DJF') %>% mutate(pred = predict(e1_n4_djf, newdata=.)),
  e1_n4_preds %>% filter(season=='MAM') %>% mutate(pred = predict(e1_n4_mam, newdata=.)),
  e1_n4_preds %>% filter(season=='JJA') %>% mutate(pred = predict(e1_n4_jja, newdata=.))) %>% 
  mutate(epoch = "AVHRR 1982-2000")
e2_n4_preds <- bind_rows(
  e2_n4_preds %>% filter(season=='SON') %>% mutate(pred = predict(e2_n4_son, newdata=.)),
  e2_n4_preds %>% filter(season=='DJF') %>% mutate(pred = predict(e2_n4_djf, newdata=.)),
  e2_n4_preds %>% filter(season=='MAM') %>% mutate(pred = predict(e2_n4_mam, newdata=.)),
  e2_n4_preds %>% filter(season=='JJA') %>% mutate(pred = predict(e2_n4_jja, newdata=.))) %>% 
  mutate(epoch = "MODIS 2001-2019")

vec_labels <- c("-50"=" -50% P:PET Anom. ",
                "0"=' 0% P:PET Anom. ',
                "50"=" +50% P:PET Anom. ")

p4_ndvi_epoch <- bind_rows(n4_preds, e1_n4_preds, e2_n4_preds) %>% 
  mutate(epoch=factor(epoch,levels=c("AVHRR 1982-2000","MODIS 2001-2019","1982-2019"), 
                      ordered=T)) %>% 
  ggplot(data=., aes(mape,pred,color=(co2), group=co2))+
  geom_line(alpha=1)+
  scale_color_viridis_c(expression(paste(CO[2]~ppm)), option='B',end=0.85)+
  scale_x_continuous(limits=c(0.08,1.5),
                     breaks=c(0,0.5,1,1.5),
                     labels = c(0,0.5,1,1.5),
                     expand=c(0,0),
                     guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
  )+
  scale_y_continuous(limits=c(0,0.9),expand=c(0,0.01))+
  labs(x=expression(paste("Mean Annual P:PET")),
       y=expression(paste(NDVI["3 mo"])))+
  facet_grid(season~epoch, labeller = labeller(pct_anom=vec_labels))+
  theme_linedraw()+
  # guides(color=guide_colorbar(title.position = 'top'))+
  theme(#panel.grid = element_blank(),
    # panel.spacing.x = unit(6, "mm"),
    axis.text = element_text(size=10),
    # axis.text.x = element_text(angle=45, vjust=-0.5),
    plot.margin = margin(t = 0.1,r = 0.4,b = 0.1,l = 0.4,'cm'),
    # legend.position = c(0.525,0.175), 
    legend.position = 'bottom',
    legend.key.width = unit(1,'cm'),
    legend.key.height = unit(0.2,'cm'),
    legend.direction = 'horizontal', 
    legend.background = element_rect(fill=NA)); p4_ndvi_epoch
# ggsave(p4_ndvi_epoch, 
#        filename = 'figures/n4_ndvi_season_by_epoch_weibull_ppet_x_co2.png',
#        width = 16, height = 16, units='cm', dpi=350, type='cairo')

# End Section ******************************************************************

vec_labels <- c("-50"=" -50% P:PET Anom. ",
                "0"=' 0% P:PET Anom. ',
                "50"=" +50% P:PET Anom. ")

library(patchwork)
(n4_preds %>% 
    mutate(epoch='all') %>% 
    ggplot(data=., aes(mape,pred,color=(co2), group=co2))+
    geom_line(alpha=0.5)+
    geom_text(data=gofs[epoch=='Merged 1982-2019'] %>% 
                mutate(r2 = paste0("R**2 ==",format(r2,digits = 2))) %>% 
        as_tibble(), 
      inherit.aes = F, 
              aes(x=1.075,y=0.2,
                  label=r2), 
              parse=TRUE, size=2.5)+
    geom_text(data=gofs[epoch=='Merged 1982-2019'] %>% 
                mutate(rmse = paste0("RMSE ==",format(rmse,digits = 2))) %>% 
        as_tibble(), inherit.aes = F, 
              aes(x=0.925,y=0.1,
                  label=rmse), 
              parse=TRUE, size=2.5)+
    scale_alpha_discrete("epoch", range=c(0.2,0.6))+
    scale_color_viridis_c(expression(paste(CO[2]~ppm)), option='B',end=0.85, 
                          limits=c(335,420))+
    scale_x_continuous(limits=c(0.08,1.5),
                       breaks=c(0,0.5,1,1.5),
                       labels = c(0,0.5,1,1.5),
                       expand=c(0,0),
                       guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
    )+
    scale_y_continuous(limits=c(0,0.9),expand=c(0,0.01))+
    labs(x=NULL, #expression(paste("Mean Annual P:PET")),
         y=expression(paste(NDVI)), 
         subtitle="Merged 1982-2019")+
    facet_grid(~season, labeller = labeller(pct_anom=vec_labels))+
    theme_linedraw()+
    # guides(color=guide_colorbar(title.position = 'top'))+
    theme(#panel.grid = element_blank(),
      # panel.spacing.x = unit(6, "mm"),
      axis.text = element_text(size=10),
      panel.grid.minor = element_blank(),
      # axis.text.x = element_text(angle=45, vjust=-0.5),
      plot.margin = margin(t = 0.1,r = 0.4,b = 0.1,l = 0.4,'cm'),
      # legend.position = c(0.525,0.175), 
      legend.position = 'none', #c(0.99,0.01),
      # legend.justification = c(0.99,0.01),
      # legend.key.width = unit(0.5,'cm'),
      # legend.key.height = unit(0.2,'cm'),
      # legend.direction = 'vertical', 
      legend.background = element_rect(fill=NA)))/(e2_n4_preds %>% 
   ggplot(data=., aes(mape,pred,color=(co2), group=co2))+
   geom_line(data=e1_n4_preds, aes(mape, pred, group=co2), col='gray70',alpha=0.05)+
   geom_line(alpha=0.25)+
     geom_text(data=gofs[epoch=='MODIS 2001-2019'] %>% 
                 as_tibble() %>% 
                 mutate(r2 = paste0("R**2 ==",format(r2,digits = 2))), inherit.aes = F, 
               aes(x=1.075,y=0.2,
                   label=r2), 
               parse=TRUE, size=2.5)+
     geom_text(data=gofs[epoch=='MODIS 2001-2019'] %>% 
                 mutate(rmse = paste0("RMSE ==",format(rmse,digits = 2))) %>% 
         as_tibble(), inherit.aes = F, 
               aes(x=0.925,y=0.1,
                   label=rmse), 
               parse=TRUE, size=2.5)+
   scale_alpha_discrete("epoch", range=c(0.2,0.6))+
   scale_color_viridis_c(expression(atop(paste(CO[2]),"(ppm)")), option='B',end=0.85, 
                         limits=c(335,420))+
   scale_x_continuous(limits=c(0.08,1.5),
                      breaks=c(0,0.5,1,1.5),
                      labels = c(0,0.5,1,1.5),
                      expand=c(0,0),
                      guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
   )+
   scale_y_continuous(limits=c(0,0.9),expand=c(0,0.01))+
   labs(x=NULL, #expression(paste("Mean Annual P:PET")),
        y=expression(paste(NDVI)), 
        subtitle="MODIS 2001-2019")+
   facet_grid(~season, labeller = labeller(pct_anom=vec_labels))+
   theme_linedraw()+
   # guides(color=guide_colorbar(title.position = 'top'))+
   theme(#panel.grid = element_blank(),
     # panel.spacing.x = unit(6, "mm"),
     axis.text = element_text(size=10),
     panel.grid.minor = element_blank(),
     # axis.text.x = element_text(angle=45, vjust=-0.5),
     plot.margin = margin(t = 0.1,r = 0.4,b = 0.1,l = 0.4,'cm'),
     # legend.position = c(0.525,0.175), 
     legend.position = 'right',
     # legend.position = c(0.99,0.01),
     # legend.justification = c(0.99,0.01),
     legend.key.width = unit(0.25,'cm'),
     legend.key.height = unit(2,'cm'),
     # legend.direction = 'horizontal', 
     legend.background = element_rect(fill=NA)))/(e1_n4_preds %>% 
  ggplot(data=., aes(mape,pred,color=(co2), group=co2))+
  geom_line(data=e2_n4_preds, aes(mape, pred, group=co2), col='gray70',alpha=0.05)+
  geom_line(alpha=0.25)+
    geom_text(data=gofs[epoch=='AVHRR 1982-2000'] %>% 
                mutate(r2 = paste0("R**2 ==",format(r2,digits = 2))) %>% 
        as_tibble(), inherit.aes = F, 
              aes(x=1.075,y=0.2,
                  label=r2), 
              parse=TRUE, size=2.5)+
    geom_text(data=gofs[epoch=='AVHRR 1982-2000'] %>% 
                mutate(rmse = paste0("RMSE ==",format(rmse,digits = 2))) %>% 
        as_tibble(), inherit.aes = F, 
              aes(x=0.925,y=0.1,
                  label=rmse), 
              parse=TRUE, size=2.5)+
  scale_alpha_discrete("epoch", range=c(0.2,0.6))+
  scale_color_viridis_c(expression(paste(CO[2]~ppm)), option='B',end=0.85, 
                        limits=c(335,420))+
  scale_x_continuous(limits=c(0.08,1.5),
                     breaks=c(0,0.5,1,1.5),
                     labels = c(0,0.5,1,1.5),
                     expand=c(0,0),
                     guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
  )+
  scale_y_continuous(limits=c(0,0.9),expand=c(0,0.01))+
  labs(x=expression(paste("Mean Annual P:PET")),
       y=expression(paste(NDVI)), 
       subtitle="AVHRR 1982-2000")+
  facet_grid(~season, labeller = labeller(pct_anom=vec_labels))+
  theme_linedraw()+
  # guides(color=guide_colorbar(title.position = 'top'))+
  theme(#panel.grid = element_blank(),
    # panel.spacing.x = unit(6, "mm"),
    axis.text = element_text(size=10),
    panel.grid.minor = element_blank(),
    # axis.text.x = element_text(angle=45, vjust=-0.5),
    plot.margin = margin(t = 0.1,r = 0.4,b = 0.1,l = 0.4,'cm'),
    # legend.position = c(0.525,0.175), 
    # legend.position = c(0.99,0.01),
    # legend.justification = c(0.99,0.01),
    legend.position = 'none',
    # legend.key.width = unit(0.5,'cm'),
    # legend.key.height = unit(0.2,'cm'),
    # legend.direction = 'horizontal', 
    legend.background = element_rect(fill=NA)))+
  plot_layout(guides='collect')+
  plot_annotation(tag_levels = 'a', 
                  tag_prefix = '(',
                  tag_suffix = ')')
# former 16*16 dims
# ggsave(filename = 'figures/Fig5_n4_ndvi_season_by_epoch_gray_overlay_weibull_ppet_x_co2.png',
#        width = 14, height = 16, units='cm', dpi=350, type='cairo')



# Map of CO2 effect
# n4_preds <- expand_grid(season=unique(train_dat$season),
#                         co2 = seq(min(dat$co2_int),max(dat$co2_int),length.out=100),
#                         mape = seq(0.05,1.5,length.out = 200), 
#                         pct_anom = c(0), 
#                         epoch = 2) %>% 
#   mutate(pe_anom_12mo = 0.01*pct_anom*mape) %>%
#   # mutate(pe_12mo = pe_anom_12mo+mape) %>% 
#   mutate(cco2 = co2-center_co2)


early_andvi <- dat %>% lazy_dt() %>% 
  filter(date>=ymd("1981-08-01") & 
           date<=ymd("1986-08-01")) %>% 
  group_by(x,y) %>% 
  summarize(early_ndvi = mean(ndvi_hyb,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()


oz_poly <- sf::read_sf("../data_general/GADM/gadm36_AUS.gpkg", 
                       layer="gadm36_AUS_1")
oz_poly <- sf::st_as_sf(oz_poly)
oz_poly <- sf::st_simplify(oz_poly, dTolerance = 1000)

colorspace::hcl_palettes() 
colorspace::hcl_palettes() %>% plot
dat[date==ymd('2001-01-01')] %>% 
  as_tibble() %>% 
  select(x,y,mape) %>% 
  merge(., early_andvi,by=c("x","y")) %>% 
  # expand_grid(., cco2 = c(max(dat$co2_int)-center_co2, 
  #                        min(dat$co2_int)-center_co2)) %>% 
  expand_grid(., cco2 = c(42.4, 
                          -34)) %>% 
  mutate(pe_anom_12mo=0, 
         epoch=1) %>% 
  mutate(pred_son = predict(n4_son,newdata=.)) %>% 
  mutate(pred_djf = predict(n4_djf,newdata=.)) %>% 
  mutate(pred_mam = predict(n4_mam,newdata=.)) %>% 
  mutate(pred_jja = predict(n4_jja,newdata=.)) %>% 
  rowwise() %>% 
  mutate(pred = mean(c(pred_son,pred_djf,pred_mam,pred_jja),na.rm=T)) %>% 
  select(x,y,pred,cco2) %>% 
  pivot_wider(id_cols=c('x','y'), names_from=cco2, values_from=pred) %>% 
  mutate(diff = 100*((`42.4`/`-34`)-1) ) %>% #pull(diff) %>% quantile(., c(0.1,0.5, 0.9))
  ggplot(data=.,aes(x,y,fill=diff))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-44,-10))+
  labs(x=NULL,y=NULL, 
       fill=expression(paste(CO[2]~Delta*NDVI)))+
  colorspace::scale_fill_binned_sequential(
    palette='Terrain 2', 
    limits=c(5,20),
    n.breaks=7, rev=T)+
  geom_tile()+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=7),
        legend.position = c(1,1), 
        legend.title.align = 1,
        legend.key.width = unit(0.33,'cm'),
        legend.justification = c(1,1))
# ggsave(filename = 'figures/Fig5_n4_ndvi_map_weibull_ppet_x_co2_modeled_effect.png',
#        width = 7, height = 16, units='cm', dpi=350, type='cairo')



# Plot Joint Map + response figures -------------------------------------------------
p1 <- (dat[date==ymd('2001-01-01')] %>% 
         as_tibble() %>% 
         select(x,y,mape) %>% 
         merge(., early_andvi,by=c("x","y")) %>% 
         expand_grid(., cco2 = c(42.4, 
                                 -34)) %>% 
         mutate(pe_anom_12mo=0, 
                epoch=1) %>% 
         mutate(pred_son = predict(n4_son,newdata=.)) %>% 
         mutate(pred_djf = predict(n4_djf,newdata=.)) %>% 
         mutate(pred_mam = predict(n4_mam,newdata=.)) %>% 
         mutate(pred_jja = predict(n4_jja,newdata=.)) %>% 
         rowwise() %>% 
         mutate(pred = mean(c(pred_son,pred_djf,pred_mam,pred_jja),na.rm=T)) %>% 
         select(x,y,pred,cco2) %>% 
         pivot_wider(id_cols=c('x','y'), names_from=cco2, values_from=pred) %>% 
         mutate(diff = 100*((`42.4`/`-34`)-1) ) %>% #pull(diff) %>% quantile(., c(0.1,0.5, 0.9))
         ggplot(data=.,aes(x,y,fill=diff))+
         geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
         geom_tile()+
         coord_equal()+
         scale_x_continuous(breaks=seq(140,154,by=5))+
         coord_sf(xlim = c(140,154),
                  ylim = c(-44,-10))+
         labs(x=NULL,y=NULL, 
              fill=expression(paste(CO[2]~Delta*NDVI~('%'))))+
         colorspace::scale_fill_binned_sequential(
           palette='Terrain 2', 
           limits=c(5,20),
           n.breaks=7, rev=T)+
         geom_tile()+
         theme(panel.background = element_rect(fill='lightblue'),
               panel.grid = element_blank(), 
               legend.title = element_text(size=7),
               legend.position = c(1,1), 
               legend.title.align = 1,
               legend.key.width = unit(0.33,'cm'),
               legend.justification = c(1,1)))


p2 <- (n4_preds %>% 
         mutate(epoch='all') %>% 
         ggplot(data=., aes(mape,pred,color=(co2), group=co2))+
         geom_line(alpha=0.5)+
         geom_text(data=gofs[epoch=='Merged 1982-2019'] %>% as_tibble() %>% 
                     mutate(r2 = paste0("R**2 ==",format(r2,digits = 2))), inherit.aes = F, 
                   aes(x=1.075,y=0.2,
                       label=r2), 
                   parse=TRUE, size=3.5)+
         geom_text(data=gofs[epoch=='Merged 1982-2019'] %>% as_tibble() %>% 
                     mutate(rmse = paste0("RMSE ==",format(rmse,digits = 2))), inherit.aes = F, 
                   aes(x=0.925,y=0.1,
                       label=rmse), 
                   parse=TRUE, size=3.5)+
         scale_alpha_discrete("epoch", range=c(0.2,0.6))+
         scale_color_viridis_c(expression(paste(CO[2]~ppm)), option='B',end=0.85, 
                               limits=c(335,420))+
         scale_x_continuous(limits=c(0.08,1.5),
                            breaks=c(0,0.5,1,1.5),
                            labels = c(0,0.5,1,1.5),
                            expand=c(0,0),
                            guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
         )+
         scale_y_continuous(limits=c(0,0.9),expand=c(0,0.01))+
         labs(x=NULL, #expression(paste("Mean Annual P:PET")),
              y=expression(paste(NDVI)), 
              subtitle="Merged 1982-2019")+
         facet_grid(~season, labeller = labeller(pct_anom=vec_labels))+
         theme_linedraw()+
         # guides(color=guide_colorbar(title.position = 'top'))+
         theme(#panel.grid = element_blank(),
           # panel.spacing.x = unit(6, "mm"),
           axis.text = element_text(size=10),
           panel.grid.minor = element_blank(),
           # axis.text.x = element_text(angle=45, vjust=-0.5),
           plot.margin = margin(t = 0.1,r = 0.4,b = 0.1,l = 0.4,'cm'),
           # legend.position = c(0.525,0.175), 
           legend.position = 'none', #c(0.99,0.01),
           # legend.justification = c(0.99,0.01),
           # legend.key.width = unit(0.5,'cm'),
           # legend.key.height = unit(0.2,'cm'),
           # legend.direction = 'vertical', 
           legend.background = element_rect(fill=NA)))


p3 <- (e2_n4_preds %>% 
         ggplot(data=., aes(mape,pred,color=(co2), group=co2))+
         geom_line(data=e1_n4_preds, aes(mape, pred, group=co2), col='gray70',alpha=0.05)+
         geom_line(alpha=0.25)+
         geom_text(data=gofs[epoch=='MODIS 2001-2019'] %>% as_tibble() %>% 
                     mutate(r2 = paste0("R**2 ==",format(r2,digits = 2))), inherit.aes = F, 
                   aes(x=1.075,y=0.2,
                       label=r2), 
                   parse=TRUE, size=3.5)+
         geom_text(data=gofs[epoch=='MODIS 2001-2019'] %>% as_tibble() %>% 
                     mutate(rmse = paste0("RMSE ==",format(rmse,digits = 2))), inherit.aes = F, 
                   aes(x=0.925,y=0.1,
                       label=rmse), 
                   parse=TRUE, size=3.5)+
         scale_alpha_discrete("epoch", range=c(0.2,0.6))+
         scale_color_viridis_c(expression(atop(paste(CO[2]),"(ppm)")), option='B',end=0.85, 
                               limits=c(335,420))+
         scale_x_continuous(limits=c(0.08,1.5),
                            breaks=c(0,0.5,1,1.5),
                            labels = c(0,0.5,1,1.5),
                            expand=c(0,0),
                            guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
         )+
         scale_y_continuous(limits=c(0,0.9),expand=c(0,0.01))+
         labs(x=NULL, #expression(paste("Mean Annual P:PET")),
              y=expression(paste(NDVI)), 
              subtitle="MODIS 2001-2019")+
         facet_grid(~season, labeller = labeller(pct_anom=vec_labels))+
         theme_linedraw()+
         theme(#panel.grid = element_blank(),
           # panel.spacing.x = unit(6, "mm"),
           axis.text = element_text(size=10),
           panel.grid.minor = element_blank(),
           # axis.text.x = element_text(angle=45, vjust=-0.5),
           plot.margin = margin(t = 0.1,r = 0.4,b = 0.1,l = 0.4,'cm'),
           # legend.position = c(0.525,0.175), 
           # legend.position = c(0.99,0.01),
           # legend.justification = c(0.99,0.01),
           legend.position = 'none',
           # legend.key.width = unit(0.5,'cm'),
           # legend.key.height = unit(0.2,'cm'),
           # legend.direction = 'horizontal', 
           legend.background = element_rect(fill=NA)))



p4 <- (e1_n4_preds %>% 
         ggplot(data=., aes(mape,pred,color=(co2), group=co2))+
         geom_line(data=e2_n4_preds, aes(mape, pred, group=co2), col='gray70',alpha=0.05)+
         geom_line(alpha=0.25)+
         geom_text(data=gofs[epoch=='AVHRR 1982-2000'] %>% as_tibble() %>% 
                     mutate(r2 = paste0("R**2 ==",format(r2,digits = 2))), inherit.aes = F, 
                   aes(x=1.075,y=0.2,
                       label=r2), 
                   parse=TRUE, size=3.5)+
         geom_text(data=gofs[epoch=='AVHRR 1982-2000'] %>% as_tibble() %>% 
                     mutate(rmse = paste0("RMSE ==",format(rmse,digits = 2))), inherit.aes = F, 
                   aes(x=0.925,y=0.1,
                       label=rmse), 
                   parse=TRUE, size=3.5)+
         scale_alpha_discrete("epoch", range=c(0.2,0.6))+
         scale_color_viridis_c(expression(paste(CO[2]~ppm)), option='B',end=0.85, 
                               limits=c(335,420))+
         scale_x_continuous(limits=c(0.08,1.5),
                            breaks=c(0,0.5,1,1.5),
                            labels = c(0,0.5,1,1.5),
                            expand=c(0,0),
                            guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
         )+
         scale_y_continuous(limits=c(0,0.9),expand=c(0,0.01))+
         labs(x=expression(paste("Mean Annual P:PET")),
              y=expression(paste(NDVI)), 
              subtitle="AVHRR 1982-2000")+
         facet_grid(~season, labeller = labeller(pct_anom=vec_labels))+
         theme_linedraw()+
         theme(
           # plot.subtitle = element_text(size=10),
           #panel.grid = element_blank(),
           # panel.spacing.x = unit(6, "mm"),
           axis.text = element_text(size=10),
           panel.grid.minor = element_blank(),
           # axis.text.x = element_text(angle=45, vjust=-0.5),
           plot.margin = margin(t = 0.1,r = 0.4,b = 0.1,l = 0.4,'cm'),
           # legend.position = c(0.525,0.175), 
           legend.position = 'bottom',
           # legend.position = c(0.99,0.01),
           # legend.justification = c(0.99,0.01),
           legend.key.width = unit(2,'cm'),
           legend.key.height = unit(0.25,'cm'),
           # legend.direction = 'horizontal', 
           legend.background = element_rect(fill=NA)))

# library(cowplot)
# p_left <- p1+draw_label(label='(a)', x=140.95,y=-11,size = 20)
# p_top <- ggdraw(p2)+draw_label(label='(b)',x=0.03,y=0.93,size=20)
# p_mid <- ggdraw(p3)+draw_label(label='(c)', x=0,y=0.96, size=20,hjust=0.1,vjust=0.1)
# p_bot <- ggdraw(p4)+draw_label(label = '(d)', x=0.035,y=0.94,size=20)
# cp_r <- cowplot::plot_grid(p_top,p_mid,p_bot,
#                            nrow = 3,
#                            rel_heights = c(1,1,1.4))
# cp_out <- cowplot::plot_grid(p_left,cp_r, 
#                              rel_widths = c(1,2), rel_heights = c(1,1),
#                              greedy = F)
# ggsave(plot=cp_out,
#        filename = "figures/Fig5_n4_ndvi_map_weibull_ppet_x_co2_modeled_effect.png",
#        width = 25, height=20, units='cm', dpi=350, type='cairo')

library(patchwork)

p_out <- (p1|(p2/p3/p4)) + 
  plot_annotation(tag_prefix = '(',tag_levels = 'a',tag_suffix = ')')&
  theme(plot.margin=margin(t=0,r=5,b=0,l=0))
p_out
ggsave(plot=p_out,
       filename = "figures/Fig5_n4_ndvi_map_weibull_ppet_x_co2_modeled_effect.png",
       width = 25, height=20, units='cm', dpi=350, device=grDevices::png)



