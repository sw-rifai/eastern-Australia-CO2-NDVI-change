library(tidyverse)
library(data.table); setDTthreads(threads = 0)
library(lubridate); 
library(dtplyr);
library(RcppArmadillo)
library(sf); library(stars)
library(patchwork); 
library(zyp)
set.seed(333)
# IMPORT ###################################################################
# load("data/gridCell_lm_ndvi_clim.Rdata") # grid cell linear regressions
oz_poly <- sf::read_sf("../data_general/GADM/gadm36_AUS.gpkg", 
                       layer="gadm36_AUS_1")
oz_poly <- st_as_sf(oz_poly)
oz_poly <- st_simplify(oz_poly, dTolerance = 0.05)

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

lvi <- lazy_dt(vi)

lvi %>% group_by(x,y) %>% summarize(ndvi_u =mean(ndvi_hyb,na.rm=TRUE)) %>% show_query()

norms_vi <- vi[,`:=`(month=month(date))] %>% 
  .[,.(ndvi_u = mean(ndvi_hyb,na.rm=TRUE), 
       ndvi_sd = sd(ndvi_hyb,na.rm=TRUE)),keyby=.(x,y,month)]

vi <- norms_vi[vi,on=.(x,y,month)] %>% 
  .[,`:=`(ndvi_anom = ndvi_hyb - ndvi_u)] %>% 
  .[,`:=`(ndvi_anom_sd = ndvi_anom/ndvi_sd)]

vi <- vi %>% select(-season)

dat <- arrow::read_parquet("/home/sami/scratch/ARD_ndvi_aclim_anoms.parquet",
                           col_select = c(
                             "date", "hydro_year", "id",
                             # "season",
                             "precip",  
                             "precip_anom",
                             "precip_anom_3mo",
                             "precip_anom_36mo",
                             "precip_anom_12mo",
                             "map", 
                             "precip_12mo",
                             # "precip_36mo",
                             # "tmax","tmax_u",
                             # "tmax_anom","tmax_anom_sd", "matmax",
                             # "tmin","tmin_anom",
                             "vpd15","vpd15_anom","vpd15_anom_sd","mavpd15",
                             # "vpd15_u",
                             "pet","mapet",
                             "pet_anom","pet_anom_3mo","pet_u","pet_sd",
                             "pet_anom_sd",
                             "pet_12mo",
                             # "pet_36mo",
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
  .[is.infinite(mape)==F] %>% 
  .[,`:=`(p_anom_12mo_frac = precip_anom_12mo/map)]

dat <- dat[,`:=`(#pe_3mo = precip_3mo/pet_3mo, 
  pe_12mo = precip_12mo/pet_12mo 
  # pe_24mo = precip_24mo/pet_24mo, 
  # pe_36mo = precip_36mo/pet_36mo, 
  # pe_48mo = precip_48mo/pet_48mo
)]
dim(dat)
dat <- merge(dat, 
             vi,
             by=c("x","y","date","month","year"), 
             all=TRUE,allow.cartesian=TRUE)
dat <- dat[order(x,y,date)][, ndvi_3mo := frollmean(ndvi_hyb,n = 3,fill = NA,align='center',na.rm=TRUE), by=.(x,y)]

rm(vi); gc(full=TRUE)

mlo <- readr::read_table("../data_general/CO2_growth_rate/co2_mm_mlo_20200405.txt", 
                         skip = 72, col_names = F) %>% 
  set_names(
    c("year","month","ddate","co2_avg","co2_int","co2_trend","ndays")
  ) %>% 
  mutate(date = ymd(paste(year,month,1))) %>% 
  select(date,co2_int,co2_trend) %>% 
  as.data.table()
dat <- merge(mlo,dat,by="date")
# dat <- dat[is.na(ndvi_3mo)==F & is.na(co2_int)==F]
center_co2 <- mean(dat$co2_int)
dat <- dat[,`:=`(cco2=co2_int-center_co2)]
gc()
dat <- dat[is.na(vc)==F]
dat <- dat[str_detect(vc,"Forests") |
             str_detect(vc, "Eucalypt") |
             str_detect(vc, "Rainforests")]
dat <- dat[x>= 140] # FILTER TO LON >= 140
# dat <- dat[ndvi_hyb>0][ndvi_anom_sd > -3.5 & ndvi_anom_sd < 3.5]
dat[,`:=`(pe_anom_12mo = pe_12mo - mape)]
dat[,`:=`(epoch = ifelse(date<ymd("2000-12-31"),'avhrr','modis'))]
dat[,`:=`(year=year(date),month=month(date))] %>%
  .[,`:=`(season = case_when(month%in%c(3:5)~'MAM',
                             month%in%c(6:8)~'JJA',
                             month%in%c(9:11)~'SON',
                             month%in%c(12,1,2)~'DJF'))]

dat <- dat %>% mutate(epoch = as_factor(epoch), 
                      season = factor(season, levels=c("SON","DJF","MAM","JJA")))
kop <- arrow::read_parquet("../data_general/Koppen_climate/BOM_Koppen_simplified7.parquet")
kop <- setDT(kop)
kop <- kop[,.(x,y,zone)]
dat <- merge(dat, kop, by=c("x","y"),all = T)
#*******************************************************************************


coords_keep <- dat %>% 
  lazy_dt() %>% 
  group_by(x,y,hydro_year) %>% 
  summarize(nobs = sum(is.na(ndvi_hyb)==F)) %>% 
  ungroup() %>% 
  as.data.table()
coords_keep <- coords_keep %>% 
  filter(nobs >= 6) %>% 
  group_by(x,y) %>% 
  summarize(nobs_annual = n()) %>% 
  ungroup()



# summarize to annual -----------------------------------------------------
dat_annual <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
  .[date>= ymd("1981-12-01") & date<= ymd("2019-08-30")] %>% 
  # .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
  .[,`:=`(epoch=ifelse(hydro_year < 2001,0,1))] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[is.na(pe_anom_12mo)==F] %>% 
  inner_join(., coords_keep,by=c('x','y')) %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982, 
          frac_p_anom = precip_anom_12mo/map, 
          frac_ppet_anom = pe_anom_12mo/mape)] %>% 
  .[,.(x,y,date,ndvi_hyb,
       hydro_year,hydro_year_c,frac_p_anom,frac_ppet_anom,epoch, 
       nobs_annual)] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE), 
       frac_p_anom = mean(frac_p_anom,na.rm=TRUE), 
       frac_ppet_anom = mean(frac_ppet_anom,na.rm=TRUE), 
       epoch = mean(epoch, na.rm=TRUE)), by=.(x,y,hydro_year_c)]
dat_annual <- dat_annual %>% mutate(epoch=as.numeric(epoch))
dat_annual %>% select(ndvi_hyb,hydro_year_c,frac_p_anom,frac_ppet_anom,epoch, x,y) %>% dim
dat_annual %>% select(ndvi_hyb,hydro_year_c,frac_p_anom,frac_ppet_anom,epoch, x,y) %>% 
  distinct() %>% dim


# Regression ---------------------------------------------------------
## RLM on annual ndvi -----------------------------------------------
tmp <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-31")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982, 
          epoch = ifelse(date<=ymd("2000-12-31"),0,1))] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE), 
       epoch=mean(epoch,na.rm=TRUE)),
    by=.(x,y,hydro_year_c)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  lazy_dt() %>% 
  select(x,y,hydro_year_c,ndvi_hyb,epoch) %>% 
  as.data.table() %>% 
  group_by(x,y) %>% 
  mutate(id = cur_group_id()) %>% 
  ungroup()
tmp <- tmp %>% distinct()  
tmp <- tmp %>% as.data.table()

# filter out pixels that only have data for one satellite epoch
vec_rlm <- tmp %>% group_by(x,y,id) %>% 
  summarize(nobs=n()) %>% 
  ungroup() %>% 
  filter(nobs > 19)
tmp <- tmp[id %in% unique(vec_rlm$id)]
rlm_ndvi_annual <-  tmp %>% 
  as.data.table() %>% 
  .[,.(beta = list(coef(MASS::rlm(ndvi_hyb~hydro_year_c+epoch)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3]), by=.(x,y)]
rlm_ndvi_annual

# Summarize annual data by epoch -------------------------------
tmp_ndvi_e1 <- dat %>% 
  .[date >= ymd("1981-12-01") & date <= ymd("2000-11-30")] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982, 
          epoch = ifelse(date<=ymd("2000-12-31"),0,1),
          frac_p = precip_12mo/map,
          frac_p_anom = precip_anom_12mo/map,
          frac_pet_anom = (pet_12mo-mapet)/mapet,
          frac_vpd_anom = vpd15_anom/mavpd15,
          frac_ppet_anom = (pe_12mo-mape)/mape)] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE), 
       co2 = mean(co2_trend,na.rm=TRUE),
       epoch=mean(epoch,na.rm=TRUE), 
       nobs = sum(is.na(ndvi_hyb)==F), 
       p_anom = mean(precip_anom_12mo,na.rm=TRUE),
       frac_p = mean(frac_p, na.rm=TRUE),
       frac_p_anom = mean(frac_p_anom,na.rm=TRUE), 
       frac_vpd_anom = mean(frac_vpd_anom,na.rm=TRUE),
       frac_pet_anom = mean(frac_pet_anom, na.rm=TRUE),
       frac_ppet_anom = mean(frac_ppet_anom,na.rm=TRUE)),
    by=.(x,y,hydro_year, hydro_year_c)] %>% 
  group_by(x,y) %>% 
  mutate(id = cur_group_id()) %>%
  ungroup() %>% 
  as.data.table()

tmp_ndvi_e2 <- dat %>% 
  .[date >= ymd("2001-01-01") & date <= ymd("2019-08-30")] %>% 
  .[,`:=`(hydro_year_c = hydro_year-2001, 
          epoch = ifelse(date<=ymd("2000-12-31"),0,1),
          frac_p = precip_12mo/map,
          frac_p_anom = precip_anom_12mo/map,
          frac_pet_anom = (pet_12mo-mapet)/mapet,
          frac_vpd_anom = vpd15_anom/mavpd15,
          frac_ppet_anom = (pe_12mo-mape)/mape)] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE), 
       co2 = mean(co2_trend,na.rm=TRUE),
       epoch=mean(epoch,na.rm=TRUE), 
       nobs = sum(is.na(ndvi_hyb)==F), 
       p_anom = mean(precip_anom_12mo,na.rm=TRUE),
       frac_p = mean(frac_p, na.rm=TRUE),
       frac_p_anom = mean(frac_p_anom,na.rm=TRUE), 
       frac_vpd_anom = mean(frac_vpd_anom,na.rm=TRUE),
       frac_pet_anom = mean(frac_pet_anom, na.rm=TRUE),
       frac_ppet_anom = mean(frac_ppet_anom,na.rm=TRUE)),
    by=.(x,y,hydro_year, hydro_year_c)] %>% 
  group_by(x,y) %>% 
  mutate(id = cur_group_id()) %>%
  ungroup() %>% 
  as.data.table()

tmp_ndvi <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-30")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982, 
          epoch = ifelse(date<=ymd("2000-12-31"),0,1),
          frac_p = precip_12mo/map,
          frac_p_anom = precip_anom_12mo/map,
          frac_pet_anom = (pet_12mo-mapet)/mapet,
          frac_vpd_anom = vpd15_anom/mavpd15,
          frac_ppet_anom = (pe_12mo-mape)/mape)] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE), 
       co2 = mean(co2_trend,na.rm=TRUE),
       epoch=mean(epoch,na.rm=TRUE), 
       nobs = sum(is.na(ndvi_hyb)==F), 
       p_anom = mean(precip_anom_12mo,na.rm=TRUE),
       frac_p_anom = mean(frac_p_anom,na.rm=TRUE), 
       frac_vpd_anom = mean(frac_vpd_anom,na.rm=TRUE),
       frac_pet_anom = mean(frac_pet_anom, na.rm=TRUE),
       frac_ppet_anom = mean(frac_ppet_anom,na.rm=TRUE)),
    by=.(x,y,hydro_year, hydro_year_c)] %>% 
  group_by(x,y) %>% 
  mutate(id = cur_group_id()) %>%
  ungroup() %>% 
  as.data.table()
tmp_annual_nobs <- tmp_ndvi %>%
  mutate(epoch=round(epoch)) %>% 
  lazy_dt() %>% 
  group_by(id,epoch) %>% 
  summarize(obs_annual = sum(nobs >= 4)) %>% 
  as.data.table()
tmp_epoch_nobs <- tmp_annual_nobs %>% group_by(id) %>% 
  summarize(obs_epoch = sum(obs_annual > 5)) %>% 
  ungroup()
  # pull(obs_epoch) %>% table
tmp_epoch_nobs %>% 
  inner_join(., tmp_ndvi %>% select(x,y,id) %>% distinct(), by='id') %>% 
  ggplot(data=.,aes(x,y,fill=obs_epoch))+geom_tile()+scale_fill_viridis_c()
vec_ids <- tmp_epoch_nobs %>% 
  filter(obs_epoch >= 1.9) %>% 
  pull(id)


rlm_ndvi_annual_co2_ppet_epoch <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,.(beta = list(coef(MASS::rlm(
    ndvi_hyb~
      co2+
      frac_ppet_anom+
      # jitter(frac_vpd_anom)+
      jitter(epoch))))),
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b4=unlist(beta)[4]
  ), by=.(x,y)]
rlm_ndvi_annual_co2_ppet_epoch


## NDVI robust regression ------------------------------------------------------

system.time(
  lt_ndvi_year <-  tmp_ndvi[id%in%vec_ids] %>% #dat_annual %>% #[nobs_annual >= 10] %>% 
    .[,.(beta = list(unname(MASS::rlm(ndvi_hyb~hydro_year_c+epoch, 
                               data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3],b3=unlist(beta)[4]
    ), by=.(x,y)]
)

system.time(
  lt_ndvi_year_noEpoch <-  tmp_ndvi[id%in%vec_ids] %>% #dat_annual %>% #[nobs_annual >= 10] %>% 
    .[,.(beta = list(unname(MASS::rlm(ndvi_hyb~hydro_year_c, 
                                      data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3],b3=unlist(beta)[4]
    ), by=.(x,y)]
)

system.time(
  lt_ndvi_p_ppet_year <-  tmp_ndvi[id%in%vec_ids] %>% #dat_annual %>% #[nobs_annual >= 10] %>% 
    .[,.(beta = list(unname(MASS::rlm(ndvi_hyb~hydro_year_c+frac_p_anom+frac_ppet_anom+epoch, 
                                      data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3],b3=unlist(beta)[4],b4=unlist(beta)[5]
    ), by=.(x,y)]
)

system.time(
  lt_ndvi_ppet_year <-  tmp_ndvi[id%in%vec_ids] %>% #dat_annual %>% #[nobs_annual >= 10] %>% 
    .[,.(beta = list(unname(MASS::rlm(ndvi_hyb~hydro_year_c+frac_ppet_anom+epoch, 
                                      data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3],b3=unlist(beta)[4]
    ), by=.(x,y)]
)


# Regression: NDVI by epoch -----------------------------------------------
system.time(
  lt_ndvi_sen_e1 <- tmp_ndvi_e1[id%in%vec_ids] %>% 
    .[,.(beta = list(unname(zyp.sen(ndvi_hyb~hydro_year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]
    ), by=.(x,y)]
)
system.time(
  lt_ndvi_sen_e2 <- tmp_ndvi_e2[id%in%vec_ids] %>% 
    .[,.(beta = list(unname(zyp.sen(ndvi_hyb~hydro_year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]
    ), by=.(x,y)]
)



## Regression VPD Thiel Sen regression ---------------------------------------
system.time(
  lt_v_sen <- dat[date>=ymd("1981-12-01")][date<=ymd("2019-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[is.na(vpd15)==F] %>% 
    .[,.(vpd15 = mean(vpd15,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(vpd15~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)

## Regression VPD Thiel Sen epoch1 ---------------------------------------
system.time(
  lt_v_sen_e1 <- dat[date>=ymd("1981-12-01")][date<=ymd("2000-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[is.na(vpd15)==F] %>% 
    .[,.(vpd15 = mean(vpd15,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(vpd15~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)

## Regression VPD Thiel Sen epoch 2 ---------------------------------------
system.time(
  lt_v_sen_e2 <- dat[date>=ymd("2000-12-01")][date<=ymd("2019-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-2000)] %>% 
    .[is.na(vpd15)==F] %>% 
    .[,.(vpd15 = mean(vpd15,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(vpd15~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)

## Regression: Precip Thiel Sen --------------------------------------------
system.time(
  lt_p_sen <- dat[date>=ymd("1981-12-01")][date<=ymd("2019-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[is.na(precip_anom_12mo)==F] %>%
    .[,.(p_tot = mean(precip_anom_12mo+map,na.rm=TRUE)), 
      by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(p_tot~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)

## Regression: Precip Thiel Sen by epoch -------------------------------------------------
system.time(
  lt_fp_sen_e1 <- tmp_ndvi_e1[id%in%vec_ids] %>% 
    .[,.(beta = list(unname(zyp.sen(frac_p~hydro_year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]
    ), by=.(x,y)]
)
system.time(
  lt_fp_sen_e2 <- tmp_ndvi_e2[id%in%vec_ids] %>% 
    .[,.(beta = list(unname(zyp.sen(frac_p~hydro_year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]
    ), by=.(x,y)]
)

## Regression: P:PET Thiel Sen -------------------------------------------------
system.time(
  lt_ppet_sen <- dat[date>=ymd("1981-12-01")][date<=ymd("2019-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[is.na(pe_12mo)==F] %>%
    .[,.(ppet = mean(pe_12mo,na.rm=TRUE)), 
      by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(ppet~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)

## Regression: P:PET Thiel Sen epoch 1 -------------------------------------------------
system.time(
  lt_ppet_sen_e1 <- dat[date>=ymd("1981-12-01")][date<=ymd("2000-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[is.na(pe_12mo)==F] %>%
    .[,.(ppet = median(pe_12mo,na.rm=TRUE)), 
      by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(ppet~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)

## Regression: P:PET Thiel Sen epoch 2 -------------------------------------------------
system.time(
  lt_ppet_sen_e2 <- dat[date>=ymd("2000-12-01")][date<=ymd("2019-08-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-2001)] %>% 
    .[is.na(pe_12mo)==F] %>%
    .[,.(ppet = median(pe_12mo,na.rm=TRUE)), 
      by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(ppet~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)



# Percent increase of VPD
lt_v_sen <- lt_v_sen[b0 > 0]
dVPD_VPD_sen <- mean(37*lt_v_sen$b1,na.rm=TRUE)/mean(lt_v_sen$b0,na.rm=TRUE)

# Percent increase of VPD
lt_ppet_sen <- lt_ppet_sen[b0 > 0]
dPPET_PPET_sen <- mean(37*lt_ppet_sen$b1,na.rm=TRUE)/mean(lt_ppet_sen$b0,na.rm=TRUE)

# Actual percent increase in NDVI
lt_ndvi_year <- lt_ndvi_year[b0 > 0]
dNDVI_NDVI_rlm <- mean(37*lt_ndvi_year$b1,na.rm=TRUE)/mean(lt_ndvi_year$b0,na.rm=TRUE)

# Percent increase in NDVI after effects of ppet were removed
lt_ndvi_ppet_year <- lt_ndvi_ppet_year[b0 > 0]
dNDVI_NDVI_rlm_ppet <- mean(37*lt_ndvi_ppet_year$b1,na.rm=TRUE)/mean(lt_ndvi_ppet_year$b0,na.rm=TRUE)

# Percent increase in Ca
dCa_Ca <- 
  diff(range(mlo[date>=ymd("1981-12-01") & date <= ymd("2019-08-30")]$co2_trend))/
  mean(mlo[date>=ymd("1981-12-01") & date <= ymd("1982-11-01")]$co2_trend)
dCa_Ca_e1 <- 
  diff(range(mlo[date>=ymd("1981-12-01") & date <= ymd("2000-11-30")]$co2_trend))/
  mean(mlo[date>=ymd("1981-12-01") & date <= ymd("1982-11-01")]$co2_trend)
dCa_Ca_e2 <- 
  diff(range(mlo[date>=ymd("2000-12-01") & date <= ymd("2019-08-30")]$co2_trend))/
  mean(mlo[date>=ymd("2000-12-01") & date <= ymd("2001-11-01")]$co2_trend)


# Expected WUE related increase (Donohue 2013)
0.5*(dCa_Ca - 0.5*dVPD_VPD_sen)

# Actual percent relative increase in NDVI
dNDVI_NDVI_rlm


dNDVI_NDVI_rlm_ppet

lt_ndvi_year_noEpoch %>% 
  filter(b0>0.05 & b0<1) %>% 
  filter(b1 > -0.05 & b1 <0.05) %>% 
  mutate(rel_pred_ndvi_2019 = (37*b1)/(b0)) %>% pull(rel_pred_ndvi_2019) %>% 
  summary

lt_ndvi_year %>% 
  mutate(rel_pred_ndvi_2019 = (37*b1+1*b2)/(b0)) %>% pull(rel_pred_ndvi_2019) %>% 
  summary
lt_ndvi_year %>% 
  mutate(rel_pred_ndvi_2019 = 100*37*b1/b0) %>%
  ggplot(data=.,aes(x,y,fill=b1))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()

# Plotting ----------------------------------------------------------------
p_vpd_sen <- lt_v_sen %>% 
  as_tibble() %>% 
  filter(between(b1,-0.1,0.1)) %>% 
  filter(b0 > 0) %>% 
  ggplot(data=.,aes(x,y,fill=100*38*b1/b0))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL)+
  scico::scale_fill_scico(expression(paste(Delta,"VPD(%)")),
                          palette ='romaO', direction=-1,
                          limits=c(-20,20), 
                          oob=scales::squish)+
  # scale_fill_viridis_c(expression(paste(Delta,"VPD(%)")),
  #   option='A', limits=c(0,20), oob=scales::squish)+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1)); p_vpd_sen


p_wue_sen <- lt_v_sen %>% 
  as_tibble() %>% 
  filter(b0 > 0) %>% 
  filter(between(b1,-0.05,0.05)) %>% 
  mutate(dVPD_VPD = b1*37/b0) %>% 
  mutate(expectation = 0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
  # pull(expectation) %>% quantile(., c(0.01,0.99))
  ggplot(data=.,aes(x,y,fill=expectation*100))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL)+
  scico::scale_fill_scico(expression(paste(Delta*NDVI[WUE~Pred.]("%"))),
                          palette = 'bamako', 
                          direction = -1,
                          limits=c(5,15), #na.value = 'red',
                          oob=scales::squish
  )+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1)); p_wue_sen

library(colorspace)
hcl_palettes("diverging", n=7, plot=T)
vec_cols <- colorspace::diverge_hcl(11, "Green-Brown",rev = T)
p_diff_sen <- inner_join({lt_ndvi_year_noEpoch %>% as_tibble() %>% 
    filter(between(b0,0.05,1)) %>% 
    mutate(b2=0) %>% 
    mutate(dNDVI = 100*(37*b1+0.5*b2)/(b0+0.5*b2)) %>% 
    select(x,y,dNDVI)},
    {lt_v_sen %>% 
        as_tibble() %>% 
        filter(between(b1,-0.1,0.1)) %>% 
        mutate(dVPD_VPD = b1*37/b0) %>% 
        mutate(expectation = 0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
        mutate(expectation=expectation*100) %>% 
        select(x,y,expectation)
    }) %>% 
  filter(between(expectation,100*0.05,100*0.15)) %>% 
  filter(between(dNDVI, 100*-0.5,100*0.5)) %>% 
  mutate(val = dNDVI-expectation) %>% #pull(val) %>% quantile(., c(0.01,0.99))
  ggplot(data=., aes(x,y,fill=dNDVI-expectation))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL)+
  scale_fill_gradient2(expression(paste(Delta,NDVI,"-",Delta,NDVI[Pred.],"(%)")), 
                       limits=c(-25,25), 
                       high=vec_cols[11], 
                       mid=vec_cols[6],
                       low=vec_cols[1])+
  # scale_fill_viridis_c(expression(paste("Expected",Delta*NDVI("%"))),
  #                      option='D', 
  #                      limits=c(0,20), 
  #                      na.value = 'red',
  #                      # oob=scales::squish
  # )+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=7),
        legend.position = c(1,1), 
        legend.title.align = 1,
        legend.key.width = unit(0.33,'cm'),
        legend.justification = c(1,1)); p_diff_sen

# scico::scico_palette_show()
# vec_cols2 <- scico::scico(n=11,palette = 'batlow',direction = -1)
vec_cols2 <- colorspace::diverge_hcl(11, "Blue-Red",rev = T)
p_violin_sen <- inner_join({lt_ndvi_year_noEpoch %>% as_tibble() %>% 
    filter(between(b0,0.05,1)) %>% 
    # mutate(b2=0) %>% 
    mutate(dNDVI = (37*b1 )/(b0)) %>% 
    select(x,y,dNDVI)},{
      lt_v_sen %>% 
        as_tibble() %>% 
        filter(between(b1,-0.05,0.05)) %>% 
        mutate(dVPD_VPD = b1*37/b0) %>% 
        mutate(expectation = 0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
        select(x,y,expectation)
    }) %>% 
  # filter(between(expectation,0.05,0.15)) %>% 
  # filter(between(dNDVI, -0.4,0.4)) %>% 
  inner_join(., kop, by=c("x","y")) %>% 
  inner_join(., {lt_ppet_sen %>% 
      as_tibble() %>% 
      # filter(between(b1,-0.1,0.1)) %>% 
      mutate(delta_precip = 100*(37*b1)/b0) %>% 
      select(x,y,delta_precip) %>% 
      inner_join(., kop, by=c("x","y")) %>% 
      group_by(zone) %>% 
      summarize(delta_precip=mean(delta_precip,na.rm=TRUE)) %>% 
      ungroup()
  }, by=c("zone")) %>% 
  mutate(zone = recode(zone, 'Desert'='Arid')) %>% 
  mutate(zone = recode(zone, 'Temperate Tas.'='Temp. Tasm.')) %>% 
  mutate(diff = 100*(dNDVI-expectation)) %>% 
  ggplot(data=.,aes(diff,zone,fill=delta_precip))+
  geom_vline(aes(xintercept=0),color='grey50',lty=3)+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), 
              trim=TRUE,color='black')+
  scale_color_viridis_d(option='B',end=0.85)+
  scale_fill_gradient2(expression(paste(Delta*P:PET[12*mo],"(%)")),
                       low = vec_cols2[1], 
                       mid= vec_cols2[6], 
                       high= vec_cols2[11], 
                       limits=c(-40,40)
  )+
  labs(y=NULL, x=expression(paste(Delta,"NDVI",-Delta*NDVI[Pred.]," (%)")))+
  scale_x_continuous(limits=c(-30,40))+
  scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
                                        .Label = c("Equatorial",
                                                   "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
                                                   "Temp. Tasm."), class = c("ordered", "factor"))))+
  theme_linedraw()+
  theme(legend.position = 'top', 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_violin_sen
p_violin_sen


library(patchwork)
ggsave((p_vpd_sen|p_wue_sen|p_diff_sen|p_violin_sen)+
         plot_layout()+
         plot_annotation(tag_levels = 'a',tag_prefix = '(',tag_suffix = ')'), 
       filename = 'figures/Fig4_map_dVpd_dNdviPred_dDifference_violin.png', 
       width = 30, height = 20, units='cm', dpi=350, type='cairo')




co2_start <- tmp_ndvi$co2 %>% min
co2_end <- tmp_ndvi$co2 %>% max
delta_co2 <- co2_end - co2_start

system.time(
  lt_ndvi_co2_p_epoch <-  tmp_ndvi[id%in%vec_ids] %>% #dat_annual %>% #[nobs_annual >= 10] %>% 
    .[,.(beta = list(unname(MASS::rlm(ndvi_hyb~I(co2-co2_start)+frac_p_anom+jitter(epoch), 
                                      data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3],b3=unlist(beta)[4]
    ), by=.(x,y)]
)

p_co2 <- lt_ndvi_co2_p_epoch %>% 
  as_tibble() %>% 
  filter(between(b0,0.05,1)) %>% 
  # mutate(b2=0) %>% 
  mutate(dNDVI = 100*(b1*delta_co2+0*b3)/(b0)) %>% 
  ggplot(data=.,aes(x,y,fill=dNDVI))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL,
       fill=expression(paste(Delta*NDVI[CO[2]]("%"))))+
  # scico::scale_fill_scico(expression(paste(Delta*NDVI[CO[2]]("%"))),
  #                         palette = 'tofino',
  #                         direction = 1,
  #                         limits=c(-5,25), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  scale_fill_continuous_divergingx(
                          palette = 'BrBG',
                          # direction = 1,
                          limits=c(-5,25), #na.value = 'red',
                          oob=scales::squish
  )+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1)); p_co2
lt_ndvi_co2_p_epoch %>% as_tibble() %>% 
  filter(between(b0,0.05,1)) %>% 
  # mutate(b2=0) %>% 
  mutate(dNDVI = (b1*delta_co2+0*b3)/(b0)) %>% 
  pull(dNDVI) %>% quantile(., c(0.05,0.5,0.95))






# hcl_palettes(plot=T)
# RColorBrewer::display.brewer.all()
# vec_cols <- colorspace::hcl_palettes(palette = 'BrBG')
# colorspace::divergingx_palettes(plot=T)
# scico::scico_palette_show()

p_box <- inner_join({lt_ndvi_co2_p_epoch %>% as_tibble() %>% 
    filter(between(b0,0.05,1)) %>% 
    # mutate(b2=0) %>% 
    mutate(dNDVI = (b1*delta_co2+0*b3)/(b0)) %>% 
    select(x,y,dNDVI)},{
      lt_v_sen %>% 
        as_tibble() %>% 
        filter(between(b1,-0.05,0.05)) %>% 
        mutate(dVPD_VPD = b1*37/b0) %>% 
        mutate(expectation10 = 0.1*(dCa_Ca - 0.5*dVPD_VPD), 
               expectation50 = 0.5*(dCa_Ca - 0.5*dVPD_VPD), 
               expectation90 = 0.9*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
        select(x,y,expectation10,expectation50,expectation90)}) %>%
  # filter(between(expectation,0.05,0.15)) %>% 
  # filter(between(dNDVI, -0.4,0.4)) %>% 
  inner_join(., kop, by=c("x","y")) %>% 
  inner_join(., {lt_ppet_sen %>% 
      as_tibble() %>% 
      # filter(between(b1,-0.1,0.1)) %>% 
      mutate(delta_precip = 100*(37*b1)/b0) %>% 
      select(x,y,delta_precip) %>% 
      inner_join(., kop, by=c("x","y")) %>% 
      group_by(zone) %>% 
      summarize(delta_precip=mean(delta_precip,na.rm=TRUE)) %>% 
      ungroup()
  }, by=c("zone")) %>% 
  mutate(zone = recode(zone, 'Desert'='Arid')) %>% 
  mutate(zone = recode(zone, 'Temperate Tas.'='Temp. Tasm.')) %>% 
  mutate(`10%` = 100*(expectation10-dNDVI), 
         `50%` = 100*(expectation50-dNDVI),
         `90%` = 100*(expectation90-dNDVI)) %>%
  select(zone,delta_precip,`10%`,`50%`,`90%`) %>% 
  gather(-zone,-delta_precip,key='key',value='value') %>% 
  mutate(key = factor(key,
                      levels = c('10%','50%','90%'))) %>% 
  ggplot(data=.,aes(value,zone,
                    fill=key))+
  geom_vline(aes(xintercept=0),color='grey50',lty=3)+
  geom_boxplot(#draw_quantiles = c(0.25,0.5,0.75), 
              # trim=TRUE,
    outlier.colour = NA,
              color='black')+
  scico::scale_fill_scico_d(expression(paste('Foliar gain from '~CO[2]~'fert.:')),
                            palette='bamako',direction = -1, begin = 0.2,end=0.95,
  guide=guide_legend(title.position = 'top'))+
  # scale_fill_manual("% WUE benefit allocated to foliar area:",
  #   values=c(vec_cols[1], vec_cols[2], vec_cols[3]))+
  # scale_fill_gradient2(expression(paste(Delta*P:PET[12*mo],"(%)")),
  #                      low = vec_cols2[1], 
  #                      mid= vec_cols2[6], 
  #                      high= vec_cols2[11], 
  #                      limits=c(-40,40)
  # )+
  labs(y=NULL, x=expression(paste(Delta*NDVI[Pred.]-Delta*NDVI[CO[2]]," (%)")))+
  scale_x_continuous(limits=c(-30,30))+
  scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
                                        .Label = c("Equatorial",
                                                   "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
                                                   "Temp. Tasm."), class = c("ordered", "factor"))))+
  theme_linedraw()+
  theme(legend.position = 'top', 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank()); p_box


# ggsave((p_vpd_sen|p_wue_sen|p_co2|p_box)+
#          plot_layout()+
#          plot_annotation(tag_levels = 'a',tag_prefix = '(',tag_suffix = ')'), 
#          filename="figures/Fig6_map_dVpd_dNdviPred_dDifference_boxplot.png",
#          width=35, height=20,units='cm',dpi=350)
library(cowplot)
p_left <- ggdraw(p_vpd_sen)+draw_label(label='(a)', x=0.07,y=0.95,size = 25)
p_mid1 <- ggdraw(p_co2)+draw_label(label='(b)',x=0.05,y=0.95,size=25)
p_mid2 <- ggdraw(p_wue_sen)+draw_label(label='(c)', x=0.05,y=0.95, size=25)
#hjust=0,vjust = 0)
p_right <- ggdraw(p_box)+draw_label(label = '(d)', x=0.05,y=0.95,size=25)

cp_r <- cowplot::plot_grid(p_left,p_mid1,p_mid2,p_right,
                           nrow = 1,
                           rel_widths = c(1,1,1,1))
cp_r
ggsave(plot=cp_r,
       filename = "figures/Fig6_map_dVpd_dNdviPred_dDifference_boxplot.png", 
       width = 35, height=20, units='cm', dpi=350, type='cairo')




## Alternative scatter plot approach to show Obs-Pred differences by 
## climate zone and change in P or P:PET 
# sen_ndvi_season_e1
# inner_join({lt_ndvi_year_noEpoch %>% as_tibble() %>% 
#     filter(between(b0,0.05,1)) %>% 
#     mutate(b2=0) %>% 
#     mutate(dNDVI = (37*b1)) %>% 
#     select(x,y,dNDVI)},
#     {lt_v_sen %>% 
#         as_tibble() %>% 
#         filter(between(b1,-0.1,0.1)) %>% 
#         mutate(dVPD_VPD = b1*37/b0) %>% 
#         mutate(expectation = 0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
#         mutate(expectation_pdiff=expectation*100) %>% 
#         select(x,y,expectation)
#     }) %>% 
#   # filter(between(expectation,100*0.05,100*0.15)) %>% 
#   # filter(between(dNDVI, 100*-0.5,100*0.5)) %>% 
#   mutate(res = dNDVI-expectation) %>% 
#   select(x,y,res,dNDVI) %>% 
#   inner_join(., lt_ppet_sen %>% 
#                as_tibble() %>% 
#                # filter(between(b1,-0.1,0.1)) %>% 
#                mutate(delta_ppet = (37*b1)/b0) %>% 
#                select(x,y,delta_ppet), 
#              by=c('x','y')) %>% 
#   inner_join(., kop, by=c('x','y')) %>% 
#   ggplot(data=.,aes(delta_ppet, res,color=zone))+
#   geom_point()+
#   geom_hline(aes(yintercept=0))+
#   geom_vline(aes(xintercept=0))+
#   scale_color_viridis_d(option='B')+
#   geom_smooth(se=F, method='lm',color='navy',
#               inherit.aes = F,aes(delta_ppet,res))+
#   geom_smooth(se=F, method=MASS::rlm,color='#cf0000',
#               inherit.aes = F,aes(delta_ppet,res))+
#   theme_linedraw()
# 
# 
# inner_join({lt_ndvi_sen_e1 %>% as_tibble() %>% 
#     filter(between(b0,0.05,1)) %>% 
#     mutate(b2=0) %>% 
#     mutate(dNDVI = (18*b1)) %>% 
#     select(x,y,dNDVI)},
#     {lt_v_sen_e1 %>% 
#         as_tibble() %>% 
#         filter(between(b1,-0.1,0.1)) %>% 
#         mutate(dVPD_VPD = b1*18/b0) %>% 
#         mutate(expectation = 0.5*(0.5*dCa_Ca - 0.5*dVPD_VPD)) %>% 
#         mutate(expectation_pdiff=expectation*100) %>% 
#         select(x,y,expectation)
#     }) %>% 
#   # filter(between(expectation,100*0.05,100*0.15)) %>% 
#   # filter(between(dNDVI, 100*-0.5,100*0.5)) %>% 
#   mutate(res = dNDVI-expectation) %>% 
#   select(x,y,res,dNDVI) %>% 
#   inner_join(., lt_ppet_sen_e1 %>% 
#                as_tibble() %>% 
#                filter(between(b0,0.1,3)) %>% 
#                filter(between(b1,-0.1,0.1)) %>%
#                mutate(delta_ppet = (18*b1)/b0) %>% 
#                select(x,y,delta_ppet), 
#              by=c('x','y')) %>% 
#   filter(delta_ppet < 1) %>% 
#   filter(between(res,-0.2,0.2)) %>% 
#   inner_join(., kop, by=c('x','y')) %>% 
#   ggplot(data=.,aes(delta_ppet, res,fill=zone,color=zone))+
#   geom_point(alpha=0.1)+
#   geom_hline(aes(yintercept=0))+
#   geom_vline(aes(xintercept=0))+
#   scale_color_viridis_d(option='B',end=0.9)+
#   scale_fill_viridis_d(option='B',end=0.9)+
#   geom_smooth(se=F, method=MASS::rlm,color='white',
#               lwd=2,
#               inherit.aes = F,aes(delta_ppet,res))+
#   geom_smooth(se=F, method=MASS::rlm,color='#cf0000',
#               inherit.aes = F,aes(delta_ppet,res))+
#   scale_y_continuous(limits=c(-0.15,0.1))+
#   labs(x=expression(paste(Delta~P:PET)), 
#        y=expression(paste(NDVI[pred]-NDVI[obs])))+
#   theme_linedraw()+
#   theme(panel.grid.minor = element_blank())
# 
# # plotting relative percent 
# p_e1 <- inner_join({lt_ndvi_sen_e1 %>% as_tibble() %>% 
#     filter(between(b0,0.05,1)) %>% 
#     mutate(dNDVI = 100*(18*b1/b0)) %>% 
#     select(x,y,dNDVI)},
#     {lt_v_sen_e1 %>% 
#         as_tibble() %>% 
#         filter(between(b1,-0.1,0.1)) %>% 
#         mutate(dVPD_VPD = 100*b1*18/b0) %>% 
#         mutate(expectation = 0.5*(dCa_Ca_e1 - 0.5*dVPD_VPD)) %>% 
#         mutate(expectation_pdiff=expectation) %>% 
#         select(x,y,expectation_pdiff)
#     }) %>% 
#   filter(between(expectation_pdiff,-50,50)) %>%
#   filter(between(dNDVI, -50,50)) %>%
#   mutate(res = dNDVI-expectation_pdiff) %>% 
#   select(x,y,res,dNDVI) %>% 
#   inner_join(., lt_fp_sen_e1 %>% 
#                as_tibble() %>% 
#                filter(between(b0,0.1,3)) %>% 
#                filter(between(b1,-0.025,0.025)) %>%
#                mutate(delta_ppet = 100*(18*b1)/b0) %>% 
#                select(x,y,delta_ppet),
#              by=c('x','y')) %>% 
#   # filter(delta_ppet < 1) %>% 
#   # filter(between(res,-0.2,0.2)) %>% 
#   inner_join(., kop, by=c('x','y')) %>% 
#   ggplot(data=.,aes(delta_ppet, res,fill=zone,color=zone))+
#   geom_point(alpha=0.1)+
#   geom_hline(aes(yintercept=0))+
#   geom_vline(aes(xintercept=0))+
#   scale_color_viridis_d(option='B',end=0.9)+
#   scale_fill_viridis_d(option='B',end=0.9)+
#   geom_smooth(se=F, method=MASS::rlm,color='white',
#               lwd=2,
#               inherit.aes = F,aes(delta_ppet,res))+
#   geom_smooth(se=F, method=MASS::rlm,color='#cf0000',
#               inherit.aes = F,aes(delta_ppet,res))+
#   # scale_y_continuous(limits=c(-0.15,0.1))+
#   labs(x=expression(paste(Delta~P:PET)), 
#        y=expression(paste(NDVI[pred]-NDVI[obs])))+
#   theme_linedraw()+
#   theme(panel.grid.minor = element_blank()); p_e1
# 
# p_e2 <- inner_join({lt_ndvi_sen_e2 %>% as_tibble() %>% 
#     filter(between(b0,0.05,1)) %>% 
#     mutate(dNDVI = 100*(18*b1/b0)) %>% 
#     select(x,y,dNDVI)},
#     {lt_v_sen_e2 %>% 
#         as_tibble() %>% 
#         filter(between(b1,-0.1,0.1)) %>% 
#         mutate(dVPD_VPD = 100*b1*18/b0) %>% 
#         mutate(expectation = 0.5*(dCa_Ca_e2 - 0.5*dVPD_VPD)) %>% 
#         mutate(expectation_pdiff=expectation) %>% 
#         select(x,y,expectation_pdiff)
#     }) %>% 
#   filter(between(expectation_pdiff,-50,50)) %>%
#   filter(between(dNDVI, -50,50)) %>%
#   mutate(res = dNDVI-expectation_pdiff) %>% 
#   select(x,y,res,dNDVI) %>% 
#   inner_join(., lt_fp_sen_e2 %>% 
#                as_tibble() %>% 
#                filter(between(b0,0.1,3)) %>% 
#                filter(between(b1,-0.025,0.025)) %>%
#                mutate(delta_ppet = 100*(18*b1)/b0) %>% 
#                select(x,y,delta_ppet),
#              by=c('x','y')) %>% 
#   # filter(delta_ppet < 1) %>% 
#   # filter(between(res,-0.2,0.2)) %>% 
#   inner_join(., kop, by=c('x','y')) %>% 
#   ggplot(data=.,aes(delta_ppet, res,fill=zone,color=zone))+
#   geom_point(alpha=0.1)+
#   geom_hline(aes(yintercept=0))+
#   geom_vline(aes(xintercept=0))+
#   scale_color_viridis_d(option='B',end=0.9)+
#   scale_fill_viridis_d(option='B',end=0.9)+
#   # geom_smooth(se=F, method=MASS::rlm,color='white',
#   #             lwd=2,
#   #             inherit.aes = F,aes(delta_ppet,res))+
#   geom_smooth(method='lm',se=F)+
#   # geom_smooth(se=F, method=MASS::rlm, #color='#cf0000',
#   #             inherit.aes = T,aes(delta_ppet,res))+
#   # scale_y_continuous(limits=c(-0.15,0.1))+
#   labs(x=expression(paste(Delta~P:PET)), 
#        y=expression(paste(NDVI[pred]-NDVI[obs])))+
#   theme_linedraw()+
#   theme(panel.grid.minor = element_blank()); p_e2
# 
# p_e1+p_e2+theme(legend.position = 'none')
# 
# # lt_ppet_sen_e1 %>% 
# #   as_tibble() %>% 
# #   filter(between(b0,0.1,3)) %>% 
# #   filter(between(b1,-0.05,0.05)) %>%
# #   mutate(delta_ppet = 100*(18*b1)/b0) %>% 
# #   select(x,y,delta_ppet) %>%   pull(delta_ppet) %>% hist
# # 
# # lt_ppet_sen_e1 %>% 
# #   as_tibble() %>% 
# #   filter(b0 > 0) %>% 
# #   filter(b1 == max(b1,na.rm=T)) %>% 
# #   pull(x)
# # dat[near(x,145.9538,tol=0.01) & near(y,-17.18028,tol=0.01)] %>% 
# #   ggplot(data=.,aes(hydro_year-1982, pe_12mo))+
# #   geom_point()+
# #   # geom_smooth(method='lm')+
# #   # geom_smooth(method=MASS::rlm,col='red')+
# #   geom_abline(aes(intercept= 0.105, slope = 0.144), color='black')+
# #   scale_y_continuous(limits=c(0,5))
# # 
# # 
# 
# dat[epoch=='avhrr'] %>% 
#   .[sample(.N,1e6)] %>% 
#   lm(ndvi_hyb~scale(co2_trend)+scale(pe_anom_12mo)+scale(precip_anom_12mo), data=.) %>% 
#   summary
# dat[epoch=='modis'] %>% 
#   .[sample(.N,1e6)] %>% 
#   lm(ndvi_hyb~scale(co2_trend)+scale(pe_anom_12mo)+scale(precip_anom_12mo), data=.) %>% 
#   summary
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# p_violin_sen <- inner_join({lt_ndvi_year_noEpoch %>% as_tibble() %>% 
#     filter(between(b0,0.05,1)) %>% 
#     mutate(b2=0) %>% 
#     mutate(dNDVI = (37*b1)/(b0)) %>% 
#     select(x,y,dNDVI)},{
#       lt_v_sen %>% 
#         as_tibble() %>% 
#         filter(between(b1,-0.05,0.05)) %>% 
#         mutate(dVPD_VPD = b1*37/b0) %>% 
#         mutate(expectation = 0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
#         select(x,y,expectation)
#     }) %>% 
#   # filter(between(expectation,0.05,0.15)) %>% 
#   # filter(between(dNDVI, -0.4,0.4)) %>% 
#   inner_join(., kop, by=c("x","y")) %>% 
#   inner_join(., {lt_ppet_sen %>% 
#       as_tibble() %>% 
#       # filter(between(b1,-0.1,0.1)) %>% 
#       mutate(delta_precip = 100*(37*b1)/b0) %>% 
#       select(x,y,delta_precip) %>% 
#       inner_join(., kop, by=c("x","y")) %>% 
#       group_by(zone) %>% 
#       summarize(delta_precip=mean(delta_precip,na.rm=TRUE)) %>% 
#       ungroup()
#   }, by=c("zone")) %>% 
#   mutate(zone = recode(zone, 'Desert'='Arid')) %>% 
#   mutate(zone = recode(zone, 'Temperate Tas.'='Temp. Tasm.')) %>% 
#   mutate(diff = 100*(dNDVI-expectation)) %>% 
#   ggplot(data=.,aes(diff,zone,fill=delta_precip))+
#   geom_vline(aes(xintercept=0),color='grey50',lty=3)+
#   geom_violin(draw_quantiles = c(0.25,0.5,0.75), 
#               trim=TRUE,color='black')+
#   scale_color_viridis_d(option='B',end=0.85)+
#   scale_fill_gradient2(expression(paste(Delta*P:PET[12*mo],"(%)")),
#                        low = vec_cols2[1], 
#                        mid= vec_cols2[6], 
#                        high= vec_cols2[11], 
#                        limits=c(-40,40)
#   )+
#   labs(y=NULL, x=expression(paste(Delta,"NDVI",-Delta*NDVI[Pred.]," (%)")))+
#   scale_x_continuous(limits=c(-30,40))+
#   scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
#                                         .Label = c("Equatorial",
#                                                    "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
#                                                    "Temp. Tasm."), class = c("ordered", "factor"))))+
#   theme_linedraw()+
#   theme(legend.position = 'top', 
#         legend.key.height = unit(0.2,'cm'),
#         panel.grid = element_blank()); p_violin_sen


tmp_ndvi %>% select(co2,frac_p_anom,frac_vpd_anom,frac_pet_anom,frac_ppet_anom) %>% drop_na() %>% cor
system.time(
  lt_ndvi_co2_p_pet_ppet_epoch <-  tmp_ndvi[id%in%vec_ids] %>% #dat_annual %>% #[nobs_annual >= 10] %>%
    .[,.(beta = list(unname(lm(ndvi_hyb~I(co2-co2_start)+
                                        frac_p_anom+
                                        frac_vpd_anom+
                                        I(frac_vpd_anom**2)+
                                        jitter(epoch),
                                      data=.SD)$coefficients))),
      by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1],
            b_co2=unlist(beta)[2]
            # b_vpd=unlist(beta)[3],
            # b_pet_anom=unlist(beta)[4],
            # b_epoch=unlist(beta)[5]
    ), by=.(x,y)]
)
lt_ndvi_co2_p_pet_ppet_epoch %>% as_tibble() %>%
  filter(between(b0,0.05,1)) %>%
  # mutate(b2=0) %>%
  mutate(dNDVI = (b_co2*delta_co2)/(b0)) %>%
  pull(dNDVI) %>% quantile(., c(0.05,0.5,0.95))

lt_ndvi_co2_p_pet_ppet_epoch %>% as_tibble() %>%
  filter(between(b0,0.05,1)) %>%
  # mutate(b2=0) %>%
  mutate(dNDVI = (b_co2*delta_co2)/(b0)) %>%
  pull(dNDVI) %>% hist(100)


system.time(
  lt_test <-  tmp_ndvi[id%in%vec_ids] %>% #dat_annual %>% #[nobs_annual >= 10] %>%
    .[,.(beta = list(unname(lm(ndvi_hyb~I(co2-co2_start)+
                                 frac_p_anom+
                                 I(frac_p_anom**2)+
                                 frac_vpd_anom+
                                 I(frac_vpd_anom**2)+
                                 frac_ppet_anom+
                                 I(frac_ppet_anom**2)+
                                 jitter(epoch),
                               data=.SD)$coefficients))),
      by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1],
            b_co2=unlist(beta)[2],
            # b_vpd=unlist(beta)[3],
            # b_pet_anom=unlist(beta)[4],
            b_epoch=unlist(beta)[9]
    ), by=.(x,y)]
)
lt_test %>% as_tibble() %>%
  # filter(between(b0,0.05,1)) %>%
  # mutate(b2=0) %>%
  mutate(dNDVI = exp(b0+b_co2*delta_co2+b_epoch)/(exp(b0)) - 1) %>%
  pull(dNDVI) %>% quantile(., c(0.05,0.5,0.95),na.rm=T)




lt_ndvi_co2_p_pet_ppet_epoch %>% 
  as_tibble() %>% 
  filter(between(b0,0.05,1)) %>% 
  # mutate(b2=0) %>% 
  mutate(dNDVI = 100*(b_co2*delta_co2)/(b0)) %>% 
  ggplot(data=.,aes(x,y,fill=dNDVI))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL,
       fill=expression(paste(Delta*NDVI[CO[2]]("%"))))+
  scico::scale_fill_scico(expression(paste(Delta*NDVI[CO[2]]("%"))),
                          palette = 'bamako',
                          direction = -1,
                          limits=c(6,12), #na.value = 'red',
                          oob=scales::squish
  )+
  # scale_fill_continuous_divergingx(
  #   palette = 'BrBG',
  #   # direction = 1,
  #   limits=c(-5,25), #na.value = 'red',
  #   oob=scales::squish
  # )+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1));




tmp_mappet <- unique(dat[,.(x,y,mape)])
ggplot(data=tmp_mappet,aes(x,y,fill=mape))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,7))

n4 <- merge(tmp_ndvi,tmp_mappet,by=c("x","y")) %>%
  .[id%in%vec_ids] %>%
  .[,`:=`(cco2 = co2-co2_start)] %>% 
  nls.multstart::nls_multstart(ndvi_hyb ~ 
                                 Asym-Drop*exp(-exp(lrc)*mape^pwr) + 
                                 B1*(frac_ppet_anom) + 
                                 B2*(cco2)+
                                 B3*(cco2*mape)+
                                 B4*(frac_vpd_anom)+
                                 B5*(frac_p_anom),
                               data = .,
                               iter = 1,
                               start_lower = c(Asym=0.0, Drop=0.6,lrc=0,pwr=0,B1=-0.5,B2=0,B3=0,B4=-0.1,B5=-0.1),
                               start_upper = c(Asym=1, Drop=1,lrc=1,pwr=2,B1=0.5,B2=0.001,B3=0.001,B4=0.1,B5=0.1),
                               # supp_errors = 'Y',
                               na.action = na.omit)
summary(n4)

yardstick::rsq_vec(truth=merge(tmp_ndvi,tmp_mappet,by=c("x","y")) %>%
                     .[id%in%vec_ids] %>% pull(ndvi_hyb), 
                   estimate=predict(n4,newdata=merge(tmp_ndvi,tmp_mappet,by=c("x","y")) %>%
                                      .[id%in%vec_ids]))

junk1 <- tmp_mappet %>% 
  as_tibble() %>%
  mutate(frac_ppet_anom=0,frac_vpd_anom=0,frac_p_anom=0,epoch=1,cco2=co2_start) %>% 
  mutate(pred1 = predict(n4,newdata=.))
junk2 <- tmp_mappet %>% 
  as_tibble() %>%
  mutate(frac_ppet_anom=0,frac_vpd_anom=0,frac_p_anom=0,epoch=1,cco2=co2_end) %>% 
  mutate(pred2 = predict(n4,newdata=.))

inner_join(junk1,junk2,by=c('x','y')) %>% 
  mutate(dNDVI = 100*(pred2-pred1)/pred1) %>% 
  # mutate(dNDVI = exp(b0+b_co2*delta_co2)/(exp(b0)) - 1) %>%
  pull(dNDVI) %>% quantile(., c(0.05,0.5,0.95),na.rm=T)

inner_join(junk1,junk2,by=c('x','y')) %>% 
  mutate(dNDVI = (pred2-pred1)/pred1) %>% pull(dNDVI) %>% summary
