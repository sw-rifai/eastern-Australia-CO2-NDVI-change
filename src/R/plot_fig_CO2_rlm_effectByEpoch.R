# !diagnostics off
#*******************************************************************************
#* Description:
#* Plot Climate stuff
#*
#*
#*
library(testthat)
library(tidyverse); library(sf)
library(data.table); setDTthreads(threads = 0)
library(lubridate); 
library(mgcv); #library(mgcViz); 
library(dtplyr); 
library(RcppArmadillo); library(patchwork)
library(stars)
library(foreach); library(doParallel)
library(zyp); library(nls.multstart)
options(pillar.sigfig = 5)
# IMPORT DATA ###################################################################

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
                             "vpd15_12mo","vpd15_anom_12mo",
                             "vpd15_u",
                             "pet","mapet","pet_anom","pet_anom_3mo","pet_u","pet_sd",
                             "pet_anom_sd", "pet_12mo","pet_36mo",
                             "pe","mape",
                             'vc','veg_class',
                             'month',
                             "x", "y", "year")) %>% 
  as.data.table() %>% 
  .[is.infinite(mape)==F]
norms_mape <- dat %>% lazy_dt() %>% 
  filter(date>=ymd("1982-01-01") & date<=ymd("2011-12-31")) %>% 
  group_by(x,y,hydro_year) %>% 
  summarize(ppet_12mo = mean(precip_12mo/pet_12mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(x,y) %>% 
  summarize(mappet = mean(ppet_12mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
dat <- merge(dat,norms_mape,by=c("x","y"))

dat <- dat[,`:=`(pe = precip/pet, 
                 pe_12mo = precip_12mo/pet_12mo)]
dat <- merge(dat, 
             vi,
             by=c("x","y","date"), 
             all=TRUE,allow.cartesian=TRUE)
dat <- dat[order(x,y,date)][, ndvi_3mo := frollmean(ndvi_hyb,n = 3,fill = NA,align='center',na.rm=TRUE), by=.(x,y)]
rm(vi); gc(full=TRUE)

# Attach season
dat[,`:=`(year=year(date),month=month(date))] %>%
  .[,`:=`(season = case_when(month%in%c(3:5)~'MAM',
                             month%in%c(6:8)~'JJA',
                             month%in%c(9:11)~'SON',
                             month%in%c(12,1,2)~'DJF'))]
dat[,`:=`(season = factor(season, levels = c('SON','DJF','MAM','JJA'),ordered = TRUE))]
dat[,`:=`(hydro_year=year(date+months(1)))]

# FILTER TO LON >= 140 !!! **********
dat <- dat[x>=140]

coords_keep <- dat %>% lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(nobs_total = sum(is.na(ndvi_hyb)==F)) %>% 
  ungroup() %>% 
  as.data.table()
dat <- merge(dat, coords_keep, by=c("x","y"))

# Load Mauna Loa CO2 record
mlo <- readr::read_table("../data_general/CO2_growth_rate/co2_mm_mlo_20200405.txt", 
                         skip = 72, col_names = F) %>% 
  set_names(
    c("year","month","ddate","co2_avg","co2_int","co2_trend","ndays")
  ) %>% 
  mutate(date = ymd(paste(year,month,1))) %>% 
  select(date,co2_int,co2_trend) %>% 
  as.data.table()
dat <- merge(mlo,dat,by="date",all = TRUE)

ldat <- dat %>% lazy_dt()
# END Load awap clim dat *****************************************************************



# Load simplified BOM Koppen climate zones --------------------------------
ref_grid <- stars::read_stars('../data_general/AVHRR_CDRv5_VI/AVHRR_SR_median_EastOz_1982_2019.tif', 
                              RasterIO = list(bands=1))
bom <- stars::read_stars("../data_general/Koppen_climate/BOM/kpngrp.txt")
bom <- st_warp(src=bom, dest=ref_grid[,,], use_gdal = F)
bom <- set_names(bom, 'koppen') %>% as_tibble()
bom <- left_join(ref_grid %>% as_tibble() %>% select(x,y), 
                 bom)

coords <- dat %>% select(x,y) %>% distinct() %>% filter(x>=140) 

g_map <- ldat %>% 
  filter(date>=ymd("1982-01-01")&date<=ymd("2010-12-31")) %>% 
  group_by(x,y) %>% 
  summarize(map = mean(precip,na.rm=TRUE)*12) %>% 
  ungroup() %>% 
  as_tibble()

kop <- left_join(coords, bom, by=c("x","y")) %>% 
  inner_join(., g_map, by=c("x","y")) %>% 
  as_tibble() %>% 
  mutate(zone = case_when(between(koppen,0,11) ~ 'Temperate', 
                          (y <= -40) ~ 'Tasmania',
                          between(koppen, 12,21)~'GD_temp', # Grassland
                          between(koppen, 22,30)~'GD_temp', # Desert
                          between(koppen, 31,34)~'Subtropical',
                          between(koppen, 35,40)~'Tropical', 
                          koppen >= 41 ~ 'Equatorial')) %>% 
  mutate(zone = ifelse(y < -40, 'Temperate Tas.', zone)) %>% #pull(zone) %>% table
  mutate(zone = ifelse(zone == "GD_temp" & map < 500, 'Arid',zone)) %>%   
  mutate(zone = ifelse(zone == "GD_temp" & map >= 500, 'Grassland',zone)) %>%   
  mutate(zone = factor(zone, levels = c("Equatorial","Tropical",
                                        "Subtropical","Grassland","Arid",
                                        "Temperate","Temperate Tas."), ordered = T))
kop <- kop %>% mutate(cz=zone)
kop <- as.data.table(kop)
#*** End Kop zone load ********************************************************




# SM Fig 14: NDVI Beta CO2 comparison by epoch ----------------------------
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
          # epoch_bin = ifelse(date<=ymd("2000-12-31"),0,1),
          frac_p = precip_12mo/map,
          frac_p_anom = precip_anom_12mo/map,
          frac_pet_anom = (pet_12mo-mapet)/mapet,
          frac_vpd_anom = vpd15_anom/mavpd15,
          frac_ppet_anom = (pe_12mo-mape)/mape)] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE), 
       ndvi_mcd = mean(ndvi_mcd,na.rm=TRUE),
       ndvi_cdr = mean(ndvi_cdr,na.rm=TRUE),
       co2 = mean(co2_trend,na.rm=TRUE),
       # epoch=mean(epoch_bin,na.rm=TRUE), 
       nobs = sum(is.na(ndvi_hyb)==F), 
       p_anom = mean(precip_anom_12mo,na.rm=TRUE),
       pet_anom = mean(pet_12mo-mapet,na.rm=TRUE),
       ppet_anom = mean(pe_12mo,na.rm=TRUE),
       vpd_anom = mean(vpd15_12mo,na.rm=TRUE),
       # epoch = mean(epoch, na.rm=TRUE), 
       map = mean(map,na.rm=TRUE), 
       mapet = mean(mapet, na.rm=TRUE), 
       mappet = mean(mappet, na.rm=TRUE), 
       mavpd15 = mean(mavpd15,na.rm=TRUE),
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
          # epoch = ifelse(date<=ymd("2000-12-31"),0,1),
          frac_p = precip_12mo/map,
          frac_p_anom = precip_anom_12mo/map,
          frac_pet_anom = (pet_12mo-mapet)/mapet,
          frac_vpd_anom = vpd15_anom/mavpd15,
          frac_ppet_anom = (pe_12mo-mape)/mape)] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE), 
       ndvi_mcd = mean(ndvi_mcd,na.rm=TRUE),
       ndvi_cdr = mean(ndvi_cdr,na.rm=TRUE),
       co2 = mean(co2_trend,na.rm=TRUE),
       nobs = sum(is.na(ndvi_hyb)==F), 
       p_anom = mean(precip_anom_12mo,na.rm=TRUE),
       pet_anom = mean(pet_12mo-mapet,na.rm=TRUE),
       ppet_anom = mean(pe_12mo,na.rm=TRUE),
       vpd_anom = mean(vpd15_12mo,na.rm=TRUE),
       # epoch = mean(epoch, na.rm=TRUE), 
       map = mean(map,na.rm=TRUE), 
       mapet = mean(mapet, na.rm=TRUE), 
       mappet = mean(mappet, na.rm=TRUE), 
       mavpd15 = mean(mavpd15,na.rm=TRUE),
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
       ndvi_mcd = mean(ndvi_mcd,na.rm=TRUE),
       ndvi_cdr = mean(ndvi_cdr,na.rm=TRUE),
       co2 = mean(co2_trend,na.rm=TRUE),
       epoch=mean(epoch,na.rm=TRUE), 
       nobs = sum(is.na(ndvi_hyb)==F), 
       p_anom = mean(precip_anom_12mo,na.rm=TRUE),
       pet_anom = mean(pet_12mo-mapet,na.rm=TRUE),
       ppet_anom = mean(pe_12mo,na.rm=TRUE),
       vpd_anom = mean(vpd15_12mo,na.rm=TRUE),
       map = mean(map,na.rm=TRUE), 
       mapet = mean(mapet, na.rm=TRUE), 
       mappet = mean(mappet, na.rm=TRUE), 
       mavpd15 = mean(mavpd15,na.rm=TRUE),
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
# tmp_epoch_nobs %>% 
#   inner_join(., tmp_ndvi %>% select(x,y,id) %>% distinct(), by='id') %>% 
#   ggplot(data=.,aes(x,y,fill=obs_epoch))+geom_tile()+scale_fill_viridis_c()
vec_ids <- tmp_epoch_nobs %>% 
  filter(obs_epoch >= 1.9) %>% 
  pull(id)

# Percent increase in NDVI after effects of ppet were removed
min_co2 <- min(tmp_ndvi[hydro_year==1982]$co2)
max_co2 <- max(tmp_ndvi[hydro_year==2019]$co2)
mid_co2 <- max(tmp_ndvi[hydro_year==2000]$co2)

rlm_ndvi_annual_co2_ppet_epoch <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_start = co2 - min_co2)] %>% 
  .[,.(beta = list(coef(MASS::rlm(
    ndvi_hyb~
      co2_start+
      frac_vpd_anom+
      frac_p_anom+
      frac_pet_anom+
      jitter(epoch)))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b3=unlist(beta)[4],
          b4=unlist(beta)[5], 
          b5=unlist(beta)[6]
  ), by=.(x,y)]


vec_id_e1 <- tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 1982:2000] %>% 
  .[is.na(ndvi_cdr)==F & is.na(frac_p_anom)==F] %>% 
  .[,.(nobs=.N),keyby=.(id)] %>% 
  .[nobs >= 10] %>% 
  pull(id)

rlm_ndvi_annual_co2_ppet_e1 <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 1982:2000] %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_start = co2 - min_co2)] %>% 
  .[id%in%vec_id_e1] %>% 
  .[,.(beta = list(coef(MASS::rlm(
    ndvi_cdr~
      co2_start+
      frac_vpd_anom+
      frac_p_anom+
      frac_pet_anom))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b3=unlist(beta)[4],
          b4=unlist(beta)[5]
  ), by=.(x,y)]


vec_id_e2 <- tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 2001:2019] %>% 
  .[is.na(ndvi_mcd)==FALSE & is.na(frac_vpd_anom)==FALSE] %>% 
  .[,.(nobs=.N),keyby=.(id)] %>% 
  .[nobs >= 10] %>% 
  pull(id)

rlm_ndvi_annual_co2_ppet_e2 <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 2001:2019] %>% 
  .[is.na(ndvi_mcd)==FALSE & is.na(frac_vpd_anom)==FALSE] %>% 
  .[,`:=`(co2_start = co2 - mid_co2)] %>% 
  .[id%in%vec_id_e2] %>% 
  .[,.(beta = list(coef(MASS::rlm(
    ndvi_hyb~
      co2_start+
      frac_vpd_anom+
      frac_p_anom+
      frac_pet_anom))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b3=unlist(beta)[4],
          b4=unlist(beta)[5]
  ), by=.(x,y)]

bind_rows(as_tibble(rlm_ndvi_annual_co2_ppet_epoch) %>% mutate(epoch='Merged 1982-2019'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e1) %>% mutate(epoch='AVHRR 1982-2000'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e2) %>% mutate(epoch='MODIS 2001-2019')) %>% 
  inner_join(., as_tibble(kop %>% select(x,y,zone)), by=c("x","y")) %>% 
  mutate(epoch=factor(epoch, ordered = T, levels=rev(c("Merged 1982-2019",
                                                       "AVHRR 1982-2000",
                                                       "MODIS 2001-2019")))) %>% 
  ggplot(data=.,aes(b1, zone, color=zone, fill=epoch))+
  geom_vline(aes(xintercept=0),color='grey')+
  geom_boxplot(outlier.colour = NA)+
  scale_color_viridis_d('Climate', option='B', end=0.85)+
  scale_fill_manual('Epoch', values=c(
    "Merged 1982-2019"='white',
    "AVHRR 1982-2000"='grey80', 
    "MODIS 2001-2019"='grey30'))+
  scale_x_continuous(limits=c(-0.0025,0.005))+
  scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
                                        .Label = c("Equatorial",
                                                   "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
                                                   "Temperate Tas."), class = c("ordered", "factor"))))+
  labs(y=NULL,
       x=expression(paste(Delta,"NDVI"~CO[2]~'ppm'**-1)))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
# ggsave(filename = 'figures/fig_5_rlm_CO2_effect_by_epoch.png',
#        width=16, height = 12, units='cm',type='cairo')
# End section ******************************************************************


# Calculate pixel level goodness of fit ----------------------------------------
fn_gof_rlm <- function(din){
  fit <- MASS::rlm(ndvi_hyb~co2_start+frac_vpd_anom+p_anom+pet_anom+jitter(epoch),data=din)
  out_r2 <- yardstick::rsq_vec(truth=din$ndvi_hyb, estimate=predict(fit,type='response'))
  out_rmse <- yardstick::rmse_vec(truth=din$ndvi_hyb, estimate=predict(fit,type='response'))
  out_df <- data.frame(r2=out_r2, rmse=out_rmse)
  return(out_df)
}
gof_rlm_ndvi_annual_co2_ppet_epoch <- tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_start = co2 - min_co2)] %>% 
  .[,fn_gof_rlm(.SD), by=.(x,y)]

gof_rlm_ndvi_annual_co2_ppet_epoch$r2 %>% summary

gof_rlm_ndvi_annual_co2_ppet_epoch$rmse %>% summary





g3_eval <- bam(ndvi_hyb~
                 te(mappet,co2,k=5)+
                 te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
                 te(map,frac_p_anom,k=5,bs='cs')+
                 te(mapet,frac_pet_anom,k=5,bs='cs')+
                 epoch,
               data=merge(tmp_ndvi[id%in%vec_ids],kop[,.(x,y,zone)],by=c('x','y')) %>% 
                 lazy_dt() %>% 
                 mutate(zone=factor(zone)) %>% 
                 as.data.table(), 
               discrete = T,select=TRUE, method='fREML')
summary(g3_eval)

yardstick::rsq_trad_vec(truth=merge(tmp_ndvi[id%in%vec_ids],kop[,.(x,y,zone)],by=c('x','y')) %>% 
                          lazy_dt() %>% 
                          mutate(zone=factor(zone)) %>% 
                          mutate(pred=predict(g3_eval,newdata=.)) %>% 
                          as.data.table() %>%
                          filter(is.na(pred)==F) %>% 
                          pull(ndvi_hyb), 
                        estimate=merge(tmp_ndvi[id%in%vec_ids],kop[,.(x,y,zone)],by=c('x','y')) %>% 
                          lazy_dt() %>% 
                          mutate(zone=factor(zone)) %>% 
                          mutate(pred=predict(g3_eval,newdata=.)) %>% 
                          as.data.table() %>%
                          filter(is.na(pred)==F) %>% 
                          pull(pred))
yardstick::rmse_vec(truth=merge(tmp_ndvi[id%in%vec_ids],kop[,.(x,y,zone)],by=c('x','y')) %>% 
                          lazy_dt() %>% 
                          mutate(zone=factor(zone)) %>% 
                          mutate(pred=predict(g3_eval,newdata=.)) %>% 
                          as.data.table() %>%
                          filter(is.na(pred)==F) %>% 
                          pull(ndvi_hyb), 
                        estimate=merge(tmp_ndvi[id%in%vec_ids],kop[,.(x,y,zone)],by=c('x','y')) %>% 
                          lazy_dt() %>% 
                          mutate(zone=factor(zone)) %>% 
                          mutate(pred=predict(g3_eval,newdata=.)) %>% 
                          as.data.table() %>%
                          filter(is.na(pred)==F) %>% 
                          pull(pred))

base_pred <- unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
  .[,`:=`(frac_vpd_anom=0, frac_p_anom=0, frac_pet_anom=0, epoch=0)]
preds <- expand_grid(base_pred, tibble(co2=c(340,370,400), tp=c('start','mid','end')))
preds <- preds %>% 
  # merge(., kop[,.(x,y,zone)] %>% mutate(zone=factor(zone)),by=c('x','y')) %>% 
  mutate(pred = predict(g3_eval,newdata=.))
preds %>% pivot_wider(names_from=c('tp','co2'),values_from='pred') %>% 
  select(x,y,start_340,mid_370,end_400) %>% 
  mutate(delta1 = (mid_370-start_340)/30, 
         delta2 = (end_400-mid_370)/30, 
         delta3 = (end_400-start_340)/60 ) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  select(delta1,delta2,delta3,zone) %>% 
  gather(-zone,key='epoch',value='delta_ndvi') %>% 
  mutate(epoch=   case_when(epoch=='delta1'~"AVHRR 1982-2000",
                            epoch=='delta2'~"MODIS 2001-2019", 
                            epoch=='delta3'~"Merged 1982-2019")) %>% 
  filter(is.na(zone)==FALSE) %>% 
  ggplot(data=.,aes(delta_ndvi, zone, color=zone, fill=epoch))+
  geom_vline(aes(xintercept=0),color='grey')+
  geom_boxplot(outlier.colour = NA)+
  scale_color_viridis_d('Climate', option='B', end=0.85)+
  scale_fill_manual('Epoch', values=c(
    "Merged 1982-2019"='white',
    "AVHRR 1982-2000"='grey80', 
    "MODIS 2001-2019"='grey30'))+
  # scale_x_continuous(limits=c(-0.0025,0.005))+
  scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
                                        .Label = c("Equatorial",
                                                   "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
                                                   "Temperate Tas."), class = c("ordered", "factor"))))+
  labs(y=NULL,
       x=expression(paste(Delta,"NDVI"~CO[2]~'ppm'**-1)))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())



deltas_gam <- preds %>% pivot_wider(names_from=c('tp','co2'),values_from='pred') %>% 
  select(x,y,start_340,mid_370,end_400) %>% 
  mutate(delta1 = (mid_370-start_340)/30, 
         delta2 = (end_400-mid_370)/30, 
         delta3 = (end_400-start_340)/60 ) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  select(delta1,delta2,delta3,zone) %>% 
  gather(-zone,key='epoch',value='delta_ndvi') %>% 
  mutate(epoch=   case_when(epoch=='delta1'~"AVHRR 1982-2000",
                            epoch=='delta2'~"MODIS 2001-2019", 
                            epoch=='delta3'~"Merged 1982-2019")) %>% 
  filter(is.na(zone)==FALSE) %>% 
  mutate(model = 'GAM') %>% 
  select(zone,epoch,model,delta_ndvi)



deltas_rlm <- bind_rows(as_tibble(rlm_ndvi_annual_co2_ppet_epoch) %>% mutate(epoch='Merged 1982-2019'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e1) %>% mutate(epoch='AVHRR 1982-2000'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e2) %>% mutate(epoch='MODIS 2001-2019')) %>% 
  inner_join(., as_tibble(kop %>% select(x,y,zone)), by=c("x","y")) %>% 
  # mutate(epoch=factor(epoch, ordered = T, levels=rev(c("Merged 1982-2019",
  #                                                      "AVHRR 1982-2000",
  #                                                      "MODIS 2001-2019")))) %>% 
  mutate(model= 'RLM') %>% 
  mutate(delta_ndvi = b1) %>% 
  select(zone,epoch,model,delta_ndvi)

bind_rows(deltas_gam,deltas_rlm) %>% 
  mutate(epoch=factor(epoch,ordered=T,level=c("Merged 1982-2019",
                                              "AVHRR 1982-2000",
                                              "MODIS 2001-2019"))) %>%
  ggplot(data=.,aes(delta_ndvi, zone, color=model, fill=epoch))+
  geom_vline(aes(xintercept=0),color='grey30',lwd=1)+
  geom_boxplot(outlier.colour = NA, varwidth = F,lwd=0.75, coef=0)+
  # stat_boxplot()+
  # scale_color_viridis_d('Climate', option='B', end=0.85)+
  scale_fill_manual('Epoch', values=c(
    "Merged 1982-2019"='white',
    "AVHRR 1982-2000"='grey80', 
    "MODIS 2001-2019"='grey30'))+
  scale_color_manual('Model', values=c("GAM"='#4287f5',"RLM"="#c74b2c"))+
  scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L),
                  .Label = c("Equatorial",
                             "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
                             "Temperate Tas."), class = c("ordered", "factor"))))+
  labs(y=NULL,
       x=expression(paste(Delta,"NDVI"~CO[2]~'ppm'**-1)))+
  coord_cartesian(xlim=c(-0.0005,0.002),clip='off')+
  theme_linedraw()+
  guides(fill=guide_legend(title.position = 'top',), 
         color=guide_legend(title.position = 'top'))+
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom', 
        legend.text = element_text(size=8))
ggsave(filename = 'figures/fig_5_rlm_CO2_effect_by_epoch.png',
       width=17, height = 16, units='cm',type='cairo')











bind_rows(as_tibble(rlm_ndvi_annual_co2_ppet_epoch) %>% mutate(epoch='Merged 1982-2019'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e1) %>% mutate(epoch='AVHRR 1982-2000'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e2) %>% mutate(epoch='MODIS 2001-2019')) %>% 
  inner_join(., as_tibble(kop %>% select(x,y,zone)), by=c("x","y")) %>% 
  mutate(epoch=factor(epoch, ordered = T, levels=rev(c("Merged 1982-2019",
                                                       "AVHRR 1982-2000",
                                                       "MODIS 2001-2019")))) %>% 
  filter(epoch!='Merged 1982-2019') %>% 
  t.test(. %>% filter()) %>% 
  lm(b1~epoch, data=.) %>% 
  plot()

  
bind_rows(as_tibble(rlm_ndvi_annual_co2_ppet_epoch) %>% mutate(epoch='Merged 1982-2019'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e1) %>% mutate(epoch='AVHRR 1982-2000'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e2) %>% mutate(epoch='MODIS 2001-2019')) %>% 
  inner_join(., as_tibble(kop %>% select(x,y,zone)), by=c("x","y")) %>% 
  mutate(epoch=factor(epoch, ordered = T, levels=rev(c("Merged 1982-2019",
                                                       "AVHRR 1982-2000",
                                                       "MODIS 2001-2019")))) %>% 
  filter(epoch!='Merged 1982-2019') %>% 
  lm(b1~epoch,data=.) %>% 
  summary()


junk <- bind_rows(as_tibble(rlm_ndvi_annual_co2_ppet_epoch) %>% mutate(epoch='Merged 1982-2019'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e1) %>% mutate(epoch='AVHRR 1982-2000'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e2) %>% mutate(epoch='MODIS 2001-2019')) %>% 
  inner_join(., as_tibble(kop %>% select(x,y,zone)), by=c("x","y")) %>% 
  mutate(epoch=factor(epoch, ordered = T, levels=rev(c("Merged 1982-2019",
                                                       "AVHRR 1982-2000",
                                                       "MODIS 2001-2019")))) %>% 
  filter(epoch!='Merged 1982-2019') %>% 
  as.data.table()
t.test(junk[epoch=='AVHRR 1982-2000']$b2, 
       junk[epoch=='MODIS 2001-2019']$b2)


junk %>% 
  ggplot(data=.,aes(b1,fill=epoch))+
  geom_density(alpha=0.25)+
  geom_vline(data=. %>% filter(epoch=="AVHRR 1982-2000"), 
             aes(xintercept=median(b1)),col='blue')+
  geom_vline(data=. %>% filter(epoch=="MODIS 2001-2019"), 
             aes(xintercept=median(b1)),col='red')+
  scale_x_continuous(limits=c(-0.005,0.005))

junk[epoch=='AVHRR 1982-2000']$b1 %>% summary
junk[epoch=='MODIS 2001-2019']$b1 %>% summary


# ndvi_hyb~
#   co2_start+
#   scale(ppet_anom)+
#   scale(p_anom)+
#   scale(pet_anom)))))

100*(1-median(junk[epoch=='MODIS 2001-2019']$b0)/median(junk[epoch=='AVHRR 1982-2000']$b0))
100*(1-median(junk[epoch=='MODIS 2001-2019']$b1)/median(junk[epoch=='AVHRR 1982-2000']$b1))
100*(1-median(junk[epoch=='MODIS 2001-2019']$b2)/median(junk[epoch=='AVHRR 1982-2000']$b2))
100*(1-median(junk[epoch=='MODIS 2001-2019']$b3)/median(junk[epoch=='AVHRR 1982-2000']$b3))



junk %>% 
  ggplot(data=.,aes(b2,fill=epoch))+
  geom_density(alpha=0.25)+
  geom_vline(data=. %>% filter(epoch=="AVHRR 1982-2000"), 
             aes(xintercept=median(b2)),col='blue')+
  geom_vline(data=. %>% filter(epoch=="MODIS 2001-2019"), 
             aes(xintercept=median(b2)),col='red')+
  scale_x_continuous(limits=c(-0.15,0.15))








tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 1982:2000] %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_start = co2 - mid_co2)] %>% 
  lm(ndvi_hyb~
       scale(co2,center = T,scale = F)+
       vpd_anom+
       p_anom+
       pet_anom, data=.) %>% 
  summary

tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 2001:2019] %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_start = co2 - mid_co2)] %>% 
   lm(ndvi_hyb~
      scale(co2,center = T,scale = F)+
      vpd_anom+
      p_anom+
      pet_anom, data=.) %>% 
  summary



b1 <- bam(ndvi_hyb~s(x,y,by=co2), 
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            as.data.table() %>% 
            .[hydro_year %in% 1982:2000], 
          select=TRUE, discrete = TRUE)
summary(b1)

b1 <- bam(ndvi_hyb~
                 te(mappet,co2_center)+
                 te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
                 te(map,frac_p_anom,k=5,bs='cs')+
                 te(mapet,frac_pet_anom,k=5,bs='cs'),
               data=tmp_ndvi_e1[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-355.2) %>% 
            as.data.table(), 
           discrete = T,select=TRUE, method='fREML')

b2 <- bam(ndvi_hyb~
            te(mappet,co2_center)+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs'),
          data=tmp_ndvi_e2[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-390.2463) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')

summary(b1)
summary(b2)



gratia::smooths(b1)
gratia::smooths(b2)
gratia::evaluate_smooth(b1,"s(x,y):co2") %>% 
  ggplot(data=.,aes(x,y,fill=est))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()

coords <- unique(tmp_ndvi[,.(x,y,mappet)])


gratia::evaluate_smooth(b1,"te(mappet,scale(co2, center = T, scale = F))")
bs <- inner_join(gratia::evaluate_smooth(b1,"te(mappet,co2_center)", newdata=coords[,`:=`(co2_center=1)]),
           gratia::evaluate_smooth(b2,"te(mappet,co2_center)", newdata=coords[,`:=`(co2_center=1)]),by=c("x","y"), 
           suffix=c("_e1","_e2"))

bs

bs %>% 
  mutate(val = est_e2-est_e1) %>% 
  pull(val) %>% 
  hist

bs %>% ggplot(data=.,aes(x,y,fill=est_e2-est_e1))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()



b3 <- bam(ndvi_hyb~
            epoch+
            te(mappet,co2_center)+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs'),
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-372.75) %>% 
            mutate(epoch = round(epoch)) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
summary(b3)


preds1 <- unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)])[,`:=`(co2_center = -20,
                          frac_vpd_anom=0,
                          frac_p_anom=0,
                          frac_pet_anom=0,
                          epoch=0)]
preds2 <- unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)])[,`:=`(co2_center = 0,
                           frac_vpd_anom=0,
                           frac_p_anom=0,
                           frac_pet_anom=0,
                           epoch=0)]

(predict(b3,newdata=preds2,type='response')-predict(b1,newdata=preds1,type='response')) -> vec_a

preds3 <- unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)])[,`:=`(co2_center = 0,
                                                                   frac_vpd_anom=0,
                                                                   frac_p_anom=0,
                                                                   frac_pet_anom=0,
                                                                   epoch=0)]
preds4 <- unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)])[,`:=`(co2_center = 20,
                                                                   frac_vpd_anom=0,
                                                                   frac_p_anom=0,
                                                                   frac_pet_anom=0,
                                                                   epoch=0)]

(predict(b3,newdata=preds4,type='response')-predict(b1,newdata=preds3,type='response')) -> vec_b


vec_b %>% hist

vec_a %>% summary
vec_b %>% summary





c1 <- bam(ndvi_cdr~
            te(mappet,co2_center)+
            te(mavpd15,scale(frac_vpd_anom,center=TRUE,scale=FALSE),k=5,bs='cs')+
            te(map,scale(frac_p_anom,center=TRUE,scale=FALSE),k=5,bs='cs')+
            te(mapet,scale(frac_pet_anom,center=TRUE,scale=FALSE),k=5,bs='cs'),
          data=tmp_ndvi_e1[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-mean(tmp_ndvi_e1$co2)) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
c2 <- bam(ndvi_mcd~
            te(mappet,co2_center)+
            te(mavpd15,scale(frac_vpd_anom,center=TRUE,scale=FALSE),k=5,bs='cs')+
            te(map,scale(frac_p_anom,center=TRUE,scale=FALSE),k=5,bs='cs')+
            te(mapet,scale(frac_pet_anom,center=TRUE,scale=FALSE),k=5,bs='cs'),
          data=tmp_ndvi_e2[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-mean(tmp_ndvi_e2$co2)) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
summary(c1)
summary(c2)

vec_a <- ((unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
  .[,`:=`(co2_center = 20,
         frac_vpd_anom=0,
         frac_p_anom=0,
         frac_pet_anom=0,
         epoch=0)] %>% predict(c1,newdata=.))-
  (unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
     .[,`:=`(co2_center = -20,
             frac_vpd_anom=0,
             frac_p_anom=0,
             frac_pet_anom=0,
             epoch=0)] %>% predict(c1,newdata=.)))

vec_b <- ((unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
             .[,`:=`(co2_center = 20,
                     frac_vpd_anom=0,
                     frac_p_anom=0,
                     frac_pet_anom=0,
                     epoch=0)] %>% predict(c2,newdata=.))-
            (unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
               .[,`:=`(co2_center = -20,
                       frac_vpd_anom=0,
                       frac_p_anom=0,
                       frac_pet_anom=0,
                       epoch=0)] %>% predict(c2,newdata=.)))

summary(vec_a)
summary(vec_b)

coords$beta_1 <- vec_a
coords$beta_2 <- vec_b




coords %>% 
    inner_join(., as_tibble(kop %>% select(x,y,zone)), by=c("x","y")) %>% 
    select(zone,beta_1,beta_2) %>% 
    gather(-zone,key='key',value='value') %>% 
  filter(is.na(zone)==F) %>% 
  ggplot(data=.,aes(x=value, 
                    y=zone, 
                    color=key))+
  geom_vline(aes(xintercept=0),color='grey')+
  geom_boxplot(outlier.colour = NA)
  scale_color_viridis_d('Climate', option='B', end=0.85)
  # scale_fill_manual('Epoch', values=c(
  #   "Merged 1982-2019"='white',
  #   "AVHRR 1982-2000"='grey80', 
  #   "MODIS 2001-2019"='grey30'))+
  # scale_x_continuous(limits=c(-0.0025,0.005))+
  # scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
  #                                       .Label = c("Equatorial",
  #                                                  "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
  #                                                  "Temperate Tas."), class = c("ordered", "factor"))))+
  # labs(y=NULL,
  #      x=expression(paste(Delta,"NDVI"~CO[2]~'ppm'**-1)))+
  # theme_linedraw()+
  # theme(panel.grid.minor = element_blank())

  
center <- function(x,na.rm=TRUE) x-mean(x,na.rm=TRUE)
########################################  
junk1 <- tmp_ndvi %>% 
  filter(hydro_year %in% 1982:2000) %>% 
  group_by(x,y,id) %>% 
  mutate_at(c("ndvi_cdr","co2", "vpd_anom","p_anom","pet_anom"),~center(.,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

junk2 <- tmp_ndvi %>% 
    filter(hydro_year %in% 2001:2019) %>% 
    group_by(x,y,id) %>% 
    mutate_at(c("ndvi_mcd","co2", "vpd_anom","p_anom","pet_anom"),~center(.,na.rm=TRUE)) %>% 
    ungroup() %>% 
    as.data.table()
########################################  

as_tibble(tmp_ndvi[id==1]) %>% 
  filter(hydro_year %in% 2001:2019) %>% 
  group_by(x,y,id) %>% 
  mutate_at(c("ndvi_mcd","co2", "vpd_anom","p_anom","pet_anom"),~center(.,na.rm=TRUE)) %>% 
  ungroup() %>% 
  select(ndvi_mcd,co2) %>% pull(ndvi_mcd) %>% plot
  
junk2[id==2]  %>% pull(ndvi_mcd) %>% plot




lm(ndvi_mcd~jitter(epoch)*co2+vpd_anom+p_anom+pet_anom,data=tmp_ndvi) %>% summary



rlm_ndvi_annual_co2_ppet_e1 <-  junk1[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  .[is.na(ndvi_cdr)==F & is.na(vpd_anom)==F] %>% 
  .[,.(beta = list(coef(lm(
    ndvi_hyb~
      scale(co2,center=T,scale=F)+
      scale(vpd_anom)+
      scale(p_anom)+
      scale(pet_anom) ))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b3=unlist(beta)[4],
          b4=unlist(beta)[5]
  ), by=.(x,y)]
rlm_ndvi_annual_co2_ppet_e2 <-  junk2[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  .[is.na(ndvi_cdr)==F & is.na(vpd_anom)==F] %>% 
  .[,.(beta = list(coef(lm(
    ndvi_mcd~
      scale(co2,center=T,scale=F)+
      scale(vpd_anom)+
      scale(p_anom)+
      scale(pet_anom) ))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b3=unlist(beta)[4],
          b4=unlist(beta)[5]
  ), by=.(x,y)]


bind_rows(#as_tibble(rlm_ndvi_annual_co2_ppet_epoch) %>% mutate(epoch='Merged 1982-2019'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e1) %>% mutate(epoch='AVHRR 1982-2000'), 
          as_tibble(rlm_ndvi_annual_co2_ppet_e2) %>% mutate(epoch='MODIS 2001-2019')) %>% 
  inner_join(., as_tibble(kop %>% select(x,y,zone)), by=c("x","y")) %>% 
  mutate(epoch=factor(epoch, ordered = T, levels=rev(c("Merged 1982-2019",
                                                       "AVHRR 1982-2000",
                                                       "MODIS 2001-2019")))) %>% 
  ggplot(data=.,aes(b1, zone, color=zone, fill=epoch))+
  geom_vline(aes(xintercept=0),color='grey')+
  geom_boxplot(outlier.colour = NA)+
  scale_color_viridis_d('Climate', option='B', end=0.85)+
  scale_fill_manual('Epoch', values=c(
    "Merged 1982-2019"='white',
    "AVHRR 1982-2000"='grey80', 
    "MODIS 2001-2019"='grey30'))+
  scale_x_continuous(limits=c(-0.0025,0.005))+
  scale_y_discrete(limits=rev(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
                                        .Label = c("Equatorial",
                                                   "Tropical", "Subtropical", "Grassland", "Arid", "Temperate",
                                                   "Temperate Tas."), class = c("ordered", "factor"))))+
  labs(y=NULL,
       x=expression(paste(Delta,"NDVI"~CO[2]~'ppm'**-1)))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())




b4 <- bam(ndvi_hyb~
            s(mappet,k=5,bs='cs')+
            epoch*co2_center+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs'),
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-372.75) %>% 
            mutate(epoch = factor(round(epoch))) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
summary(b4)


predict(b4, newdata=tmp_ndvi[id==1][20,][,`:=`(epoch=1,co2_center=20)])-
  predict(b4, newdata=tmp_ndvi[id==1][20,][,`:=`(epoch=1,co2_center=0)])

predict(b4, newdata=tmp_ndvi[id==1][20,][,`:=`(epoch=0,co2_center=0)])-
  predict(b4, newdata=tmp_ndvi[id==1][20,][,`:=`(epoch=0,co2_center=-20)])




b5 <- bam(ndvi_hyb~
            te(mappet,co2_center,bs='gp')+
            epoch+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs'),
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-370.32) %>% 
            mutate(epoch = factor(round(epoch))) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
summary(b5)
plot(b5,scheme=1,scale = 0)


bad1 <- bam(ndvi_hyb~
            s(mappet,k=5,bs='cs')+
            epoch+
            s(co2_center,k=5)+
            te(map,frac_p_anom,k=5,bs='cs'),
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            lazy_dt() %>% 
            mutate(co2_center = co2-372.75) %>% 
            mutate(epoch = factor(round(epoch))) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
summary(bad1)
plot(bad1,scheme=1,scale = 0)


gratia::evaluate_smooth(object = b5, 'te(mappet') %>% 
  filter(mappet<=2) %>% 
  # filter(between(mappet,0.9,1.1)) %>% 
  mutate(mappet = cut_interval(mappet,n=5)) %>% 
  ggplot(data=.,aes(co2_center, est,color=cut_interval(co2_center,2)))+
  geom_smooth(method='lm')+
  facet_wrap(~mappet,ncol = 1,scales = 'free')



b5 <- bam(ndvi_hyb~
            s(x,y,fx = TRUE)+
            te(mappet,co2_center,bs='gp')+
            epoch+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs'),
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-370.32) %>% 
            mutate(epoch = factor(round(epoch))) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')

gratia::evaluate_smooth(object = b5, 'te(mappet') %>% 
  filter(mappet<=2) %>% 
  # filter(between(mappet,0.9,1.1)) %>% 
  mutate(mappet = cut_interval(mappet,n=5)) %>% 
  ggplot(data=.,aes(co2_center, est,color=cut_interval(co2_center,2)))+
  geom_smooth(method='lm')+
  facet_wrap(~mappet,ncol = 1,scales = 'free')



b6 <- bam(ndvi_hyb~
            s(x,y,fx = TRUE)+
            te(mappet,co2_center,k=2,bs='gp')+
            epoch+
            te(mavpd15,frac_vpd_anom,k=2,bs='cs')+
            te(map,frac_p_anom,k=2,bs='cs')+
            te(mapet,frac_pet_anom,k=2,bs='cs'),
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-370.32) %>% 
            mutate(epoch = factor(round(epoch))) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
summary(b6)
gratia::evaluate_smooth(object = b6, 'te(mappet') %>% 
  filter(mappet<=2) %>% 
  # filter(between(mappet,0.9,1.1)) %>% 
  mutate(mappet = cut_interval(mappet,n=5)) %>% 
  ggplot(data=.,aes(co2_center, est,color=cut_interval(co2_center,2)))+
  geom_smooth(method='lm')+
  facet_wrap(~mappet,ncol = 1,scales = 'free')





# rlm_ndvi_annual_co2_ppet_epoch <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
#   as.data.table() %>% 
#   .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
#   .[,`:=`(co2_center = co2 - mid_co2)] %>% 
#   .[,.(beta = list(coef(MASS::rlm(
#     ndvi_hyb~
#       co2+
#       vpd_anom+
#       p_anom+
#       pet_anom+
#       jitter(epoch)))))
#     ,
#     by=.(x,y)] %>%
#   .[,`:=`(b0=unlist(beta)[1], 
#           b1=unlist(beta)[2], 
#           b2=unlist(beta)[3], 
#           b3=unlist(beta)[4],
#           b4=unlist(beta)[5], 
#           b5=unlist(beta)[6]
#   ), by=.(x,y)]

b7 <- bam(ndvi_hyb~
            vpd_anom+
            ti(x,y,by=co2)+
            ti(x,y,by=vpd_anom)+
            ti(x,y,by=p_anom)+
            ti(x,y,by=pet_anom)+
            epoch,
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-370.32) %>% 
            mutate(epoch = factor(round(epoch))) %>% 
            as.data.table(), 
          discrete = T,select=TRUE, method='fREML')
summary(b7)
plot(b7)
gratia::smooths(b7)
gratia::evaluate_smooth(b7,'ti(x,y):co2',newdata=unique(rlm_ndvi_annual_co2_ppet_epoch[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., rlm_ndvi_annual_co2_ppet_epoch,by=c('x','y')) %>% 
  ggplot(data=.,aes(x,y,fill=est))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2(limits=c(-1,1))

gratia::evaluate_smooth(b7,'ti(x,y):co2',newdata=unique(rlm_ndvi_annual_co2_ppet_epoch[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., rlm_ndvi_annual_co2_ppet_epoch,by=c('x','y')) %>% 
  ggplot(data=.,aes(b1,est))+
  geom_point()+
  geom_smooth()

gratia::evaluate_smooth(b7,'ti(x,y):vpd_anom',newdata=unique(rlm_ndvi_annual_co2_ppet_epoch[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., rlm_ndvi_annual_co2_ppet_epoch,by=c('x','y')) %>% 
  ggplot(data=.,aes(b2,est))+
  geom_point()+
  geom_smooth(method='lm')

gratia::evaluate_smooth(b7,'ti(x,y):p_anom',newdata=unique(rlm_ndvi_annual_co2_ppet_epoch[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., rlm_ndvi_annual_co2_ppet_epoch,by=c('x','y')) %>% 
  ggplot(data=.,aes(b3,est))+
  geom_point()+
  geom_smooth(method='lm')

gratia::evaluate_smooth(b7,'ti(x,y):pet_anom',newdata=unique(rlm_ndvi_annual_co2_ppet_epoch[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., rlm_ndvi_annual_co2_ppet_epoch,by=c('x','y')) %>% 
  ggplot(data=.,aes(b4,est))+
  geom_point()+
  geom_smooth(method='lm')






test <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>%
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>%
  .[,`:=`(co2_center = co2 - mid_co2)] %>%
  .[,.(beta = list(coef(lm(
    ndvi_hyb~
      vpd_anom+
      p_anom+
      jitter(epoch)))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1],
          b1=unlist(beta)[2],
          b2=unlist(beta)[3]
          # b3=unlist(beta)[4],
          # b4=unlist(beta)[5],
          # b5=unlist(beta)[6]
  ), by=.(x,y)]

btest <- bam(ndvi_hyb~
            ti(x,y)+
            ti(x,y,vpd_anom)+
            ti(x,y,p_anom)+
            epoch,
          data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
            mutate(co2_center = co2-370.32) %>% 
            mutate(epoch = factor(round(epoch))) %>% 
            as.data.table(), 
          discrete = T,select=T, method='fREML')
summary(btest)

gratia::evaluate_smooth(btest,gratia::smooths(btest)[1])
gratia::evaluate_smooth(btest,'ti(x,y)',newdata=unique(test[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., test,by=c('x','y')) %>% 
  ggplot(data=.,aes(b0,est))+
  geom_point()+
  geom_smooth(method='lm')

gratia::evaluate_smooth(btest,'ti(x,y,vpd_anom',newdata=unique(test[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., test,by=c('x','y')) %>% 
  ggplot(data=.,aes(b1,est))+
  geom_point()+
  geom_smooth(method='lm')

gratia::evaluate_smooth(btest,'te(x,y):p_anom',newdata=unique(test[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., test,by=c('x','y')) %>% 
  ggplot(data=.,aes(b2,est))+
  geom_point()+
  geom_smooth(method='lm')

gratia::evaluate_smooth(btest,'te(x,y):vpd_anom',newdata=unique(test[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., test,by=c('x','y')) %>% 
  ggplot(data=.,aes(x,y,fill=b1))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()
gratia::evaluate_smooth(btest,'te(x,y):p_anom',newdata=unique(test[,.(x,y)])) %>% 
  select(x,y,est) %>% 
  inner_join(., test,by=c('x','y')) %>% 
  ggplot(data=.,aes(x,y,fill=b2))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()



btest <- bam(ndvi_hyb~
               te(map,by=p_anom)+
               epoch,
             data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
               as_tibble() %>% 
               filter(map<=2500) %>% 
               mutate(epoch = factor(round(epoch))) %>% 
               mutate(M = cbind(x,y)),             
             discrete = T,select=TRUE, method='fREML')
plot(btest)
summary(btest)
gratia::evaluate_smooth(btest,gratia::smooths(btest))


bam(ndvi_hyb~
      p_anom+
      epoch,
    data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
      as_tibble() %>% 
      filter(map>=1500 & map<=2500) %>% 
      mutate(epoch = factor(round(epoch)))) %>% summary

bam(ndvi_hyb~
      p_anom:map+
      epoch,
    data={tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
      as_tibble() %>% 
      filter(map<2000) %>% 
      mutate(map = cut_interval(map,10)) %>% 
      mutate(epoch = factor(round(epoch)))}) %>% summary




b8 <- bam(ndvi_hyb~
            te(x,y,fx = T)+
            te(mappet,co2,k=5,bs='cs')+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs')+
               epoch,
             data=tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
              lazy_dt() %>% 
            filter(mappet<=2) %>% 
           mutate(co2_center = co2-370.32) %>% 
               mutate(epoch = factor(round(epoch))) %>% 
               as.data.table(), 
             discrete = T,select=T, method='fREML')
summary(b8)
plot(b8,scale=0,scheme = 2)


library(mgcViz)
plot(sm(getViz(b8),1))+l_fitRaster()+l_fitContour()+
  scale_fill_gradient2()

library(gratia)
evaluate_smooth(b8,"te(mappet") %>% 
  mutate(mappet_d = cut_interval)
  ggplot(data=.,aes(mappet,est,color=co2))+
  geom_point()+
  scale_color_viridis_c()

vec_e1 <- predict(b8,newdata=unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
      .[,`:=`(epoch=1,frac_vpd_anom=0,frac_p_anom=0,
              frac_pet_anom=0, co2=370)])-
  predict(b8,newdata=unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
            .[,`:=`(epoch=1,frac_vpd_anom=0,frac_p_anom=0,
                    frac_pet_anom=0, co2=330)])  
vec_e2 <- predict(b8,newdata=unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
                    .[,`:=`(epoch=1,frac_vpd_anom=0.05,frac_p_anom=0,
                            frac_pet_anom=0, co2=400)])-
  predict(b8,newdata=unique(tmp_ndvi[,.(x,y,mappet,mavpd15,map,mapet)]) %>% 
            .[,`:=`(epoch=1,frac_vpd_anom=-0.05,frac_p_anom=0,
                    frac_pet_anom=0, co2=370)])  

summary(vec_e1)
summary(vec_e2)  

plot(vec_e2~vec_e1); abline(0,1,col='red')
  
tmp_ndvi %>% sa
mple_n(1)
predict(b8, 
        tibble(mavpd15=3.14, 
               mappet=0.8,map=1055,mapet=1365,frac_vpd_anom=0,frac_p_anom=0,frac_pet_anom=0, 
               epoch=0,
       co2=370))-
  predict(b8, 
          tibble(mavpd15=3.14, 
                 mappet=0.8,map=1055,mapet=1365,frac_vpd_anom=0,frac_p_anom=0,frac_pet_anom=0, 
                 epoch=1,
                 co2=340))

predict(b8, 
        tibble(mavpd15=3.14, 
               mappet=0.8,map=1055,mapet=1365,frac_vpd_anom=0,frac_p_anom=0,frac_pet_anom=0, 
               epoch=0,
               co2=400, hydro_year = 2018))-
  predict(b8, 
          tibble(mavpd15=3.14, 
                 mappet=0.8,map=1055,mapet=1365,frac_vpd_anom=0,frac_p_anom=0,frac_pet_anom=0, 
                 epoch=1,
                 co2=370))

       



junk <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_center = co2 - mid_co2)] %>% 
  .[,.(beta = list(coef(MASS::rlm(
    ndvi_hyb~
      co2*hydro_year+
      vpd_anom+
      p_anom+
      pet_anom+
      jitter(epoch)))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b3=unlist(beta)[4],
          b4=unlist(beta)[5], 
          b5=unlist(beta)[6]
  ), by=.(x,y)]

fit <- lm(ndvi_hyb~co2*hydro_year_c*mappet+
     frac_vpd_anom+
     frac_p_anom+
     frac_pet_anom+
     epoch,data=tmp_ndvi)
tmp_ndvi$hydro_year_c %>% summary
411*37*2.007e-05 - 9.1e-4
0.02777


low <- predict(fit, 
  newdata=expand_grid(frac_vpd_anom=0,frac_p_anom=0,frac_pet_anom=0,epoch=0,co2=340,hydro_year_c=0,mappet=1))
mid <- predict(fit, 
                newdata=expand_grid(frac_vpd_anom=0,frac_p_anom=0,frac_pet_anom=0,epoch=0,co2=371,hydro_year_c=18,mappet=1))
high <- predict(fit, 
        newdata=expand_grid(frac_vpd_anom=0,frac_p_anom=0,frac_pet_anom=0,epoch=0,co2=410,hydro_year_c=37,mappet=1))


high-mid
mid-low
