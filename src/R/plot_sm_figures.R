# !diagnostics off
#*******************************************************************************
#* Description:
#* Plot Climate stuff
#*
#*
#*
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
dt <- as.data.table
tb <- as_tibble
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
arrow::write_parquet(kop, sink='../data_general/Koppen_climate/BOM_Koppen_simplified7.parquet')
#*** End Kop zone load ********************************************************

# MODIS Vegetation continuous cover data ----------------------------------------
mod_tree <- stars::read_stars("../data_general/Oz_misc_data/MOD44BPercent_Tree_Cover_5000m_East_Oz_noMask_2000_2019.tif") %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-01-01"),ymd("2019-01-01"),by="1 year"), 
                    names = 'date') %>% 
  set_names(c("tree_cover")) %>% 
  as_tibble() %>% 
  as.data.table()

mod_nontree <- stars::read_stars("../data_general/Oz_misc_data/MOD44BPercent_NonTree_Vegetation_5000m_East_Oz_noMask_2000_2019.tif") %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-01-01"),ymd("2019-01-01"),by="1 year"), 
                    names = 'date') %>% 
  set_names(c("nontree_cover")) %>% 
  as_tibble() %>% 
  as.data.table()

mod_nonveg <- stars::read_stars("../data_general/Oz_misc_data/MOD44BPercent_NonVegetated_5000m_East_Oz_noMask_2000_2019.tif") %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-01-01"),ymd("2019-01-01"),by="1 year"), 
                    names = 'date') %>% 
  set_names(c("nonveg_cover")) %>% 
  as_tibble() %>% 
  as.data.table()

mod <- merge(mod_tree, mod_nontree, by = c("x","y","date"))
mod <- merge(mod, mod_nonveg, by = c("x","y","date"))
rm(mod_tree, mod_nontree, mod_nonveg); gc()

# add the NVIS vegetation classes
# base <- stars::read_stars('../data_general/AVHRR_CDRv5_VI/AVHRR_SR_median_EastOz_1982_2019.tif', 
#                           RasterIO = list(bands=1))
base <- stars::read_stars("../data_general/Oz_misc_data/MOD44BPercent_Tree_Cover_5000m_East_Oz_noMask_2000_2019.tif", 
                          RasterIO = list(bands=1))
nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif")
nvis2 <- st_warp(src=nvis, dest=base[,,], use_gdal = T)
names(nvis2) <- "veg_class"
nvis <- nvis2 %>% as_tibble() %>% as.data.table()
codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
                         fwf_widths(c(2,100)), skip = 1) %>% 
  set_names(c("veg_class","veg_class_descrip")) %>% 
  mutate(vc = as.factor(veg_class_descrip)) 
vc <- left_join(nvis, codes, by='veg_class')
mod <- vc[mod, on=.(x,y)]
mod[,`:=`(year=year(date))]

# mod <- svi[mod, on=.(x,y,year)]
# end section ******************************************************************
# END DATA IMPORT SECTION ******************************************************



# Calc NDVI linear change by season ***********************************************************
library(RcppArmadillo)
system.time(
  lt_ndvi_season <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
    .[date>= ymd("1981-09-01") & date<= ymd("2019-08-30")] %>% 
    .[nobs_total > 200] %>% 
    .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
    .[is.na(val)==F] %>% 
    .[,.(b1 = fastLm(X = cbind(1,hydro_year-2000.5), y=val, data=.SD)$coefficients[2]), 
      by=.(x,y,season)]
)
system.time(
  lt_ndvi_season_wEpoch <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
    .[date>= ymd("1981-09-01") & date<= ymd("2019-08-30")] %>% 
    .[nobs_total > 200] %>% 
    .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
    .[,`:=`(epoch=ifelse(hydro_year < 2001,0,1))] %>% 
    .[is.na(val)==F] %>% 
    .[,.(b1 = fastLm(X = cbind(1,hydro_year-2000.5,epoch), y=val, data=.SD)$coefficients[2]), 
      by=.(x,y,season)]
)
# END **************************************************************************



# REGRESSIONS *** ------------------------------------------------------------------
# Linear change in MODIS VCF  ---------------------------------------------------
library(RcppArmadillo)
system.time(
  lt_tree_sen <- mod[year<=2018][,`:=`(year_c = year-2009.5)] %>% 
    .[,.(beta = list(unname(zyp.sen(tree_cover~year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]  
)
system.time(
  lt_nontree_sen <- mod[year<=2018][,`:=`(year_c = year-2009.5)] %>% 
    .[,.(beta = list(unname(zyp.sen(nontree_cover~year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]  
)
system.time(
  lt_nonveg_sen <- mod[year<=2018][,`:=`(year_c = year-2009.5)] %>% 
    .[,.(beta = list(unname(zyp.sen(nonveg_cover~year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]  
)



vcf_sen <- inner_join(lt_nontree_sen %>% rename(grass_u=b0, grass_b=b1) %>% select(-beta), 
                      lt_nonveg_sen %>% rename(nonveg_u=b0, nonveg_b=b1) %>% select(-beta))
vcf_sen <- lt_tree_sen %>% lazy_dt() %>% 
  select(-beta) %>% 
  rename(tree_u=b0, 
         tree_b=b1) %>% 
  as.data.table() %>% 
  inner_join(., vcf_sen)

# NDVI 38 year seasonal trend -----------------------------------------------------------
# OLS
system.time(
  lt_ndvi_season <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
    .[nobs_total > 200] %>% 
    .[ndvi_hyb > 0] %>% 
    .[date>= ymd("1981-09-01") & date<= ymd("2019-08-30")] %>% 
    .[,.(val = mean(ndvi_hyb, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
    .[is.na(val)==F] %>% 
    .[,.(b1 = fastLm(X = cbind(1,hydro_year-2000.5), y=val, data=.SD)$coefficients[2]), 
      by=.(x,y,season)]
)


# RLM 
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
# END RLM NDVI ANNUAL **********************************************************


# NDVI trend using Theil-Sen --------------------------------------------------------
sen_ndvi_annual <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-31")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE)),by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

sen_ndvi_annual %>% 
  filter(b0>0) %>% 
  mutate(val = 100*b1*38/b0) %>% 
  pull(val) %>% na.omit() %>% 
  quantile(., c(0.05,0.5,0.95))

sen_ndvi_season <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-31")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y,season)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y,season)]
sen_ndvi_annual_e1 <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2000-12-31")] %>% 
  .[hydro_year %in% c(1982:2000)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE)),by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]
sen_ndvi_season_e1 <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2000-12-31")] %>% 
  .[hydro_year %in% c(1982:2000)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y,season)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y,season)]
sen_ndvi_annual_e2 <- dat %>% 
  .[date >= ymd("2001-01-01") & date <= ymd("2019-08-30")] %>% 
  .[hydro_year %in% c(2001:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-2001)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE)),by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

sen_ndvi_season_e2 <- dat %>% 
  .[date >= ymd("2001-01-01") & date <= ymd("2019-08-30")] %>% 
  .[hydro_year %in% c(2001:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-2001)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y,season)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y,season)]

sen_ndvi_season_e1$epoch <- "AVHRR NDVI 1982-2000"
sen_ndvi_season_e2$epoch <- "MODIS NDVI 2001-2019"










sen_ndvi_annual_e1$b1 %>% na.omit() %>% quantile(., c(0.05,0.5,0.95))
sen_ndvi_annual_e2$b1 %>% na.omit() %>% quantile(., c(0.05,0.5,0.95))



# P:PET annual mean trend w/Thiel Sen ------------------------------------------
sen_ppet_annual <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-31")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(pe_12mo)==F] %>% 
  .[,.(ppet = mean(pe_12mo,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ppet~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

sen_ppet_annual_e1 <- dat %>% 
  .[date >= ymd("1981-12-01") & date <= ymd("2000-11-30")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(pe_12mo)==F] %>% 
  .[,.(ppet = mean(pe_12mo,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ppet~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

sen_ppet_annual_e2 <- dat %>% 
  .[date >= ymd("2000-12-01") & date <= ymd("2019-11-30")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(pe_12mo)==F] %>% 
  .[,.(ppet = mean(pe_12mo,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ppet~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

# PET annual mean trend w/Thiel Sen ------------------------------------------
sen_pet_annual <- dat %>% 
  .[date >= ymd("1981-12-01") & date <= ymd("2019-11-30")] %>% # need full year for total
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(pet)==F] %>% 
  .[,.(pet_tot = sum(pet,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(pet_tot~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

# Precip annual mean trend w/Thiel Sen ------------------------------------------
sen_p_annual <- dat %>% 
  .[date >= ymd("1981-12-01") & date <= ymd("2019-11-30")] %>% # need full year for total
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(precip)==F] %>% 
  .[,.(p_tot = sum(precip,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(p_tot~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

# NDVI 38 year annual trend of ----------------------------------------------
lt_ndvi_annual <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
  .[date>= ymd("1981-09-01") & date<= ymd("2019-08-30")] %>% 
  .[nobs_total > 200] %>% 
  # .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
  .[,`:=`(epoch=ifelse(hydro_year < 2001,0,1))] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(beta = list(unname(fastLm(
    X = cbind(1,hydro_year-1982,epoch), 
    y=ndvi_hyb, data=.SD)$coefficients))), 
    by=.(x,y)] %>% 
  .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3]#,b3=unlist(beta)[4]
  ), by=.(x,y)]

lt_ndvi_annual <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
  .[nobs_total > 200] %>% 
  .[date>= ymd("1981-09-01") & date<= ymd("2019-08-30")] %>% 
  # .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[, .(val=mean(ndvi_hyb,na.rm=TRUE)), 
    keyby=.(x,y,hydro_year)] %>% 
  .[,`:=`(epoch=ifelse(hydro_year < 2001,0,1))] %>% 
  .[,.(beta = list(unname(fastLm(
    X = cbind(1,hydro_year-1982,epoch), 
    y=val, data=.SD)$coefficients))), 
    by=.(x,y)] %>% 
  .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3]), by=.(x,y)]

lt_ndvi_annual_v2 <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
  .[nobs_total > 200] %>% 
  .[date>= ymd("1981-09-01") & date<= ymd("2019-08-30")] %>% 
  # .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[, .(val=mean(ndvi_hyb,na.rm=TRUE)), 
    keyby=.(x,y,hydro_year)] %>% 
  .[,`:=`(epoch=ifelse(hydro_year < 2001,0,1))] %>% 
  .[,.(beta = list(unname(fastLm(
    X = cbind(1,hydro_year-1982), 
    y=val, data=.SD)$coefficients))), 
    by=.(x,y)] %>% 
  .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]

# end section ******************************************************************

# Calc NDVI linear change by satellite epoch -----------------------------------------------------------
c_year <- mean(seq(ymd("1982-01-01"),ymd("2000-12-31"),by='1 month')) %>% decimal_date();
lt_ndvi_season_p1 <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
  .[nobs_total > 200] %>% 
  .[ndvi_hyb > 0] %>% 
  .[date>= ymd("1981-01-01") & date<= ymd("2000-12-31")] %>% 
  .[,.(val = mean(ndvi_hyb, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
  .[is.na(val)==F] %>% 
  .[,.(b1 = fastLm(X = cbind(1,hydro_year-c_year), y=val, data=.SD)$coefficients[2]), 
    by=.(x,y,season)]

c_year <- mean(seq(ymd("2001-07-01"),ymd("2019-08-30"),by='1 month')) %>% decimal_date()
lt_ndvi_season_p2 <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
  .[nobs_total > 200] %>% 
  .[ndvi_hyb > 0] %>% 
  .[date>= ymd("2001-07-01") & date<= ymd("2019-08-30")] %>% 
  .[,.(val = mean(ndvi_hyb, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
  .[is.na(val)==F] %>% 
  .[,.(b1 = fastLm(X = cbind(1,hydro_year-c_year), y=val, data=.SD)$coefficients[2]), 
    by=.(x,y,season)]

lt_ndvi_season_p1$epoch <- "AVHRR NDVI 1982-2000"
lt_ndvi_season_p2$epoch <- "MODIS NDVI 2001-2019"




#QUESTIONS *** ---------------------------------------------------------------------
#Q: % of pixels exper. greater aridity? ------------------------------
# overall 1982-2019
sen_ppet_annual %>% 
  as_tibble() %>% 
  filter(b0>0) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(arid_frac = 100*neg/nobs)

# 1982-2000 
sen_ppet_annual_e1 %>% 
  as_tibble() %>% 
  filter(b0>0) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(arid_frac = 100*neg/nobs)

# 2001-2019
sen_ppet_annual_e2 %>% 
  as_tibble() %>% 
  filter(b0>0) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(arid_frac = 100*neg/nobs)

# How much did relative P:PET change? 
sen_ppet_annual %>% 
  as_tibble() %>% 
  filter(b0>0) %>% 
  mutate(rel_change = 100*(b1*37)/b0) %>% 
  pull(rel_change) %>% 
  quantile(., c(0.05,0.5,0.95))

sen_ppet_annual_e1 %>% 
  as_tibble() %>% 
  filter(b0>0) %>% 
  mutate(rel_change = 100*(b1*18)/b0) %>% 
  pull(rel_change) %>% 
  quantile(., c(0.05,0.5,0.95))

sen_ppet_annual_e2 %>% 
  as_tibble() %>% 
  filter(b0>0) %>% 
  mutate(rel_change = 100*(b1*18)/b0) %>% 
  pull(rel_change) %>% 
  quantile(., c(0.05,0.5,0.95))



# How many pixels experienced less rainfall?
sen_p_annual %>% 
  as_tibble() %>% 
  filter(b0>0) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(frac = 100*neg/nobs)


# What % of pixels experienced greater PET?
sen_pet_annual %>% 
 as_tibble() %>% 
  filter(b0>0) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(frac = 100*pos/nobs)

sen_pet_annual %>% 
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(b0>100) %>% 
  mutate(rel_change = 100*(b1*37)/b0) %>% 
  pull(rel_change) %>% 
  quantile(., c(0.05,0.5,0.95))


# what % of pixels experienced greening by satellite epoch?
sen_ndvi_annual_e1 %>% # AVHRR epoch
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(green_frac = 100*pos/nobs)
sen_ndvi_annual_e2 %>% # MODIS epoch
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(green_frac = 100*pos/nobs)


#Q: What percent of pixels experienced greening? --------------------------
# 1982-2019 using rlm and accounting for the change in sensor
rlm_ndvi_annual %>% 
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(green_frac = 100*pos/nobs)

# 1982-2019 using Thiel Sen, although this doesn't work because it 
# does not account for the sensor difference
sen_ndvi_annual %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(green_frac = 100*pos/nobs)

# AHVRR epoch
sen_ndvi_annual_e1 %>% 
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(green_frac = 100*pos/nobs)

# MODIS epoch 
sen_ndvi_annual_e2 %>% 
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  summarize(pos = sum(b1>0), 
            neg = sum(b1<0), 
            nobs = sum(is.na(b1)==F)) %>% 
  mutate(green_frac = 100*pos/nobs)


#Q: How much relative greening? --------------------------
rlm_ndvi_annual %>% 
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.05,1)) %>% 
  filter(between(b1,-0.1,0.1)) %>% 
  mutate(val = 100*(b1*37)/b0) %>% 
  pull(val) %>% 
  quantile(., c(0.05,0.5,0.95))
sen_ndvi_annual_e1 %>% 
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  mutate(val = 100*(b1*18)/b0) %>% 
  pull(val) %>% 
  quantile(., c(0.05,0.5,0.95))
sen_ndvi_annual_e2 %>% 
  as_tibble() %>% 
  filter(is.na(b0)==F) %>% 
  filter(between(b0,0.1,1)) %>% 
  mutate(val = 100*(b1*18)/b0) %>% 
  pull(val) %>% 
  quantile(., c(0.05,0.5,0.95))


# WUE model comparison with NDVI observations -----------------------------
# RLM 
tmp_ndvi <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-30")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982, 
          epoch = ifelse(date<=ymd("2000-12-31"),0,1),
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
    by=.(x,y,hydro_year_c)] %>% 
  group_by(x,y) %>% 
  mutate(id = cur_group_id()) %>%
  ungroup() %>% 
  as.data.table()
tmp_annual_nobs <- tmp_ndvi %>% 
  lazy_dt() %>% 
  group_by(id) %>% 
  summarize(obs_annual = sum(nobs >= 6)) %>% 
  as.data.table()
vec_ids <- tmp_annual_nobs %>% 
  filter(obs_annual >= 25) %>% 
  pull(id)

# for(i in 1:length(vec_ids)){ # Debugging
#   tmp_ndvi[id%in%vec_ids[i]] %>% 
#     .[,.(beta = list(coef(MASS::rlm(ndvi_hyb~co2+frac_p_anom+frac_ppet_anom+epoch)))),by=.(x,y)] %>%
#     .[,`:=`(b0=unlist(beta)[1], 
#             b1=unlist(beta)[2], 
#             b2=unlist(beta)[3], 
#             b4=unlist(beta)[4], 
#             b5=unlist(beta)[5]), by=.(x,y)]
#   # print(i)
# }
# tmp_ndvi[id%in%vec_ids[i]] 

tmp_ndvi[id%in%sample(vec_ids,50)] %>% 
  ggplot(data=.,aes(frac_ppet_anom,ndvi_hyb,group=id))+
  geom_smooth(method='lm', formula=y~poly(x,2))

tmp_ndvi[id%in%sample(vec_ids,50)] %>% 
  ggplot(data=.,aes(frac_p_anom,ndvi_hyb,group=id))+
  geom_point()+
  geom_smooth(method='lm',se=F)+
  geom_abline(aes(intercept=0,slope=1),col='red')

rlm_ndvi_annual_p_co2_epoch <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
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
          b2=unlist(beta)[3] 
          # b4=unlist(beta)[4]
          ), by=.(x,y)]
rlm_ndvi_annual_p_co2_epoch


lt_v_sen <- dat[date>=ymd("1981-12-01")][date<=ymd("2019-08-30")] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(vpd15)==F] %>% 
  .[,.(vpd15 = mean(vpd15,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(vpd15~hydro_year_c)))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]

# Percent increase in Ca
dCa_Ca <- 
  diff(range(mlo[date>=ymd("1981-12-01") & date <= ymd("2019-09-30")]$co2_trend))/
  mean(mlo[date>=ymd("1982-01-01") & date <= ymd("1983-01-01")]$co2_trend)

# Percent increase of VPD
dVPD_VPD_sen <- mean(37*lt_v_sen$b1,na.rm=TRUE)/mean(lt_v_sen$b0,na.rm=TRUE)

# Percent increase in NDVI attributable to CO2 
inner_join(tmp_ndvi, rlm_ndvi_annual_p_co2_epoch, by=c("x","y")) %>% 
  inner_join(., lt_v_sen[,.(x,y,b0,b1)], by=c("x","y"), suffix=c("","_vpd")) %>% 
  lazy_dt() %>% 
  filter(hydro_year_c==0) %>% 
  mutate(rel_ndvi_co2 = 37*b1/ndvi_hyb, 
         dVPD_VPD = 37*b1_vpd/b0_vpd) %>% 
  mutate(e_wue = 0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
  as_tibble() %>% 
  select(rel_ndvi_co2, e_wue) %>% 
  filter(is.na(rel_ndvi_co2)==F) %>% 
  summarize_all(., quantile, c(0.05,0.5,0.95))
  # ggplot(data=.,aes(e_wue, rel_ndvi_co2))+
  # geom_point(size=0.1,alpha=0.1)+
  # geom_abline(aes(intercept=0,slope=1),col='red')+
  # geom_smooth(method='lm')+
  # scale_x_continuous(limits=c(0.05,0.15))+
  # scale_y_continuous(limits=c(-0.2,0.2))+
  # scale_color_viridis_c(option='b')

#%>% 
  # pull(rel) %>% 
  # na.omit() %>% 
  # quantile(., c(0.05,0.5,0.95))
# (38*rlm_ndvi_annual_p_co2_epoch$b1)/(rlm_ndvi_annual_p_co2_epoch$b0)

# Percent increase in NDVI
dNDVI_NDVI_sen <- mean(38*lt_ndvi_year$b1,na.rm=TRUE)/mean(lt_ndvi_year$b0,na.rm=TRUE)
mean(37*rlm_ndvi_annual$b1)/mean(rlm_ndvi_annual$b0)
median(37*rlm_ndvi_annual$b1)/median(rlm_ndvi_annual$b0)

# Expected WUE related increase (Donohue 2013)
0.5*(dCa_Ca - 0.5*dVPD_VPD_sen)

# Actual percent relative increase in NDVI
dNDVI_NDVI_sen



# SUPPLEMENTAL FIGURES ---------------------------------------------------------
# SM Fig 1: MA P:PET & NVIS -----------------------------------------------------

p_vc <- dat[date==ymd("2000-01-01")] %>%
  .[veg_class %in% c(1:3,5:9,11,12,13,14)] %>%
  ggplot(data=., aes(x,y, fill=str_wrap(vc, 20)))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray70',color='gray10')+
  geom_tile()+
  labs(x=NULL,y=NULL)+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  # guides(fill=guide_colorbar(barwidth = 0.5))+
  # scale_fill_brewer(palette = 'Paired', direction = -1)+
  scico::scale_fill_scico_d("NVIS 5.1 Major Veg. Class",
                            direction=-1, palette = 'batlow',end=0.95)+
  theme(panel.background = element_rect(fill = '#99A3C4'), 
        panel.grid = element_blank(), 
        legend.position = 'right',
        # legend.text = element_text(size=5),
        legend.title = element_text(size=8),
        # legend.key.width = unit(0.25,'cm'),
        strip.text = element_text(face='bold'),
        axis.text = element_blank(), 
        axis.ticks = element_blank()); p_vc
# ggsave(p_vc, 
#        filename = 'figures/map_eastOz_NVIS.png', 
#        width=10, height = 15, units='cm', dpi = 'retina', type='cairo')

# Mean Annual P:PET ****
mape <- dat %>% lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(mape = mean(pe_12mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble()
p_mappet <- mape %>% 
  ggplot(data=., aes(x,y,fill=mape))+
  geom_sf(inherit.aes = F, data=oz_poly,
          fill='gray70',color='gray10')+
  geom_tile()+
  coord_sf(xlim = c(140,155),
           ylim = c(-45,-10), 
           expand = FALSE)+
  scale_x_continuous(breaks=seq(140,155,by=5))+
  scale_fill_viridis_c(option='B',direction = -1,end=0.95, 
                       limits=c(0,2),
                       oob=scales::squish, 
                       breaks=c(0,0.5,1,1.5,2),
                       labels=c("0","0.5","1","1.5","2+"))+
  # scico::scale_fill_scico_d(end=0.9,direction = 1)+
  labs(x=NULL,y=NULL)+
  theme_linedraw()+
  guides(fill=guide_colorbar(title='Mean Annual P:PET',
                             title.position = 'top'))+
  theme(panel.background = element_rect(fill = '#99A3C4'),
        # legend.position = c(0.7,0.925),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.direction = 'horizontal',
        legend.background = element_rect(fill='white'),
        panel.grid = element_blank())
p_mappet+p_vc+plot_annotation(tag_levels = 'a',tag_prefix = '(',tag_suffix = ')')
ggsave(
  filename = 'figures/SM_Fig1_map_eastOz_MAPPET_NVIS.png', 
  width=20, height = 20, units='cm', dpi = 'retina', type='cairo')



# SM Fig 2: Long-term change P, PET, PE ---------------------------------------------
# .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y)] %>%

lt_P_season_sen <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-31")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(precip)==F] %>% 
  .[,.(val = mean(precip,na.rm=TRUE)), by=.(x,y,season,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(val~hydro_year_c)))),by=.(x,y,season)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y,season)]

lt_PET_season_sen <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2019-08-31")] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(pet)==F] %>% 
  .[,.(val = mean(pet,na.rm=TRUE)), by=.(x,y,season,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(val~hydro_year_c)))),by=.(x,y,season)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y,season)]

p_P <- lt_P_season_sen %>% 
  ggplot(data=., aes(x,y,fill=b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray70',color='gray10')+
  geom_tile()+
  scale_fill_gradient2(expression(atop("Monthly P  ",(mm~yr**-1))),
                       limits=c(-3,3),
                       oob=scales::squish, 
                       na.value = 'gray')+
  labs(x=NULL,y=NULL)+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  facet_wrap(~season, nrow = 1)+
  theme(panel.background = element_rect(fill = '#99A3C4'), 
        panel.grid = element_blank(), 
        legend.position = 'right',
        # legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        legend.key.width = unit(0.25,'cm'),
        strip.text = element_text(face='bold'),
        axis.text = element_blank(), 
        axis.ticks = element_blank()); p_P
p_PET <- lt_PET_season_sen %>% 
  ggplot(data=., aes(x,y,fill=b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray70',color='gray10')+
  geom_tile()+
  scale_fill_gradient2(expression(atop("Monthly PET",(mm~yr**-1))),
                       limits=c(-1,1),
                       high=RColorBrewer::brewer.pal(7,'PuOr')[1], 
                       mid=RColorBrewer::brewer.pal(7,'PuOr')[4], 
                       low=RColorBrewer::brewer.pal(7,'PuOr')[7],
                       oob=scales::squish, 
                       na.value = 'gray')+
  labs(x=NULL,y=NULL)+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  facet_wrap(~season, nrow = 1)+
  # guides(fill=guide_colorbar(barwidth = 0.5))+
  theme(panel.background = element_rect(fill = '#99A3C4'), 
        panel.grid = element_blank(), 
        legend.position = 'right',
        # legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        legend.key.width = unit(0.25,'cm'),
        strip.text = element_text(face='bold'),
        axis.text = element_blank(), 
        axis.ticks = element_blank()); p_PET

library(cowplot)
p_out <- cowplot::plot_grid(p_P,p_PET, ncol = 1,labels='auto')
p_out
ggsave(p_out, 
       filename = "figures/SM_Fig2_Map_EastOz_lt_TheilSen_seasonal_P_PET.png", 
       type='cairo', dpi=300,
       width=11.75, height=13, units='cm')

#*******************************************************************************
# SM Fig 3: Long-term NDVI by season Thiel-Sen regression --------------------------------------------------------------
#*******************************************************************************
vec_col <- RColorBrewer::brewer.pal(n=7, name='BrBG')
quantile(na.omit(sen_ndvi_season$b1),c(0.01,0.99))
sen_ndvi_season %>% 
  ggplot(data=., aes(x,y,fill=b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray60',color='black')+
  geom_tile()+
  scale_fill_gradient2(expression(paste(Delta*NDVI~yr^-1)),
                       high=vec_col[7], mid=vec_col[4], low=vec_col[1],
                       limits=c(-0.005,0.005),
                       breaks=c(-0.005,-0.0025,0,0.0025,0.005),
                       labels=c("<-0.0005",-0.0025,0,0.0025,">0.005"),
                       oob=scales::squish,
                       na.value='gray')+
  labs(x=NULL,y=NULL)+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  facet_wrap(~season, nrow = 1)+
  guides(fill = 
           guide_colourbar(label = T)) +
  theme(panel.background = element_rect(fill = '#99A3C4'), 
        panel.grid = element_blank(), 
        legend.position = 'bottom',
        legend.key.width = unit(1.5,'cm'),
        strip.text = element_text(face='bold'),
        axis.text = element_blank(), 
        axis.ticks = element_blank())
ggsave(filename = 
         paste0("figures/SM_fig3_ndvi_f_of_timeANDsensor_seasonal_longterm_ThielSen_Trend_1982_2019_",Sys.Date(),".png"), 
       dpi=350, width=15,height=12,units='cm',type='cairo')
# END **************************************************************************


# #*******************************************************************************
# # SM Fig 4: NDVI Linear Model time series by Koppen climate zone ------------------------------------
# #*******************************************************************************
# library(mgcv)
# vec_ids <- unique(dat[,.(id,cz)]) %>% .[is.na(id)==F & is.na(cz)==F]
# vec_ids <- vec_ids[,.SD[sample(.N, min(10000,.N))],by=cz]
# o <- dat[is.na(season)==F] %>%
#   .[id %in% vec_ids$id] %>% 
#   .[date >= ymd('1981-09-01')] %>% 
#   .[date <= ymd('2019-08-30')] %>% 
#   .[,.(x,y,date,season,cz,id,ndvi_hyb,ndvi_3mo,ndvi_mcd)]
# vec_cols <- viridis::viridis(10, begin = 0.1,end=0.9)
# factor(o$season[1], levels=c("SON","DJF","MAM","JJA"),ordered = T)
# lut_kop <- c("Equatorial" = "Equat.",
#              "Tropical" = "Trop.", 
#              "Subtropical" = "Subtr.", 
#              "Grassland" = "Grass.",
#              "Arid" = "Arid",
#              "Temperate" = "Temp.",
#              "Temperate Tas." = "Tasm.")
# p_rlm <- o %>%   
#   # mutate(season=factor(o$season, levels=c("SON","DJF","MAM","JJA"),ordered = T)) %>% 
#   ggplot(data=., aes(date, ndvi_hyb))+
#   geom_smooth(method=MASS::rlm, color='black',se=F)+
#   geom_smooth(method=MASS::rlm,color=vec_cols[2],se=F,
#               data=o[sample(.N,1e6)][date %between% c("1981-01-01","1991-01-01")])+
#   geom_smooth(method=MASS::rlm,color=vec_cols[3],se=F,
#               data=o[sample(.N,1e6)][date %between% c("1986-01-01","1996-01-01")])+
#   geom_smooth(method=MASS::rlm,color=vec_cols[4],se=F,
#               data=o[sample(.N,1e6)][date %between% c("1991-01-01","2001-01-01")])+
#   geom_smooth(method=MASS::rlm,color=vec_cols[5],se=F,
#               data=o[sample(.N,1e6)][date %between% c("1996-01-01","2006-01-01")])+
#   geom_smooth(method=MASS::rlm,color=vec_cols[6],se=F,
#               data=o[sample(.N,1e6)][date %between% c("2001-01-01","2011-01-01")])+
#   geom_smooth(method=MASS::rlm,color=vec_cols[7],se=F,
#               data=o[sample(.N,1e6)][date %between% c("2006-01-01","2016-01-01")])+
#   geom_smooth(method=MASS::rlm,color=vec_cols[8],se=F,
#               data=o[sample(.N,1e6)][date %between% c("2011-01-01","2019-09-01")])+
#   geom_smooth(method=MASS::rlm,color='red',se=F,
#               data=o[date %between% c("2001-01-01","2019-08-30")],
#               aes(date, ndvi_mcd))+
#   geom_smooth(method=MASS::rlm,color=scales::muted('red'),se=F,
#               data=o[date %between% c("1982-01-01","2000-12-01")],
#               aes(date, ndvi_hyb))+
#   scale_x_date(expand=c(0,0))+
#   labs(x=NULL, y="NDVI")+
#   facet_grid(cz~season, scales = 'free_y', 
#              # labeller = label_wrap_gen(width=10, multi_line = TRUE)
#              labeller = labeller(cz = lut_kop)
#   )+
#   theme_linedraw()+
#   theme(panel.grid = element_blank(), 
#         strip.text = element_text(face='bold'))
# ggsave(p_rlm, filename = "figures/SM_fig4_ndvi_lin_rlm_trend_10yr_segs_by_Koppen.png", 
#        width=20, height=15, units='cm', dpi=350, type='cairo')
# # END **************************************************************************


# ******************************************************************************
# SM Fig 5: NDVI~CO2+PPETanom Linear Model by discretized P:PET range ---------------------------------------------------------------
# ******************************************************************************
library(broom); library(broom.mixed)

rlm_rm_na <- function(...){
  out <- MASS::rlm(...)
}

dat[is.na(season)==FALSE][mape<2][ndvi_anom_sd>-3.5&ndvi_anom_sd<3.5] %>% 
  as_tibble() %>% 
  mutate(pe_anom_12mo = pe_12mo - mape) %>% 
  mutate(epoch = ifelse(date < ymd("2000-12-31"),'avhrr','modis')) %>% 
  mutate(mape_d = cut_width(mape, 0.125)) %>% 
  group_by(mape_d) %>% 
  mutate(mape_dc= mean(mape,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(season = factor(season, levels=c("SON","DJF","MAM","JJA"),ordered = T)) %>% 
  nest(data = c(-mape_d,-mape_dc,-season)) %>% View
  mutate(fit = map(data, 
                   ~MASS::rlm(ndvi_3mo~scale(co2_trend)+scale(I(pe_anom_12mo/mape))+vc+epoch, 
                       # ~lm(ndvi_3mo~scale(co2_trend)+scale(I(pe_anom_12mo/mape)), 
                       data=.x)), 
         tidied = map(fit, tidy)) %>% 
  unnest(tidied) %>% 
  filter(str_detect(term,"vc")==F) %>%
  filter(str_detect(term,"Intercept")==F) %>%
  filter(str_detect(term,"epoch")==F) %>%
  # filter(str_detect(term, 'co2') | 
  #          str_detect(term, 'tmax_anom_3mo')|
  #          str_detect(term, 'pe_anom_12mo')) %>% 
  ggplot(data=., aes(mape_dc, estimate))+
  geom_col()+
  geom_errorbar(aes(ymin=estimate-2*std.error, 
                    ymax=estimate+2*std.error),
                width=0.05)+
  # scale_fill_viridis_c(option='B',direction=1,begin = 0.1,end=0.9, 
  #                       limits=c(0,0.05),oob=scales::squish)+
  # scale_color_viridis_c(option='B',direction=1,begin = 0.1,end=0.9, 
  #                       limits=c(0,0.05),oob=scales::squish)+
  # scale_x_discrete(guide = guide_axis(n.dodge = 1,angle = 90))+
  scale_x_continuous(limits=c(0,2.1), 
                     # breaks=c(0.25,0.5,0.75,1,1.25,1.5,1.75),
                     expand=c(0,0.05))+
  labs(x='Mean Annual P:PET Range', 
       title='NDVI Linear Model Effects (1982-2019)', 
       subtitle = 'N = 14.3e6   NDVI~CO2+P:PET_anom+Veg.Class.+Sensor')+
  facet_grid(term~season, scales = 'free_y', 
             labeller = labeller(
               term=c("scale(co2_trend)"="CO2",
                      "scale(co2_trend):scale(I(pe_anom_12mo/mape))"="CO2 x P:PET_anom",
                      "scale(I(pe_anom_12mo/mape))"="P:PET anom"))
  )+
  theme_linedraw()
ggsave(filename = 'figures/SM_fig5_big_linearModel_by_MAPPET_range.png', 
       height=12, width=16, units='cm')
# END **************************************************************************


# ******************************************************************************
# SM Fig 6: Logistic Asymptotic model ----------------------------------------------
# ******************************************************************************
library(nls.multstart)
# fn_xmid <- function(dat){
#   fit <- nls_multstart(ndvi_3mo ~ Asym/(1+exp((xmid-mape)/scal))+
#                          beta1*(precip_12mo-map)/(map)+
#                          beta2*(pe_12mo-mape)/(mape)+
#                          beta3*(co2_trend),
#                        data = dat,
#                        iter = 10,
#                        start_lower = c(Asym=0.5, xmid=0, scal=0, beta1=-0.1 ,beta2=-0.1,beta3=-0.1),
#                        start_upper = c(Asym=1, xmid=0.5, scal=1, beta1=0.1 ,beta2=0.1,beta3=0.1),
#                        supp_errors = 'Y',
#                        na.action = na.omit,
#                        lower= c(0.65,  0,   0, -0.1, -0.2,-0.1),
#                        upper= c(1,    0.5, 2,  0.1,  0.2,0.1))
#   bfit <- broom::tidy(fit)
#   bfit$hydro_year <- unique(dat$hydro_year);
#   bfit$season <- unique(dat$season)
#   return(bfit)
# }
fn_xmid <- function(dat){
  fit <- nls_multstart(ndvi_12mo ~ Asym/(1+exp((xmid-pe_12mo)/scal)),
                       data = dat,
                       iter = 10,
                       start_lower = c(Asym=0.5, xmid=0, scal=0),
                       start_upper = c(Asym=1, xmid=0.5, scal=1),
                       supp_errors = 'Y',
                       na.action = na.omit,
                       lower= c(0.65,  0,   0),
                       upper= c(1,    0.5, 2))
  bfit <- broom::tidy(fit)
  bfit$hydro_year <- unique(dat$hydro_year);
  # bfit$season <- unique(dat$season)
  return(bfit)
}
xdat <- dat[is.na(season)==FALSE][mape<1.5][ndvi_anom_sd>-3.5&ndvi_anom_sd<3.5] %>% 
    .[hydro_year %in% c(1982:2019)] %>% 
    .[date <= ymd("2019-08-30")] %>% 
    # .[sample(.N,1e6)] %>%
    .[,.(ndvi_12mo = mean(ndvi_hyb,na.rm=TRUE), 
            pe_12mo = mean(pe_12mo, na.rm=TRUE), 
            co2_12mo = mean(co2_trend,na.rm=TRUE)), 
      by=.(x,y,hydro_year)] %>% 
    .[,fn_xmid(.SD),by=.(hydro_year)]
gc()

mlo_annual <- dat[is.na(season)==FALSE][mape<1.5][ndvi_anom_sd>-3.5&ndvi_anom_sd<3.5] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[sample(.N,10000)] %>% 
  .[date <= ymd("2019-08-30")] %>% 
  # .[sample(.N,1e6)] %>%
  .[,.(co2_12mo = mean(co2_trend,na.rm=TRUE)), 
    by=.(hydro_year)] 
  
fn_logis <- function(x,Asym,xmid,scal) {
  # input <- seq(0.01,2,length.out = 100)
  Asym/(1+exp((xmid-x)/scal))}
wdat <- xdat %>%
  as_tibble() %>%
  # filter(p.value < 0.1) %>%
  select(#season,
         hydro_year,term,estimate) %>%
  spread(key = term, value=estimate) %>% 
  inner_join(., mlo_annual, by=c("hydro_year"))

p_asym <- wdat %>% 
  ggplot(data=.,aes(hydro_year, Asym,color=xmid))+
  geom_point()+
  geom_smooth(method='lm',se=F,color='black')+
  scale_color_viridis_c(expression(italic(m)), 
                        option='C',end=0.85)+
  labs(x=NULL,y='NDVI Asymptote')+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size=10), 
        legend.position = c(0.99,0.01),
        legend.justification = c(1,0),
        legend.key.width = unit(1,'cm'),
        legend.direction = 'horizontal'); p_asym

p_logistic <- expand_grid(wdat, 
                          x=seq(0.01,2,length.out = 100)) %>% 
  rowwise() %>% 
  mutate(p1 = fn_logis(x,Asym,xmid,scal)) %>% 
  mutate(pred = p1) %>% 
  # mutate(pred = fn_logis(x,Asym = max(c(Asym,Asym2,Asym3)),
  #                        xmid = (xmid+xmid2+xmid3)/3 ,
  #                        scal = (scal+scal2+scal3)/3 )) %>% 
  ggplot(data=., aes(x,pred,color=(co2_12mo), group=co2_12mo))+
  geom_line(alpha=0.5)+
  # geom_vline(aes(xintercept=p50, color=hydro_year), 
  #            data=wdat %>% 
  #        mutate(p50 = (log((Asym - R0)/Asym) + 0.693147180559945)*exp(-lrc)))+
  scale_color_viridis_c(expression(paste(CO[2]~(ppm))), option='B',end=0.85)+
  scale_x_continuous(limits=c(0,2), expand=c(0,0))+
  scale_y_continuous(limits=c(0,1), expand=c(0,0))+
  labs(x=expression(paste(paste(sum(Precip[t], "1 mo", "12 mo")," / ", sum(PET[t], "1 mo", "12 mo")))), 
       y=expression(paste(NDVI["annual"])))+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size=10), 
        legend.position = c(0.7,0.1), 
        legend.key.width = unit(1,'cm'),
        legend.direction = 'horizontal'); p_logistic

p_out <- p_logistic+p_asym+patchwork::plot_layout(ncol=2)
p_out
ggsave(p_out, filename = "figures/SM_fig6_ndvi3mo_PE12mo_logistic_Asym_xmid.png", 
       dpi=300, type='cairo', 
       width = 30, height=20, units='cm')
# END # ******************************************************************************

#*******************************************************************************
# SM Fig: 7 Richards Function NDVI~CO2+P:PET ---------------------------------------------
#*******************************************************************************
center_co2 <- mlo[date>=ymd("1981-12-01")&date<=ymd("2019-08-30")]$co2_trend %>% mean
dat[,`:=`(cco2=co2_trend-center_co2)]
dat[,`:=`(epoch = ifelse(date<=ymd("2000-12-31"),0,1))]
train_son <- dat[season=='SON'][mape <= 1.5][ndvi_anom_sd >= -3.5 & ndvi_anom_sd <=3.5] %>% 
  .[is.na(ndvi_3mo)==F] %>% 
  .[sample(.N,1000000)]
train_djf <- dat[season=='DJF'][mape <= 1.5][ndvi_anom_sd >= -3.5 & ndvi_anom_sd <=3.5] %>% 
  .[is.na(ndvi_3mo)==F] %>% 
  .[sample(.N,1000000)]
train_mam <- dat[season=='MAM'][mape <= 1.5][ndvi_anom_sd >= -3.5 & ndvi_anom_sd <=3.5] %>% 
  .[is.na(ndvi_3mo)==F] %>% 
  .[sample(.N,1000000)]
train_jja <- dat[season=='JJA'][mape <= 1.5][ndvi_anom_sd >= -3.5 & ndvi_anom_sd <=3.5] %>% 
  .[is.na(ndvi_3mo)==F] %>% 
  .[sample(.N,1000000)]

ric_x3_son <- nls_multstart(ndvi_3mo ~ 
              (Asym+Asym2*cco2+Asym3*I((pe_12mo-mape)/mape)) * 
              (1+exp(((xmid+xmid2*cco2+xmid3*I((pe_12mo-mape)/mape)) - mape)/
                       (scal+scal2*cco2+scal3*I((pe_12mo-mape)/mape))))^
              (-exp(-(lpow+lpow2*cco2+lpow3*I((pe_12mo-mape)/mape)))) + 
              B1*as.numeric(epoch), 
            data=train_son, 
            iter = 1, 
            supp_errors = 'N',
            control=nls.control(maxiter=100),
            start_lower = c(Asym=0.7561,Asym2=5e-4,Asym3=0,
                            xmid=0.27,xmid2=3e-4,xmid3=0,
                            scal=0.28,scal2=0,scal3=0,
                            lpow=-0.28,lpow2=0, lpow3=0, 
                            B3=0), 
            start_upper = c(Asym=0.7561,Asym2=5e-4,Asym3=0.0001,
                            xmid=0.27,xmid2=3e-4,xmid3=0.0001,
                            scal=0.28,scal2=0,scal3=0,
                            lpow=-0.28,lpow2=0, lpow3=0, 
                            B3=0))
summary(ric_x3_son)
yardstick::rsq_trad_vec(truth = train_son$ndvi_3mo, 
                        estimate = predict(ric_x3_son, newdata=train_son))
yardstick::rmse_vec(truth = train_son$ndvi_3mo, 
                    estimate = predict(ric_x3_son,newdata=train_son))


ric_x3_djf <- nls_multstart(ndvi_3mo ~ 
                              (Asym+Asym2*cco2+Asym3*I((pe_12mo-mape)/mape)) * 
                              (1+exp(((xmid+xmid2*cco2+xmid3*I((pe_12mo-mape)/mape)) - mape)/
                                       (scal+scal2*cco2+scal3*I((pe_12mo-mape)/mape))))^
                              (-exp(-(lpow+lpow2*cco2+lpow3*I((pe_12mo-mape)/mape)))) + 
                              B1*as.numeric(epoch), 
                            data=train_djf, 
                            iter = 1, 
                            supp_errors = 'N',
                            control=nls.control(maxiter=100),
                            start_lower = c(Asym=0.7561,Asym2=5e-4,Asym3=0,
                                            xmid=0.27,xmid2=3e-4,xmid3=0,
                                            scal=0.28,scal2=0,scal3=0,
                                            lpow=-0.28,lpow2=0, lpow3=0, 
                                            B3=0), 
                            start_upper = c(Asym=0.7561,Asym2=5e-4,Asym3=0.0001,
                                            xmid=0.27,xmid2=3e-4,xmid3=0.0001,
                                            scal=0.28,scal2=0,scal3=0,
                                            lpow=-0.28,lpow2=0, lpow3=0, 
                                            B3=0))
ric_x3_djf
yardstick::rsq_trad_vec(truth = train_djf$ndvi_3mo, 
                        estimate = predict(ric_x3_djf,newdata=train_djf))
yardstick::rmse_vec(truth = train_djf$ndvi_3mo,
                    estimate = predict(ric_x3_djf,newdata=train_djf))

ric_x3_mam <- nls_multstart(ndvi_3mo ~ 
                              (Asym+Asym2*cco2+Asym3*I((pe_12mo-mape)/mape)) * 
                              (1+exp(((xmid+xmid2*cco2+xmid3*I((pe_12mo-mape)/mape)) - mape)/
                                       (scal+scal2*cco2+scal3*I((pe_12mo-mape)/mape))))^
                              (-exp(-(lpow+lpow2*cco2+lpow3*I((pe_12mo-mape)/mape)))) + 
                              B1*as.numeric(epoch), 
                            data=train_mam, 
                            iter = 1, 
                            supp_errors = 'N',
                            control=nls.control(maxiter=100),
                            start_lower = c(Asym=0.7561,Asym2=5e-4,Asym3=0,
                                            xmid=0.27,xmid2=3e-4,xmid3=0,
                                            scal=0.28,scal2=0,scal3=0,
                                            lpow=-0.28,lpow2=0, lpow3=0, 
                                            B3=0), 
                            start_upper = c(Asym=0.7561,Asym2=5e-4,Asym3=0.0001,
                                            xmid=0.27,xmid2=3e-4,xmid3=0.0001,
                                            scal=0.28,scal2=0,scal3=0,
                                            lpow=-0.28,lpow2=0, lpow3=0, 
                                            B3=0))
ric_x3_mam
yardstick::rsq_trad_vec(truth = train_mam$ndvi_3mo, 
                        estimate = predict(ric_x3_mam, newdata = train_mam))
yardstick::rmse_vec(truth = train_mam$ndvi_3mo, 
                    estimate = predict(ric_x3_mam,newdata=train_mam))


ric_x3_jja <- nls_multstart(ndvi_3mo ~ 
                              (Asym+Asym2*cco2+Asym3*I((pe_12mo-mape)/mape)) * 
                              (1+exp(((xmid+xmid2*cco2+xmid3*I((pe_12mo-mape)/mape)) - mape)/
                                       (scal+scal2*cco2+scal3*I((pe_12mo-mape)/mape))))^
                              (-exp(-(lpow+lpow2*cco2+lpow3*I((pe_12mo-mape)/mape)))) + 
                              B1*as.numeric(epoch), 
                            data=train_jja, 
                            iter = 1, 
                            supp_errors = 'N',
                            control=nls.control(maxiter=100),
                            start_lower = c(Asym=0.7561,Asym2=5e-4,Asym3=0,
                                            xmid=0.27,xmid2=3e-4,xmid3=0,
                                            scal=0.28,scal2=0,scal3=0,
                                            lpow=-0.28,lpow2=0, lpow3=0, 
                                            B3=0), 
                            start_upper = c(Asym=0.7561,Asym2=5e-4,Asym3=0.0001,
                                            xmid=0.27,xmid2=3e-4,xmid3=0.0001,
                                            scal=0.28,scal2=0,scal3=0,
                                            lpow=-0.28,lpow2=0, lpow3=0, 
                                            B3=0))
ric_x3_jja
yardstick::rsq_trad_vec(truth = train_jja$ndvi_3mo, 
                        estimate = predict(ric_x3_jja, newdata = train_jja))
yardstick::rmse_vec(truth = train_jja$ndvi_3mo, 
                    estimate = predict(ric_x3_jja, newdata = train_jja))


pred_son <- expand_grid(#season=unique(tmp$season), 
  co2 = seq(min(dat$co2_trend,na.rm=T),
            max(dat$co2_trend,na.rm=T),length.out=50),
  epoch=0,
  pe_12mo = seq(0.01,1.5,length.out = 200)) %>% 
  mutate(cco2 = co2-center_co2) %>% 
  mutate(mape=pe_12mo) %>% 
  mutate(pred = predict(ric_x3_son, newdata=.), 
         season=train_son$season[1])
pred_djf <- expand_grid(#season=unique(tmp$season), 
  co2 = seq(min(dat$co2_trend,na.rm=T),
            max(dat$co2_trend,na.rm=T),length.out=50),
  epoch=0,
  pe_12mo = seq(0.01,1.5,length.out = 200)) %>% 
  mutate(cco2 = co2-center_co2) %>% 
  mutate(mape=pe_12mo) %>% 
  mutate(pred = predict(ric_x3_djf, newdata=.), 
         season = train_djf$season[1])
pred_mam <- expand_grid(#season=unique(tmp$season), 
  co2 = seq(min(dat$co2_trend,na.rm=T),
            max(dat$co2_trend,na.rm=T),length.out=50),
  epoch=0,
  pe_12mo = seq(0.01,1.5,length.out = 200)) %>% 
  mutate(cco2 = co2-center_co2) %>% 
  mutate(mape=pe_12mo) %>% 
  mutate(pred = predict(ric_x3_mam, newdata=.), 
         season = train_mam$season[1])
pred_jja <- expand_grid(#season=unique(tmp$season), 
  co2 = seq(min(dat$co2_trend,na.rm=T),
            max(dat$co2_trend,na.rm=T),length.out=50),
  epoch=0,
  pe_12mo = seq(0.01,1.5,length.out = 200)) %>% 
  mutate(cco2 = co2-center_co2) %>% 
  mutate(mape=pe_12mo) %>% 
  mutate(pred = predict(ric_x3_jja, newdata=.), 
         season=train_jja$season[1])

bind_rows(pred_son, pred_djf, pred_mam, pred_jja) %>% 
  mutate(season = factor(season,levels=c("SON","DJF","MAM","JJA"), 
                         ordered = TRUE)) %>% 
  ggplot(data=., aes(pe_12mo,pred,color=(co2), group=co2))+
  geom_line(alpha=1)+
  # geom_vline(aes(xintercept=p50, color=hydro_year), 
  #            data=wdat %>% 
  #        mutate(p50 = (log((Asym - R0)/Asym) + 0.693147180559945)*exp(-lrc)))+
  scale_color_viridis_c(expression(paste(CO[2]~ppm)), option='B',end=0.85)+
  scale_x_continuous(limits=c(0,1.55),
                     breaks=c(0,0.5,1,1.5),
                     labels = c(0,0.5,1,1.5),
                     expand=c(0,0)
                     # guide = guide_axis(n.dodge=1, angle=0,check.overlap = TRUE)
  )+
  scale_y_continuous(limits=c(0,0.9), expand=c(0,0))+
  labs(x="Mean Annual P:PET",
       # x=expression(paste(paste(sum(Precip[t], "1 mo", "12 mo")," / ", sum(PET[t], "1 mo", "12 mo")))), 
       y=expression(paste(NDVI)))+
  facet_wrap(~season,nrow=2)+
  theme_linedraw()+
  guides(color=guide_colorbar(title.position = 'top'))+
  theme(panel.grid.minor = element_blank(),
    # panel.spacing.x = unit(6, "mm"),
    axis.text = element_text(size=10),
    # axis.text.x = element_text(angle=45, vjust=-0.5),
    legend.position = c(0.85,0.1), 
    legend.key.width = unit(0.65,'cm'),
    legend.direction = 'horizontal', 
    legend.background = element_rect(fill=NA))
ggsave("figures/SM_fig7_ndvi3mo_PE12mo_richard_x3_nlsFit_bySeason.png", 
       width=15, height=12, units='cm', dpi=350, type='cairo')
# END SECTION *******************************************************************



# END # ******************************************************************************

#*******************************************************************************
# SM Fig 8: quantile GAM regression ---------------------------------------
#*******************************************************************************
library(qgam)
center_co2 <- mlo[date>=ymd("1981-11-01") & date<=ymd("2019-08-30")]$co2_trend %>% mean
dat <- dat[,`:=`(cco2=co2_trend-center_co2)]

train_dat <- dat[is.na(season)==FALSE][mape<1.5][ndvi_anom_sd>-3.5&ndvi_anom_sd<3.5] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[date <= ymd("2019-08-30")] %>% 
  # .[sample(.N,1e6)] %>%
  .[,.(ndvi_12mo = mean(ndvi_hyb,na.rm=TRUE), 
       pe_anom_12mo = mean(pe_12mo-mape, na.rm=TRUE),
       mape = mean(mape,na.rm=TRUE),
       cco2 = mean(cco2,na.rm=TRUE),
       co2_trend = mean(co2_trend,na.rm=TRUE)), 
    by=.(x,y,hydro_year)] %>% 
  .[,`:=`(epoch = ifelse(hydro_year<=2000,'avhrr','modis'))] %>% 
  .[,`:=`(epoch=factor(epoch))]

quSeq <- c(0.10, 0.5, 0.975)
fit <- mqgam(ndvi_12mo ~ 
               ti(cco2,k=3,bs='cs')+
               ti(mape,bs='cs',k=3) +
               ti(co2_trend,mape,k=3,bs='cs')+
               s(I(pe_anom_12mo/mape),k=3,bs='cs')+
               # s(I(pe_anom_12mo/mape),bs='cs',k=5) + 
               epoch, 
             data=train_dat[sample(.N, 100000)], 
             multicore=TRUE,
             ncores=10,
             qu = quSeq)
# qdo(fit,0.975,plot, scale=0,scheme=2)
# qdo(fit,quSeq[1],summary, scale=0)
# qdo(fit,quSeq[2],summary, scale=0)
# qdo(fit,quSeq[3],summary, scale=0)

# 0.74, 0.78, 0.71

zz <- expand_grid(season=unique(train_dat$season),
                  co2=c(340,380,420),
                  # co2 = seq(min(dat$co2_trend),max(dat$co2_trend),length.out=3),
                  mape = seq(0.05,1.5,length.out = 100), 
                  epoch=train_dat$epoch[1],
                  pct_anom = c(-66,0,33), 
                  iq = quSeq) %>% 
  mutate(co2_trend=co2) %>% 
  mutate(pe_anom_12mo = 0.01*pct_anom*mape) %>%
  mutate(pe_12mo = pe_anom_12mo+mape) %>%
  mutate(cco2 = co2-center_co2)

zz <- zz %>% 
  mutate(pred_high = qdo(fit, qu=quSeq[3], predict, newdata=.)) %>% 
  mutate(pred_med = qdo(fit, qu=quSeq[2], predict, newdata=.)) %>% 
  mutate(pred_low = qdo(fit, qu=quSeq[1], predict, newdata=.))

d_mean <- train_dat %>% summarize(value=mean(ndvi_12mo),
                                  mape=mean(mape))


p_q <- zz %>% 
  select(co2,mape,pct_anom,pred_high, pred_med, pred_low) %>% 
  gather(-co2,-mape,-pct_anom, key='key',value='value') %>% 
  mutate(key=factor(key,levels = c('pred_high','pred_med','pred_low'),ordered = T)) %>% 
  ggplot(data=., aes(mape, value,color=co2,group=co2))+
  # geom_hline(aes(yintercept=c(0.1)),color='red',lty=3)+
  # geom_hline(aes(yintercept=c(0.9)),color='red',lty=3)+
  # geom_vline(aes(xintercept=0.25),color='blue',lty=3)+
  # geom_vline(aes(xintercept=0.5),color='blue',lty=3)+
  # geom_vline(aes(xintercept=0.75),color='blue',lty=3)+
  # geom_vline(aes(xintercept=1),color='blue',lty=3)+
  geom_line()+
  # scale_color_discrete_sequential('ag_GrnYl',rev = TRUE)+
  scale_color_viridis_c(expression(paste(CO[2]~ppm)),
                        end=0.8,option='B',direction = 1)+
  scale_y_continuous(limits=c(0,1),expand=c(0,0.01))+
  scale_x_continuous(limits=c(0,1.5),
                     expand=c(0,0.1), 
                     breaks=seq(0,1.5,by=0.25),
                     minor_breaks = NULL)+
  labs(x="Mean Annual P:PET", 
       y=expression(NDVI[annual]))+
  facet_grid(key~pct_anom,
             labeller = labeller(key=c("pred_high"="Quantile 97.5%", 
                                       "pred_med"="Quantile 50%",
                                       "pred_low"="Quantile 10%"), 
                                 pct_anom=c("-66"="-66% P:PET anom.", 
                                            "0"="0% P:PET anom.",
                                            "33"="33% P:PET anom.")))+
  theme_linedraw()+
  theme(panel.background = element_blank(), 
        panel.grid.minor.y = element_blank()
        #panel.grid = element_blank()
  ); p_q
ggsave(p_q, filename = "figures/SM_fig8_quantile_reg_ndvi_ppet_x_co2.png", 
       width=20, height=14, units='cm', dpi=350, type='cairo')

# END # ******************************************************************************

#*******************************************************************************
# SM Fig: 9 Vegetation Index comparison of Weibull Function NDVI, EVI2, NIRv ~ CO2+P:PET 
#*******************************************************************************
# END # ******************************************************************************


#*******************************************************************************
# SM Fig:  Zonal VPD Trends  ---------------------------------------------
#*******************************************************************************
aa <- dat %>%
  lazy_dt() %>% 
  filter(hydro_year %in% 1982:2019) %>% 
  group_by(x,y,hydro_year) %>% 
  summarize(
    vpd = mean(vpd15_12mo, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(is.na(vpd)==F) %>% 
  as_tibble()
aa <- aa %>% 
  inner_join(., kop %>% select(x,y,zone), by=c("x","y"))
aa_mavpd <- aa %>% filter(hydro_year %in% c(1982:2010)) %>% 
  group_by(x,y) %>% 
  summarize(mavpd = mean(vpd,na.rm=TRUE)) %>% 
  ungroup()
p_vpd <- aa %>% 
  inner_join(aa_mavpd,by=c("x","y")) %>%
  mutate(vpd_pct = 100*vpd/mavpd - 100) %>% 
  sample_frac(0.5) %>% 
  rename(`Climate Zone` = zone) %>% 
  mutate(epoch = ifelse(hydro_year<=2000,'AVHRR 1982-2000','MODIS 2001-2019')) %>% 
  ggplot(data=., aes(hydro_year, vpd_pct,
                     color=`Climate Zone`,
                     group=paste(epoch,`Climate Zone`)))+
  # geom_point(alpha=0.05,color='gray')+
  geom_smooth(method=MASS::rlm)+
  # geom_smooth(inherit.aes = F, 
  #             aes(hydro_year, vpd_pct,
  #                    color=`Climate Zone`), 
  #             method=MASS::rlm, lty=3)+
  scale_x_continuous(expand=c(0,0), breaks = c(1982,1990,2000,2010,2019))+
  scale_y_continuous(expand=c(0,0), labels = scales::format_format(3))+
  scale_color_viridis_d(option='B',end=0.95)+
  # facet_wrap(~`Climate Zone`,scales = 'free',labeller = label_value, 
  #            ncol = 2)+
  labs(x=NULL, y="% Change of Annual VPD")+
  theme_linedraw()+
  theme(strip.text = element_text(face='bold'), 
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=7), 
        legend.position = c(0.01,0.99), 
        legend.justification = c(0.01,0.99)); p_vpd
ggsave(p_vpd, filename = "figures/SM_fig_zonal_vpd_trend_by_epoch.png", 
       width=12, height=10, units='cm', dpi=350, type='cairo')
# END # ******************************************************************************


#*******************************************************************************
# SM Fig:  Zonal PRECIP Trends  ---------------------------------------------
#*******************************************************************************
aa <- dat %>%
  lazy_dt() %>% 
  filter(hydro_year %in% 1982:2019) %>% 
  group_by(x,y,hydro_year) %>% 
  summarize(
    precip = mean(precip_12mo, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(is.na(precip)==F) %>% 
  as_tibble()
aa <- aa %>% 
  inner_join(., kop %>% select(x,y,zone), by=c("x","y"))
aa_map <- aa %>% filter(hydro_year %in% c(1982:2010)) %>% 
  group_by(x,y) %>% 
  summarize(map = mean(precip,na.rm=TRUE)) %>% 
  ungroup()
p_precip <- aa %>% 
  inner_join(aa_map,by=c("x","y")) %>%
  mutate(p_pct = 100*precip/map - 100) %>% 
  sample_frac(0.5) %>% 
  rename(`Climate Zone` = zone) %>% 
  mutate(epoch = ifelse(hydro_year<=2000,'AVHRR 1982-2000','MODIS 2001-2019')) %>% 
  ggplot(data=., aes(hydro_year, p_pct,
                     color=`Climate Zone`,
                     group=paste(epoch,`Climate Zone`)))+
  # geom_point(alpha=0.05,color='gray')+
  geom_smooth(method=MASS::rlm)+
  # geom_smooth(inherit.aes = F, 
  #             aes(hydro_year, vpd_pct,
  #                    color=`Climate Zone`), 
  #             method=MASS::rlm, lty=3)+
  scale_x_continuous(expand=c(0,0), breaks = c(1982,1990,2000,2010,2019))+
  scale_y_continuous(expand=c(0,0), labels = scales::format_format(3))+
  scale_color_viridis_d(option='B',end=0.95)+
  # facet_wrap(~`Climate Zone`,scales = 'free',labeller = label_value, 
  #            ncol = 2)+
  labs(x=NULL, y="% Change of Annual Precip.")+
  theme_linedraw()+
  theme(strip.text = element_text(face='bold'), 
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=7), 
        legend.position = c(0.01,0.99), 
        legend.justification = c(0.01,0.99)); p_precip
ggsave(p_precip, filename = "figures/SM_fig_zonal_precip_trend_by_epoch.png", 
       width=12, height=10, units='cm', dpi=350, type='cairo')
# END # ******************************************************************************

#*******************************************************************************
# SM Fig:  Zonal PET Trends  ---------------------------------------------
#*******************************************************************************
aa <- dat %>%
  lazy_dt() %>% 
  filter(hydro_year %in% 1982:2019) %>% 
  group_by(x,y,hydro_year) %>% 
  summarize(
    pet = mean(pet_12mo, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(is.na(pet)==F) %>% 
  as_tibble()
aa <- aa %>% 
  inner_join(., kop %>% select(x,y,zone), by=c("x","y"))
aa_mapet <- aa %>% filter(hydro_year %in% c(1982:2010)) %>% 
  group_by(x,y) %>% 
  summarize(mapet = mean(pet,na.rm=TRUE)) %>% 
  ungroup()
p_pet <- aa %>% 
  inner_join(aa_mapet,by=c("x","y")) %>%
  mutate(pet_pct = 100*pet/mapet - 100) %>% 
  # sample_frac(0.5) %>% 
  rename(`Climate Zone` = zone) %>% 
  mutate(epoch = ifelse(hydro_year<=2000,'AVHRR 1982-2000','MODIS 2001-2019')) %>% 
  ggplot(data=., aes(hydro_year, pet_pct,
                     color=`Climate Zone`,
                     group=paste(epoch,`Climate Zone`)))+
  geom_smooth(method=MASS::rlm)+
  scale_x_continuous(expand=c(0,0), breaks = c(1982,1990,2000,2010,2019))+
  scale_y_continuous(expand=c(0,0), labels = scales::format_format(3))+
  scale_color_viridis_d(option='B',end=0.95)+
  labs(x=NULL, y="% Change of Annual PET")+
  theme_linedraw()+
  theme(strip.text = element_text(face='bold'), 
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=7), 
        legend.position = c(0.01,0.99), 
        legend.justification = c(0.01,0.99)); p_pet
ggsave(p_pet, filename = "figures/SM_fig_zonal_pet_trend_by_epoch.png", 
       width=12, height=10, units='cm', dpi=350, type='cairo')
# END # ******************************************************************************

# SM Fig Merge VPD P PET relative change figures -----------------------------------
p_epoch_joint <- p_vpd/p_precip/p_pet+plot_layout(guides='collect')
p_epoch_joint
ggsave(p_epoch_joint, 
       filename = "figures/SM_fig_zonal_vpd_precip_pet_trend_by_epoch.png", 
       width=20, height=20, units='cm', dpi=350, type='cairo')
# END # **********************************************************************

#*******************************************************************************
# SM Fig:  Zonal P:PET Trends  ---------------------------------------------
#*******************************************************************************
aa <- dat %>%
  lazy_dt() %>% 
  filter(hydro_year %in% 1982:2019) %>% 
  group_by(x,y,hydro_year) %>% 
  summarize(
    ppet = mean(pe_12mo, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(is.na(ppet)==F) %>% 
  as_tibble()
aa <- aa %>% 
  inner_join(., kop %>% select(x,y,cz), by=c("x","y"))
aa_mappet <- aa %>% filter(hydro_year %in% c(1982:2010)) %>% 
  group_by(x,y) %>% 
  summarize(mappet = mean(ppet,na.rm=TRUE)) %>% 
  ungroup()
p_ppet <- aa %>% 
  inner_join(aa_mappet,by=c("x","y")) %>%
  mutate(ppet_pct = 100*ppet/mappet - 100) %>% 
  sample_frac(0.5) %>% 
  rename(`Climate Zone` = cz) %>% 
  mutate(epoch = ifelse(hydro_year<=2000,'AVHRR 1982-2000','MODIS 2001-2019')) %>% 
  ggplot(data=., aes(hydro_year, ppet,
                     color=`Climate Zone`,
                     group=paste(epoch,`Climate Zone`)))+
  # geom_point(alpha=0.05,color='gray')+
  geom_smooth(method=MASS::rlm)+
  # geom_smooth(inherit.aes = F, 
  #             aes(hydro_year, vpd_pct,
  #                    color=`Climate Zone`), 
  #             method=MASS::rlm, lty=3)+
  scale_x_continuous(expand=c(0,0), breaks = c(1982,1990,2000,2010,2019))+
  scale_y_continuous(expand=c(0,0), labels = scales::format_format(3))+
  scale_color_viridis_d(option='B')+
  # facet_wrap(~`Climate Zone`,scales = 'free',labeller = label_value, 
  #            ncol = 2)+
  labs(x=NULL, y="% Change of Annual P:PET")+
  facet_wrap(~`Climate Zone`,ncol = 1,scales = 'free')+
  theme_linedraw()+
  theme(strip.text = element_text(face='bold'), 
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=7), 
        legend.position = c(0.01,0.99), 
        legend.justification = c(0.01,0.99)); p_ppet
ggsave(p_ppet, filename = "figures/SM_fig_zonal_P:PET_trend_by_epoch.png", 
       width=12, height=10, units='cm', dpi=350, type='cairo')
# END # ******************************************************************************


# #*******************************************************************************
# # SM Fig:  Zonal P:PET Trends  ---------------------------------------------
# #*******************************************************************************
# aa <- dat %>%
#   lazy_dt() %>% 
#   filter(hydro_year %in% 1982:2019) %>% 
#   group_by(x,y,hydro_year) %>% 
#   summarize(
#     precip = mean(precip_12mo, na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   filter(is.na(precip)==F) %>% 
#   as_tibble()
# aa <- aa %>% 
#   inner_join(., kop %>% select(x,y,cz), by=c("x","y"))
# aa_maprecip <- aa %>% filter(hydro_year %in% c(1982:2010)) %>% 
#   group_by(x,y) %>% 
#   summarize(maprecip = mean(precip,na.rm=TRUE)) %>% 
#   ungroup()
# p_precip <- aa %>% 
#   inner_join(aa_maprecip,by=c("x","y")) %>%
#   mutate(precip_pct = 100*precip/maprecip - 100) %>% 
#   sample_frac(0.5) %>% 
#   rename(`Climate Zone` = cz) %>% 
#   mutate(epoch = ifelse(hydro_year<=2000,'AVHRR 1982-2000','MODIS 2001-2019')) %>% 
#   ggplot(data=., aes(hydro_year, precip,
#                      color=`Climate Zone`,
#                      group=paste(epoch,`Climate Zone`)))+
#   # geom_point(alpha=0.05,color='gray')+
#   geom_smooth(method=MASS::rlm)+
#   geom_smooth(method='lm',lty=3)+
#   geom_smooth()+
#   # geom_smooth(inherit.aes = F, 
#   #             aes(hydro_year, vpd_pct,
#   #                    color=`Climate Zone`), 
#   #             method=MASS::rlm, lty=3)+
#   scale_x_continuous(expand=c(0,0), breaks = c(1982,1990,2000,2010,2019))+
#   scale_y_continuous(expand=c(0,0), labels = scales::format_format(3))+
#   scale_color_viridis_d(option='B')+
#   # facet_wrap(~`Climate Zone`,scales = 'free',labeller = label_value, 
#   #            ncol = 2)+
#   labs(x=NULL, y="Annual Precip")+
#   facet_wrap(~`Climate Zone`,ncol = 1,scales = 'free')+
#   theme_linedraw()+
#   theme(strip.text = element_text(face='bold'), 
#         panel.grid = element_blank(), 
#         legend.position = 'none',# c(0.01,0.99), 
#         # legend.justification = c(0.01,0.99), 
#         axis.text.x = element_text(size=7)); p_precip
# ggsave(p_precip, filename = "figures/SM_fig_zonal_Precip_trend_by_epoch.png", 
#        width=12, height=10, units='cm', dpi=350, type='cairo')
# # END # ******************************************************************************




#*******************************************************************************
# SM Fig: 10 Burn Area Trends  ---------------------------------------------
#*******************************************************************************
library(sf)
oz_poly <- sf::read_sf("../data_general/GADM/gadm36_AUS.gpkg", 
                       layer="gadm36_AUS_1")
oz_poly <- st_as_sf(oz_poly)
oz_poly <- st_simplify(oz_poly, dTolerance = 0.05)

ba <- stars::read_stars("../data_general/MCD64/MCD64_Oz/MCD64A1_BurnArea_reducedToAVHRRres_2001_2019_.tif", 
                        proxy = F) %>% 
  stars::st_set_dimensions(., 3, values=seq(ymd("2001-01-01"),ymd("2019-12-01"),by='1 month')) %>% 
  as_tibble() %>% 
  purrr::set_names(c('x','y','date','ba'))
ba <- ba %>% inner_join(., kop, by=c("x","y"))
ba <- ba %>% filter(date <= ymd("2018-12-31"))

ba <- ba %>% mutate(year=year(date))
ba$year %>% mean
system.time(
  lt_ba <- ba %>% as.data.table() %>% .[,`:=`(year_c = year-2001)] %>%
    .[,`:=`(ba_km2 = ba/(1000**2))] %>% 
    .[,.(beta = list(lm(ba_km2~year_c, data=.SD)$coefficients)), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]  
)

p_mab <- lt_ba %>% 
  ggplot(data=.,aes(x,y,fill=b0+9*b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray40',color='black')+
  geom_tile()+
  coord_sf(xlim = c(140,155),
           ylim = c(-44,-10), 
           expand = FALSE)+
  # scale_fill_gradient2()+
  scale_x_continuous(breaks=seq(140,155,by=5))+
  scale_fill_viridis_c(expression(paste((km**2~yr**-1))),
                       limits=c(0,1),oob=scales::squish,
                       breaks=c(0,0.25,0.5,0.75,1),
                       labels=c(0,0.25,0.5,0.75,">1"),
                       option='B', begin = 0.1, trans='identity')+
  labs(x=NULL,y=NULL,
       title='Mean Burn Area: 2001-2018')+
  guides(fill = guide_colorbar(title.position = 'top'))+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(),
        legend.title = element_text(size=8),
        legend.key.width = unit(0.2,'cm'),
        legend.position = c(1,1), 
        legend.justification = c(1,1)); p_mab

# Absolute trend in BA
p_at <- lt_ba %>% 
  ggplot(data=.,aes(x,y,fill=b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray30',color='black')+
  geom_tile()+
  coord_sf(xlim = c(140,155),
           ylim = c(-44,-10), expand = FALSE)+
  scale_x_continuous(breaks=seq(140,155,by=5))+
  khroma::scale_fill_sunset(expression(paste("(",km**2~yr**-1,")")), 
                            limits=c(-0.1,0.1),oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title='Absolute Trend: 2001-2018')+
  guides(fill = guide_colorbar(title.position = 'top'))+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.key.width = unit(0.2,'cm'),
        legend.position = c(1,1), 
        legend.justification = c(1,1)); p_at

# Burn area trend by Koppen climate zone
d_tmp <- inner_join({ba %>% mutate(hydro_year=year(date+months(1))) %>% 
    group_by(hydro_year,zone) %>% 
    summarize(ba_km2 = sum(ba/1e6,na.rm=TRUE)) %>% 
    ungroup()}, 
    {ba %>% filter(date<=ymd("2005-12-31")) %>% group_by(zone) %>% 
        summarize(mab=sum(ba/(1e6),na.rm=TRUE)/5) %>% 
        ungroup()}, 
    by=c("zone")) %>% 
  mutate(val = ba_km2/mab)
p_ts <- ba %>% 
  # mutate(hydro_year=year(date+months(1))) %>% 
  mutate(year=year(date)) %>% 
  group_by(year,zone) %>% 
  summarize(ba_km2 = sum(ba/1e6,na.rm=TRUE)) %>% 
  ungroup() %>% 
  inner_join(.,   {ba %>% filter(date<=ymd("2005-12-31")) %>% group_by(zone) %>% 
      summarize(mab=sum(ba/(1e6),na.rm=TRUE)/5) %>% 
      ungroup()}, 
      by=c("zone")) %>% 
  mutate(val = ba_km2/mab) %>% 
  ggplot(data=.,aes(year, 100*val,color=zone))+
  geom_segment(aes(x=2001,xend=2005,y=50,yend=50), 
               arrow = arrow(ends = 'both',length=unit(0.1,'cm')), 
               inherit.aes = F)+
  geom_text(aes(x=2003,y=57),
            inherit.aes = F,
            size=4,
            label='reference period')+
  geom_smooth(method=MASS::rlm,se=F)+
  scale_color_viridis_d("",
                        option='B',end=0.9, direction=-1,
                        limits=rev(c("Equatorial","Tropical", 
                                     "Subtropical", "Grassland","Arid", 
                                     "Temperate","Temperate Tas.")), 
                        labels=str_wrap(rev(c("Equat.","Trop.", 
                                              "Subtrop.", "Grassland", "Arid", 
                                              "Temperate","Temp. Tasm.")),width = 10))+
  scale_x_continuous(expand=c(0.03,0.03), limits=c(2001,2018), 
                     breaks=c(2001,2005,2010,2015,2018))+
  labs(y=expression(paste('Burn Area Change (%)')), 
       x=NULL)+
  theme_linedraw()+
  theme(legend.position = c(0.01,0.01), 
        legend.justification = c(0.01,0.01), 
        panel.grid.minor = element_blank(),
        legend.direction = 'horizontal', 
        legend.background = element_blank()); p_ts


p_out <- p_mab|p_at+plot_layout(guides='keep')
ggsave(p_out/p_ts+plot_layout(heights=c(2,0.75))+plot_annotation(tag_levels = 'a'), 
       filename = 'figures/SM_fig10_mcd64_burnArea_mean_trend_v2.png', 
       width=16, height = 20, units='cm', dpi=350, type='cairo')
# END # ******************************************************************************

#*******************************************************************************
# SM Fig: 11 Map of MODIS Cover Change Trends ---------------------------------------------
#*******************************************************************************

lt_veg <- vcf_sen %>% rename(delta_treecover = tree_b, 
                   delta_nontree_veg = grass_b, 
                   delta_nonveg=nonveg_b,
                   nontree_veg_u = grass_u)

map_theme <- theme(panel.background = element_rect(fill = '#99A3C4'), 
                   panel.grid = element_blank(), 
                   legend.position = 'bottom', 
                   axis.text = element_blank(), 
                   axis.ticks = element_blank())

map_theme2 <- theme(panel.background = element_rect(fill = '#99A3C4'),
                    panel.grid = element_blank(),
                    legend.position = 'bottom',
                    axis.text = element_blank(),
                    axis.ticks = element_blank())

p_change <- lt_veg %>% 
  inner_join(., vc) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  mutate(change = 
           case_when(delta_treecover>=0.1 & delta_nontree_veg >= 0.1 ~ 'tree & grass incr.', 
                     delta_treecover<0.1 & delta_nontree_veg >= 0.1 ~ 'grass incr. & tree decr.',
                     delta_treecover>=0.1 & delta_nontree_veg < 0.1 ~ 'tree incr. & grass decr.', 
                     delta_nonveg > 0.1 ~ 'bare incr.', 
                     (between(delta_treecover,-0.1,0.1)==T) &
                       (between(delta_nontree_veg,-0.1,0.1)==T)~ "min change",
                     is.na(delta_treecover)==TRUE ~ NA_character_)) %>% #pull(change) %>% table
  filter(is.na(change)==F) %>% 
  ggplot(data=., aes(x,y,fill=change))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray40',color='black')+
  geom_tile()+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  # scale_fill_viridis_d()+
  # scale_fill_brewer(type='qual')+
  scale_fill_manual("", 
                    values = c('tree & grass incr.' = 'blue', 
                               'tree incr. & grass decr.' = '#228833', 
                               'grass incr. & tree decr.' = '#F0E442', 
                               'bare incr.' = '#D55E00', 
                               'min change' = '#AAAAAA'))+
  labs(x=NULL,y=NULL)+
  theme_linedraw()+
  guides(fill=guide_legend(ncol=1))+
  theme(panel.background = element_rect(fill = '#99A3C4'), 
        panel.grid = element_blank(), 
        # legend.position = c(0.7,0.88),
        legend.position = 'right',
        legend.direction = 'vertical',
        axis.text = element_blank(), 
        axis.ticks = element_blank()); p_change

# Delta treecover 
vec_col3 <- RColorBrewer::brewer.pal(n=7,'BrBG')
p_treecover <- lt_veg %>% 
  inner_join(., vc) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  filter(is.na(delta_treecover)==F) %>% 
  ggplot(data=., aes(x,y,fill=delta_treecover))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray40',color='black')+
  geom_tile()+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL,title='Tree Cover')+
  scale_fill_gradient2(expression(paste('%'~yr**-1)), 
                       limits=c(-1,1),oob=scales::squish, 
                       high=vec_col3[7], 
                       mid=vec_col3[4],
                       low=vec_col3[1])+
  theme_linedraw()+
  guides(fill=guide_colorbar(title.position = 'top'))+
  theme(panel.background = element_rect(fill = '#99A3C4'), 
        panel.grid = element_blank(),
        text = element_text(size=7),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key.width = unit(0.1,'cm'),
        legend.key.height = unit(0.4,'cm'),
        axis.text = element_blank(), 
        axis.ticks = element_blank()); p_treecover

# Delta grasscover
vec_col3 <- RColorBrewer::brewer.pal(n=7,'BrBG')
p_nontree_veg_cover <- lt_veg %>% 
  inner_join(., vc) %>% 
  inner_join(., kop,by=c("x","y")) %>% 
  filter(is.na(delta_nontree_veg)==F) %>% 
  ggplot(data=., aes(x,y,fill=delta_nontree_veg))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray40',color='black')+
  geom_tile()+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL,title='Non-tree Veg.')+
  scale_fill_gradient2(expression(paste("%"~yr**-1)), 
                       limits=c(-1,1),oob=scales::squish, 
                       high=vec_col3[7], 
                       mid=vec_col3[4],
                       low=vec_col3[1])+
  theme_linedraw()+
  guides(fill=guide_colorbar(title.position = 'top'))+
  theme(panel.background = element_rect(fill = '#99A3C4'), 
        panel.grid = element_blank(),
        text = element_text(size=7),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.title = element_text(),
        legend.key.width = unit(0.1,'cm'),
        legend.key.height = unit(0.4,'cm'),
        axis.text = element_blank(), 
        axis.ticks = element_blank()); p_nontree_veg_cover
p_join <- p_treecover+p_nontree_veg_cover+p_change+plot_layout(ncol=3)
p_join
ggsave(p_join, 
       filename = "figures/SM_fig11_EOz_map_delta_treecover_nontreeVegcover.png", 
       width=17, height = 12, units='cm', dpi=350, type='cairo')


# END # ******************************************************************************

#*******************************************************************************
# SM Fig: 12 NDVI GAM P:PET density plot ---------------------------------------------
#*******************************************************************************
library(viridisLite)
library(mgcv)
dat[pe_12mo <= 2][ndvi_mcd>0][sample(.N, 1e6)] %>% 
  ggplot(data=., aes(pe_12mo,ndvi_mcd))+
  # ggpointdensity::geom_pointdensity(size=0.5)+
  stat_density_2d(geom='raster',
                  aes(fill=after_stat(density)),
                  contour = F)+
  geom_smooth(se=F,color='gray80', method='gam',
              method.args=list(select=TRUE))+
  # scale_color_viridis_c()+
  scale_fill_gradientn(colors=c('black', viridis(99,end = 0.9)), 
                       trans='identity')+
  # khroma::scale_fill_land()+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(x=expression(paste(paste(sum(Precip[t], "1 mo", "12 mo")," / ", sum(PET[t], "1 mo", "12 mo")))), 
       y=expression(paste(NDVI)))+
  guides(fill=guide_colorbar(title.position = 'top', 
                             title = 'density'))+
  theme(legend.direction = 'horizontal',
        legend.position = c(0.75,0.15), 
        legend.background = element_rect(fill='black'), 
        legend.text = element_text(color='gray90'), 
        legend.title = element_text(color='gray90'), 
        legend.key.width = unit(1,'cm'))
ggsave(filename = 'figures/SM_fig12_mcd43_ndvi_ppet_2d_density.png',
       width=16, height = 10, units='cm',type='cairo')



# SM Fig 13: NDVI GAM P:VPD density plot ---------------------------------------------
library(viridisLite)
library(mgcv)
p_13 <- dat[mape <= 1.5][ndvi_mcd>0][sample(.N, 1e6)] %>% 
  ggplot(data=., aes(vpd15_12mo-mavpd15,ndvi_3mo))+
  # ggpointdensity::geom_pointdensity(size=0.5)+
  stat_density_2d(geom='raster',
                  aes(fill=after_stat(density)),
                  contour = F)+
  geom_smooth(se=F,color='white', method='gam',
              formula=y~s(x,bs='cs'),
              method.args=list(select=TRUE, 
                               method='REML'), 
              lwd=1.5)+
  geom_smooth(se=F,color='#cf0000', method='gam',
              formula=y~s(x,bs='cs'),
              method.args=list(select=TRUE, 
                               method='REML'), 
              lwd=1)+
  # scale_color_viridis_c()+
  scale_fill_gradientn(colors=c('black', inferno(99)), 
                       trans='identity')+
  # khroma::scale_fill_land()+
  scale_x_continuous(expand=c(0,0), limits=c(-0.5,0.5))+
  scale_y_continuous(expand=c(0,0))+
  labs(x=expression(paste(VPD~anom["12 mo"])),
       y=expression(paste(NDVI)))+
  guides(fill=guide_colorbar(title.position = 'top', 
                             title = 'density'))+
  theme(legend.direction = 'horizontal',
        legend.position = c(0.025,0),
        legend.justification = c(0,0),
        legend.background = element_rect(fill='black'), 
        legend.text = element_text(color='gray90'), 
        legend.title = element_text(color='gray90'), 
        legend.key.width = unit(0.5,'cm'))
ggsave(p_13, filename = 'figures/SM_fig_13_vpdanom12mo_ndvi_2d_density.png',
       width=16, height = 10, units='cm',type='cairo')
# END ********


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
       pet_anom = mean(pet_12mo-mapet,na.rm=TRUE),
       ppet_anom = mean(pe_12mo,na.rm=TRUE),
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
       pet_anom = mean(pet_12mo-mapet,na.rm=TRUE),
       ppet_anom = mean(pe_12mo,na.rm=TRUE),
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
       pet_anom = mean(pet_12mo-mapet,na.rm=TRUE),
       ppet_anom = mean(pe_12mo,na.rm=TRUE),
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
      scale(ppet_anom)+
      scale(p_anom)+
      scale(pet_anom)+
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


rlm_ndvi_annual_co2_ppet_e1 <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 1982:2000] %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_start = co2 - min_co2)] %>% 
  .[,.(beta = list(coef(MASS::rlm(
    ndvi_hyb~
      co2_start+
      scale(ppet_anom)+
      scale(p_anom)+
      scale(pet_anom)))))
    ,
    by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2], 
          b2=unlist(beta)[3], 
          b3=unlist(beta)[4],
          b4=unlist(beta)[5]
  ), by=.(x,y)]

rlm_ndvi_annual_co2_ppet_e2 <-  tmp_ndvi[id%in%vec_ids] %>% # filter out pixels that only have data for one satellite epoch
  as.data.table() %>% 
  .[hydro_year %in% 2001:2019] %>% 
  .[is.na(ndvi_hyb)==F & is.na(frac_p_anom)==F] %>% 
  .[,`:=`(co2_start = co2 - mid_co2)] %>% 
  .[,.(beta = list(coef(MASS::rlm(
    ndvi_hyb~
      co2_start+
      scale(ppet_anom)+
      scale(p_anom)+
      scale(pet_anom)))))
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
ggsave(filename = 'figures/SM_fig_14_rlm_CO2_effect_by_epoch.png',
       width=16, height = 12, units='cm',type='cairo')
# End section ******************************************************************


# SM Fig : Compare ERA5, AWAP, CRU PET ------------------------------------------
library(stars)
base <- stars::read_stars("../data_general/MCD43/MCD43A4_NDVI_5000m_EastOz_mmean_maskFireDefor_2001_2019.tif", 
                          RasterIO=list(bands=1))
cru <- stars::read_stars("../data_general/clim_grid/cru/cru_eoz_1901_2017.nc")
cru <- stars::st_warp(cru,base)
cru <- cru %>% st_set_dimensions(.,3,values=seq(ymd("1901-01-01"),ymd("2017-12-01"),by='1 month'),names = 'date')
dcru <- cru[,,,c(80*12):c(117*12)] %>% as_tibble()
dcru <- dcru %>% set_names(c("x","y","date","pet_cru")) %>% 
  filter(is.na(pet_cru)==F)
dcru <- dcru %>% units::drop_units()

e5pet <- stars::read_ncdf("../data_general/clim_grid/era5-land/Oz/Oz/Oz_era5-land_pet_1981_2019.nc")
names(e5pet) <- "pet_e5"
st_crs(e5pet) <- st_crs(4326)
eoz_box <- st_bbox(c(xmin = min(nvis$x),
                     ymin = min(nvis$y),
                     xmax = max(nvis$x),
                     ymax = max(nvis$y)), 
                   crs = st_crs(4326))
e5pet <- st_crop(e5pet, eoz_box)
e5pet <- e5pet %>% units::drop_units()
e5pet <- stars::st_warp(e5pet,base)
e5pet <- e5pet %>% as_tibble() %>% rename(date=time)

dawap <- dat %>% select(x,y,date,pet) %>% tb %>% rename(pet_awap = pet)
d <- inner_join(dcru, dawap, by=c("x","y","date"))
d <- inner_join(d, e5pet, by=c("x","y","date"))
d <- inner_join(d,kop %>% select(x,y,zone) %>% tb,by=c("x","y"))

# Correlation between AWAP PET and CRU PET
d %>% 
  filter(date >= ymd("1990-01-01")) %>% 
  filter(is.na(pet_cru)==F & is.na(pet_awap)==F) %>% 
  select(pet_cru,pet_awap) %>% 
  cor

# Correlation between AWAP PET and CRU PET
d %>% 
  filter(date < ymd("1990-01-01")) %>% 
  filter(is.na(pet_cru)==F & is.na(pet_awap)==F) %>% 
  select(pet_cru,pet_awap) %>% 
  cor

d %>% 
  mutate(period = ifelse(date<ymd("1990-01-01"),
                         'calib-ERA5 1981-1989',
                         'AWAP 1990-1999')) %>% 
  filter(date<ymd("2000-01-01")) %>% 
  sample_n(0.1e6) %>% 
  ggplot(data=.,aes(pet_cru,pet_awap))+
  geom_point(alpha=0.1,size=0.1)+
  geom_smooth(method='lm',color='blue',lwd=1)+
  geom_abline(color='#cf0000')+
  labs(x='CRU PET mm/month', 
       y='PET mm/month')+
  facet_grid(zone~period)+
  ggpubr::stat_regline_equation(label.y = 230, 
                                label.x = 0)+
  ggpubr::stat_cor(label.y = 180, 
                   label.x=0)+
  theme_linedraw()
ggsave(filename = "figures/compare_PET_awap_era5_cru.png", 
       width=16,height = 20,units='cm')

# d %>% filter(date>=ymd("1990-01-01")&date<=ymd("2010-12-31")) %>% 
#   select(date,pet_cru, pet_awap, pet_e5) %>% 
#   filter(is.na(pet_e5)==F) %>% 
#   select(-date) %>% 
#   cor

# End section ******************************************************************




