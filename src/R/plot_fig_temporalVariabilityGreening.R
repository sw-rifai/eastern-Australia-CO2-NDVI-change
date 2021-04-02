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

# mod <- svi[mod, on=.(x,y,year)]
# end section ******************************************************************
# END DATA IMPORT SECTION ******************************************************


#*******************************************************************************
# SM Fig 4: NDVI Linear Model time series by Koppen climate zone ------------------------------------
#*******************************************************************************
library(mgcv)
dat <- merge(dat,kop,by=c("x","y"))
vec_ids <- unique(dat[,.(id,cz)]) %>% .[is.na(id)==F & is.na(cz)==F]
vec_ids <- vec_ids[,.SD[sample(.N, min(10000,.N))],by=cz]
o <- dat[is.na(season)==F] %>%
  .[id %in% vec_ids$id] %>% 
  .[date >= ymd('1981-09-01')] %>% 
  .[date <= ymd('2019-08-30')] %>% 
  .[,.(x,y,date,season,cz,id,ndvi_hyb,ndvi_3mo,ndvi_mcd)]
vec_cols <- colorspace::sequential_hcl(n=12,palette = 'Light Grays',rev=TRUE)[2:10]
vec_reds <- colorspace::sequential_hcl(n=10,palette = 'Greens3',rev = TRUE)[c(5,8)]
# colorspace::hcl_palettes()
colorspace::demoplot(colorspace::sequential_hcl(n=10,palette = 'Greens 3'))
colorspace::demoplot(vec_reds)


factor(o$season[1], levels=c("SON","DJF","MAM","JJA"),ordered = T)
lut_kop <- c("Equatorial" = "Equat.",
             "Tropical" = "Trop.", 
             "Subtropical" = "Subtr.", 
             "Grassland" = "Grass.",
             "Arid" = "Arid",
             "Temperate" = "Temp.",
             "Temperate Tas." = "Tasm.")
p_rlm <- o %>%   
  # mutate(season=factor(o$season, levels=c("SON","DJF","MAM","JJA"),ordered = T)) %>% 
  ggplot(data=., aes(date, ndvi_hyb))+
  geom_smooth(method=MASS::rlm,color=vec_cols[2],se=F,
              data=o[sample(.N,1e6)][date %between% c("1981-01-01","1991-01-01")])+
  geom_smooth(method=MASS::rlm,color=vec_cols[3],se=F,
              data=o[sample(.N,1e6)][date %between% c("1986-01-01","1996-01-01")])+
  geom_smooth(method=MASS::rlm,color=vec_cols[4],se=F,
              data=o[sample(.N,1e6)][date %between% c("1991-01-01","2001-01-01")])+
  geom_smooth(method=MASS::rlm,color=vec_cols[5],se=F,
              data=o[sample(.N,1e6)][date %between% c("1996-01-01","2006-01-01")])+
  geom_smooth(method=MASS::rlm,color=vec_cols[6],se=F,
              data=o[sample(.N,1e6)][date %between% c("2001-01-01","2011-01-01")])+
  geom_smooth(method=MASS::rlm,color=vec_cols[7],se=F,
              data=o[sample(.N,1e6)][date %between% c("2006-01-01","2016-01-01")])+
  geom_smooth(method=MASS::rlm,color=vec_cols[8],se=F,
              data=o[sample(.N,1e6)][date %between% c("2011-01-01","2019-09-01")])+
  geom_smooth(method=MASS::rlm,color=vec_reds[2],se=F,
              data=o[date %between% c("2001-01-01","2019-08-30")],
              aes(date, ndvi_mcd))+
  geom_smooth(method=MASS::rlm,color=vec_reds[1],se=F,
              data=o[date %between% c("1982-01-01","2000-12-01")],
              aes(date, ndvi_hyb))+
  geom_smooth(method=MASS::rlm, color='black',se=F,lwd=0.75)+
  scale_x_date(expand=c(0,0))+
  labs(x=NULL, y="NDVI")+
  facet_grid(cz~season, scales = 'free_y', 
             # labeller = label_wrap_gen(width=10, multi_line = TRUE)
             labeller = labeller(cz = lut_kop)
  )+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.text = element_text(face='bold'))
ggsave(p_rlm, filename = "figures/fig3_ndvi_lin_rlm_trend_10yr_segs_by_Koppen.png", 
       width=20, height=15, units='cm', dpi=350, type='cairo')
# END **************************************************************************



# # 7 Koppen Climate Zones & P:PET Trend & NDVI & VCF distributions --------------------------------------------------
# lut_kop <- c("Equatorial" = "Equat.",
#              "Tropical" = "Trop.", 
#              "Subtropical" = "Subtr.", 
#              "Grassland" = "Grass.", 
#              "Arid" = "Arid",
#              "Temperate" = "Temp.",
#              "Temperate Tas." = "Tasm.")
# 
# p_kop <- kop %>% 
#   # mutate(zone = fct_recode(zone, zone="Desert",zone="Arid")) %>% 
#   filter(is.na(zone)==F) %>% 
#   ggplot(data=., aes(x,y,fill=zone))+
#   geom_sf(inherit.aes = F, data=oz_poly,
#           fill='gray40',color='gray10')+
#   geom_tile()+
#   coord_sf(xlim = c(140,155.5),
#            ylim = c(-45,-10), 
#            expand = FALSE)+
#   scale_x_continuous(breaks=seq(140,155,by=5))+
#   scale_fill_viridis_d(option='B',direction = 1,end=0.95, label = lut_kop,na.translate=F)+
#   # scico::scale_fill_scico_d(end=0.9,direction = 1)+
#   labs(x=NULL,y=NULL)+
#   theme_linedraw()+
#   guides(fill=guide_legend(title='Zone', 
#                            title.position = 'top'))+
#   theme(#legend.position = c(0.75,0.85),
#     legend.position = c(1,1), 
#     legend.justification = c(1.01,1.005),
#     legend.box.background = element_rect(fill="#FFFFFF11",color=NA,linetype = NULL),
#     legend.text = element_text(size=8),
#     legend.title = element_text(size=8.5),
#     legend.spacing = unit(10,'points'),
#     legend.key.width = unit(0.3,'cm'),
#     legend.direction = 'vertical',
#     panel.grid = element_blank(), 
#     panel.background = element_rect(fill='lightblue')); p_kop
# 
# ggsave(plot=p_kop,
#        filename = "figures/map_7KoppenZones.png", 
#        width = 7.25, height=15, units='cm', dpi=350, type='cairo')
# 
# p_kop$theme$legend.key.size
# 
# library(magick)
# p_left <- image_read("figures/fig3_ndvi_lin_rlm_trend_10yr_segs_by_Koppen.png")
# p_left <- magick::image_annotate(image = p_left,text = '(a)',location="+0+15",size = 100)
# 
# p_right <- magick::image_read("figures/map_7KoppenZones.png")
# p_right <- magick::image_annotate(image = p_right,text = '(b)',location="+0+35",size = 100)
# p_out <- image_append(c(p_left,p_right),stack = F)
# p_out
# magick::image_write(p_out, path="figures/Fig3_map-Koppen-ndvi_lin_rlm_trend_10yr_segs_by_Koppen.png")
# 
# 
# 
