# !diagnostics off
#*************************************************************************
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
options(pillar.sigfig = 5)
# IMPORT DATA ##############################################################

# load("data/gridCell_lm_ndvi_clim.Rdata") # grid cell linear regressions
oz_poly <- sf::read_sf("../data_general/GADM/gadm36_AUS.gpkg", 
                       layer="gadm36_AUS_1")
oz_poly <- st_as_sf(oz_poly)
oz_poly <- st_simplify(oz_poly, dTolerance = 0.05)

mcf <- stars::read_ncdf("../data_general/Oz_misc_data/csiro_FC.v310.MCD43A4_0p5_2001_2019.nc")
mcf[,,,1] %>% as_tibble() %>% filter(is.na(soil)==F)
eoz <- stars::read_stars("../data_general/Oz_misc_data/MOD44BPercent_Tree_Cover_5000m_East_Oz_noMask_2000_2019.tif",
                         RasterIO = list(bands=1))
mcf <- stars::st_warp(src=mcf, dest=eoz, use_gdal = F)
mcf <- mcf %>% as.data.table() %>% lazy_dt() %>% 
  rename(date=time) %>%
  mutate(date=as.Date(date)) %>% 
  filter(is.na(npv)==F) %>% 
  as.data.table()
# mcf[,`:=`(year=year(date),month=month(date))] %>% 
#   .[,`:=`(season = case_when(month%in%c(3:5)~'MAM',
#                              month%in%c(6:8)~'JJA',
#                              month%in%c(9:11)~'SON',
#                              month%in%c(12,1,2)~'DJF'))]
# mcf[,`:=`(season = factor(season, levels = c('SON','DJF','MAM','JJA'),ordered = TRUE))]



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

vi <- mcf[,.(x,y,date,soil,gv,npv)][vi,on=.(x,y,date)]

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
vc <- left_join(nvis, codes, by='veg_class') %>% 
  as.data.table()
mod <- vc[mod, on=.(x,y)]
mod[,`:=`(year=year(date))]



# Regressions -------------------------------------------------------------
library(zyp)
## NDVI trend using Theil-Sen ---------------------------------------------
sen_ndvi_season_e1 <- dat %>% 
  .[date >= ymd("1981-11-01") & date <= ymd("2000-12-31")] %>% 
  .[hydro_year %in% c(1982:2000)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE)), by=.(x,y,season,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y,season)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y,season)]
sen_ndvi_season_e2 <- dat %>% 
  .[date >= ymd("2001-01-01") & date <= ymd("2019-08-30")] %>% 
  .[hydro_year %in% c(2001:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-2001)] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE)), by=.(x,y,season,hydro_year_c)] %>% 
  .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y,season)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y,season)]
sen_ndvi_season_e1$epoch <- "AVHRR NDVI 1982-2000"
sen_ndvi_season_e2$epoch <- "MODIS NDVI 2001-2019"

# MODIS VCF trends using Theil Sen
system.time(
  lt_tree_sen <- inner_join(mod[year<=2018][,`:=`(year_c = year-2009.5)],
                            kop %>% select(x,y,zone), by=c("x","y")) %>%
    as.data.table() %>% 
    .[,.(beta = list(unname(zyp.sen(tree_cover~year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]  
)
system.time(
  lt_nontree_sen <- inner_join(mod[year<=2018][,`:=`(year_c = year-2009.5)],
                               kop %>% select(x,y,zone), by=c("x","y")) %>% 
    as.data.table() %>% 
    .[,.(beta = list(unname(zyp.sen(nontree_cover~year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]  
)
system.time(
  lt_nonveg_sen <- inner_join(mod[year<=2018][,`:=`(year_c = year-2009.5)],
                              kop %>% select(x,y,zone), by=c("x","y")) %>% 
    as.data.table() %>% 
    .[,.(beta = list(unname(zyp.sen(nonveg_cover~year_c, 
                                    data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2]), by=.(x,y)]  
)

vcf_sen <- inner_join(lt_nontree_sen %>% rename(grass_u=b0, grass_b=b1) %>% select(-beta), 
                      lt_nonveg_sen %>% rename(nonveg_u=b0, nonveg_b=b1) %>% select(-beta)) %>% as.data.table()

vcf_sen <- lt_tree_sen %>% lazy_dt() %>% 
  select(-beta) %>% 
  rename(tree_u=b0, 
         tree_b=b1) %>% 
  as.data.table() %>% 
  inner_join(., vcf_sen) %>% 
  as.data.table()


# 7 Koppen Climate Zones & P:PET Trend & NDVI & VCF distributions --------------------------------------------------
lut_kop <- c("Equatorial" = "Equat.",
             "Tropical" = "Trop.", 
             "Subtropical" = "Subtr.", 
             "Grassland" = "Grass.", 
             "Arid" = "Arid",
             "Temperate" = "Temp.",
             "Temperate Tas." = "Tasm.")

p_kop <- kop %>% 
  # mutate(zone = fct_recode(zone, zone="Desert",zone="Arid")) %>% 
  filter(is.na(zone)==F) %>% 
  as.data.table() %>% 
  ggplot(data=., aes(x,y,fill=zone))+
  geom_sf(inherit.aes = F, data=oz_poly,
          fill='gray40',color='gray10')+
  geom_tile()+
  coord_sf(xlim = c(140,155),
           ylim = c(-45,-10))+
  scale_x_continuous(breaks=seq(140,155,by=5))+
  scale_fill_viridis_d(option='B',direction = 1,end=0.95, label = lut_kop,na.translate=F)+
  # scico::scale_fill_scico_d(end=0.9,direction = 1)+
  labs(x=NULL,y=NULL)+
  theme_linedraw()+
  guides(fill=guide_legend(title='Climate Zone', 
                           title.position = 'top'))+
  theme(#legend.position = c(0.75,0.85),
    legend.position = c(1,1), 
    legend.justification = c(1.01,1.01),
    legend.text = element_text(size=17),
    legend.title = element_text(size=18),
    legend.direction = 'vertical',
    panel.grid = element_blank(), 
    panel.background = element_rect(fill='lightblue')); p_kop

# Zonal P:PET trend
aa <- ldat %>% 
  filter(hydro_year %in% 1982:2019) %>% 
  group_by(x,y,hydro_year) %>% 
  summarize(
    ppet = mean(pe_12mo, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(is.na(ppet)==F) %>% 
  as_tibble()
aa <- aa %>% 
  inner_join(., kop %>% select(x,y,cz) %>% as.data.table, 
    by=c("x","y")) %>% 
    as.data.table() 
aa_mappet <- aa %>% filter(hydro_year %in% c(1982:2010)) %>% 
  group_by(x,y) %>% 
  summarize(mappet = mean(ppet,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
p_ppet <- aa %>% 
  inner_join(aa_mappet,by=c("x","y")) %>% 
  mutate(ppet_pct = 100*ppet/mappet - 100) %>% 
  # sample_frac(0.5) %>% 
  filter(is.na(cz)==FALSE) %>% 
  rename(`Climate Zone` = cz) %>% 
  as.data.table() %>% 
  ggplot(data=., aes(hydro_year, ppet_pct,color=`Climate Zone`,group=`Climate Zone`))+
  # geom_point(alpha=0.05,color='gray')+
  geom_smooth(method=MASS::rlm)+
  scale_x_continuous(expand=c(0,0), breaks = c(1982,1990,2000,2010,2019), limits=c(1982,2019.75))+
  scale_y_continuous(expand=c(0,0), labels = scales::format_format(3))+
  scale_color_viridis_d(option='B',end=0.9)+
  # facet_wrap(~`Climate Zone`,scales = 'free',labeller = label_value, 
  #            ncol = 2)+
  labs(x=NULL, y="% Change of Annual P:PET")+
  theme_linedraw()+
  theme(strip.text = element_text(face='bold'), 
        panel.grid = element_blank(), 
        axis.text.x = element_text(size=10), 
        legend.position = 'none'); p_ppet

aa %>% 
  inner_join(aa_mappet,by=c("x","y")) %>% 
  mutate(ppet_pct = 100*ppet/mappet - 100) %>% 
  # sample_frac(0.5) %>% 
  filter(is.na(cz)==FALSE) %>% 
  # rename(`Climate Zone` = cz) %>% 
  as.data.table() %>% 
  split(.$cz) %>% 
  map(~MASS::rlm(ppet_pct~hydro_year, data=.x)) %>% 
  map(summary)


sen_ndvi_season_e1$epoch <- "AVHRR NDVI 1982-2000"
sen_ndvi_season_e2$epoch <- "MODIS NDVI 2001-2019"
j_sen <- bind_rows(sen_ndvi_season_e1, sen_ndvi_season_e2) %>% 
  as.data.table()
p_ndvi_sen <- inner_join(kop,j_sen,by=c("x","y")) %>% 
  filter(is.na(cz)==FALSE) %>% 
  rename(`Zone` = cz) %>% 
  filter(between(b1,-0.006,0.006)) %>% 
  mutate(season = factor(season, levels=c("SON","DJF","MAM","JJA"),ordered = T)) %>% 
  as.data.table() %>% 
  ggplot(data=., aes(b1,after_stat(density),
                     fill=as_factor(`Zone`),
                     color=epoch, #as_factor(`Zone`),
                     alpha=epoch))+
  # geom_histogram(bins = 30, position = 'identity')+
  geom_density(position='identity')+
  # geom_freqpoly()+
  geom_vline(aes(xintercept=0),col='#cf0000',lwd=1,alpha=0.5)+
  geom_vline(aes(xintercept=0),col='red',lwd=0.5,alpha=0.5)+
  # scico::scale_fill_scico_d()+
  scale_alpha_discrete(range=c(0.2,0.8))+
  scale_fill_viridis_d(option='B',direction = 1,end=0.925)+
  # scale_color_viridis_d(option='B',direction = 1,end=0.925)+
  scale_color_manual(values=c("black",NA),
                     breaks=c("AVHRR NDVI 1982-2000","MODIS NDVI 2001-2019"))+
  scale_x_continuous(expand=c(0,0), 
                     limits=c(-0.005,0.006),
                     breaks=c(-0.003,0,0.003))+
  scale_y_continuous(expand=c(0,0),position = 'left')+
  facet_grid(`Zone`~season, scales = 'free_y', 
             # labeller = label_wrap_gen(width=10, multi_line = TRUE)
             labeller = labeller(`Zone` = lut_kop)
  )+
  labs(x=expression(paste(Delta~NDVI~yr**-1)))+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.text = element_text(face='bold'), 
        legend.position = 'none', 
        axis.text = element_text(size=8)); p_ndvi_sen


library(ggridges)
p_vcf <- vcf_sen %>% 
  as_tibble() %>% 
  filter(is.na(tree_u)==F) %>% 
  inner_join(., kop %>% select(x,y,cz) %>% as_tibble(), by=c("x","y")) %>% 
  select(cz, nonveg_b, grass_b, tree_b) %>%
  rename(`Non. Veg.`=nonveg_b, 
         `Non. Tree Veg.`=grass_b, 
         `Tree Veg.`=tree_b) %>% 
  gather(-cz, key = 'measure', value='estimate') %>% 
  mutate(measure = as_factor(measure)) %>% 
  filter(is.na(estimate)==F) %>% 
  filter(is.na(cz)==F) %>% 
  ggplot(data=., aes(x=estimate,
                     y=cz,
                     fill=cz,
                     after_stat(scaled)))+
  ggridges::stat_density_ridges(
    quantile_lines = TRUE, 
    quantiles = 2,
    rel_min_height=0.02, 
    bandwidth = 0.1,
    alpha=0.5, 
    scale=1, 
    lwd=0.5,
    color='black')+
  scale_fill_viridis_d(option='B',direction = 1,end=0.9,na.translate=F)+
  geom_vline(aes(xintercept=0),
    color='red',lwd=1,alpha=0.5)+
  scale_x_continuous(limits=c(-1,1))+
  scale_y_discrete(expand=c(0,0),
                   limits=rev(c("Equatorial","Tropical",
                                "Subtropical", "Grassland","Arid",
                                "Temperate","Temperate Tas.")),
                   labels=str_wrap(rev(c("Equatorial","Tropical",
                                         "Subtropical", "Grassland", "Arid",
                                         "Temperate","Temperate Tasmania")),width = 10),
                   position = 'left')+
  labs(x=expression(paste(Cover~Trend~('%'~yr**-1))), 
       y=NULL)+
  facet_wrap(~ measure)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.text = element_text(face='bold'),
        legend.position = 'none'); p_vcf




# library(cowplot)
# p_left <- p_kop+draw_label(label='(a)', x=140.75,y=-11,size = 25)
# p_top <- ggdraw(p_ppet)+draw_label(label='(b)',x=0.03,y=0.93,size=25)
# p_mid <- ggdraw(p_ndvi_sen)+draw_label(label='(c)', x=0,y=0.96, size=25,hjust=0.1,vjust=0.1)
#                                        #hjust=0,vjust = 0)
# p_bot <- ggdraw(p_vcf)+draw_label(label = '(d)', x=0.035,y=0.94,size=25)
# 
# cp_r <- cowplot::plot_grid(p_top,p_mid,p_bot,
#                            nrow = 3,
#                            rel_heights = c(2,3,3))
# ggsave(plot=cowplot::plot_grid(p_left,cp_r),
#        filename = "figures/Fig8_map_7KoppenZones_PPET_change_ThielSen_VCF.png", 
#        width = 25, height=30, units='cm', dpi=350, type='cairo')
# ggsave(plot=cowplot::plot_grid(p_left,cp_r),
#        filename = "doc/submission_pnas_1/Fig5_map_7KoppenZones_PPET_change_ThielSen_VCF.pdf", 
#        width = 25, height=30, units='cm', dpi=350)
# ggsave(plot=cowplot::plot_grid(p_left,cp_r),
#        filename = "doc/submission_pnas_1/Fig5_map_7KoppenZones_PPET_change_ThielSen_VCF.tiff", 
#        width = 25, height=30, units='cm', dpi=350)


library(patchwork)

p_out <- p_ppet/p_ndvi_sen/p_vcf+
    plot_annotation(tag_prefix = '(', tag_levels = 'a',tag_suffix = ')')+
  plot_layout(byrow = T, 
    heights = c(0.5,1,0.75))
ggsave(
       filename = "figures/Fig8_PPET_change_ThielSen_VCF.png", 
       width = 19, height=26, units='cm', dpi=350, device = grDevices::png)


