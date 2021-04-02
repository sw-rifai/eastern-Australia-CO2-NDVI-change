# !diagnostics off
#*******************************************************************************
#* Description:
#* Plot Climate stuff
#*
#*
#*
library(tidyverse); 
library(data.table); setDTthreads(threads = 0)
library(lubridate); 
library(dtplyr); 
library(RcppArmadillo); library(patchwork)
library(sf); library(zyp)
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

ldat <- dat %>% lazy_dt()
# END Load awap clim dat *****************************************************************


# Load simplified BOM Koppen climate zones --------------------------------
library(stars)
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
# arrow::write_parquet(kop, sink='../data_general/Koppen_climate/BOM_Koppen_simplified7.parquet')
#*** End Kop zone load ********************************************************



# long-term trend of annual P:PET ----------------------------------------------
system.time(
  lt_PPET_annual <- dat[hydro_year %in% c(1982:2019)] %>% 
    .[,.(val = mean(pe_12mo, na.rm=TRUE)), by=.(x,y,hydro_year)] %>% 
    .[is.na(val)==F] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[,.(beta = list(coef(zyp.sen(val~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)
lt_PPET_annual[is.na(b0)==F]
sum(lt_PPET_annual$b1<0)/length(lt_PPET_annual$b1) # fraction of grid cells with negative trend
quantile(lt_PPET_annual$b1,c(0.01,0.99))




# Join NDVI linear change by epoch & P:PET change -------------------------------
lt_PPET_annual$b1 %>% quantile(., c(0.01,0.99))
p_left <- lt_PPET_annual %>%
  ggplot(data=., aes(x,y,fill=b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray70',color='gray10')+
  geom_tile()+
  scale_fill_gradient2(expression(paste(Delta~frac(P,PET)~(yr**-1))),
                       limits=c(-0.008,0.008),
                       breaks=c(-0.008,-0.004,0,0.004,0.008), 
                       labels=c("< -0.008",-0.004,0,0.004,"> 0.008"),
                       high=scico::scico(11,palette = 'vik')[1],
                       mid="#ffffe3",
                       # mid=scico::scico(11,palette = 'vik')[6],
                       low=scico::scico(11,palette = 'vik')[11],
                       oob=scales::squish,
                       na.value = 'gray')+
  labs(x=NULL,y=NULL)+
  coord_sf(xlim = c(140,155),
           ylim = c(-45,-10), expand = FALSE)+
  scale_x_continuous(breaks=seq(140,155,by=5))+
  guides(fill = guide_colorbar(title.position = 'top'))+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(),
        legend.position = c(0.825,0.88),
        # legend.position = 'bottom',
        legend.title = element_text(size=8),
        legend.text.align = 0,
        legend.key.width = unit(0.2, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        strip.text = element_text(face='bold')
        # axis.text = element_blank(),
        # axis.ticks = element_blank()
        ); p_left


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
  ggplot(data=., aes(x,y,fill=zone))+
  geom_sf(inherit.aes = F, data=oz_poly,
          fill='gray40',color='gray10')+
  geom_tile()+
  coord_sf(xlim = c(140,155.5),
           ylim = c(-45,-10), 
           expand = FALSE)+
  scale_x_continuous(breaks=seq(140,155,by=5))+
  scale_fill_viridis_d(option='B',direction = 1,end=0.95, label = lut_kop,na.translate=F)+
  # scico::scale_fill_scico_d(end=0.9,direction = 1)+
  labs(x=NULL,y=NULL)+
  theme_linedraw()+
  guides(fill=guide_legend(title='Zone', 
                           title.position = 'top'))+
  theme(#legend.position = c(0.75,0.85),
    legend.position = c(1,1), 
    legend.justification = c(1.01,1.005),
    legend.box.background = element_rect(fill="#FFFFFF11",color=NA,linetype = NULL),
    legend.text = element_text(size=8),
    legend.title = element_text(size=8.5),
    legend.spacing = unit(10,'points'),
    legend.key.width = unit(0.3,'cm'),
    legend.direction = 'vertical',
    panel.grid = element_blank(), 
    panel.background = element_rect(fill='lightblue')); p_kop



# p_left+p_right+plot_layout(widths = c(1.25,2), heights = c(1.1,1),
#                            guides = 'keep')+plot_annotation(tag_levels = 'A')
# ggsave(p_left+p_right+
#          plot_layout(widths = c(1.25,2), guides = 'keep')+
#          plot_annotation(tag_levels = 'A'),
#        filename = "figures/Fig2_PPET_ndvi_seasonal_TheilSen_trend_by_epoch.png", 
#        dpi=350, width=15,height=15,units='cm',type='cairo')

ggsave(cowplot::plot_grid(p_left,p_kop,ncol=2,
                          labels=c('(a)','(b)'),
                          label_x = c(-0.035,-0.05),
                          label_y = c(0.985,0.985),
                          rel_widths = c(1,1)),
       filename = "figures/Fig_PPET-TheilSen_map-Koppen.png", 
       dpi=350, width=15,height=15.5,units='cm',type='cairo')
# ggsave(cowplot::plot_grid(p_left,p_right,ncol=2,labels=c('(a)','(b)'), 
#                           rel_widths = c(1,2)),
#        filename = "doc/submission_pnas_1/Fig2_PPET_ndvi_seasonal_TheilSen_trend_by_epoch.pdf", 
#        dpi=350, width=15,height=15.5,units='cm')






# # Calc Thiel Sen Trends ---------------------------------------------------
# vec_col <- RColorBrewer::brewer.pal(n=7, name='BrBG')
# p_right <- bind_rows(sen_ndvi_season_e1, sen_ndvi_season_e2) %>% 
#   ggplot(data=., aes(x,y,fill=b1))+
#   geom_sf(inherit.aes = F, data=oz_poly,fill='gray70',color='gray10')+
#   geom_tile()+
#   scale_fill_gradient2(expression(paste(Delta*NDVI~yr^-1)),
#                        high=vec_col[7], 
#                        # mid=vec_col[4],
#                        mid="#ffffe3",
#                        low=vec_col[1],
#                        limits=c(-0.006,0.006),
#                        breaks=c(-0.006,0,0.006),
#                        labels=c("< -0.006",0,"> 0.006"),
#                        oob=scales::squish,
#                        na.value='gray')+
#   labs(x=NULL,y=NULL)+
#   coord_sf(xlim = c(140,154),
#            ylim = c(-45,-10), expand = FALSE)+
#   facet_grid(epoch~season)+
#   guides(fill = guide_colourbar(label = T)) +
#   theme(panel.background = element_rect(fill = '#99A3C4'), 
#         panel.grid = element_blank(), 
#         # legend.position = c(0.475, 0.01),
#         legend.position = 'bottom',
#         legend.direction = 'horizontal',
#         legend.title = element_text(size=8),
#         legend.margin=margin(0,0,0,0),
#         # legend.box.margin=margin(-10,-10,-10,-10),
#         # legend.key.width = unit(1.5,'cm'),
#         strip.text = element_text(face='bold'),
#         axis.text = element_blank(), 
#         axis.ticks = element_blank()); p_right
# 
# library(zyp)
# sen_ndvi_season_e1 <- dat %>% 
#   .[date >= ymd("1981-11-01") & date <= ymd("2000-12-31")] %>% 
#   .[hydro_year %in% c(1982:2000)] %>% 
#   .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
#   .[is.na(ndvi_hyb)==F] %>% 
#   .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE)), by=.(x,y,season,hydro_year_c)] %>% 
#   .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y,season)] %>%
#   .[,`:=`(b0=unlist(beta)[1], 
#           b1=unlist(beta)[2]), by=.(x,y,season)]
# sen_ndvi_season_e2 <- dat %>%
#   .[date >= ymd("2001-01-01") & date <= ymd("2019-08-30")] %>%
#   .[hydro_year %in% c(2001:2019)] %>%
#   .[,`:=`(hydro_year_c = hydro_year-2001)] %>%
#   .[is.na(ndvi_hyb)==F] %>%
#   .[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE)), by=.(x,y,season,hydro_year_c)] %>% 
#   .[,.(beta = list(coef(zyp.sen(ndvi_hyb~hydro_year_c)))),by=.(x,y,season)] %>%
#   .[,`:=`(b0=unlist(beta)[1],
#           b1=unlist(beta)[2]), by=.(x,y,season)]
# 
# sen_ndvi_season_e1$epoch <- "AVHRR NDVI 1982-2000"
# sen_ndvi_season_e2$epoch <- "MODIS NDVI 2001-2019"
