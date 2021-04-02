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

ldat <- dat %>% lazy_dt()
# END Load awap clim dat *****************************************************************


# Calc Thiel Sen Trends by epoch ---------------------------------------------------
library(zyp)
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



dat_annual <- dat %>% 
  .[date >= ymd("1981-08-01") & date <= ymd("2019-08-31")] %>% 
  .[,.(ndvi = mean(ndvi_hyb,na.rm=TRUE)), by=.(x,y,hydro_year)] %>% 
  .[hydro_year %in% c(1982:2019)] %>% 
  .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
  .[,`:=`(epoch = ifelse(hydro_year<=2000, 0, 1))] %>% 
  .[is.na(ndvi)==F]
dat_nobs <- dat_annual[,.(nobs=sum(is.na(ndvi)==F)),by=.(x,y)]
dat_annual <- merge(dat_annual,dat_nobs,by=c('x','y'))
dat_annual <- dat_annual[nobs>18] # ONLY locations with more than 18 observed year (~50% of the ts)

rlm_ndvi_annual <-  dat_annual[,.(beta = list(coef(MASS::rlm(ndvi~hydro_year_c+jitter(epoch))))),by=.(x,y)] %>%
  .[,`:=`(b0=unlist(beta)[1], 
          b1=unlist(beta)[2]), by=.(x,y)]


# Calc Thiel Sen Trends ---------------------------------------------------
vec_col <- RColorBrewer::brewer.pal(n=7, name='BrBG')
p_ndvi <- bind_rows(sen_ndvi_season_e1, sen_ndvi_season_e2) %>% 
  ggplot(data=., aes(x,y,fill=b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray60',color='gray10')+
  geom_tile()+
  scale_fill_gradient2(expression(paste(Delta*NDVI~yr^-1)),
                       high=vec_col[7], 
                       # mid=vec_col[4],
                       mid="#ffffe3",
                       low=vec_col[1],
                       limits=c(-0.006,0.006),
                       breaks=c(-0.006,-0.003,0,0.003,0.006),
                       labels=c("< -0.006","  -0.003","  0","  0.003","> 0.006"),
                       oob=scales::squish,
                       na.value='gray')+
  labs(x=NULL,y=NULL)+
  coord_sf(xlim = c(140,154),
           ylim = c(-44,-10), expand = FALSE)+
  facet_grid(epoch~season)+
  guides(fill = guide_colourbar(label = T)) +
  theme(#panel.background = element_rect(fill = '#99A3C4'), 
        panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        # legend.position = c(0.475, 0.01),
        legend.position = 'none',
        legend.direction = 'vertical',
        legend.title = element_text(size=8),
        legend.margin=margin(0,0,0,0),
        # legend.box.margin=margin(-10,-10,-10,-10),
        legend.key.width = unit(0.35,'cm'),
        legend.key.height = unit(1,'cm'),
        strip.text = element_text(face='bold'),
        axis.text = element_blank(), 
        axis.ticks = element_blank()); p_ndvi


colorspace::hcl_palettes(plot=T)
p_annual <- rlm_ndvi_annual[b0>0.1] %>% 
  ggplot(data=., aes(x,y,fill=b1))+
  geom_sf(inherit.aes = F, data=oz_poly,fill='gray60',color='gray10')+
  geom_tile()+
  # scale_fill_steps2(expression(paste(Delta*NDVI~yr^-1)),
  #                      high=vec_col[7],
  #                      mid="#ffffe3",
  #                      low=vec_col[1],
  #                      n.breaks=9,
  #                      nice.breaks=FALSE,
  #                      limits=c(-0.005,0.005),
  #                      # breaks=c(-0.006,-0.003,0,0.003,0.006),
  #                      # labels=c("< -0.006","  -0.003","  0","  0.003","> 0.006"),
  #                      oob=scales::squish,
  #                      na.value='gray')+
  # scale_fill_binned(n.breaks=7, 
  #                   steps=RColorBrewer::brewer.pal(n=7, name='BrBG'))+
  # colorspace::scale_fill_binned_diverging(palette='Green-Brown',rev=T)+
  scale_fill_gradient2(expression(paste(Delta*NDVI~yr^-1)),
                       high=vec_col[7],
                       # mid=vec_col[4],
                       mid="#ffffe3",
                       low=vec_col[1],
                       limits=c(-0.006,0.006),
                       breaks=c(-0.006,-0.003,0,0.003,0.006),
                       labels=c("< -0.006","  -0.003","  0","  0.003","> 0.006"),
                       oob=scales::squish,
                       na.value='gray')+
  labs(x=NULL,y=NULL)+
  coord_sf(xlim = c(140,154),
           ylim = c(-44,-10), expand = FALSE)+
  scale_x_continuous(breaks=c(140,145,150))+
  guides(fill = guide_colourbar(label = T)) +
  theme(#panel.background = element_rect(fill = '#99A3C4'), 
    panel.background = element_rect(fill='lightblue'),
    panel.grid = element_blank(), 
    # legend.position = c(0.475, 0.01),
    legend.position = c(1,1),
    legend.justification = c(1,1),
    legend.direction = 'vertical',
    legend.title = element_text(size=9),
    legend.margin=margin(5,5,5,5),
    # legend.box.margin=margin(-10,-10,-10,-10),
    legend.key.width = unit(0.2,'cm'),
    legend.key.height = unit(0.75,'cm'),
    legend.background = element_rect(fill="#FFFFFF99"),
    strip.text = element_text(face='bold')
    # axis.text = element_blank(), 
    # axis.ticks = element_blank()
    ); p_annual

library(cowplot)
p_left <- ggdraw(p_annual)+draw_label(label='(a)', x=0.075,y=0.95,size = 15,fontface='bold')
p_right <- ggdraw(p_ndvi)+draw_label(label='(b)', x=0.055,y=0.95,size = 15,fontface='bold')

p_joint <- cowplot::plot_grid(p_left,p_right,
                           nrow = 1,ncol = 2,
                           rel_widths = c(1,1.85))
p_joint
ggsave(plot=p_joint,
       filename = "figures/Fig_ndvi-annual-rlm_ndvi-seasonal-TheilSen-trend-by-epoch.png",
       width = 20, height=19, units='cm', dpi=350, type='cairo')

