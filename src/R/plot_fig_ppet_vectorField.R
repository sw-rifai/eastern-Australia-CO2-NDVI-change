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
options(pillar.sigfig = 5)
# IMPORT DATA ###################################################################

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


# Arrow Vector Field P:PET & NDVI Shift WITH INSET{MAPPET 0-0.75} -------------------------------------------------------------
set.seed(321)
epoch1 <- ldat %>% 
  filter(ndvi_anom_sd > -3.5) %>% 
  filter(ndvi_anom_sd < 3.5) %>% 
  filter(hydro_year %in% c(1982:1986)) %>% 
  group_by(x,y) %>% 
  summarize(ppet = mean(pe_12mo,na.rm=TRUE), 
            ndvi_u = mean(ndvi_hyb,na.rm=TRUE), 
            nobs = n(), 
            mape = mean(mape)) %>%
  ungroup() %>% 
  as_tibble()
epoch2 <- ldat %>% 
  filter(date<=ymd("2019-08-30")) %>% 
  filter(ndvi_anom_sd > -3.5) %>% 
  filter(ndvi_anom_sd < 3.5) %>% 
  filter(hydro_year %in% c(2015:2019)) %>% 
  group_by(x,y) %>% 
  summarize(ppet = mean(pe_12mo,na.rm=TRUE), 
            ndvi_u = mean(ndvi_3mo,na.rm=TRUE), 
            nobs = n()) %>% 
  ungroup() %>% 
  as_tibble()

o <- inner_join(epoch1,epoch2,by=c("x","y"),suffix=c("_1","_2")) %>% 
  # filter(ppet_1 <= 1.5) %>% 
  # filter(mape<=2) %>% 
  filter(ndvi_u_1 > 0 & ndvi_u_2 > 0) %>% 
  mutate(delta_x = ppet_2 - ppet_1) %>% 
  mutate(id = cur_group_rows())
o <- o %>% filter(nobs_1 >= 25 & 
                    nobs_2 >= 25)

set.seed(321)
vec_ids <- sample(o$id, 1000)
p_vector <- o %>% 
  filter(id %in% vec_ids) %>% 
  ggplot(data=., aes(ppet_1, ndvi_u_1,color=delta_x))+
  geom_point(data=epoch2 %>% filter(ndvi_u>0), 
             aes(ppet, ndvi_u), color=NA)+
  geom_segment(aes(xend=ppet_2,yend=ndvi_u_2),
               arrow=arrow(length=unit(0.1,'cm')), 
               lwd=0.25)+
  scale_color_gradientn(colours = c("#cf0000",
                                    "#b87b7b",
                                    "#7b81b8","navy"), 
                        limits=c(-0.15,0.15), oob=scales::squish)+
  labs(x='Annual P:PET',y='NDVI')+
  scale_x_continuous(limits=c(0,2.25),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme_linedraw()+
  guides(color=guide_colorbar(title.position = 'top', 
                              title = expression(paste(Delta~P*":"*PET,~(mm~mm**-1)))))+
  theme(legend.direction = 'vertical',
        legend.justification = c(0, 1), 
        legend.position = c(0, 1),
        legend.background = element_rect(fill=NA),
        legend.key.width = unit(0.2,units = 'cm'),
        legend.key.height = unit(0.7,units = 'cm'),
        # legend.position = c(0.05,0.75),
        axis.text = element_text(size=12),
        panel.grid = element_blank()); p_vector
p_out <- ggExtra::ggMarginal(p_vector, type='histogram')
p_out

p_vector_inset <- o %>% 
  filter(id %in% vec_ids) %>%
  # filter(ppet_1 < 0.6) %>%
  filter(mape < 0.4) %>% 
  # sample_frac(0.02) %>%
  # sample_n(10) %>% 
  mutate(ord = abs(delta_x)) %>% 
  arrange(ord) %>% 
  ggplot(data=., aes(ppet_1, ndvi_u_2,color=delta_x))+
  geom_segment(aes(xend=ppet_2,yend=ndvi_u_1),
               arrow=arrow(length=unit(0.1,'cm')), 
               lwd=0.4)+
  scale_color_gradientn(colours = c("#cf0000","#b87b7b",
                                    "#7b81b8","navy"), 
                        limits=c(-0.15,0.15), oob=scales::squish)+
  labs(x=NULL,y=NULL)+
  # labs(x='P:PET',y='NDVI')+
  scale_x_continuous(limits=c(0.1,0.4),
                     expand=c(0,0.01))+
  scale_y_continuous(expand=c(0,0.01), 
                     limits=c(0.2,0.5)
  )+
  coord_equal()+
  theme_linedraw()+
  guides(color=guide_colorbar(title.position = 'top' 
                              # title = expression(paste(Delta~P*":"*PET))
  ))+
  theme(legend.direction = 'horizontal',
        legend.position = 'none',# c(0.75,0.15), 
        panel.grid = element_blank(),
        panel.background = element_rect(fill=NA),
        axis.text = element_text(size=10)); p_vector_inset

ggsave(p_out, 
       filename = 'figures/ndvi_ppet_1981_1985_shift_2015_2019_wMarginDistribution_wInset_p1.png',
       width=16, height = 12, units='cm',type='cairo', dpi=350)
ggsave(p_vector_inset, 
       filename = 'figures/ndvi_ppet_1981_1985_shift_2015_2019_wMarginDistribution_wInset_p2.png',
       width=6.5, height = 6.5, units='cm',type='cairo')

library(magick)
p0 <- magick::image_read("figures/diagram_PPET_CO2_response.png")
p1 <- magick::image_read('figures/ndvi_ppet_1981_1985_shift_2015_2019_wMarginDistribution_wInset_p1.png')
p2 <- magick::image_read('figures/ndvi_ppet_1981_1985_shift_2015_2019_wMarginDistribution_wInset_p2.png')

p1_2 <- image_composite(p1,p2,offset="+1040+700")

p_0_1_2 <- magick::image_append(c(p0,p1_2),
                                stack=TRUE)
magick::image_write(p1_2, path="figures/Fig1_ndvi_ppet_vectorPlot_wInset.png")
magick::image_write(p1_2, format = 'tiff', 
                    path="doc/submission_pnas_1/Fig1_ndvi_ppet_vectorPlot_wInset.tiff")
