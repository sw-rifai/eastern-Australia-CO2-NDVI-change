library(tidyverse)
library(data.table); setDTthreads(threads = 0)
library(lubridate); 
library(dtplyr);
library(RcppArmadillo)
library(sf); library(stars)
library(patchwork); library(colorspace) 
library(zyp); 
library(mgcv)
set.seed(333)
# IMPORT ###################################################################
oz_poly <- sf::read_sf("../data_general/GADM/gadm36_AUS.gpkg", 
                       layer="gadm36_AUS_1")
oz_poly <- st_as_sf(oz_poly)
oz_poly <- st_simplify(oz_poly, dTolerance = 1000)

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
                             "precip", "precip_12mo","precip_anom", "precip_anom_12mo","map",
                             "vpd15","vpd15_12mo","vpd15_anom", "vpd15_anom_12mo","mavpd15",
                             "pet","pet_anom","pet_anom_12mo","pet_12mo", "mapet",
                             "pe","pe_12mo","pe_anom_12mo", 
                             'vc','veg_class',
                             'month',
                             "x", "y", "year")) %>% 
  as.data.table()

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
setDT(dat)
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
  ungroup() %>% 
  as.data.table()

# summarize to annual -----------------------------------------------------
dat_annual <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
  .[date>= ymd("1981-11-01") & date<= ymd("2019-08-30")] %>% 
  # .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
  .[,`:=`(epoch=ifelse(hydro_year < 2001,0,1))] %>% 
  .[is.na(ndvi_hyb)==F] %>% 
  .[is.na(pe_anom_12mo)==F] %>% 
  merge(., coords_keep,by=c('x','y'))

dat_annual <- dat_annual[,`:=`(hydro_year_c = hydro_year-1982, 
                               frac_p_anom = precip_anom_12mo/map, 
                               frac_ppet_anom = pe_anom_12mo/mappet, 
                               frac_pet_anom = pet_anom_12mo/mapet, 
                               frac_vpd_anom = vpd15_anom_12mo/mavpd15)]

dat_annual <- dat_annual %>% 
  select(x,y,date,ndvi_hyb,co2_trend,
         hydro_year,hydro_year_c,
         frac_p_anom,frac_ppet_anom,frac_pet_anom,frac_vpd_anom,
         map,mappet,mapet,mavpd15,
         epoch,nobs_annual) %>% 
  as.data.table()

dat_annual <- dat_annual[nobs_annual>=10]
dat_annual <- dat_annual[,.(ndvi_hyb = mean(ndvi_hyb,na.rm=TRUE),
                            co2 = mean(co2_trend,na.rm=TRUE),
                            frac_p_anom = mean(frac_p_anom,na.rm=TRUE), 
                            frac_pet_anom = mean(frac_pet_anom,na.rm=TRUE),
                            frac_ppet_anom = mean(frac_ppet_anom,na.rm=TRUE),
                            frac_vpd_anom = mean(frac_vpd_anom,na.rm=TRUE),
                            epoch = mean(epoch, na.rm=TRUE), 
                            map = mean(map,na.rm=TRUE), 
                            mapet = mean(mapet, na.rm=TRUE), 
                            mappet = mean(mappet, na.rm=TRUE), 
                            mavpd15 = mean(mavpd15,na.rm=TRUE)), 
                         by=.(x,y,hydro_year_c)]
dat_annual <- dat_annual[is.infinite(frac_ppet_anom)==F]



# Fit GAM of CO2 effect and climate  --------------------------------------
dat_train <- dat_annual[sample(.N, floor(nrow(dat_annual)/2))]
dat_test <- fsetdiff(dat_annual,dat_train)

g3_eval <- bam(ndvi_hyb~
            te(mappet,co2,k=5)+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs')+
            # te(mappet,frac_ppet_anom,k=5,bs='cs')+
            epoch,
          data=dat_train, 
          discrete = T,select=TRUE, method='fREML')
yardstick::rsq_trad_vec(truth=dat_test$ndvi_hyb, estimate=predict(g3_eval,newdata=dat_test))
yardstick::rmse_vec(truth=dat_test$ndvi_hyb, estimate=predict(g3_eval,newdata=dat_test))

g3_full <- bam(ndvi_hyb~
            te(mappet,co2,k=5)+
            te(mavpd15,frac_vpd_anom,k=5,bs='cs')+
            te(map,frac_p_anom,k=5,bs='cs')+
            te(mapet,frac_pet_anom,k=5,bs='cs')+
            # te(mappet,frac_ppet_anom,k=5,bs='cs')+
            epoch,
          data=dat_annual, 
          discrete = T,select=TRUE, method='fREML')




## Regression NDVI robust linear model ---------------------------------------
system.time(
  lt_ndvi_year <-  dat_annual[is.na(ndvi_hyb)==FALSE] %>% 
    .[,.(beta = list(unname(MASS::rlm(ndvi_hyb~hydro_year_c+jitter(epoch), 
                                      data=.SD)$coefficients))), 
      by=.(x,y)] %>% 
    .[,`:=`(b0=unlist(beta)[1], b1=unlist(beta)[2],b2=unlist(beta)[3],b3=unlist(beta)[4]
    ), by=.(x,y)]
)

# % increase in annual NDVI ------
lt_ndvi_year %>% as_tibble() %>% 
  mutate(val = 100*(37*b1)/(b0)) %>% 
  pull(val) %>% 
  quantile(., c(0.01,0.05,0.25,0.5,0.75,0.95,0.99))

# fraction of area experiencing greening ------
sum(lt_ndvi_year$b1>0)/nrow(lt_ndvi_year)


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

pred_vpd_e1 <- merge(unique(dat_annual[,.(x,y,mavpd15)]), 
  lt_v_sen[,.(x,y,b0)], 
  by=c("x","y")) %>% 
  .[,`:=`(e='e1', 
          frac_vpd_anom = (b0-mavpd15)/mavpd15)]
pred_vpd_e2 <- lt_v_sen[,`:=`(vpd15_0 = b0,
                              vpd15_1 = b0+16*b1,
                              vpd15_2 = b0+37*b1)] %>% 
  merge(unique(dat_annual[,.(x,y,mavpd15)]), ., by=c("x","y")) %>% 
  .[,`:=`(frac_vpd_anom = (vpd15_2-mavpd15)/mavpd15, 
          e='e2')]


system.time(
  lt_p_sen <- dat[date>=ymd("1981-12-01")][date<=ymd("2019-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[is.na(vpd15)==F] %>% 
    .[,.(precip_12mo = mean(precip_12mo,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(precip_12mo~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)
pred_p_e1 <- merge(unique(dat_annual[,.(x,y,map)]), 
  lt_p_sen[,.(x,y,b0,b1)], 
  by=c("x","y")) %>% 
  .[,`:=`(e='e1', 
          frac_p_anom = (b0-map)/map)]
pred_p_e2 <- merge(unique(dat_annual[,.(x,y,map)]), 
  lt_p_sen[,.(x,y,b0,b1)], 
  by=c("x","y")) %>%
  .[,`:=`(
          frac_p_anom = (b0+37*b1 - map)/map, 
          e = 'e2')]

system.time(
  lt_pet_sen <- dat[date>=ymd("1981-12-01")][date<=ymd("2019-11-30")] %>% 
    .[,`:=`(hydro_year_c = hydro_year-1982)] %>% 
    .[is.na(vpd15)==F] %>% 
    .[,.(pet_12mo = mean(pet_12mo,na.rm=TRUE)), by=.(x,y,hydro_year_c)] %>% 
    .[,.(beta = list(coef(zyp.sen(pet_12mo~hydro_year_c)))),by=.(x,y)] %>%
    .[,`:=`(b0=unlist(beta)[1], 
            b1=unlist(beta)[2]), by=.(x,y)]
)

pred_pet_e1 <- merge(unique(dat_annual[,.(x,y,mapet)]), 
  lt_pet_sen[,.(x,y,b0,b1)], 
  by=c("x","y")) %>% 
  .[,`:=`(e='e1', 
          frac_pet_anom = (b0-mapet)/mapet)]
pred_pet_e2 <- merge(unique(dat_annual[,.(x,y,mapet)]), 
  lt_pet_sen[,.(x,y,b0,b1)], 
  by=c("x","y")) %>% 
  .[,`:=`(
          frac_pet_anom = (b0+37*b1-mapet)/mapet, 
          e = 'e2')]
pred_pet_e1$frac_pet_anom %>% hist
pred_pet_e2$frac_pet_anom %>% hist


d_preds <- inner_join(
 bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
          tibble(pred_vpd_e2) %>% select(x,y,e,frac_vpd_anom)), 
  bind_rows(
  pred_p_e1 %>% select(x,y,e,frac_p_anom) %>% as_tibble(),
  pred_p_e2 %>% select(x,y,e,frac_p_anom) %>% as_tibble()), 
  by=c("x","y","e"))

d_preds <- inner_join(d_preds, 
  bind_rows(
  pred_pet_e1 %>% select(x,y,e,frac_pet_anom) %>% as_tibble(),
  pred_pet_e2 %>% select(x,y,e,frac_pet_anom) %>% as_tibble()), 
  by=c("x","y","e"))

d_preds %>% 
    inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  mutate(co2=340.8, 
    epoch=0) %>% 
  mutate(pred = predict(g3_full, newdata=.)) %>% 
  pivot_wider(.,id_cols = c('x','y'), 
    names_from=c('e'), values_from=c('pred')) %>% 
  filter(e1>0.1) %>% 
  mutate(dNDVI = 100*((e2-e1)/e1))

summary(g3_full)
dat_annual %>% 
  ggplot(data=.,aes(frac_pet_anom))+
  geom_density()
d_preds %>% 
  ggplot(data=.,aes(frac_p_anom,fill=e))+
  geom_density()
d_preds %>% 
  ggplot(data=.,aes(frac_vpd_anom,fill=e))+
  geom_density()
d_preds %>% 
  ggplot(data=.,aes(frac_pet_anom,fill=e))+
  geom_density()

d_no_co2 <- d_preds %>% 
    inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  mutate(co2=340.8, 
    epoch=0) %>% 
  mutate(pred_no_co2 = predict(g3_full, newdata=.)) %>% 
    pivot_wider(.,id_cols = c('x','y'), 
    names_from=c('e'), 
    values_from=c('pred_no_co2')) %>% 
  as.data.table() %>% 
  filter(e1>0.1) %>% 
  mutate(dNDVI_no_co2 = 100*((e2-e1)/e1)) %>% 
  as.data.table()

d_co2 <- d_preds %>% 
    inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  mutate(co2=case_when(
    e=='e1'~340.8, 
    e=='e2'~411.9), 
    epoch=case_when(e=='e1'~0,
                    e=='e2'~1)) %>% 
  mutate(pred_co2 = predict(g3_full, newdata=.)) %>% 
    pivot_wider(.,id_cols = c('x','y'), 
    names_from=c('e'), 
    values_from=c('pred_co2')) %>% 
  as.data.table() %>% 
  filter(e1>0.1) %>% 
  mutate(dNDVI_co2 = 100*((e2-e1)/e1)) %>% 
  as.data.table()

vec_cols <- pals::brewer.brbg(n=9)
merge(d_no_co2[,.(x,y,dNDVI_no_co2)], 
  d_co2[,.(x,y,dNDVI_co2)], 
  by=c("x","y")) %>% 
  .[,`:=`(difference = dNDVI_co2 - dNDVI_no_co2)] %>% 
  select(x,y,dNDVI_no_co2,dNDVI_co2,difference) %>%
  as_tibble() %>% 
  pivot_longer(cols=c("dNDVI_no_co2","dNDVI_co2","difference")) %>% 
  filter(name!='difference') %>% 
  mutate(name = case_when(
    name=="dNDVI_no_co2"~"No CO2 change",
    name=="dNDVI_co2"~"With CO2 change"
  )) %>% 
  ggplot(data=.,aes(x,y,fill=value))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  scale_fill_gradient2(limits=c(-20,20), 
    low=first(vec_cols),
    mid=c("#FFFED4"), 
    high=last(vec_cols),
    oob=scales::squish)+
  # scale_fill_viridis_b(n.breaks=12, option='H', direction = -1)+
  # scale_fill_binned_diverging(
  #   palette='Berlin',rev=T,n.breaks=10, limits=c(-6,1),oob=scales::squish)+
  # scico::scale_fill_scico("%",
  #                         palette = 'bamako', 
  #                         direction = -1,
  #                         limits=c(5,15), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  labs(x=NULL,y=NULL, 
    title=expression(paste("Counterfactual Relative Difference in NDVI: 1982-2019")), 
    fill=expression(paste(Delta*NDVI[GAM~Pred.],('%'))))+
  facet_wrap(~name,nrow = 1)+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank());
scale_factor <- 0.5
ggsave(
  filename = "figures/Fig_map_counterfactual_no_co2_comparison.png",
  device=grDevices::png,
  width = 32.5*scale_factor,
  height=35*scale_factor, 
  units='cm',
  dpi=350)





# SCRAP -------------------------------------------------
d_preds %>% 
    inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  mutate(co2=340.8, 
    epoch=0) %>% 
  mutate(pred = predict(g3_full, newdata=.)) %>% 
  pivot_wider(.,id_cols = c('x','y'), 
    names_from=c('e'), 
    values_from=c('pred','frac_p_anom','frac_vpd_anom','frac_pet_anom')) %>% 
  filter(pred_e1>0.1) %>% 
  mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
  mutate(dP = 100*((frac_p_anom_e2-frac_p_anom_e1)/frac_p_anom_e1),
         dPET = 100*((frac_pet_anom_e2-frac_pet_anom_e1)/frac_pet_anom_e1),
         dVPD = 100*((frac_vpd_anom_e2-frac_vpd_anom_e1)/frac_vpd_anom_e1)) %>% 
  select(x,y,dNDVI,dP,dPET,dVPD) %>% 
  pivot_longer(cols=c("dNDVI","dP","dPET","dVPD")) %>% 
  ggplot(data=.,aes(x,y,fill=value))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  scale_fill_gradient2(limits=c(-20,20), 
    mid=c("#FFFED4"), 
    oob=scales::squish)+
  # scale_fill_viridis_b(n.breaks=12, option='H', direction = -1)+
  # scale_fill_binned_diverging(
  #   palette='Berlin',rev=T,n.breaks=10, limits=c(-6,1),oob=scales::squish)+
  # scico::scale_fill_scico("%",
  #                         palette = 'bamako', 
  #                         direction = -1,
  #                         limits=c(5,15), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  labs(x=NULL,y=NULL, 
    title=expression(paste(No~CO[2]~VPD~Effect)), 
    fill=expression(paste(Delta*NDVI[GAM~Pred.],('%'))))+
  facet_wrap(~name,nrow = 1)+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank());






# GAM attributed CO2 effect on relative dNDVI
bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
          tibble(pred_vpd_e2) %>% select(x,y,e,frac_vpd_anom)) %>% 
  inner_join(., tibble(co2=c(340.8,411.9),e=c('e1','e2')), 
             by=c('e')) %>% 
  mutate(frac_p_anom=0, 
         frac_ppet_anom=0,
         frac_pet_anom=0,
         # frac_vpd_anom=0,
         epoch=0) %>%
  inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  mutate(pred = predict(g3_full,newdata=.,type='response')) %>%   
  filter(is.na(pred)==F) %>% 
  select(x,y,e,frac_vpd_anom,pred,co2) 
%>% 
  pivot_wider(., names_from=c('e'), values_from=c('pred','co2','frac_vpd_anom')) %>% 
  mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
  pull(dNDVI) %>% summary


pred_p_e1
# GAM attributed CO2 effect on relative dNDVI
bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
          tibble(pred_vpd_e2) %>% select(x,y,e,frac_vpd_anom)) %>% 
  inner_join(., tibble(co2=c(340.8,340.8),e=c('e1','e2')), 
             by=c('e')) %>% 
  mutate(frac_p_anom=0, 
         frac_ppet_anom=0,
         frac_pet_anom=0,
         # frac_vpd_anom=0,
         epoch=0) %>%
  inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  mutate(pred = predict(g3_full,newdata=.,type='response')) %>%   
  filter(is.na(pred)==F) %>% 
  select(x,y,e,frac_vpd_anom,pred,co2) %>% 
  pivot_wider(., names_from=c('e'), values_from=c('pred','co2','frac_vpd_anom')) %>% 
  mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
  pull(dNDVI) %>% summary

# PLOTTING --------------------------------------------
## VPD effect w/no change in CO2 -----
p_gam_pred_noCO2_wVPD <- bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
          tibble(pred_vpd_e2) %>% select(x,y,e,frac_vpd_anom)) %>% 
  inner_join(., tibble(co2=c(340.8,340.8),e=c('e1','e2')), 
             by=c('e')) %>% 
  mutate(frac_p_anom=0, 
         frac_ppet_anom=0,
         frac_pet_anom=0,
         # frac_vpd_anom=0,
         epoch=0) %>%
  inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  mutate(pred = predict(g3_full,newdata=.,type='response')) %>%   
  filter(is.na(pred)==F) %>% 
  select(x,y,e,frac_vpd_anom,pred,co2) %>% 
  pivot_wider(., names_from=c('e'), values_from=c('pred','co2','frac_vpd_anom')) %>% 
  mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
  mutate(d_vpd = frac_vpd_anom_e2 - frac_vpd_anom_e1) %>%
  # pull(dNDVI)
  ggplot(data=.,aes(x,y,fill=dNDVI))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  scale_fill_viridis_b(n.breaks=12, option='B', direction = -1)+
  # scale_fill_binned_diverging(
  #   palette='Berlin',rev=T,n.breaks=10, limits=c(-6,1),oob=scales::squish)+
  # scico::scale_fill_scico("%",
  #                         palette = 'bamako', 
  #                         direction = -1,
  #                         limits=c(5,15), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  labs(x=NULL,y=NULL, 
    title=expression(paste(No~CO[2]~VPD~Effect)), 
    fill=expression(paste(Delta*NDVI[GAM~Pred.],('%'))))+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()); p_gam_pred_noCO2_wVPD


p_gam_pred_noCO2_P <- bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
          tibble(pred_vpd_e2) %>% select(x,y,e,frac_vpd_anom)) %>% 
  inner_join(., tibble(co2=c(340.8,340.8),e=c('e1','e2')), 
             by=c('e')) %>% 
  mutate(frac_p_anom=0, 
         frac_ppet_anom=0,
         frac_pet_anom=0,
         # frac_vpd_anom=0,
         epoch=0) %>%
  inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  mutate(pred = predict(g3_full,newdata=.,type='response')) %>%   
  filter(is.na(pred)==F) %>% 
  select(x,y,e,frac_vpd_anom,pred,co2) %>% 
  pivot_wider(., names_from=c('e'), values_from=c('pred','co2','frac_vpd_anom')) %>% 
  mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
  mutate(d_vpd = frac_vpd_anom_e2 - frac_vpd_anom_e1) %>%
  # pull(dNDVI)
  ggplot(data=.,aes(x,y,fill=dNDVI))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  scale_fill_viridis_b(n.breaks=12, option='B', direction = -1)+
  # scale_fill_binned_diverging(
  #   palette='Berlin',rev=T,n.breaks=10, limits=c(-6,1),oob=scales::squish)+
  # scico::scale_fill_scico("%",
  #                         palette = 'bamako', 
  #                         direction = -1,
  #                         limits=c(5,15), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  labs(x=NULL,y=NULL, 
    title=expression(paste(No~CO[2]~VPD~Effect)), 
    fill=expression(paste(Delta*NDVI[GAM~Pred.],('%'))))+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()); p_gam_pred_noCO2_P





