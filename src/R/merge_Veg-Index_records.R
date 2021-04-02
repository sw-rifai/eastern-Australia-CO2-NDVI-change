library(stars); library(tidyverse); library(data.table); library(lubridate)
library(dtplyr, warn.conflicts = FALSE)
setDTthreads(threads=0)
library(mgcv)
library(RcppArmadillo)
#*******************************************************************************
kop <- arrow::read_parquet(file = '../data_general/Koppen_climate/BOM_Koppen_simplified7.parquet')


# Import NOAA AHVRR CDR --------------------------------------------------------
tmp <- stars::read_stars("../data_general/AVHRR_CDRv5_VI/AVHRR_CDRv5_ndvi_median_nobs_stdDev_5000m_EastOz_1981_2019.tif",
                         proxy = F)
tmp_median <- tmp %>% 
  slice('band', seq(1,by=3,length.out = 461)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("1981-08-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names(c('ndvi_cdr'))
tmp_nobs <- tmp %>% 
  slice('band', seq(2,by=3,length.out = 461)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("1981-08-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names(c('nobs'))
tmp_sd <- tmp %>% 
  slice('band', seq(3,by=3,length.out = 461)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("1981-08-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names(c('sd'))
tmp_sza <- stars::read_stars("../data_general/AVHRR_CDRv5_VI/AVHRR_CDRv5_sza_timeOfDay_median_5000m_EastOz_1981_2019.tif", 
                             proxy=F) %>% 
  slice('band', seq(1,by=2,length.out = 461)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("1981-08-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names("sza")
tmp_tod <- stars::read_stars("../data_general/AVHRR_CDRv5_VI/AVHRR_CDRv5_sza_timeOfDay_median_5000m_EastOz_1981_2019.tif", 
                             proxy=F) %>% 
  slice('band', seq(2,by=2,length.out = 461)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("1981-08-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names("tod")

tmp <- c(tmp_median, tmp_nobs, tmp_sd, tmp_sza, tmp_tod)
d_cdr <- tmp %>% as.data.frame() %>% as.data.table()
d_cdr <- merge(d_cdr, as.data.table(kop))
d_cdr <- d_cdr[is.na(zone)==F & is.na(nobs)==F]
d_cdr <- d_cdr %>% mutate(month=month(date))

rm(tmp_median, tmp_nobs, tmp_sd, tmp_sza, tmp_tod)
gc()

d_cdr_norms <- d_cdr %>% 
  lazy_dt() %>%
  filter(date <= ymd("2017-01-01")) %>% 
  mutate(month=month(date)) %>% 
  group_by(x,y,month) %>% 
  summarize(ndvi_u = mean(ndvi_cdr,na.rm=TRUE), 
            sza_u = mean(sza,na.rm=TRUE), 
            tod_u = mean(tod, na.rm=TRUE), 
            ndvi_sd = sd(ndvi_cdr, na.rm=TRUE), 
            sza_sd = sd(sza,na.rm=TRUE), 
            tod_sd = sd(tod,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
d_cdr <- merge(d_cdr_norms, d_cdr, by=c("x","y","month"))
d_cdr <- d_cdr %>% 
  lazy_dt() %>% 
  mutate(ndvi_anom = ndvi_cdr - ndvi_u, 
         sza_anom = sza - sza_u, 
         tod_anom = tod - tod_u) %>% 
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_sd, 
         sza_anom_sd = sza_anom/sza_sd, 
         tod_anom_sd = tod_anom/tod_sd) %>% 
  as.data.table()


# Filtering AVHRR CDR -----------------------------------------------------
d_cdr2 <- d_cdr %>% 
  lazy_dt() %>% 
  filter(nobs >= 3) %>% 
  filter(ndvi_cdr >= 0.1) %>% 
  mutate(cv_cdr = sd/ndvi_cdr) %>% 
  filter(cv_cdr < 0.25) %>% 
  filter(between(ndvi_anom_sd, -3.5,3.5)) %>% 
  filter(between(sza_anom_sd, -3.5,3.5)) %>% 
  filter(between(tod_anom_sd, -3.5,3.5)) %>% 
  as.data.table()
rm(d_cdr); gc()

# Import MCD43 w/fire & deforestation mask --------------------------------------------------------
tmp <- stars::read_stars("../data_general/MCD43/MCD43A4_ndvi_median_count_stdDev_5000m_EastOz_mMean_maskFireDefor_redRes_2001-01-01_to_2019-12-31.tif",
                         proxy = F) 
tmp_median <- tmp %>% 
  slice('band', seq(1,by=3,length.out = 228)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names(c('ndvi_mcd'))
tmp_nobs <- tmp %>% 
  slice('band', seq(2,by=3,length.out = 228)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names(c('nobs_mcd'))
tmp_sd <- tmp %>% 
  slice('band', seq(3,by=3,length.out = 228)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),ymd("2019-12-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names(c('sd_mcd'))
tmp <- c(tmp_median, tmp_nobs, tmp_sd)
d_mcd <- tmp %>% as.data.frame() %>% as.data.table()

# Filter masked MCD43
d_mcd <- merge(d_mcd, as.data.table(kop))
d_mcd <- d_mcd[is.na(zone)==F & is.na(nobs_mcd)==F]
d_mcd <- d_mcd[nobs_mcd>=3]
d_mcd <- d_mcd[sd_mcd/ndvi_mcd <= 0.25]
d_mcd <- d_mcd[,`:=`(month=month(date))]


# Import MCD43 no mask --------------------------------------------------------
tmp_nm <- stars::read_stars("../data_general/MCD43/MCD43A4_NDVI_5000m_EastOz_mMedian_noMasking_2001-01-01_to_2020-07-30.tif",
                            proxy = F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),ymd("2020-07-01"),by="1 month"), 
                    names = 'date') %>% 
  set_names(c('ndvi_mcd_nm'))
d_nm <- as_tibble(tmp_nm) %>% as.data.table()
rm(tmp_nm); gc()


# Calibrate AVHRR CDR to approximate MCD43 NDVI -----------------------------------------
# This is being done to counteract the effects of changes in solar zenith angle 
# and time of day.
dc2 <- merge(d_mcd[,.(x,y,date,ndvi_mcd)],
             d_cdr2[,.(x,y,date,ndvi_cdr,sza,tod)], all=TRUE, 
             by=c("x","y","date"))
dc2 <- merge(dc2, d_nm[,.(x,y,date,ndvi_mcd_nm)], all=TRUE)
dc2 <- dc2[,`:=`(month=month(date))]
dc2_train <- dc2[is.na(ndvi_mcd)==F][between(date,ymd("2001-01-01"),ymd("2016-12-31"))==T][sample(.N, 1e6)]
dc2_test <- dc2[is.na(ndvi_mcd)==F][between(date,ymd("2001-01-01"),ymd("2016-12-31"))==T][sample(.N, 1e6)]

mc10 <- bam(ndvi_mcd ~ s(ndvi_cdr,bs='cs')+
              s(month,bs='cs')+
              s(sza,bs='cs')+
              s(tod,bs='cs')+
              te(x,y), data=dc2_train[ndvi_mcd_nm>0.1], 
            discrete=T, 
            select=T)
summary(mc10)
summary(predict(mc10, newdata=dc2_test[is.na(ndvi_cdr)==F]))

# dc3 <- merge(as.data.table(d_nm)[,.(x,y,date,ndvi_mcd_nm)],
#              d_cdr2[,.(x,y,date,ndvi_cdr,sza,tod)], all=TRUE, 
#              by=c("x","y","date"))
# dc3 <- dc3[,`:=`(month=month(date))]
# dc3_train <- dc3[is.na(ndvi_mcd_nm)==F][between(date,ymd("2001-01-01"),ymd("2016-12-31"))==T][sample(.N, 1e6)]
# dc3_test <- dc3[is.na(ndvi_mcd_nm)==F][between(date,ymd("2001-01-01"),ymd("2016-12-31"))==T][sample(.N, 1e6)]


mc10_nm <- bam(ndvi_mcd_nm ~ s(ndvi_cdr,bs='cs')+
                 s(month,bs='cs')+
                 s(sza,bs='cs')+
                 s(tod,bs='cs')+
                 te(x,y), data=dc2_train[ndvi_mcd_nm>0.1], 
               discrete=T, 
               select=T)
summary(mc10_nm)

# junk <- left_join(as_tibble(tmp_nm),d_mcd, by=c("x","y","date")) 
# junk %>% select(ndvi_mcd, ndvi_mcd_nm) %>% summary


# Apply Calibration prediction in parallel --------------------------------
library(foreach); library(doParallel)
n_cores <- 4
cl <- makeCluster(n_cores)
registerDoParallel(cl)
gc()
dc2 <- dc2 %>% mutate(year=year(date), month=month(date))
# vec_years <- unique(dc2$year) %>% na.omit() %>% sort()
vec_years <- 1981:2019

out <- foreach(i = 1:length(vec_years), 
               .packages = c("mgcv","data.table","tidyverse"),
               .combine=rbind) %dopar% {
                 out <- dc2[year==vec_years[i]][is.na(ndvi_cdr)==F] %>% 
                   mutate(ndvi_mcd_pred = predict(mc10, 
                                                  newdata=., 
                                                  newdata.guaranteed = T,
                                                  discrete = TRUE), 
                          ndvi_mcd_nm_pred = predict(mc10_nm, 
                                                  newdata=., 
                                                  newdata.guaranteed = T,
                                                  discrete = TRUE))
                 out
               } 


# Join datasets -----------------------------------------------------------

coords_keep <- as.data.table(kop)[,.(x,y)]
d_nm <- as.data.table(d_nm)[,.(x,y,date,ndvi_mcd_nm)][is.na(ndvi_mcd_nm)==F]
d_nm <- merge(d_nm,coords_keep,by=c("x","y"))

d_export <- merge(d_mcd[,.(x,y,date,ndvi_mcd,nobs_mcd,sd_mcd,month)], 
      d_nm, by=c("x","y","date"), all=TRUE)
d_export <- merge(d_export,
  out[, .(x, y, date, ndvi_cdr, sza, tod, year, ndvi_mcd_pred, 
     ndvi_mcd_nm_pred)], by=c("x",'y','date'),all=TRUE)

d_export <- d_export %>% 
  mutate(month=month(date), 
       season = case_when(month%in%c(1,2,12)~'DJF',
                          month%in%c(3,4,5)~'MAM',
                          month%in%c(6,7,8)~'JJA',
                          month%in%c(9,10,11)~'SON')) %>% 
  mutate(season = factor(season,levels = c("SON","DJF","MAM","JJA"), 
                         ordered=TRUE)) %>% 
  .[,`:=`(year=year(date))]

# Filter bad values
d_export <- d_export %>% 
  lazy_dt() %>% 
  mutate(ndvi_mcd_pred = ifelse(between(ndvi_mcd_pred,0,1),ndvi_mcd_pred, NA), 
         ndvi_mcd_nm_pred = ifelse(between(ndvi_mcd_nm_pred,0,1),ndvi_mcd_nm_pred, NA), 
         ndvi_mcd = ifelse(between(ndvi_mcd,0,1),ndvi_mcd,NA)) %>% 
  as.data.table()

arrow::write_parquet(d_export, 
       sink =paste0("../data_general/MCD43/MCD43_AVHRR_NDVI_hybrid_",Sys.Date(),".parquet"), 
         compression='snappy')


# Comparison figures ------------------------------------------------------
d_export <- merge(d_export, as.data.table(kop %>% select(x,y,zone)), by=c('x','y'))
coords_keep <- d_export[is.na(ndvi_mcd)==F] %>% 
  lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(nobs_mcd_2 = sum(is.na(ndvi_mcd)==F)) %>% 
  ungroup() %>% 
  as.data.table()


# Plot Annual Intercept & Trend Differences --------------------------------------
creg1 <- merge(coords_keep, d_export, by=c('x','y')) %>% 
  .[date>= ymd("2001-01-01") & date<= ymd("2016-12-31")] %>% 
  .[nobs_mcd_2 >= 100] %>% 
  .[,`:=`(year=year(date))] %>% 
  .[, .(val_cdr=mean(ndvi_cdr,na.rm=TRUE), 
        val_mcd=mean(ndvi_mcd,na.rm=TRUE), 
        val_mcd_pred = mean(ndvi_mcd_pred,na.rm=TRUE), 
        val_mcd_nm = mean(ndvi_mcd_nm,na.rm=TRUE)), 
    keyby=.(x,y,year)] %>% 
  .[,.(beta_cdr = list(unname(fastLm(
    X = cbind(1,year-2001), 
    y=val_cdr, data=.SD)$coefficients)), 
    beta_mcd = list(unname(fastLm(
      X = cbind(1,year-2001), 
      y=val_mcd, data=.SD)$coefficients)), 
    beta_cdr_pred = list(unname(fastLm(
      X = cbind(1,year-2001), 
      y=val_mcd_pred, data=.SD)$coefficients)), 
    beta_mcd_nm = list(unname(fastLm(
      X = cbind(1,year-2001), 
      y=val_mcd_nm, data=.SD)$coefficients))), 
    by=.(x,y)] %>% 
  .[,`:=`(b0_cdr=unlist(beta_cdr)[1], b1_cdr=unlist(beta_cdr)[2], 
          b0_mcd=unlist(beta_mcd)[1], b1_mcd=unlist(beta_mcd)[2],
          b0_mcd_nm=unlist(beta_mcd_nm)[1], b1_mcd_nm=unlist(beta_mcd_nm)[2],
          b0_mcd_pred=unlist(beta_cdr_pred)[1], b1_mcd_pred=unlist(beta_cdr_pred)[2]), by=.(x,y)]

creg1 %>% 
  select(b0_mcd, b0_mcd_pred, b0_cdr, b0_mcd_nm) %>% 
  gather(-b0_mcd, key='key',value='value') %>% 
  filter(b0_mcd >= 0.1) %>% 
  ggplot(data=.,aes(b0_mcd, value, color=key))+
  geom_point(size=0.25,alpha=0.1)+
  geom_abline(aes(intercept=0,slope=1),col='gray',lwd=2)+
  geom_smooth(method='lm')+
  labs(x='Intercept MCD NDVI', 
       y='Intercept', 
       title='2001-2016')+
  scale_color_viridis_d("",
                        option='B',end=0.8, 
                        labels=c("b0_cdr"="Intercept AVHRR CDR", 
                                 "b0_mcd_pred"="Intercept AVHRR CDR calibrated to MCD", 
                                 "b0_mcd_nm"="Intercept MCD43 w/no masking"))+
  scale_y_continuous(limits=c(0,1))+
  theme_linedraw()+
  theme(legend.position = c(0.99,0.01), 
        legend.justification = c(0.99,0.01))
ggsave(filename = "figures/compare_lm_intercept_mcd_avhrr_recalibAvhrr.png", 
       width=15, height=12, units='cm')

creg1 %>% 
  select(b1_mcd, b1_mcd_pred, b1_cdr, b1_mcd_nm) %>% 
  gather(-b1_mcd, key='key',value='value') %>% 
  filter(between(b1_mcd,-0.01,0.01)) %>% 
  filter(between(value,-0.01,0.01)) %>% 
  ggplot(data=.,aes(b1_mcd, value, color=key))+
  geom_point(size=0.25,alpha=0.1)+
  geom_abline(aes(intercept=0,slope=1),col='gray',lwd=2)+
  geom_smooth(method='lm')+
  labs(x='Intercept MCD NDVI', 
       y='Intercept', 
       title='2001-2016')+
  scale_color_viridis_d("",
                        option='B',end=0.8, 
                        labels=c("b1_cdr"="Intercept AVHRR CDR", 
                                 "b1_mcd_pred"="Intercept AVHRR CDR calibrated to MCD", 
                                 "b1_mcd_nm"="Intercept MCD43 w/no masking"))+
  theme_linedraw()+
  theme(legend.position = c(0.99,0.01), 
        legend.justification = c(0.99,0.01))
ggsave(filename = "figures/compare_lm_trend_mcd_avhrr_recalibAvhrr.png", 
       width=15, height=12, units='cm')


creg1
