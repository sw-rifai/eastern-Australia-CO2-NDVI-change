library(raster); library(rasterVis)
library(arrow)
library(sf); library(stars)
library(tidyverse); 
library(data.table); library(lubridate);
library(dtplyr)
setDTthreads(threads=8)


#*******************************************************************************
# Get NVIS Coords ---------------------------------------------------------
#*******************************************************************************
# base <- stars::read_stars("../data_general/AVHRR_CDRv5_VI/AVHRR_NIRV_monmean_EastOz_1982_2019.tif",
#                           RasterIO = list(bands=1))
nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif") %>% 
  set_names('veg_class')

# nvis2 <- st_warp(src=nvis, dest=base[,,], use_gdal = T)
# names(nvis2) <- "veg_class"
# nvis <- nvis2 %>% as_tibble() %>% as.data.table()
codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt",
                         fwf_widths(c(2,100)), skip = 1) %>%
  set_names(c("veg_class","veg_class_descrip")) %>%
  mutate(vc = as.factor(veg_class_descrip))
# nvis <- inner_join(nvis, codes, by='veg_class')
# nvis <- nvis %>% filter(veg_class <= 15) # !!! only forests and woodlands !!!
# rm(base); #
base <- arrow::read_parquet("../data_general/MCD43/MCD43_AVHRR_NDVI_hybrid_2020-10-12.parquet") %>% 
  as.data.table()
base_coords <- unique(base[,.(x,y)])
base_coords <- st_as_sf(base_coords,crs=st_crs(4326),coords=c("x","y"))
nvis <- st_extract(nvis, pts = base_coords)
nvis <- bind_cols(nvis %>% st_coordinates() %>% as.data.table,nvis$veg_class) %>% set_names(c("x","y","veg_class"))
nvis <- merge(nvis,codes,by='veg_class')
nvis <- nvis[veg_class<=15]
# nvis <- base %>% 
#   group_by(x,y) %>% 
#   summarize(veg_class = first(veg_class), 
#             vc = first(vc)) %>% 
#   ungroup() %>% 
#   as.data.table()
rm(base); gc()
rm(base_coords)
#*******************************************************************************
#* END SECTION
#*******************************************************************************

#*******************************************************************************
# Extract Tmax Threshold Frequencies
#*******************************************************************************
attmax <- c(
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt35.nc") %>% 
    set_names("t35"), 
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt36.nc") %>% 
    set_names("t36"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt37.nc") %>% 
    set_names("t37"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt38.nc") %>% 
    set_names("t38"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt39.nc") %>% 
    set_names("t39"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt40.nc") %>% 
    set_names("t40"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt41.nc") %>% 
    set_names("t41"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt42.nc") %>% 
    set_names("t42"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt43.nc") %>% 
    set_names("t43"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt44.nc") %>% 
    set_names("t44"),
  stars::read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmax_gt45.nc") %>% 
    set_names("t45")
)

attmax <- as_tibble(attmax) %>% as.data.table()
attmax <- attmax[is.na(t35)==F]

coords_awap <- attmax %>% select(lon,lat) %>% distinct()
coords_awap <- coords_awap %>% rename(x=lon,y=lat)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))

coords_vi <- nvis %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])

# vis check that vi and clim coords are close
coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
coords_dict %>% head
coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

attmax <- attmax %>% rename(x=lon,y=lat)

# complicated way of doing full join
attmax <- merge(attmax,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
attmax <- attmax[is.na(x_vi)==F]
attmax %>% head

round(attmax$x[1:5])==round(attmax$x_vi[1:5])
round(attmax$y[1:5])==round(attmax$y_vi[1:5])

# visual check
attmax[time%in%c(ymd("1990-01-01",tz='UTC'),
                ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=t41))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))
gc()

#*******************************************************************************
#* END SECTION
#*******************************************************************************

#*******************************************************************************
# Extract AWAP tmax grid cells for east Oz ----------------------------------------------
#*******************************************************************************
atmax <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/tmax/AWAP_monthly_tmax_1970_2019.nc")
names(atmax) <- "tmax"
st_crs(atmax) <- st_crs(4326)

eoz_box <- st_bbox(c(xmin = min(nvis$x),
                     ymin = min(nvis$y),
                     xmax = max(nvis$x),
                     ymax = max(nvis$y)), 
                   crs = st_crs(4326))
atmax <- st_crop(atmax, eoz_box)
atmax <- atmax %>% as_tibble() %>% as.data.table()
atmax <- atmax %>% units::drop_units()
coords_awap <- atmax %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))

coords_vi <- nvis %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])



# vis check that vi and clim coords are close
coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
coords_dict %>% head
coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

atmax <- atmax %>% rename(x=longitude,y=latitude)

# complicated way of doing full join
atmax <- merge(atmax,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
atmax <- atmax[is.na(x_vi)==F]
atmax %>% head

round(atmax$x[1:5])==round(atmax$x_vi[1:5])
round(atmax$y[1:5])==round(atmax$y_vi[1:5])

# visual check
atmax[time%in%c(ymd("1990-01-01",tz='UTC'),
                ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=tmax))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))
gc()
#*******************************************************************************
# END SECTION
#*******************************************************************************

#*******************************************************************************
# Extract AWAP tmin grid cells for east Oz -------------------------------------
#*******************************************************************************
atmin <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/tmin/AWAP_monthly_tmin_1970_2019.nc")
names(atmin) <- "tmin"
st_crs(atmin) <- st_crs(4326)

eoz_box <- st_bbox(c(xmin = min(nvis$x),
                     ymin = min(nvis$y),
                     xmax = max(nvis$x),
                     ymax = max(nvis$y)), 
                   crs = st_crs(4326))
atmin <- st_crop(atmin, eoz_box)
atmin <- atmin %>% as_tibble() %>% as.data.table()
atmin <- atmin %>% units::drop_units()
coords_awap <- atmin %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))

coords_vi <- nvis %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])



# vis check that vi and clim coords are close
# coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
# coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
# coords_dict %>% head
coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

atmin <- atmin %>% rename(x=longitude,y=latitude)

# complicated way of doing full join
atmin <- merge(atmin,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
atmin <- atmin[is.na(x_vi)==F]
atmin %>% head

round(atmin$x[1:5])==round(atmin$x_vi[1:5])
round(atmin$y[1:5])==round(atmin$y_vi[1:5])

# visual check
atmin[time%in%c(ymd("1990-01-01",tz='UTC'),
                ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=tmin))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))
gc()
#*******************************************************************************
# END SECTION
#*******************************************************************************

#*******************************************************************************
# Extract AWAP VP3pm (-> VPD) grid cells for east Oz ----------------------------------------------
#*******************************************************************************
#' Calculates saturation vapour pressure
#' @return saturation vapour pressure
calc_esat <- function(airtemp){
  #Tair in degrees C
  
  #From Jones (1992), Plants and microclimate: A quantitative approach 
  #to environmental plant physiology, p110
  esat <- 613.75 * exp(17.502 * airtemp / (240.97+airtemp))
  
  return(esat)
}
avp9 <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/vph09/AWAP_monthly_vph09_1970_2019.nc")
avp15 <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/vph15/AWAP_monthly_vph15_1970_2019.nc")
names(avp9) <- "vp9"
names(avp15) <- "vp15"
st_crs(avp9) <- st_crs(4326)
st_crs(avp15) <- st_crs(4326)

avp <- c(avp9,avp15)


eoz_box <- st_bbox(c(xmin = min(nvis$x),
                     ymin = min(nvis$y),
                     xmax = max(nvis$x),
                     ymax = max(nvis$y)), 
                   crs = st_crs(4326))
avp <- st_crop(avp, eoz_box)
avp <- avp %>% as_tibble() %>% as.data.table() %>% units::drop_units()
coords_awap <- avp %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))

coords_vi <- unique(nvis[,.(x,y)])
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])

coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!

avp <- avp %>% rename(x=longitude,y=latitude)

# complicated way of doing full join
avp <- merge(avp,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
avp <- avp[is.na(x_vi)==F]


avp <- merge(avp, atmax, by=c("x","y","x_vi","y_vi","time"))

avp <- avp %>% lazy_dt() %>% 
  mutate(vpd15 = 0.01*(calc_esat(tmax)/10 - vp15)) %>% 
  as.data.table()

#*******************************************************************************
#* END SECTION
#*******************************************************************************
 

#*******************************************************************************
# Extract AWAP precip grid cells for east Oz ----------------------------------------------
#*******************************************************************************
aprecip <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/rain/AWAP_monthly_rain_1970_2019.nc")
names(aprecip) <- "precip"
st_crs(aprecip) <- st_crs(4326)

eoz_box <- st_bbox(c(xmin = min(nvis$x),
                     ymin = min(nvis$y),
                     xmax = max(nvis$x),
                     ymax = max(nvis$y)), 
                   crs = st_crs(4326))
aprecip <- st_crop(aprecip, eoz_box)
aprecip <- aprecip %>% as_tibble() %>% as.data.table()
aprecip <- aprecip %>% units::drop_units()
coords_awap <- aprecip %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))

coords_vi <- nvis %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])



# vis check that vi and clim coords are close
coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
coords_dict %>% head
coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

aprecip <- aprecip %>% rename(x=longitude,y=latitude)

# complicated way of doing full join
aprecip <- merge(aprecip,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
aprecip <- aprecip[is.na(x_vi)==F]
aprecip %>% head

round(aprecip$x[1:5])==round(aprecip$x_vi[1:5])
round(aprecip$y[1:5])==round(aprecip$y_vi[1:5])

# visual check
aprecip[time%in%c(ymd("1990-01-01",tz='UTC'),
               ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=precip))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))
gc()
#*******************************************************************************
# END SECTION
#*******************************************************************************

#*******************************************************************************
# Extract AWAP pet grid cells for east Oz ----------------------------------------------
#*******************************************************************************
apet <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/pet/AWAP_monthly_PriestleyTaylor_PET_1990_2019.nc")
names(apet) <- "pet"
st_crs(apet) <- st_crs(4326)

eoz_box <- st_bbox(c(xmin = min(nvis$x),
                     ymin = min(nvis$y),
                     xmax = max(nvis$x),
                     ymax = max(nvis$y)), 
                   crs = st_crs(4326))
apet <- st_crop(apet, eoz_box)
apet <- apet %>% as_tibble() %>% as.data.table()
apet <- apet %>% units::drop_units()
coords_awap <- apet %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))

coords_vi <- nvis %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])



# vis check that vi and clim coords are close
coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
coords_dict %>% head
coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

apet <- apet %>% rename(x=longitude,y=latitude)

# complicated way of doing full join
apet <- merge(apet,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
apet <- apet[is.na(x_vi)==F]
apet %>% head

round(apet$x[1:5])==round(apet$x_vi[1:5])
round(apet$y[1:5])==round(apet$y_vi[1:5])

# visual check
apet[time%in%c(ymd("1990-01-01",tz='UTC'),
              ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=pet))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))

#*******************************************************************************
# END SECTION
#*******************************************************************************





#*******************************************************************************
# Extract ERA5 pet grid cells for east Oz ----------------------------------------------
#*******************************************************************************
e5pet <- stars::read_ncdf("../data_general/clim_grid/era5-land/Oz/Oz/Oz_era5-land_pet_1981_2019.nc")
names(e5pet) <- "pet"
st_crs(e5pet) <- st_crs(4326)

eoz_box <- st_bbox(c(xmin = min(nvis$x),
                     ymin = min(nvis$y),
                     xmax = max(nvis$x),
                     ymax = max(nvis$y)), 
                   crs = st_crs(4326))
e5pet <- st_crop(e5pet, eoz_box)
e5pet <- e5pet %>% as_tibble() %>% as.data.table()
e5pet <- e5pet %>% units::drop_units()
coords_awap <- e5pet %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))

coords_vi <- nvis %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])



# vis check that vi and clim coords are close
coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
coords_dict %>% head
coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

e5pet <- e5pet %>% rename(x=longitude,y=latitude)

# complicated way of doing full join
e5pet <- merge(e5pet,coords_dict,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
e5pet <- e5pet[is.na(x_vi)==F]
e5pet %>% head

round(e5pet$x[1:5])==round(e5pet$x_vi[1:5])
round(e5pet$y[1:5])==round(e5pet$y_vi[1:5])

# visual check
e5pet[time%in%c(ymd("1990-01-01",tz='UTC'),
                ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=pet))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))


# ATTEMPTING TO RECALIBRATE ERA5-LAND PET TO AWAP PET
tmp_apet <- apet %>% select(x_vi,y_vi,time,pet) %>% 
  mutate(time=as.Date(time))
tmp_e5pet <- e5pet %>% select(x_vi,y_vi,time,pet) %>% 
  mutate(time=as.Date(time))
tmp_pet <- tmp_apet[tmp_e5pet,on=.(x_vi,y_vi,time)]

system.time(
fit_pet <- tmp_pet %>% 
  filter(is.na(pet)==F) %>% 
  filter(is.na(i.pet)==F) %>% 
  sample_frac(0.35) %>% 
  group_by(x_vi,y_vi) %>% 
  summarize(beta0 = coef(lm(pet~i.pet+I(i.pet**2)))[1], 
            beta1 = coef(lm(pet~i.pet+I(i.pet**2)))[2],
            beta2 = coef(lm(pet~i.pet+I(i.pet**2)))[3]) %>% 
  ungroup()
)

setDT(fit_pet)
tmp_pet <- fit_pet[tmp_pet,on=.(x_vi,y_vi)]
tmp_pet <- tmp_pet[, `:=`(pet_pred = beta0 + beta1*i.pet + beta2*i.pet**2)]
  
# Visual Checks ***************
# tmp_pet %>% 
#   filter(is.na(pet_pred)==F) %>% 
#   sample_frac(0.005) %>% 
#   ggplot(data=., aes(pet_pred,pet))+
#   ggpointdensity::geom_pointdensity(alpha=0.1)+
#   scale_color_viridis_c()+
#   geom_smooth(method='lm')+
#   geom_abline(aes(intercept=0,slope=1),color='red')

# tmp_pet %>% 
#   filter(time==ymd("1995-11-01")) %>% 
#   ggplot(data=., aes(x_vi,y_vi,fill=pet_pred-pet))+
#   geom_tile()+
#   coord_equal()+
#   scale_fill_gradient2()
gc()
#*******************************************************************************
# END SECTION
#*******************************************************************************

#*******************************************************************************
# Join climate data to maintain date structure ---------
#*******************************************************************************
tmp_pet[time>=ymd("1981-01-01") & time<=ymd("1989-12-31")] %>% 
  .[,.(x_vi,y_vi,time,pet_pred)] %>% 
  .[,.(x_vi, y_vi,time,pet=pet_pred)] %>% 
  .[is.na(pet)==F]

jpet <- rbindlist(list(
  tmp_pet[time>=ymd("1981-01-01") & time<=ymd("1989-12-31")] %>% 
    .[,.(x_vi,y_vi,time,pet_pred)] %>% 
    .[,.(x_vi, y_vi,time,pet=pet_pred)],
  apet[,.(x_vi,y_vi,time,pet)][,`:=`(time=as.Date(time))]
))

jpet <- jpet[is.na(pet)==F]
gc(reset = T, full = T)

# atmax <- atmax[,.(x_vi,y_vi,time,tmax)][,`:=`(time=as.Date(time))]
avp <- avp[,.(x_vi,y_vi,time,tmax,vp9,vp15,vpd15)][,`:=`(time=as.Date(time))]
attmax <- attmax[,`:=`(time=as.Date(time))]
atmin <- atmin[,`:=`(time=as.Date(time))]
aprecip <- aprecip[,.(x_vi,y_vi,time,precip)][,`:=`(time=as.Date(time))]
tmp_clim <- merge(avp,jpet,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)
tmp_clim <- tmp_clim[is.na(pet)==F][order(time,x_vi,y_vi)]
tmp_clim <- merge(tmp_clim,aprecip,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)
tmp_clim <- merge(tmp_clim,attmax,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)
tmp_clim <- merge(tmp_clim,atmin,by=c("x_vi","y_vi","time"),all=TRUE,allow.cartesian=TRUE)

gc(reset = T, full=T)
#*******************************************************************************
# END SECTION
#*******************************************************************************


#*******************************************************************************
# NVIS and NDVI (forest & woodlands) ------------------------------------------------
#*******************************************************************************
# base <- arrow::read_parquet("../data_general/MCD43/MCD64_AVHRR_NDVI_hybrid_2020-05-18.parquet") %>% 
#   as.data.table()


# base <- base %>% lazy_dt() %>% 
#   mutate(x = round(x, 2) ,
#          y = round(y, 2)) %>% 
#   as.data.table()
# names(base) <- "nirv"
# vec_dates <- seq(ymd("1982-01-01"),ymd("2019-12-01"),by="1 month")
# base <- st_set_dimensions(base, 3, values=vec_dates, names='date')
# nvis <- stars::read_stars("../data_general/NVIS/nvis51_majorVegClass_0p05.tif")
# nvis2 <- st_warp(src=nvis, dest=base[,,,1], use_gdal = T)
# names(nvis2) <- "veg_class"
# nvis <- nvis2 %>% as_tibble() %>% as.data.table()
# codes <- readr::read_fwf("../data_general/NVIS/nvis51_majorVegClass_codes.txt", 
#                          fwf_widths(c(2,100)), skip = 1) %>% 
#   set_names(c("veg_class","veg_class_descrip")) %>% 
#   mutate(vc = as.factor(veg_class_descrip))
# nvis <- inner_join(nvis, codes, by='veg_class')
# nvis <- nvis %>% filter(veg_class <= 15) # !!! only forests and woodlands !!!
# nvis <- nvis %>% select(-veg_class_descrip)
# base <- base %>% as_tibble() %>% as.data.table()
# base_vc <- base <- base[nvis,on=.(x,y)]
# base <- base_vc; rm(base_vc)

#*******************************************************************************
# END SECTION
#*******************************************************************************


#*******************************************************************************
# JOIN ALL THE PIECES -----------------------------------------------------
#*******************************************************************************
tmp_clim <- tmp_clim %>% rename(date=time)
# base <- base %>% rename(x_vi=x,y_vi=y)

# tmp_clim <- merge(tmp_clim, 
#               base,
#              by=c("x_vi","y_vi","date"), 
#              all=TRUE,allow.cartesian=TRUE)
tmp_clim %>% head
tmp_clim <- tmp_clim[is.na(pet)==F][order(date,x_vi,y_vi)]

arrow::write_parquet(tmp_clim, 
                     sink=paste0("../data_general/Oz_misc_data/ARD_ndvi_aclim_",Sys.Date(),".parquet"),
                     compression='gzip',compression_level = 9)
#*******************************************************************************
#*
#*******************************************************************************
rm(fit_pet, tmp_pet, tmp_apet, apet);
rm(aprecip, atmax,e5pet,jpet);
rm(tmp_e5pet)
gc(reset = T, full = T)
