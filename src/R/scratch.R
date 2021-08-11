p_vpd_sen <- lt_v_sen %>% 
  as_tibble() %>% 
  filter(between(b1,-0.1,0.1)) %>% 
  filter(b0 > 0) %>% 
  ggplot(data=.,aes(x,y,fill=100*38*b1/b0))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL,
    title=expression(paste(Delta,"VPD")), 
    fill="%")+
  # scico::scale_fill_scico(
  #                         palette ='romaO', direction=-1,
  #                         limits=c(-20,20), 
  #                         oob=scales::squish)+
  scale_fill_viridis_b('%',
    option='F', limits=c(0,20), oob=scales::squish, n.breaks=8)+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1)); p_vpd_sen


p_wue_sen <- lt_v_sen %>% 
  as_tibble() %>% 
  filter(b0 > 0) %>% 
  filter(between(b1,-0.05,0.05)) %>% 
  mutate(dVPD_VPD = b1*37/b0) %>% 
  mutate(expectation = 0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
  # pull(expectation) %>% quantile(., c(0.01,0.99))
  ggplot(data=.,aes(x,y,fill=expectation*100))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  scale_fill_binned_sequential(
    palette='ag_GrnYl',rev=T,n.breaks=10, limits=c(5,15),oob=scales::squish)+
  labs(x=NULL,y=NULL, 
    title=expression(paste(Delta*NDVI[WUE~Pred.])), 
    fill='%')+
  # scico::scale_fill_scico("%",
  #                         palette = 'bamako', 
  #                         direction = -1,
  #                         limits=c(5,15), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()); p_wue_sen



p_gam_pred <- bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
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
  select(x,y,e,frac_vpd_anom,pred,co2) %>% 
  pivot_wider(., names_from=c('e'), values_from=c('pred','co2','frac_vpd_anom')) %>% 
  mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
  mutate(d_vpd = frac_vpd_anom_e2 - frac_vpd_anom_e1) %>%
  ggplot(data=.,aes(x,y,fill=dNDVI))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  scale_fill_binned_sequential(
    palette='ag_GrnYl',rev=T,n.breaks=10, limits=c(5,15),oob=scales::squish)+
  # scico::scale_fill_scico("%",
  #                         palette = 'bamako', 
  #                         direction = -1,
  #                         limits=c(5,15), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  labs(x=NULL,y=NULL, 
    title=expression(paste(Delta*NDVI[GAM~Pred.])), 
    fill='%')+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()); p_gam_pred

p_res <- inner_join({lt_v_sen %>% 
  as_tibble() %>% 
  filter(b0 > 0) %>% 
  filter(between(b1,-0.05,0.05)) %>% 
  mutate(dVPD_VPD = b1*37/b0) %>% 
  mutate(expectation = 100*0.5*(dCa_Ca - 0.5*dVPD_VPD)) %>% 
    select(x,y,expectation)}, 
  {bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
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
      select(x,y,e,frac_vpd_anom,pred,co2) %>% 
      pivot_wider(., names_from=c('e'), values_from=c('pred','co2','frac_vpd_anom')) %>% 
      mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
      mutate(d_vpd = frac_vpd_anom_e2 - frac_vpd_anom_e1) %>% 
      select(x,y,dNDVI)}, by=c('x','y')) %>% 
 ggplot(data=.,aes(x,y,fill=expectation-dNDVI))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  scale_x_continuous(breaks=seq(140,154,by=5))+
  coord_sf(xlim = c(140,154),
           ylim = c(-45,-10), expand = FALSE)+
  labs(x=NULL,y=NULL,
    title=expression(paste(Delta~Pred[WUE]-Delta~Pred[GAM])), 
    fill='%')+
  scale_fill_binned_diverging(palette='Blue-Red',rev=F,n.breaks=10,
                              # cmax=90,
                              p1=1.5,p2=1.5)+
  # scico::scale_fill_scico(expression(paste(Delta*NDVI[WUE~Pred.]("%")-Delta*NDVI[GAM~Pred.]("%"))),
  #                         palette = 'roma', 
  #                         direction = -1,
  #                         limits=c(-15,15), #na.value = 'red',
  #                         oob=scales::squish
  # )+
  guides(fill=guide_colorbar())+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.title = element_text(size=8),
        legend.position = c(1,1), 
        legend.justification = c(1,1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()); p_res


vec_cols <- scico::scico(n=4, palette = 'bamako', 
                         direction = -1, end=0.9, begin=0.2)
p_box <- bind_rows(tibble(pred_vpd_e1) %>% select(x,y,e,frac_vpd_anom), 
          tibble(pred_vpd_e2) %>% select(x,y,e,frac_vpd_anom)) %>% 
  inner_join(., tibble(co2=c(340.8,411.9),e=c('e1','e2')), 
             by=c('e')) %>% 
  mutate(frac_p_anom=0, 
         frac_ppet_anom=0,
         frac_pet_anom=0,
         epoch=0) %>%
  inner_join(., tibble(unique(dat_annual[,.(x,y,map,mapet,mappet,mavpd15)])), by=c('x','y')) %>% 
  inner_join(., kop,by=c('x','y')) %>% 
  mutate(pred = predict(g3_full,newdata=.,type='response')) %>%   
  filter(is.na(pred)==F) %>% 
  select(x,y,e,frac_vpd_anom,pred,co2) %>% 
  pivot_wider(., names_from=c('e'), values_from=c('pred','co2','frac_vpd_anom')) %>% 
  mutate(dNDVI = 100*((pred_e2-pred_e1)/pred_e1)) %>% 
  inner_join(., {
    lt_v_sen %>% 
      as_tibble() %>% 
      filter(b0 > 0) %>% 
      filter(between(b1,-0.05,0.05)) %>% 
      mutate(dVPD_VPD = b1*37/b0) %>% 
      mutate(wue25 = 100*(0.25*(dCa_Ca - 0.5*dVPD_VPD)),
             wue50 = 100*(0.5*(dCa_Ca - 0.5*dVPD_VPD)),
             wue75 = 100*(0.75*(dCa_Ca - 0.5*dVPD_VPD)), 
             wue100 = 100*(1*(dCa_Ca - 0.5*dVPD_VPD))) %>% 
      select(x,y,wue25,wue50,wue75,wue100)
  },by=c("x","y")) %>% 
  inner_join(., kop, by=c('x','y')) %>% 
  mutate(res25 = dNDVI-wue25,
         res50 = dNDVI-wue50,
         res75 = dNDVI-wue75
  ) %>%
  select(zone,
         wue25,
         wue50,
         wue75,
         wue100,
         dNDVI
  ) %>%
  gather(-zone,key='allocation',
         value='res') %>% #ggplot(data=.,aes(value,y=zone,fill=key))+geom_boxplot()
  mutate(allocation = factor(allocation,
                             levels=c("wue25","wue50","wue75","wue100","dNDVI"),
                             labels = c('10%','50%','90%','100%','GAM. Estimate'), 
                             ordered = TRUE)) %>%
  ggplot(data=.,aes(y=res,
                    x=zone,
                    fill=allocation))+
  geom_hline(aes(yintercept=0),color='grey50',lty=3)+
  geom_boxplot(#draw_quantiles = c(0.25,0.5,0.75), 
    # trim=TRUE,
    outlier.colour = NA,
    color='black')+
  # scico::scale_fill_scico_d(
  #   palette='bamako',direction = -1, begin = 0.2,end=0.95,
  #   guide=guide_legend(title.position = 'top'))+
  scale_fill_manual(#"% WUE benefit allocated to foliar area:",
    values=c(vec_cols[1], vec_cols[2], vec_cols[3], vec_cols[4],"grey"))+
  labs(x=NULL, 
       y=expression(paste(Predicted~Delta*NDVI," (%)")), 
       fill=expression(paste(CO[2]*' x '*WUE,' gain allocated to foliar area')))+
  scale_x_discrete(limits=(structure(c(1L,2L,3L,4L,5L,6L,7L),# c(5L, 4L, 6L, 2L, 1L, 3L, 7L), 
                                     .Label = c("Equatorial",
                                                "Tropical", 
                                                "Subtropical", 
                                                "Grassland", 
                                                "Arid", 
                                                "Temperate",
                                                "Temperate Tas."), 
                                     class = c("ordered", "factor"))))+
  coord_cartesian(ylim = c(0,22))+
  # coord_flip()+
  theme_linedraw()+
  guides(fill=guide_legend(title.position='left'))+
  theme(legend.position = 'top', 
        legend.justification = c(1,1),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        # legend.key.height = unit(0.2,'cm'),
        axis.text = element_text(size=15),
        panel.grid = element_blank());
p_box

# library(cowplot)
# p_left <- ggdraw(p_vpd_sen)+draw_label(label='(a)', x=0.07,y=0.985,size = 25)
# p_mid1 <- ggdraw(p_gam_pred)+draw_label(label='(b)',x=0.05,y=0.985,size=25)
# p_mid2 <- ggdraw(p_wue_sen)+draw_label(label='(c)', x=0.05,y=0.985, size=25)
# p_right <- ggdraw(p_res)+draw_label(label = '(d)', x=0.05,y=0.985,size=25)
# p_bottom <- ggdraw(p_box)+draw_label(label = '(e)', x=0.015,y=0.95,size=25)
# cp_r <- cowplot::plot_grid(p_left,p_mid1,p_mid2,p_right,
#                            nrow = 1,
#                            rel_widths = c(1,1,1,1))
# cp_r
# ggsave(plot=cp_r/p_bottom+patchwork::plot_layout(heights = c(1,0.5)),
#        # filename='figures/delete.png',
#        filename = "figures/Fig6_map_dVpd_gamCO2Pred_wueCO2Pred_dDifferenceBoxplot.png",
#   width = 30,
#   height=30, 
#   units='cm', 
#   dpi=350,
#   type='cairo')

p_out <- (p_vpd_sen+theme(plot.margin=margin(t=0,r = -20,b=0,l=-10))|
    p_gam_pred+theme(plot.margin=margin(t=0,r = -20,b=0,l=-20))|
    p_wue_sen+theme(plot.margin=margin(t=0,r = -20,b=0,l=-20))|
    p_res+theme(plot.margin=margin(t=0,r = -10,b=0,l=-20)))/
  p_box+
  plot_annotation(tag_levels = 'a',
  tag_prefix = '(',
  tag_suffix=')', 
  theme=theme(plot.margin=margin(t=0,r = -10,b=0,l=-10))
    )+
  plot_layout(heights=c(1,0.4))&theme(plot.margin = margin(0,5,0,5))
ggsave(p_out,
  filename = "figures/Fig7_map_dVpd_gamCO2Pred_wueCO2Pred_dDifferenceBoxplot.png",
  device=grDevices::png,
  width = 32.5,
  height=35, 
  units='cm',
  dpi=350)

library(ggplot2)
p <- ggplot(data=data.frame(x=rnorm(10),y=rnorm(10)), aes(x,y))+
  geom_point()+
  labs(title=expression(paste(Delta~delta~'Delta')), 
    x=expression(paste(Delta~delta~'Delta'))); p
ggsave(p,
  filename = "figures/Fig6_map_dVpd_gamCO2Pred_wueCO2Pred_dDifferenceBoxplot.png", 
  device=grDevices::png)

grDevices::png
png(filename = "figures/Fig6_map_dVpd_gamCO2Pred_wueCO2Pred_dDifferenceBoxplot.png")
plot(rnorm(10),rnorm(10), 
  xlab=expression(paste(Delta~delta~'Delta')))
dev.off()

Cairo::CairoPNG(filename = "figures/Fig6_map_dVpd_gamCO2Pred_wueCO2Pred_dDifferenceBoxplot.png")
plot(rnorm(10),rnorm(10), 
  xlab=expression(paste(Delta~delta~'Delta')))
dev.off()



# Found pkg-config cflags and libs!
# Using PKG_CFLAGS=-I/usr/include/freetype2 -I/usr/include/libpng16 -I/usr/include/x86_64-linux-gnu
# Using PKG_LIBS=-lfreetype -lpng16 -lz -ltiff -ljpeg