devtools::load_all(here::here())
devtools::install()
library(weatheRcues)
library(tidyverse)
library(sf)
library(ggsci)

################################
#example with all picea bies 
#subset data of interest here 
Picea.seed = initial.data.mastree %>%
  filter(!Variable == "flower" & VarType == "C" & !Unit=="index" & !Variable == "pollen") %>% 
  group_by(Species, VarType) %>% 
  mutate(nMax = n()) %>% 
  filter(str_detect(Species, 'Picea abies')) %>% 
  filter(!str_detect(Country, 'United States of America|Russian Federation (the)')) %>% 
  filter(Year > 1952 & Year < 2023) %>% 
  mutate(sitenewname= paste0(Alpha_Number, "_",Site_number, "_", Species_code)) %>% 
  group_by(sitenewname) %>% 
  mutate(log.seed = log(1+Value),
         scale.seed = scale(Value),
         n = n()) %>% 
  filter(n > 19) %>%
  mutate() %>% 
  mutate(Date = paste0( "15/06/",Year)) %>% 
  mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
         plotname.lon.lat = paste0("longitude=",Longitude, "_","latitude=", Latitude)) %>% 
  ungroup() %>% 
  as.data.frame()

#df for collection method
methods.collection.picea = Picea.seed %>% dplyr::select(sitenewname, Country, Collection_method, Length, Longitude, Latitude) %>% distinct()

#for after 
picea.site.all = unique(Picea.seed$plotname.lon.lat)
climate.picea.path <- here::here('climate_dailyEOBS/')


########################################################################################################################
########################################################################################################################
########################################################################################################################
#TEST with climwin
########################################################################################################################
run.climwin <- T
#take some time ! alternative it will load the outputs 
if(run.climwin==T){
  #for parallel to take less time usually ... 
  library(furrr)
  plan(multisession)
  devtools::load_all(here::here())
  #for package in progress
  #
  #devtools::install()
  library(weatheRcues)
  library(dplyr)
  #formulanull = as.formula('log.seed ~ 1')
  
  #run loop for each sites 
  statistics_absolute_climwin_picea <- future_map_dfr(
    picea.site.all, 
    ~ {
      
      #format climate each site 
      climate_data <- format_climate_data(site = .x ,
                                          path = climate.picea.path, 
                                          scale.climate = T) 
      #formulanull = as.formula('log.seed ~ 1')
      #print(paste("Processing site:", .x, "with", nrow(bio_data), "rows"))
      #run climwin each sites
      climwin_site_days(
        climate_data = climate_data,
        data = data.frame(Picea.seed %>% dplyr::filter(plotname.lon.lat == .x) %>% tidyr::drop_na("log.seed")),  
        site.name = .x, 
        #formulanull = formulanull,
        range = c(600, 0),
        cinterval = 'day',
        refday = c(01, 11),  
        optionwindows = 'absolute',  
        climate_var = 'TMEAN'  
      )
    }#,
    #.options = furrr_options(seed = TRUE, globals = c('formulanull', "climate.picea.path", "Picea.seed"),
    #                         packages = c("weatheRcues", "dplyr", "climwin", "broom"))
  )
  
  qs::qsave(statistics_absolute_climwin_picea, 
            here::here('outputs/statistics_absolute_climwin_PICEA.qs'))
}else{statistics_absolute_climwin_quercus = qs::qread('outputs/statistics_absolute_climwin_PICEA.qs')}


########################################################################################################################
########################################################################################################################
########################################################################################################################
#TEST with CSP
########################################################################################################################

results.moving.site.picea = FULL.moving.climate.analysis(seed.data.all = Picea.seed,
                                                         climate.path = climate.picea.path,
                                                         refday = 305,
                                                         lastdays = 600,
                                                         rollwin = 1)
Results_daily.picea = results.moving.site.picea %>%
  arrange(sitenewname) %>% 
  group_by(sitenewname) %>%
  group_split()

name_daily.picea =   results.moving.site.picea %>% 
  arrange(sitenewname) %>% 
  group_by(sitenewname) %>%
  group_keys()

# Set names of the list based on group keys
names(Results_daily.picea) <- apply(name_daily.picea, 1, paste)

statistics_csp_method.picea = map_dfr(
  1:length(Results_daily.picea), 
  ~CSP_function_site(
    Results_daily.picea[[.]], 
    unique(Results_daily.picea[[.]]$plotname.lon.lat),
    seed.data = Picea.seed %>% arrange(sitenewname), 
    climate.path = climate.picea.path,
    refday = 305,
    lastdays = 600,
    rollwin = 1
  ))


########################################################################################################################
########################################################################################################################
########################################################################################################################
#Final figures 
########################################################################################################################

#make map 
plot_locations.picea<- st_as_sf(Picea.seed %>% 
                                  dplyr::select(Longitude, Latitude, plotname.lon.lat, Collection_method) %>% distinct(), coords = c("Longitude", "Latitude")) %>% 
  st_set_crs(4326)

library("rworldmap")
mapBase <- getMap(resolution = "high") %>%   
  st_as_sf() %>% 
  st_make_valid()%>%  st_crop(mapBase, xmin = -12, xmax = 35, ymin = 35, ymax = 70)

Maps.picea <- ggplot(data=mapBase$geometry)+ 
  geom_sf(fill = "grey90", colour = "black", alpha = .8, linewidth = .25)+ #plot map of France
  xlab(" ")+ ylab(" ")+ #white or aliceblue
  geom_point(data = plot_locations.picea %>%  mutate(x = unlist(map(geometry,1)),
                                                     y = unlist(map(geometry,2)))%>% 
               as.data.frame(), aes(x = x , y = y, col = Collection_method, fill = Collection_method), 
             stroke = 1, alpha = .8,
             shape = 21, size = 3)+
  scale_size_continuous(
    breaks = c(0.1,0.3,0.8),
    range = c(0,7)
  )+
  coord_sf(xlim = c(-10, 25), ylim = c(40, 60))+
  scale_y_continuous(breaks = c(40, 50, 60, 70))+
  scale_x_continuous(breaks = c(-10, 0, 10, 20))+
  ylab("Latitude")+xlab("Longitude")+
  ggpubr::theme_pubr()+
  scale_color_futurama()+
  scale_fill_futurama()+
  guides(col = "none", fill=guide_legend(title=NULL, override.aes = list(size=8)))+
  theme(legend.position = 'none')
Maps.picea

#make moving window plot comparison methods csp vs climwin  
picea2method = statistics_csp_method.picea %>%
  group_by(sitenewname) %>% 
  dplyr::slice(which.max(r2))   %>% 
  mutate(method = 'climate sensitivity profile') %>% 
  ungroup() %>% 
  bind_rows(statistics_absolute_climwin_picea %>% mutate(method = 'sliding moving window')) %>% 
  left_join(methods.collection.picea) %>% 
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot()+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname, col = Collection_method), size = 2) +
  geom_point(aes(y = window.open, x = sitenewname, col = Collection_method, fill = Collection_method), size = 1.5, shape = 22) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed')+
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_futurama()+
  scale_fill_futurama()+
  ggpubr::theme_cleveland()+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  ylab('Days reversed')+
  facet_grid(method~.)+
  theme(axis.text.y = element_blank())

#extract R2 from models 
R2picea = statistics_csp_method.picea %>%
  group_by(sitenewname) %>% 
  dplyr::slice(which.max(r2)) %>% 
  dplyr::select(sitenewname, r2) %>% 
  mutate(method = 'climate sensivitiy profil') %>% 
  bind_rows(statistics_absolute_climwin_picea %>% 
              dplyr::select(sitenewname, R2) %>% 
              rename(r2 = R2) %>% 
              mutate(method = 'sliding moving window')) %>% 
  ggplot(aes(x = r2, fill = method)) +
  geom_density(alpha = 0.42) +
  labs(x = "R-squared (r2)", y = "Density") +
  ggpubr::theme_pubclean()+
  theme(legend.position = c(.9, .8))+
  theme(legend.title = element_blank())+
  scale_color_viridis_d(option = 'magma')+
  scale_fill_viridis_d(option = 'magma')

#combine all the map
Figure.picea = cowplot::plot_grid((Maps.picea/R2picea), picea2method, labels = 'auto')
cowplot::save_plot(here::here("figures/Figure.picea.png"), 
                   Figure.picea, 
                   ncol = 1.8, nrow = 1.5, dpi = 300, bg = "white")
