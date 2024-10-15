#https://community.rstudio.com/t/unable-to-install-packages-from-github/124372
#Sys.unsetenv("GITHUB_PAT")
library(tidyverse)
library(climwin)
library(broom)
library(mgcv)
library(here)
#for maps 
library(sf)
library("rworldmap")
library("rworldxtra")
require(rgdal)
library(ggspatial)
library(ggsci)


functions <- list.files(here("fun"), full.names = T) %>%
  purrr::map(source)



new.folder = FALSE 
if(new.folder == TRUE){
  dir.create(here('climate_dailyEOBS'))
  dir.create(here('figures'))
}

#access climate zip file 
#https://osf.io/287ms/ and add it in climate_dailyEOBS folder

nfile = check_folder_and_contents(here('climate_dailyEOBS'), file_pattern = '.qs')
if(nfile$folder_exists == FALSE){stop(print('missing folder climate_dailyEOBS'))}
if(length(nfile$files_present) == 0){stop(print('missing file to add in the folder of interest name climate_dailyEOBS'))}
#SHOULD HAVE NO WARNING HERE

# I am generating calendar 
#double check values
calendar = goingbackpastdayscalendar(refday = 305, #double check later, because climwin start at 1 or 0
                                     lastdays = 600, 
                                     yearback = 2) %>% 
  filter(YEAR == 1950) %>% 
  dplyr::select(MONTHab, DOY, days.reversed) %>% 
  mutate(datefake  = as.Date(DOY, origin = "1948-01-01"),
         day.month = format(datefake,"%m-%d"), 
         year = case_when(
           days.reversed < 366 ~ 1,
           days.reversed >= 366 & days.reversed < 730 ~ 2,
           days.reversed >= 730 & days.reversed < 1095 ~ 3,
           TRUE ~ 4))

initial.data.mastree <- read.csv('/Users/vjourne/Library/CloudStorage/Dropbox/Mastree/MASTREEplus_2024-06-26_V2.csv',stringsAsFactors = F)

#load data from JJ Foest, last V2 version 
url.mastreev2 = 'https://raw.githubusercontent.com/JJFoest/MASTREEplus/main/Data/MASTREEplus_2024-06-26_V2.csv'

initial.data.mastree <- data.table::fread(url.mastreev2)

#now use the filtering 
#keep cont, remove index, flower and pollen 
Fagus.seed = initial.data.mastree %>%
  filter(!Variable == "flower" & VarType == "C" & !Unit=="index" & !Variable == "pollen") %>% 
  group_by(Species, VarType) %>% 
  mutate(nMax = n()) %>% 
  filter(Species == "Fagus sylvatica") %>% 
  #filter(Unit == "seeds/m2" | Unit == "seeds/individual") %>% 
  filter(Year > 1952 & Year < 2021) %>% 
  mutate(sitenewname= paste0(Alpha_Number, "_",Site_number, "_", Species_code)) %>% 
  group_by(sitenewname) %>% 
  mutate(log.seed = log(1+Value),
         scale.seed = scale(Value),
         n = n()) %>% 
  filter(n > 20) %>% #initially 14
  mutate() %>% 
  mutate(Date = paste0( "15/06/",Year)) %>% 
  mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
         plotname.lon.lat = paste0("longitude=",Longitude, "_","latitude=", Latitude)) %>% 
  group_by(plotname.lon.lat) %>%   
  #sample_frac(0.5) %>%        
  ungroup() %>% 
  as.data.frame()

#get the mothod collection
methods.collection.mv2 = Fagus.seed %>% dplyr::select(sitenewname, Country, Collection_method, Length) %>% distinct()


seed.production.plot = Fagus.seed %>% 
  group_by(plotname.lon.lat) %>% 
  complete(Year = full_seq(Year, 1)) %>% 
  ggplot(aes(x=Year, y = log.seed, group = plotname.lon.lat))+
  geom_line(na.rm = FALSE, alpha = .1)+
  ggpubr::theme_pubr()+
  ylab('Seed production (log)')


plot_locations<- st_as_sf(Fagus.seed %>% 
                            dplyr::select(Longitude, Latitude, plotname.lon.lat, Collection_method) %>% distinct(), coords = c("Longitude", "Latitude")) %>% 
  st_set_crs(4326)
#faguseuforgen = st_read(dsn = "/Users/vjourne/Documents//Projets_annexes/phenology ventoux/transfer_4516897_files_eef1955b/Chorological data for the main European woody species/chorological_maps_dataset_20221025/chorological_maps_dataset/Fagus sylvatica/shapefiles/Fagus_sylvatica_sylvatica_plg_clip.shp", stringsAsFactors = F)

mapBase <- getMap(resolution = "high") %>%   
  st_as_sf() %>% 
  st_make_valid()%>%  st_crop(mapBase, xmin = -12, xmax = 35, ymin = 35, ymax = 70)

MapsMastree <- ggplot(data=mapBase$geometry)+ 
  geom_sf(fill = "grey", colour = "grey")+ #plot map of France
  xlab(" ")+ ylab(" ")+ #white or aliceblue
  #geom_sf(data = faguseuforgen, col = "#FFDB6D", fill = "#FFDB6D", alpha = 0.55, linewidth = .35)+
  geom_point(data = plot_locations %>%  mutate(x = unlist(map(geometry,1)),
                                               y = unlist(map(geometry,2)))%>% 
               as.data.frame(), aes(x = x , y = y, col = Collection_method, fill = Collection_method), 
             stroke = 1, alpha = .8,
             shape = 21, size = 4.3)+ #size = rmse,
  scale_size_continuous(
    breaks = c(0.1,0.3,0.8),
    range = c(0,7)
  )+
  coord_sf(xlim = c(-10, 25), ylim = c(40, 60))+
  scale_y_continuous(breaks = c(40, 50, 60, 70))+
  scale_x_continuous(breaks = c(-10, 0, 10, 20))+
  ylab("Latitude")+xlab("Longitude")+
  ggpubr::theme_pubr()+
  #theme(legend.position = c(.1,.8))+
  #scale_color_aaas()+
  #scale_fill_aaas()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  guides(col = "none", fill=guide_legend(title=NULL, override.aes = list(size=8)))

library(patchwork)

cowplot::save_plot(here("figures/data.des.png"),MapsMastree+seed.production.plot, 
                   ncol = 2, nrow = 1.5, dpi = 300)

#some option 
#here for climwin 
optionwindows = 'absolute'
referenceclimwin = c(01, 11)
#going back to one year before + 6 month more or less
#range = c((round(365/2)+365), 0)
range = c(600, 0)

#climate path 
#temperature.path = here::here('~/Documents/GITprojects/breakdownmast/climatedailyERA5')
#climate.beech.path <- "/Users/vjourne/Library/CloudStorage/Dropbox/MastreeForemast/climateSiteDaily"
climate.beech.path <- here::here('climate_dailyEOBS/')

#beech site unique name 
#beech.site.all = unique(Fagus.seed$sitenewname)
beech.site.all = unique(Fagus.seed$plotname.lon.lat)
#they have same length, each sites are unique here! 

# Define a function to process each site
#I want the function to work only for one site
#option to run in parallel to make it faster (faster than map_dfr, but need double checking because if not well defined = issues)
run.climwin <- F
#take some time ! 
if(run.climwin==T){
  library(furrr)#make it in parallele, https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
  plan(multisession)#background R sessions (on current machine) , the computer is flying (alt cluster)
  
  statistics_absolute_climwin <- future_map_dfr(
    beech.site.all, 
    ~ {
      # format the climate, .x is the current site
      climate_data <- format_climate_data(site = .x ,
                                          path = climate.beech.path, 
                                          scale.climate = T)  
      # Run climwin per site
      climwin_site_days(
        climate_data = climate_data,
        data = Fagus.seed %>% filter(plotname.lon.lat == .x),  
        site.name = .x,  
        range = range,
        cinterval = 'day',
        refday = c(01, 11),  
        optionwindows = 'absolute',  
        climate_var = 'TMEAN'  
      )
    }
  )
  
qs::qsave(statistics_absolute_climwin, 
            here('statistics_absolute_climwin.qs'))}else{
  statistics_absolute_climwin = qs::qread('statistics_absolute_climwin.qs')
            }


mean(statistics_absolute_climwin$WindowOpen)
mean(statistics_absolute_climwin$WindowClose)
#159-427


climwin.dd = statistics_absolute_climwin %>% 
  dplyr::select(sitenewname, WindowOpen, WindowClose) %>% 
  pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed)) 



climwin.all.sites = climwin.dd %>% dplyr::select(sitenewname, name, days.reversed) %>% 
  pivot_wider(names_from = "name", values_from = "days.reversed") %>% 
  left_join(methods.collection.mv2) %>% 
  #arrange(Collection_method) %>% 
  mutate(sitenewname = fct_reorder(sitenewname, Collection_method)) %>%
ggplot()+
  geom_segment(aes(y = WindowOpen, yend = WindowClose, x = sitenewname, col = Collection_method), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")

cowplot::save_plot(here("figures/climwin_allsites_methods.png"),climwin.all.sites, ncol = 1.5, nrow = 1.5, dpi = 300)

  
quibble2(statistics_absolute_climwin$WindowOpen, q = c(0.25, 0.5, 0.75))
quibble2(statistics_absolute_climwin$WindowClose, q = c(0.25, 0.5, 0.75))



#it is still useful for me because here I will also work with simple correlation later 
results.moving.site = FULL.moving.climate.analysis(seed.data = Fagus.seed,
                             climate.path = climate.beech.path,
                             refday = 305,
                             lastdays = max(range),
                             rollwin = 1)


ggplot(results.moving.site)+
  geom_bar(aes(y = r.squared, x = days.reversed), size = 2, stat="identity" ) +
  geom_vline(xintercept = 133, color = 'red')+
  geom_vline(xintercept = 498, color = 'red')



#create 61 list - 1 per site 
Results_daily = results.moving.site %>%
  arrange(plotname.lon.lat) %>% 
  group_by(plotname.lon.lat) %>%
  group_split()

name_daily =   results.moving.site %>% 
  arrange(plotname.lon.lat) %>% 
  group_by(plotname.lon.lat) %>%
  group_keys()

# Set names of the list based on group keys
names(Results_daily) <- apply(name_daily, 1, paste)

rollwin = 1
lastdays = max(range)
myform.fin = formula('log.seed ~ mean.temperature')
refday = 305

#map dfr is enough, quite fast if not optimize gam model 
statistics_csp_method = map_dfr(
  1:length(Results_daily), 
  ~CSP_function_site(
    Results_daily[[.]], 
    unique(Results_daily[[.]]$plotname.lon.lat),
    seed.data = Fagus.seed %>% arrange(plotname.lon.lat), 
    climate.path = climate.beech.path,
    refday = 305,
    lastdays = 600,
    rollwin = 1
  ))
qs::qsave(statistics_csp_method, 
          here('statistics_csp.qs'))




output_fit_summary.best.csp = statistics_csp_method %>%
  group_by(sitenewname) %>% 
  dplyr::slice(which.max(r2))


csp.all.sites = output_fit_summary.best.csp  %>% 
  left_join(methods.collection.mv2) %>% 
  ungroup() %>% 
  arrange(Collection_method) %>% 
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot()+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname, col = Collection_method), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")
csp.all.sites

cowplot::save_plot("csp_allsites_methods.png",csp.all.sites, ncol = 1.5, nrow = 1.5, dpi = 300)


quibble2(output_fit_summary.best.csp$window.open, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.best.csp$window.close, q = c(0.25, 0.5, 0.75))


#no subset 
ggplot(statistics_csp_method)+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,547)

#change name here between closing and opning I mixed them 
quibble2(statistics_csp_method$window.open, q = c(0.25, 0.5, 0.75))
quibble2(statistics_csp_method$window.close, q = c(0.25, 0.5, 0.75))



#################################################################
#############################################
#####
#now test last method P spline 
#####
#############################################

site = unique(Fagus.seed$plotname.lon.lat)#[30]
#refday = 305
statistics_psr_method = map_dfr(
  site, 
  ~ PSR_function_site(
    site = .x,
    climate.path = climate.beech.path,
    seed.data = Fagus.seed,  # Use .x correctly here
    tot_days = max(range),
    matrice = c(3,1),
    knots = NULL
  )
)



ggplot(statistics_psr_method %>% 
         group_by(sitenewname) %>% 
         slice(which.max(r2)))+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
  geom_point(aes(y = window.open, x = sitenewname), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,600)

output_fit_summary.psr.best = statistics_psr_method  %>% 
  group_by(sitenewname) %>% 
  slice(which.max(r2))
psr.all.sites =  output_fit_summary.psr.best%>%  
  left_join(methods.collection.mv2) %>% 
  ungroup() %>% 
  arrange(Collection_method) %>% 
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot()+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname, col = Collection_method), size = 2) +
  geom_point(aes(y = window.open, x = sitenewname, col = Collection_method), size = 2, shape = 15) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")
psr.all.sites

cowplot::save_plot(here("figures/psr_allsites_methods.png"),psr.all.sites, ncol = 1.5, nrow = 1.5, dpi = 300)


quibble2(output_fit_summary.psr.best$window.close, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.psr.best$window.open, q = c(0.25, 0.5, 0.75))

#####################################################################################
######################################################################################################
#####################################################################################
####################################################################
#SIMPLE METHOD TO IDENTIFY WEATHER CUES 
#here I am using this only for fagus...
statistics_basic_method = map_dfr(
  1:length(Results_daily), 
  ~basiccues_function_site(
    Results_daily[[.]], 
    unique(Results_daily[[.]]$plotname.lon.lat),
    Fagus.seed = Fagus.seed, 
    climate.path = climate.beech.path,
    refday = 305,
    lastdays = 600,
    rollwin = 1
  ))


output_fit_summary.basic.best = statistics_basic_method  %>% 
  group_by(sitenewname) %>% 
  slice(which.max(r2))%>% 
  ungroup()

basic.all.sites = output_fit_summary.basic.best %>% 
  left_join(methods.collection.mv2) %>% 
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot()+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname, col = Collection_method), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")
basic.all.sites
cowplot::save_plot(here("figures/basic_allsites_methods.png"),basic.all.sites, ncol = 1.5, nrow = 1.5, dpi = 300)

quibble2(output_fit_summary.basic.best$window.open, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.basic.best$window.close, q = c(0.25, 0.5, 0.75))

#put a long lag 
windows.avg = output_fit_summary.basic.best %>% 
  dplyr::select(sitenewname, window.open,  window.close) %>% 
  mutate(method = 'signal.process') %>% 
  bind_rows(output_fit_summary.psr.best%>% 
              dplyr::select(sitenewname, window.open, window.close) %>% 
              mutate(method = 'p.spline.reg')) %>% 
  bind_rows(statistics_absolute_climwin %>% 
              dplyr::select(sitenewname, WindowOpen, WindowClose) %>% rename(window.open = WindowOpen,
                                                                             window.close = WindowClose) %>% 
              mutate(method = 'climwin.abs')) %>% 
  bind_rows(output_fit_summary.best.csp %>% 
              dplyr::select(sitenewname, window.open , window.close) %>% 
              mutate(method = 'climate.sens.profil') )%>%
  group_by(method) %>%
  summarise(wind.open_median = median(window.open, na.rm = TRUE),
            wind.close_median = median(window.close, na.rm = TRUE),
            wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
            wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
            wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
            wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE),
            wind.open_mean = mean(window.open),
            wind.close_mean = mean(window.close),
            wind.open_se = se(window.open),
            wind.close_se = se(window.close))


median.windows.plot = windows.avg %>% 
  dplyr::select(method, wind.open_median:wind.open_q75) %>% 
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>% 
  separate_wider_delim(typewind, delim = '_', names=c('windows.type', 'metric')) %>% 
  pivot_wider(names_from = 'metric', values_from = value) %>% 
ggplot(aes(x = method, group = windows.type, col = windows.type)) +
  geom_point(aes(y = median), size = 2,
             position = position_dodge(width = 0.2))+
  geom_pointrange(mapping=aes(y=median, 
                              ymin=q25, 
                              ymax=q75), 
                  position = position_dodge(width = 0.2))+
  coord_flip()+
  xlab('')+ylab('Days reversed')+
  scale_color_brewer(palette = "Paired")+
  scale_fill_brewer(palette = "Paired")+
  ggpubr::theme_pubr()+
  theme(legend.position = c(.8,.8),
        legend.title = element_blank())+
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')

#GENERIC PLOT R2 and AIC
bestr2 = output_fit_summary.basic.best %>% 
  dplyr::select(sitenewname, r2, AIC) %>% 
  mutate(method = 'signal.process') %>% 
  bind_rows(output_fit_summary.psr.best%>% 
              dplyr::select(sitenewname, r2, AIC) %>% 
              mutate(method = 'p.spline.reg')) %>% 
  bind_rows(statistics_absolute_climwin %>% 
              dplyr::select(sitenewname, R2, AIC) %>% rename(r2 = R2) %>% 
              mutate(method = 'climwin.abs')) %>% 
  bind_rows(output_fit_summary.best.csp %>% 
              dplyr::select(sitenewname, r2, AIC) %>% 
              mutate(method = 'climate.sens.profil') )

mean_r2 <- bestr2 %>%
  group_by(method) %>%
  summarise(mean_r2 = median(r2, na.rm = TRUE))

averagedensr2method = bestr2 %>% 
  ggplot(aes(x = r2, fill = method)) +
  geom_density(alpha = 0.4) +
  labs(x = "R-squared (r2)", y = "Density") +
  geom_vline(data = mean_r2, aes(xintercept = mean_r2, color = method),
             linetype = "dashed", size = 1) +
  ggpubr::theme_pubclean()+
  theme(legend.position = c(.9, .8))

averagedensr2method
cowplot::save_plot(here("figures/averager2.method.png"),averagedensr2method+median.windows.plot, ncol = 1.8, nrow = 1.4, dpi = 300)


####################################################################################
##############################################################################
#STOP HERE FOR NOW

####################################################################################
##############################################################################
#######################
#now test the frac dataset 
run_sampling_fraction_all_methods = function(fraction = 0.5){
  
  Fagus.seed.subseting50 = initial.data.mastree %>%
    filter(!Variable == "flower" & VarType == "C" & !Unit=="index" & !Variable == "pollen") %>% 
    group_by(Species, VarType) %>% 
    mutate(nMax = n()) %>% 
    filter(Species == "Fagus sylvatica") %>% 
    #filter(Unit == "seeds/m2" | Unit == "seeds/individual") %>% 
    filter(Year > 1952 & Year < 2021) %>% 
    mutate(sitenewname= paste0(Alpha_Number, "_",Site_number, "_", Species_code)) %>% 
    group_by(sitenewname) %>% 
    mutate(log.seed = log(1+Value),
           scale.seed = scale(Value),
           n = n()) %>% 
    filter(n > 20) %>% #initially 14
    mutate() %>% 
    mutate(Date = paste0( "15/06/",Year)) %>% 
    mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
           plotname.lon.lat = paste0("longitude=",Longitude, "_","latitude=", Latitude)) %>% 
    group_by(plotname.lon.lat) %>%   
    sample_frac(fraction) %>%        
    ungroup() %>% 
    as.data.frame()
  
  #run climwin 
  library(furrr)#make it in parallele, https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
  plan(multisession)#background R sessions (on current machine) , the computer is flying (alt cluster)
  
  statistics_absolute_climwin_50 = future_map_dfr(
    beech.site.all, 
    ~ climwin_site_days(
      climate.beech.path = climate.beech.path,
      data = Fagus.seed.subseting50 %>% 
        filter(plotname.lon.lat == .x),  # Use .x correctly here
      site.name = .x,
      range = range,
      cinterval = 'day',
      refday = c(01, 11),
      optionwindows = 'absolute',
      climate_var = 'TMEAN'
    )
  )
  
  results.moving.site_50 = FULL.moving.climate.analysis(Fagus.seed = Fagus.seed.subseting50,
                                                        climate.beech.path = climate.beech.path,
                                                        refday = 305,
                                                        lastdays = max(range),
                                                        rollwin = 1)
  
  Results_daily_50 = results.moving.site_50 %>%
    arrange(plotname.lon.lat) %>% 
    group_by(plotname.lon.lat) %>%
    group_split()
  
  name_daily_50 =   results.moving.site_50 %>% 
    arrange(plotname.lon.lat) %>% 
    group_by(plotname.lon.lat) %>%
    group_keys()
  
  # Set names of the list based on group keys - run csp
  names(Results_daily_50) <- apply(name_daily_50, 1, paste)
  
  statistics_csp_method_50 = map_dfr(
    1:length(Results_daily_50), 
    ~CSP_function_site(
      Results_daily_50[[.]], 
      unique(Results_daily_50[[.]]$plotname.lon.lat),
      Fagus.seed = Fagus.seed.subseting50, 
      climate.beech.path = climate.beech.path,
      refday = 305,
      lastdays = 600,
      rollwin = 1
    ))
  
  #run psr 
  statistics_psr_method_50 = map_dfr(
    site, 
    ~ PSR_function_site(
      site = .x,
      Fagus.seed = Fagus.seed.subseting50,  # Use .x correctly here
      tot_days = max(range),
      lastdays = 600,
      matrice = c(3,1),
      knots = NULL
    )
  )
  
  #run basic 
  statistics_basic_method_50 = map_dfr(
    1:length(Results_daily_50), 
    ~basiccues_function_site(
      Results_daily_50[[.]], 
      unique(Results_daily_50[[.]]$plotname.lon.lat),
      Fagus.seed = Fagus.seed.subseting50, 
      climate.beech.path = climate.beech.path,
      refday = 305,
      lastdays = 600,
      rollwin = 1
    ))
  
  list.all = list(statistics_absolute_climwin_50,
                  statistics_csp_method_50,
                  statistics_psr_method_50,
                  statistics_basic_method_50)
  return(list.all)
}


sample30percent = run_sampling_fraction_all_methods(fraction = 0.3)
sample50percent = run_sampling_fraction_all_methods(fraction = 0.5)
sample70percent = run_sampling_fraction_all_methods(fraction = 0.7)
#sample90percent = run_sampling_fraction_all_methods(fraction = 0.9)





climwin.sensi = sample30percent[[1]] %>% 
  dplyr::select(sitenewname, WindowOpen, WindowClose) %>% 
  pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed)) %>% dplyr::select(sitenewname, name, days.reversed) %>% 
  pivot_wider(names_from = "name", values_from = "days.reversed") %>% 
  mutate(sample.size = '0.3') %>% 
  bind_rows(
    sample50percent[[1]] %>% 
  dplyr::select(sitenewname, WindowOpen, WindowClose) %>% 
  pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed)) %>% dplyr::select(sitenewname, name, days.reversed) %>% 
  pivot_wider(names_from = "name", values_from = "days.reversed") %>% 
  mutate(sample.size = '0.5')) %>% 
  bind_rows(
    sample70percent[[1]] %>% 
      dplyr::select(sitenewname, WindowOpen, WindowClose) %>% 
      pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
      left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed)) %>% dplyr::select(sitenewname, name, days.reversed) %>% 
      pivot_wider(names_from = "name", values_from = "days.reversed") %>% 
      mutate(sample.size = '0.7')) %>% 
  rename(window.open = WindowOpen,window.close = WindowClose) %>% 
  mutate(method = 'climwin.abs')

csp.sensi = sample30percent[[2]]  %>% 
  group_by(sitenewname) %>% 
  slice(which.max(r2))%>% 
  ungroup()%>% 
  dplyr::select(sitenewname, window.open , window.close)%>% 
  mutate(sample.size = '0.3') %>% 
  bind_rows(sample50percent[[2]]  %>% 
              group_by(sitenewname) %>% 
              slice(which.max(r2))%>% 
              ungroup()%>% 
              dplyr::select(sitenewname, window.open , window.close)%>% 
              mutate(sample.size = '0.5') ) %>% 
  bind_rows(sample70percent[[2]]  %>% 
              group_by(sitenewname) %>% 
              slice(which.max(r2))%>% 
              ungroup()%>% 
              dplyr::select(sitenewname, window.open , window.close)%>% 
              mutate(sample.size = '0.7') )%>% 
  mutate(method = 'csp')

psr.sensi = sample30percent[[3]]  %>% 
  group_by(sitenewname) %>% 
  slice(which.max(r2))%>% 
  ungroup()%>% 
  dplyr::select(sitenewname, window.open , window.close)%>% 
  mutate(sample.size = '0.3') %>% 
  bind_rows(sample50percent[[3]]  %>% 
              group_by(sitenewname) %>% 
              slice(which.max(r2))%>% 
              ungroup()%>% 
              dplyr::select(sitenewname, window.open , window.close)%>% 
              mutate(sample.size = '0.5') ) %>% 
  bind_rows(sample70percent[[3]]  %>% 
              group_by(sitenewname) %>% 
              slice(which.max(r2))%>% 
              ungroup()%>% 
              dplyr::select(sitenewname, window.open , window.close)%>% 
              mutate(sample.size = '0.7') )%>% 
  mutate(method = 'psr')

basic.sensi = sample30percent[[4]]  %>% 
  group_by(sitenewname) %>% 
  slice(which.max(r2))%>% 
  ungroup()%>% 
  dplyr::select(sitenewname, window.open , window.close)%>% 
  mutate(sample.size = '0.3') %>% 
  bind_rows(sample50percent[[4]]  %>% 
              group_by(sitenewname) %>% 
              slice(which.max(r2))%>% 
              ungroup()%>% 
              dplyr::select(sitenewname, window.open , window.close)%>% 
              mutate(sample.size = '0.5') ) %>% 
  bind_rows(sample70percent[[4]]  %>% 
              group_by(sitenewname) %>% 
              slice(which.max(r2))%>% 
              ungroup()%>% 
              dplyr::select(sitenewname, window.open , window.close)%>% 
              mutate(sample.size = '0.7') )%>% 
  mutate(method = 'basic')

windows.avg.sensitivity = basic.sensi%>%
  bind_rows(psr.sensi) %>% 
  bind_rows(csp.sensi) %>% 
  bind_rows(climwin.sensi) %>% 
  group_by(method,sample.size) %>%
  summarise(wind.open_median = median(window.open, na.rm = TRUE),
            wind.close_median = median(window.close, na.rm = TRUE),
            wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
            wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
            wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
            wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE),
            wind.open_mean = mean(window.open),
            wind.close_mean = mean(window.close),
            wind.open_se = se(window.open),
            wind.close_se = se(window.close))


sensi.data.plot = windows.avg.sensitivity %>% 
  dplyr::select(method, sample.size, wind.open_median:wind.open_q75) %>% 
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>% 
  separate_wider_delim(typewind, delim = '_', names=c('windows.type', 'metric')) %>% 
  pivot_wider(names_from = 'metric', values_from = value) %>% 
  mutate(sample.size = as.numeric(sample.size),
         method_sample = paste(method, sample.size, sep = " (n=") %>% paste0(")")) %>% 
  ggplot(aes(x = method_sample, group = windows.type, col = windows.type)) +
  geom_point(aes(y = median), size = 2,
             position = position_dodge(width = 0.5)) +
  geom_pointrange(mapping = aes(y = median, ymin = q25, ymax = q75), 
                  position = position_dodge(width = 0.5)) +
  coord_flip() +
  xlab('Method (Sample Size)') + 
  ylab('Days reversed') +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  ggpubr::theme_pubr() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')

sensi.data.plot
cowplot::save_plot(here("figures/sensitivity.method.png"),sensi.data.plot, ncol = 1.4, nrow = 1.4, dpi = 300)



################################
#example with all quercus robur petraea and spp (hybrid)
Quercus.seed = initial.data.mastree %>%
  filter(!Variable == "flower" & VarType == "C" & !Unit=="index" & !Variable == "pollen") %>% 
  group_by(Species, VarType) %>% 
  mutate(nMax = n()) %>% 
  filter(str_detect(Species, 'Quercus robur|Quercus petraea|Quercus spp.')) %>% 
  filter(!str_detect(Country, 'United States of America|Russian Federation (the)')) %>% 
  #filter(Unit == "seeds/m2" | Unit == "seeds/individual") %>% 
  filter(Year > 1952 & Year < 2021) %>% 
  mutate(sitenewname= paste0(Alpha_Number, "_",Site_number, "_", Species_code)) %>% 
  group_by(sitenewname) %>% 
  mutate(log.seed = log(1+Value),
         scale.seed = scale(Value),
         n = n()) %>% 
  filter(n > 20) %>% #initially 14
  mutate() %>% 
  mutate(Date = paste0( "15/06/",Year)) %>% 
  mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
         plotname.lon.lat = paste0("longitude=",Longitude, "_","latitude=", Latitude)) %>% 
  group_by(plotname.lon.lat) %>%   
  #sample_frac(0.5) %>%        
  ungroup() %>% 
  as.data.frame()

quercus.site.all = unique(Quercus.seed$plotname.lon.lat)
#they have same length, each sites are unique here! 

# Define a function to process each site
#I want the function to work only for one site
#option to run in parallel to make it faster
run.climwin <- T
#take some time ! 
if(run.climwin==T){
  library(furrr)#make it in parallele, https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
  plan(multisession)#background R sessions (on current machine) , the computer is flying (alt cluster)
  
  statistics_absolute_climwin_quercus = future_map_dfr(
    quercus.site.all, 
    ~ climwin_site_days(
      climate.path = climate.beech.path,
      data = Quercus.seed %>% 
        filter(plotname.lon.lat == .x),  # Use .x correctly here
      site.name = .x,
      range = range,
      cinterval = 'day',
      refday = c(01, 11),
      optionwindows = 'absolute',
      climate_var = 'TMEAN'
    )
  )
  qs::qsave(statistics_absolute_climwin_quercus, 
            here('statistics_absolute_climwin_QUERCUS.qs'))}else{
              statistics_absolute_climwin_quercus = qs::qread('statistics_absolute_climwin_QUERCUS.qs')
            }


format.quercus.win = statistics_absolute_climwin_quercus %>% 
  dplyr::select(sitenewname, estimate, WindowOpen, WindowClose) %>% 
  pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed)) 


methods.collection.quercus = Quercus.seed %>% dplyr::select(sitenewname, Country, Collection_method, Length, Longitude, Latitude) %>% distinct()


format.quercus.win %>% dplyr::select(sitenewname, estimate, name, days.reversed) %>% 
  pivot_wider(names_from = "name", values_from = "days.reversed") %>% 
  mutate(sign = ifelse(estimate >0 , '+', '-')) %>% 
  left_join(methods.collection.quercus) %>% 
  arrange(Latitude) %>% 
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot()+
  geom_segment(aes(y = WindowOpen, yend = WindowClose, x = sitenewname), size = 2) +
 # geom_point(aes(y = WindowOpen, x = sitenewname, col = Collection_method), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")

plot_locations.quercus<- st_as_sf(Quercus.seed %>% 
                            dplyr::select(Longitude, Latitude, plotname.lon.lat, Collection_method) %>% distinct(), coords = c("Longitude", "Latitude")) %>% 
  st_set_crs(4326)

Maps.quercus <- ggplot(data=mapBase$geometry)+ 
  geom_sf(fill = "grey", colour = "grey")+ #plot map of France
  xlab(" ")+ ylab(" ")+ #white or aliceblue
  #geom_sf(data = faguseuforgen, col = "#FFDB6D", fill = "#FFDB6D", alpha = 0.55, linewidth = .35)+
  geom_point(data = plot_locations.quercus %>%  mutate(x = unlist(map(geometry,1)),
                                               y = unlist(map(geometry,2)))%>% 
               as.data.frame(), aes(x = x , y = y, col = Collection_method, fill = Collection_method), 
             stroke = 1, alpha = .8,
             shape = 21, size = 4.3)+ #size = rmse,
  scale_size_continuous(
    breaks = c(0.1,0.3,0.8),
    range = c(0,7)
  )+
  coord_sf(xlim = c(-10, 25), ylim = c(40, 60))+
  scale_y_continuous(breaks = c(40, 50, 60, 70))+
  scale_x_continuous(breaks = c(-10, 0, 10, 20))+
  ylab("Latitude")+xlab("Longitude")+
  ggpubr::theme_pubr()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  guides(col = "none", fill=guide_legend(title=NULL, override.aes = list(size=8)))
Maps.quercus

#test csp 
results.moving.site.quercus = FULL.moving.climate.analysis(Fagus.seed = Quercus.seed,
                                                   climate.path = climate.beech.path,
                                                   refday = 305,
                                                   lastdays = max(range),
                                                   rollwin = 1)
Results_daily.quercus = results.moving.site.quercus %>%
  arrange(sitenewname) %>% 
  group_by(sitenewname) %>%
  group_split()

name_daily.quercus =   results.moving.site.quercus %>% 
  arrange(sitenewname) %>% 
  group_by(sitenewname) %>%
  group_keys()

# Set names of the list based on group keys
names(Results_daily.quercus) <- apply(name_daily.quercus, 1, paste)

statistics_csp_method.quercus = map_dfr(
  1:length(Results_daily.quercus), 
  ~CSP_function_site(
    Results_daily.quercus[[.]], 
    unique(Results_daily.quercus[[.]]$sitenewname),
    Fagus.seed = Quercus.seed%>% 
      arrange(sitenewname), 
    climate.path = climate.beech.path,
    refday = 305,
    lastdays = 600,
    rollwin = 1
  ))

statistics_csp_method.quercus %>%
  group_by(sitenewname) %>% 
  dplyr::slice(which.max(r2))  %>% 
  left_join(methods.collection.quercus) %>% 
  ungroup() %>% 
  #arrange(Collection_method) %>% 
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot()+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")

statistics_csp_method.quercus %>%
  group_by(sitenewname) %>% 
  dplyr::slice(which.max(r2)) %>% 
  dplyr::select(sitenewname, r2) %>% 
  mutate(method = 'CSP') %>% 
  bind_rows(statistics_absolute_climwin_quercus %>% 
              dplyr::select(sitenewname, R2) %>% 
              rename(r2 = R2) %>% 
              mutate(method = 'climwin')) %>% 
  ggplot(aes(x = r2, fill = method)) +
  geom_density(alpha = 0.4) +
  labs(x = "R-squared (r2)", y = "Density") +
  ggpubr::theme_pubclean()+
  theme(legend.position = c(.9, .8))


########################################################################
####################################################################################
########################################################################
############################################################
######MAKE simulation study case with climwin 

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(lubridate)


# Simulate climate data: assume daily data for 20 years
set.seed(123)  
startyear = 1940
years <- startyear:2020
days_per_year <- 365
climate_data <- data.frame(
  date = seq.Date(from = as.Date(paste0(startyear,"-01-01")), by = "day", length.out = length(years) * days_per_year),
  temp = runif(length(years) * days_per_year, min = 10, max = 30)  # random temperatures
)
  
climate_data <- climate_data %>%
  mutate(year = format(date, "%Y"),
         day_of_year = yday(date),
         year = as.numeric(year))

# Select ~ one week in June 
june_week <- climate_data %>%
  filter(day_of_year >= 150 & day_of_year <= 160)

#then make it strongly correlated to the climate window
test.simulation <- function(simulation.bio.data = june_week,
                           seed_production,
                           sample.fraction = 1,
                           correlation = 0.2){
  seed_production <- simulation.bio.data %>%
    group_by(year) %>%
    summarise(mean_temp_june_week = mean(temp)) %>%
    mutate(random_noise = rnorm(n(), mean = 0, sd = sd(mean_temp_june_week) * sqrt((1 - correlation^2) / correlation^2))) %>%
    mutate(log.seed = log(round(100 + mean_temp_june_week * 5 + random_noise))) %>%
    filter(year > startyear)%>% #to select one year less 
    mutate(year = as.numeric(as.character(year))) %>% 
    mutate(Date = paste0( "15/06/",year)) %>% #need a date, even if absolute window
    mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y")) %>% #make it as in the climate date
    sample_frac(sample.fraction)  
  
  print(cor(seed_production$mean_temp_june_week, seed_production$log.seed))
  
  
  simulation1 = climwin_site_days(climate_data,
                                  data = seed_production,
                                  site.name = 'simulation', #character will be added in the final form file
                                  range = c(600, 0),
                                  cinterval = 'day',
                                  refday = c(01, 11),
                                  optionwindows = 'absolute',
                                  climate_var = 'temp',
                                  formulanull = formula(log.seed ~ 1)) 
  return(simulation1 %>% mutate(sample.fraction.used = sample.fraction))
}

seq.sim = seq(0.3,0.9,0.1)
test.simulation.results = NULL
for(i in 1:length(seq.sim)){
  simulation1 = test.simulation(simulation.bio.data = june_week,
                     seed_production,
                     sample.fraction = seq.sim[i])
  test.simulation.results = rbind(test.simulation.results, simulation1)
}

test.simulation.results %>% 
  dplyr::select(sample.fraction.used, WindowOpen, WindowClose) %>% 
  pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed))







#climwin need to adjust for 0 - 1, row start at 0 and not 1 
calendar = goingbackpastdayscalendar(refday = 304,
                                     lastdays = 600, 
                                     yearback = 2) %>% 
  filter(YEAR == 1950) %>% 
  dplyr::select(MONTHab, DOY, days.reversed) %>% 
  mutate(datefake  = as.Date(DOY, origin = "1948-01-01"),
         day.month = format(datefake,"%m-%d"), 
         year = case_when(
           days.reversed < 366 ~ 1,
           days.reversed >= 366 & days.reversed < 730 ~ 2,
           days.reversed >= 730 & days.reversed < 1095 ~ 3,
           TRUE ~ 4))

climwin_output_simulation$combos %>% 
  dplyr::select(WindowOpen, WindowClose) %>% 
  pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed)) 


plot.mysim(seed_production)
plot.mysim = function(seed_production){
  # Visualize seed production vs. June week temperature
  ggplot(seed_production, aes(x = mean_temp_june_week, y = log.seed)) +
    geom_point(color = "blue", size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    ggtitle("Seed Production (Perfect Correlation) vs June Week Temperature") +
    xlab("Mean Temperature (June Week)") +
    ylab("Seed Production (Count Data)") +
    theme_minimal()
  
  # Histogram of the seed production counts to visualize the distribution
  #ggplot(seed_production, aes(x = exp(log.seed))) +
  #  geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) 

  print(cor(seed_production$mean_temp_june_week, seed_production$log.seed))
}


climwin_output_simulation <- climwin::slidingwin(
  xvar = list(temperature.degree = climate_data$temp),
  cdate = climate_data$date,
  bdate = seed_production$Date2,
  baseline = lm(seed_production~1, data = seed_production),#i Needed to specify the formula here, if not it is not working properly
  cinterval = 'day',
  range = c(400, 0),
  refday = c(01, 11),
  type = 'absolute',
  stat = "mean",
  cmissing = 'method2',
  func = "lin"
)
climwin_output_simulation

plotall(dataset = climwin_output_simulation[[1]]$Dataset,
        bestmodel = climwin_output_simulation[[1]]$BestModel, 
        bestmodeldata = climwin_output_simulation[[1]]$BestModelData)
