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


source('utils.R')
source('method_CSP.R')
source('method_climwin.R')
source('method_signal.R')

functions <- list.files(here("fun"), full.names = T) %>%
  purrr::map(source)



new.folder = FALSE 
if(new.folder == TRUE){
  dir.create(here('climate_dailyEOBS'))
  dir.create(here('figures'))
}

nfile = check_folder_and_contents(here('climate_dailyEOBS'), file_pattern = '.qs')
if(nfile$folder_exists == FALSE){stop(print('missing folder climate_dailyEOBS'))}
if(length(nfile$files_present) == 0){stop(print('missing file to add in the folder of interest name climate_dailyEOBS'))}
#SHOULD HAVE NO WARNING HERE

# I am generating calendar 
calendar = goingbackpastdayscalendar(refday = 305,
                                     lastdays = max(range), 
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

initial.data.mastree <- read.csv('/Users/vjourne/Library/CloudStorage/Dropbox/Mastree/MASTREEplus_2022-02-03_V1.csv',stringsAsFactors = F)

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
  sample_frac(0.5) %>%        
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
#option to run in parallel to make it faster
run.climwin <- T
#take some time ! 
if(run.climwin==T){
  library(furrr)#make it in parallele, https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
  plan(multisession)#background R sessions (on current machine) , the computer is flying (alt cluster)

  statistics_absolute_climwin = future_map_dfr(
    beech.site.all, 
    ~ climwin_site_days(
      climate.beech.path = climate.beech.path,
      data = Fagus.seed %>% 
        filter(plotname.lon.lat == .x),  # Use .x correctly here
      site.name = .x,
      range = range,
      cinterval = 'day',
      refday = c(01, 11),
      optionwindows = 'absolute',
      climate_var = 'TMEAN'
    )
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
results.moving.site = FULL.moving.climate.analysis(Fagus.seed = Fagus.seed,
                             climate.beech.path = climate.beech.path,
                             refday = 305,
                             lastdays = max(range),
                             rollwin = 1)


ggplot(results.moving.site)+
  geom_bar(aes(y = r.squared, x = days.reversed), size = 2, stat="identity" ) +
  geom_vline(xintercept = 133, color = 'red')+
  geom_vline(xintercept = 498, color = 'red')



#create 61 list - 1 per site 
Results_CSP = results.moving.site %>%
  arrange(plotname.lon.lat) %>% 
  group_by(plotname.lon.lat) %>%
  group_split()

nameCSP =   results.moving.site %>% 
  arrange(plotname.lon.lat) %>% 
  group_by(plotname.lon.lat) %>%
  group_keys()

# Set names of the list based on group keys
names(Results_CSP) <- apply(nameCSP, 1, paste)

rollwin = 1
lastdays = max(range)
myform.fin = formula('log.seed ~ mean.temperature')
refday = 305

statistics_csp_method = map_dfr(
  1:length(Results_CSP), 
  ~CSP_function_site(
    Results_CSP[[.]], 
    unique(Results_CSP[[.]]$plotname.lon.lat),
    Fagus.seed = Fagus.seed, 
    climate.beech.path = climate.beech.path,
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
    Fagus.seed = Fagus.seed,  # Use .x correctly here
    tot_days = max(range),
    lastdays = 600,
    matrice = c(3,1),
    knots = NULL
  )
)



ggplot(statistics_psr_method %>% 
         group_by(sitenewname) %>% 
         slice(which.max(r2)))+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
  geom_point(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
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
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname, col = Collection_method), size = 2, shape = 15) +
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





output_fit_summary.basic = NULL
#Results_CSP

siteforsub = site[1]


Results_CSPsub = Results_CSP[[i]]




statistics_basic_method = map_dfr(
  1:length(Results_CSP), 
  ~basiccues_function_site(
    Results_CSP[[.]], 
    unique(Results_CSP[[.]]$plotname.lon.lat),
    Fagus.seed = Fagus.seed, 
    climate.beech.path = climate.beech.path,
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
cowplot::save_plot("basic_allsites_methods.png",basic.all.sites, ncol = 1.5, nrow = 1.5, dpi = 300)

quibble2(output_fit_summary.basic.best$window.open, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.basic.best$window.close, q = c(0.25, 0.5, 0.75))

#put a long lag 

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
  theme(legend.position = c(.8, .8))

averagedensr2method
cowplot::save_plot("averager2.method.png",averagedensr2method, ncol = 1, nrow = 1.4, dpi = 300)


#now simple correlaation and extract window size based on the correlation 
# sub.correlation.10days <- results.moving.site %>%
#   group_by(sitenewname) %>%
#   arrange(days.reversed) %>%
#   mutate(
#     max_corr_date = days.reversed[which.max(abs(correlation))]
#   ) %>%
#   filter(days.reversed >= max_corr_date - 10 & days.reversed <= max_corr_date + 10) %>%
#   ungroup() %>%
#   group_by(sitenewname) %>%
#   filter(days.reversed == min(days.reversed) | days.reversed == max(days.reversed)) %>%
#   ungroup()
# 
# #plot with simple correlation, extract top highest correlation , and then 10+- days before highest correlation
# ggplot(sub.correlation.10days %>% dplyr::select(sitenewname, days.reversed) %>% 
#          group_by(sitenewname) %>% 
#          mutate(windows = 1:2) %>% 
#          mutate(windows = ifelse(windows == 1, 'close', 'open')) %>% 
#          pivot_wider(names_from = "windows", values_from = "days.reversed"))+
#   geom_segment(aes(y = open, yend = close, x = sitenewname), size = 2) +
#   coord_flip()+
#   geom_hline(yintercept = 133, color = 'red')+
#   geom_hline(yintercept = 498, color = 'red')+
#   ylim(0,547)
# 
# stat.correlation = sub.correlation.10days %>% dplyr::select(sitenewname, days.reversed) %>% 
#   group_by(sitenewname) %>% 
#   mutate(windows = 1:2) %>% 
#   mutate(windows = ifelse(windows == 1, 'close', 'open')) %>% 
#   pivot_wider(names_from = "windows", values_from = "days.reversed")
# quibble2(stat.correlation$open, q = c(0.25, 0.5, 0.75))
# quibble2(stat.correlation$close, q = c(0.25, 0.5, 0.75))


# 
# quantile_threshold <- com.both %>%
#   group_by(Partition) %>%
#   summarize(combined_score_threshold = quantile(combined_score, 0.9, na.rm = TRUE))
# 
# 
# days_combinationr2slope_threshold = com.both %>%
#   left_join(quantile_threshold, by = "Partition") %>%  # Join the threshold back to the original data
#   filter(combined_score >= combined_score_threshold) %>%    # Keep only rows where predicted slope meets/exceeds the threshold
#   group_by(Partition) %>%
#   summarize(min.day = min(day, na.rm = TRUE),
#             max.day = max(day, na.rm = TRUE)) 
# 
# com.both <- com.both %>%
#   group_by(Partition) %>% 
#   mutate(cluster = kmeans(combined_score, centers = 5)$cluster)
# 
# (kmeans(com.both$combined_score, centers = 2))
# 
# days_combinationr2slope_threshold <- com.both %>%
#   left_join(quantile_threshold, by = "Partition") %>%  
#   filter(combined_score >= combined_score_threshold) %>%  # Filter based on the combined score threshold
#   group_by(Partition, cluster) %>%  # Group by Partition and cluster
#   summarize(min.day = min(day, na.rm = TRUE),
#             max.day = max(day, na.rm = TRUE),
#             .groups = 'drop')

#testData <- temporary[-partitions[[i]], ]
#rf_model_slope[[i]] <- randomForest(slope ~ day, data = trainData[[i]])
#rf_model_r2[[i]] <- randomForest(r_s ~ day, data = trainData[[i]])
#predict both for slope and r2 - use my combine to combine the different random forest models 
# rf.all.slope <- my_combine(rf_model_slope[[1]], rf_model_slope[[2]], rf_model_slope[[3]], rf_model_slope[[4]] , rf_model_slope[[5]],
#                      rf_model_slope[[6]], rf_model_slope[[7]], rf_model_slope[[8]], rf_model_slope[[9]] , rf_model_slope[[10]])
#   
# pred_slope = predict(rf.all.slope, newdata = temporary)
#   
# rf.all.slope <- my_combine(rf_model_r2[[1]], rf_model_r2[[2]], rf_model_r2[[3]], rf_model_r2[[4]] , rf_model_r2[[5]],
#                              rf_model_r2[[6]], rf_model_r2[[7]], rf_model_r2[[8]], rf_model_r2[[9]] , rf_model_r2[[10]])
# pred_r2 = predict(rf.all.slope, newdata = temporary)
# 
# predictions_slope <- data.frame(day = temporary$day, 
#                                      actual_slope = temporary$slope, 
#                                      predicted_slope = pred_slope)
# 
# predictions_r2 <- data.frame(day = temporary$day, 
#                                   actual_r2 = temporary$r_s, 
#                                   predicted_r2 = pred_r2)

# Combine predictions from all partitions for slope
#all_predictions_slope <- bind_rows(predictions_slope, .id = "Partition")%>% as_tibble()

# Combine predictions from all partitions for r2
#all_predictions_r2 <- bind_rows(predictions_r2, .id = "Partition") %>% as_tibble()
#com = list(predictions_slope, predictions_r2)

# Reformat the climate data
# climate_csv <- read_csv(climate_beech_unique, col_types = cols(...1 = col_skip())) %>%
#   mutate(Tmean = base::scale(Tmean),
#          Prp = base::scale(Prp),
#          across(c(Tmean, Prp), as.vector),
#          date = as.Date(paste0(day,'/',month,'/', year), format = "%d/%m/%Y"),
#          yday = lubridate::yday(date))

#now for me try another alternaive - I want to identify huge peak of correlation 
#the idea here is that we have seasonal changes of the correlation
#i just want to identify the breakpoints instead of using a psecifc 10 days b-a windows 
#install.packages("bfast", repos="http://R-Forge.R-project.org")
# library(bfast)
# set_default_options()    # use modifications, same as set_fast_options()
# set_fallback_options()   # use default implementation
# NDVIa <- as.ts(zoo(som$NDVI.a, som$Time))
# f <- function() bfastmonitor(NDVIa, start = c(2010, 13)) 
# 
# set_fallback_options()
# x = f() 
# system.time(replicate(100, f()))
# 
# set_fast_options()
# y = f()
# system.time(replicate(100, f()))
# 
# par(mfrow = c(1,2))
# plot(x) ; plot(y)
# 
# summary(x$model)
# x$mefp
# x$breakpoint
# x$magnitude
# x$history
# x$monitor
# outdir = "/Users/vjourne/Documents/GITprojects/weather_cues_method"
# rmarkdown::render(system.file("reports/report.test.Rmd",package = "bfast"),output_file = paste(outdir,"/report.test.html",sep=""))
# 
# test = results.moving.site %>% 
#   dplyr::filter(sitenewname == '413_1_FAGSYL') %>% 
#   mutate(year = ifelse(days.reversed<365, 0, 1)) %>% 
#   group_by(year) %>% 
#   mutate(days = ifelse(year == 0, days.reversed, row_number())) %>% 
#   ungroup() %>% 
#   mutate(time = as.double(paste0(year, '.', sprintf("%03d", days))),
#          timebis = as.double(paste0(0, '.', sprintf("%03d", days.reversed)))) %>% 
#   as.data.frame()
# 
# suba <- as.ts(zoo(test$correlation, test$time))
# f <- function() bfastmonitor(suba, 
#                              #start = c(1.110),
#                              start = c(0, 364)) 
# set_fast_options()
# y = f()
# system.time(replicate(100, f()))
# plot(y)
#combine different random forest model based on different train dataset 

# optimization.ml <- function(data,
#                             num_partitions = 10,
#                             partition_size = 0.8) {
#   
#   # Train Random Forest models
#   #https://stackoverflow.com/questions/19170130/combining-random-forests-built-with-different-training-sets-in-r
#   #https://stackoverflow.com/questions/70842819/is-there-an-r-function-in-the-ranger-package-that-can-combine-two-random-forest
#   library(caret) #for parititon 
#   library(ranger) #for rf (it takes weights in it eventually )
#   
#   set.seed(9999)
#   #create x partition 
#   partitions <- list()
#   #do not specify list T in partitioon 
#   for (i in 1:num_partitions) {
#     trainIndex <- createDataPartition(data$day, p = partition_size, times = 1, list = F)
#     partitions[[i]] <- trainIndex
#   }
#   # Train the model on each partition and make predictions
#   trainData = list()
#   #testData = list()
#   #rf_model_slope = list()
#   rf_model = list()
#   predictions.ranger.data=NULL
#   
#   for (i in 1:num_partitions) {
#     # Subset data into training and testing sets
#     trainData[[i]] <- data[partitions[[i]], ]
#     rf_model[[i]] <- ranger::ranger(
#       #formula = slope ~ day, #this one max slope and take into account r2 as weights 
#       formula = r_s ~ day, #i ends up with this, just because it maximise r2 
#       data = trainData[[i]],
#       #case.weights = trainData[[i]]$r_s,    # Include R2 as weights
#       importance = "impurity"    
#     )
#     #do not forget to make prediction for the full dataset ! 
#     predictions.ranger.data.temp = data.frame(sitenewname = unique(data$sitenewname),
#                                               day = data$day, 
#                                               actual_slope = data$slope, 
#                                               actual_r2 = data$r_s,
#                                               predicted_slope = predict(rf_model[[i]], data = data)$predictions,
#                                               iter = i)
#     predictions.ranger.data = rbind(predictions.ranger.data, predictions.ranger.data.temp)
#     
#   }
#   
#   #use ranger to make prediction on the full sampel   
#   #have same actual r2 and day and slope 
#   test.pred = predictions.ranger.data %>% 
#     group_by(sitenewname, day, actual_r2, actual_slope) %>% 
#     summarise(pred.avg = mean(predicted_slope))
#   
#   return(test.pred)
# }
# 
# test.pred = optimization.ml(temporary,
#                             num_partitions = 20,
#                             partition_size = 0.8)

# ggplot(test.pred, aes(x = day)) +
#   geom_line(aes(y = actual_r2, color = "Actual Slope"), size = 1, linetype = "dashed") +
#   geom_line(aes(y = pred.avg), size = 1) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.position = 'none')
# 
# 
# 
# find.windows.optimal.valleys = function(data.prediction, #the one use for windows identification 
#                                         rolling.temperature.data,#the data with rolling 
#                                         quantilehigh = 0.9,
#                                         minpeakdistance = 15,
#                                         plot = FALSE){
#   
#   if(unique(data.prediction$sitenewname) != unique(rolling.temperature.data$sitenewname)){
#     stop("Error: Not the same name between the prediction and the formatted data :(")
#   }
#   
#   library(pracma)
#   #use r pracma
#   data_for_peaks <- data.prediction %>%
#     ungroup() %>% 
#     dplyr::select(day, pred.avg)
#   #find peaks with pracma fuction 
#   peaks <- findpeaks(abs(data_for_peaks$pred.avg), 
#                      minpeakheight = quantile(abs(data_for_peaks$pred.avg), quantilehigh),
#                      minpeakdistance = minpeakdistance) 
#   peaks_df <- data_for_peaks[peaks[,2],] %>%
#     mutate(peak_id = row_number())%>% 
#     mutate(type = 'peak') 
#   WinCloValley = data_for_peaks[peaks[,3],] %>%
#     mutate(id = row_number())%>% 
#     rename(windowsclose = 1) %>% 
#     dplyr::select(windowsclose, id)
#   WinOpValley = data_for_peaks[peaks[,4],] %>%
#     mutate(id = row_number()) %>% 
#     rename(windowsopen = 1,) %>% 
#     dplyr::select(windowsopen, id)
#   
#   window_ranges_df = left_join(WinCloValley, WinOpValley) %>% 
#     rename(windows.sequences.number = id)
#   
#   output_fit_summary.new = NULL
#   for (z in 1:nrow(window_ranges_df)) {
#     windowsopen <- window_ranges_df$windowsopen[z]
#     windowsclose <- window_ranges_df$windowsclose[z]
#     window_number <- window_ranges_df$windows.sequences.number[z]
#     
#     # Filter the rolling temperature data according to the current window range
#     #does not change anything if sign different, just important the date smaller and bigger than
#     climate_windows_best <- rolling.temperature.data %>%
#       dplyr::filter(days.reversed <= windowsopen & days.reversed >= windowsclose) %>%
#       group_by(sitenewname, year) %>%
#       summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = TRUE)) %>%
#       ungroup() %>%
#       mutate(window_number = window_number)  # Add the window number for identification
#     
#     #need to change year colnames :( 
#     fit.best.model = bio_data %>% 
#       rename(year = Year) %>% 
#       left_join(climate_windows_best) %>%
#       lm(myform.fin, data = .)
#     
#     intercept = tidy(fit.best.model)$estimate[1]
#     intercept.se = tidy(fit.best.model)$std.error[1]
#     estimate.model = tidy(fit.best.model)$estimate[2]
#     estimate.se.model = tidy(fit.best.model)$std.error[2]
#     pvalue.model = tidy(fit.best.model)$p.value[2]
#     r2 = glance(fit.best.model)$r.squared
#     AIC = glance(fit.best.model)$AIC
#     nobs = glance(fit.best.model)$nobs
#     nsequence.id = z
#     
#     output_fit_summary.temp.new <- data.frame(sitenewname = unique(bio_data$sitenewname),
#                                               reference.day = refday,
#                                               windows.size = rollwin,
#                                               window.open = windowsopen,
#                                               window.close = windowsclose,
#                                               intercept = intercept,
#                                               intercept.se = intercept.se,
#                                               estimate = estimate.model,
#                                               estimate.se = estimate.se.model,
#                                               pvalue = pvalue.model,
#                                               r2 = r2,
#                                               AIC = AIC,
#                                               nobs = nobs,
#                                               nsequence.id = nsequence.id)
#     output_fit_summary.new = rbind(output_fit_summary.new, output_fit_summary.temp.new)
#     
#   }
#   if(plot == T){
#     plott = ggplot(data_for_peaks, aes(x = day)) +
#       geom_line(aes(y = pred.avg), size = 1) +
#       theme_minimal() +
#       theme(plot.title = element_text(hjust = 0.5),
#             legend.position = 'none')+
#       geom_vline(xintercept = c(peaks_df$day), col = 'red')+
#       geom_vline(xintercept = c(WinCloValley$windowsclose), col = 'green')+
#       geom_vline(xintercept = c(WinOpValley$windowsopen), col = 'blue')
#     print(plott)
#   }
#   return(output_fit_summary.new)
# }
# 
# 
# 
# ttt = find.windows.optimal.valleys(test.pred, #the one use for windows identification 
#                                    rolling.temperature.data,#the data with rolling 
#                                    quantilehigh = 0.9,
#                                    minpeakdistance = 15, 
#                                    plot = TRUE)
# 
# # 

# for(i in 1:length(Results_CSP)){
#   #subset the data, 
#   #scale seed production and do logit transformation (will use this logit for climwin later)
#   tible.sitelevel = Fagus.seed %>% #site = bio_data 
#     dplyr::filter(plotname.lon.lat == unique(Results_CSP[[i]]$plotname.lon.lat))
#   
#   
#   #now make subset and run regression
#   list_slope <- as.list(Results_CSP[[i]]$estimate)
#   list_rs <- as.list(Results_CSP[[i]]$r.squared)
#   
#   day = seq(rollwin,(lastdays-1),1)
#   slope = list_slope
#   r_s = list_rs
#   #k = nrow(site)
#   temporary <- data.frame(slope = unlist(slope), 
#                           day = day, 
#                           r_s = unlist(r_s))
#   #use the function to optimize k 
#   results <- optimize_and_fit_gam(temporary, optim.k = F, plots = F, k = (nrow(tible.sitelevel)-1) ) #(12+6) #I will specify just number of month 
#   #one year and 6 month for 1.5 years 
#   #get days windows 
#   #need to add rollwin because it start at 1... 
#   days = get_predictions_windows(slope_gam = results[[1]], 
#                                  rs_gam = results[[2]], 
#                                  temporary)$days
#   
#   sequences_days = extract_consecutive_sequences(days, keep_all = TRUE)
#   
#   
#   
#   climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = unique(Results_CSP[[i]]$plotname.lon.lat))
#   
#   
#   climate_csv <- qs::qread(climate_beech_unique) %>%
#     as.data.frame() %>%
#     mutate(DATEB = as.Date(DATEB, format = "%m/%d/%y")) %>%
#     mutate(date = foo(DATEB, 1949)) %>%
#     mutate(yday = lubridate::yday(date)) %>%
#     mutate(year = as.numeric(str_sub(as.character(date),1, 4)) ) %>%
#     mutate(TMEAN = as.numeric(TMEAN),
#            TMAX = as.numeric(TMAX),
#            TMIN = as.numeric(TMIN),
#            PRP = as.numeric(PRP)) %>%
#     mutate(across(c(TMEAN, TMAX, TMIN, PRP), scale))%>% 
#     mutate(across(c(TMEAN, TMAX, TMIN, PRP), as.vector))
#   
#   # Define the year period
#   yearneed <- 2
#   yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
#   
#   
#   # Apply the function across all years in yearperiod and combine results
#   #here the map dfr will basically do same as aplly by runing the function over all the time period we want 
#   rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
#                                       climate = climate_csv, 
#                                       yearneed = yearneed, 
#                                       refday = refday, 
#                                       lastdays = lastdays, 
#                                       rollwin = rollwin, 
#                                       variablemoving = 'TMEAN')
#   
#   window_ranges_df <- save_window_ranges(sequences_days) %>% 
#     mutate(windows.sequences.number = 1:nrow(.))
#   
#   
#   # Loop through each window range in window_ranges_df
#   #i need to do this because sometimes the algo identify more sequences of days 
#   # Replacing the for loop with map_dfr for optimization
#   
#   output_fit_summary.temp <- map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = tible.sitelevel, 
#                                                                                            rolling.temperature.data = rolling.temperature.data))
#   output_fit_summary = rbind(output_fit_summary, output_fit_summary.temp)
# }
# 

# t = NULL
# for(i in 1:length(Results_CSP)){#
#   #Results_CSP
#   #i=20
#   
#   Results_CSPsub = Results_CSP[[i]]
#   siteneame.forsub = unique(Results_CSP[[i]]$plotname.lon.lat)
#   
#   data.sub.fagus = Fagus.seed %>% #site = bio_data 
#     dplyr::filter(plotname.lon.lat == siteneame.forsub)
#   
#   ttt = runing_csp_site(Results_CSPsub = Results_CSPsub,
#                         data = data.sub.fagus,
#                         siteneame.forsub = siteneame.forsub,
#                         climate.beech.path = climate.beech.path,
#                         refday = 305,
#                         lastdays = max(range),
#                         rollwin = 1)
#   t = rbind(t, ttt)
# }
# 
# for(p in 1:length(unique(Fagus.seed$sitenewname))){
#   site = unique(Fagus.seed$plotname.lon.lat)[p]
#   bio_data = Fagus.seed %>% 
#     dplyr::filter(plotname.lon.lat == site )
#   
#   #bio_data$seeds <- bio_data$log.seed
#   #climate 
#   climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern =site)
#   
#   
#   climate_csv <- qs::qread(climate_beech_unique) %>%
#     as.data.frame() %>%
#     mutate(DATEB = as.Date(DATEB, format = "%m/%d/%y")) %>%
#     mutate(date = foo(DATEB, 1949)) %>%
#     mutate(yday = lubridate::yday(date)) %>%
#     mutate(year = as.numeric(str_sub(as.character(date),1, 4)) ) %>%
#     mutate(TMEAN = as.numeric(TMEAN),
#            TMAX = as.numeric(TMAX),
#            TMIN = as.numeric(TMIN),
#            PRP = as.numeric(PRP)) %>%
#     mutate(across(c(TMEAN, TMAX, TMIN, PRP), scale))%>% 
#     mutate(across(c(TMEAN, TMAX, TMIN, PRP), as.vector))
#   
#   # need climate data to be arranged with year as row
#   # need to reduce climate dataframe to only year, yday and temp
#   climate2 <- data.frame(year = climate_csv$year, yday = climate_csv$yday, temp = climate_csv[,covariates.of.interest])
#   tempmat <- climate2 %>% spread(, key = yday, value = temp)
#   tempmat <- tempmat[,-1]
#   #number years monitoring seeds 
#   ny<-length(bio_data$Year)
#   nt<-tot_days-1
#   ## Formatting data
#   #based on Simmonds et al code
#   index.matrix=matrix(1:nt,ny,nt,byrow=TRUE)
#   # Define the year period
#   yearneed <- 2
#   #will fiter the year period here to the year needeed 
#   yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
#   
#   # Apply the function across all years in yearperiod and combine results
#   rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
#                                       climate = climate_csv, 
#                                       yearneed = yearneed, 
#                                       refday = 305, 
#                                       lastdays = tot_days, 
#                                       rollwin = 1, 
#                                       variablemoving = covariates.of.interest)
#   
#   #merge data seed to moving climate
#   tible.sitelevel = bio_data %>% #site = bio_data 
#     rename(year = Year) %>% 
#     left_join(rolling.temperature.data) %>% 
#     drop_na(!!sym('rolling_avg_tmean'))
#   
#   climate2 <- data.frame(year = tible.sitelevel$year, 
#                          yday = tible.sitelevel$days.reversed, 
#                          temp = tible.sitelevel$rolling_avg_tmean)
#   covariate.matrix = climate2 %>%
#     spread(key = yday, value = temp) %>%
#     dplyr::select(-year) %>%                    
#     as.matrix()    
#   covariate.matrix <- unname(as.matrix(covariate.matrix))
#   
#   #second order, but from https://link.springer.com/article/10.1007/s00484-007-0141-4
#   #it is possible to adjust between 1-3 (instead of 2) like c(1,1) instead of c(2,1)
#   #here I specify to 0 instead of 1 as in SImmonds et al, meaning that I have no penaties on slope
#   #if using 1 instead of 0, the model is not able to identify a cues 
#   #and I think it is normal and ok to be more flexible because we might have multiple peaks with time 
#   #So if I want to specify K = 30, for example it is working 
#   #not if knots are too low 
#   #which make sense when looking https://link.springer.com/article/10.1007/s00484-011-0472-z 
#   if(is.null(knots) ){
#     K = ny-1
#     model<-mgcv::gam(log.seed~s(index.matrix,k=K,m=matrice,bs="ps",by=covariate.matrix), 
#                      data = bio_data, method="GCV.Cp")
#     summary(model)
#   }else{
#     model<-gam(log.seed~s(index.matrix,k=K,m=c(2,1),bs="ps",by=covariate.matrix), 
#                data = bio_data, method="GCV.Cp")
#     summary(model)}
#   
#   
#   plotted <- plot(model, ylab = c('partial effect'), xlab = c('days prior tree growth'))
#   coefs <- data.frame(fit = plotted[[1]]$fit)
#   #"The most important days of temperature influence on phenology in the year may be identified by as those 
#   #with partial coefficients greater than or less than zero by more than two times their standard error."
#   #I corrected here, rounding values create issues 
#   #marker <- c(which(coefs$fit > round((1.96*sd(coefs$fit))+mean(coefs$fit),2)), 
#   #            which(coefs$fit < round((1.96*sd(coefs$fit))-mean(coefs$fit),2)))
#   # marker <- c(
#   #    which(coefs$fit > (1.96 * sd(coefs$fit) + mean(coefs$fit))), 
#   #    which(coefs$fit < (1.96 * sd(coefs$fit) - mean(coefs$fit)))
#   #  )
#   
#   #wihtout rounding now, will provide different output, and adjust by what is in the main study . 
#   #In the case of Simmods et al, was certainly ok, but not here, since obs are not days////
#   upper_limit <- mean(coefs$fit) + (1.96 * sd(coefs$fit))
#   lower_limit <- mean(coefs$fit) - (1.96 * sd(coefs$fit))
#   
#   # Find the markers (days) where the fit exceeds the threshold
#   marker <- which(coefs$fit > upper_limit | coefs$fit < lower_limit)
#   
#   if(length(marker) == 0){
#     output_fit_summary.psr.temp = data.frame(sitenewname = unique(bio_data$sitenewname),
#                                              plotname.lon.lat = unique(bio_data$plotname.lon.lat),
#                                              reference.day = refday,
#                                              windows.size = NA, 
#                                              window.open = NA, window.close = NA, intercept = NA,
#                                              intercept.se = NA, estimate = NA, estimate.se = NA,
#                                              pvalue = NA, r2 = NA, AIC = NA, nobs = ny,
#                                              nsequence.id = NA)
#   }else{
#     window <- round(plotted[[1]]$x[find_concurrent_period(marker, coefs)])
#     
#     #here does not work because not consecutive 
#     #extract_consecutive_sequences(window, keep_all = T) 
#     sequences_days = extract_sequences_auto(window, tolerance = tolerancedays) 
#     #i guess now will do the same shit as the other methods 
#     window_ranges_df <- save_window_ranges(sequences_days) %>% 
#       mutate(windows.sequences.number = 1:nrow(.))
#     
#     output_fit_summary.psr.temp <- map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = bio_data, 
#                                                                                                  window_ranges_df = window_ranges_df,
#                                                                                                  rolling.temperature.data = rolling.temperature.data,
#                                                                                                  myform.fin = formula('log.seed ~ mean.temperature')))
#     
#   }
#   
#   output_fit_summary.psr = rbind(output_fit_summary.psr, output_fit_summary.psr.temp)
#   
#   
# }
# 

# for(i in 1:length(Results_CSP)){
#   
#   site.name = unique(Results_CSP[[i]]$plotname.lon.lat)
#   
#   #first reformat data 
#   bio_data <- Fagus.seed %>%
#     filter(plotname.lon.lat == site.name) %>%
#     as.data.frame() 
#   
#   # Climate load and format
#   climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = site.name)
#   # Reformat the climate data
#   climate_csv <- qs::qread(climate_beech_unique) %>%
#     as.data.frame() %>%
#     mutate(DATEB = as.Date(DATEB, format = "%m/%d/%y")) %>%
#     mutate(date = foo(DATEB, 1949)) %>%
#     mutate(yday = lubridate::yday(date)) %>%
#     mutate(year = as.numeric(str_sub(as.character(date),1, 4)) ) %>%
#     mutate(TMEAN = as.numeric(TMEAN),
#            TMAX = as.numeric(TMAX),
#            TMIN = as.numeric(TMIN),
#            PRP = as.numeric(PRP)) %>%
#     mutate(across(c(TMEAN, TMAX, TMIN, PRP), scale))%>% 
#     mutate(across(c(TMEAN, TMAX, TMIN, PRP), as.vector))
#   
#   # Define the year period
#   yearneed <- 2
#   yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
#   
#   # Apply the function across all years in yearperiod and combine results
#   rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
#                                       climate = climate_csv, 
#                                       yearneed = yearneed, 
#                                       refday = 305, 
#                                       lastdays = lastdays, 
#                                       rollwin = 1, 
#                                       variablemoving = 'TMEAN')
#   
#   #now load data related to lm ~ days ~ provide r2 and slope 
#   
#   
#   result <- ThresholdingAlgo(y,lag,threshold,influence)
#   
#   #now made them together
#   formatzek = data.frame(day = seq(rollwin,(lastdays-1),1),
#                          signal = replace_0((result$signals), threshold = tolerancedays)) %>% 
#     dplyr::filter(signal !=0)
#   
#   window.basic = as.vector(formatzek$day)
#   
#   #tolerance of 20 days gaps between 0 
#   sequences_days = extract_sequences_auto(window.basic, tolerance = tolerancedays) 
#   #i guess now will do the same shit as the other methods 
#   window_ranges_df <- save_window_ranges(sequences_days) %>% 
#     mutate(windows.sequences.number = 1:nrow(.))
#   
#   output_fit_summary.temp.basic <- map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = bio_data, 
#                                                                                                  window_ranges_df = window_ranges_df,
#                                                                                                  rolling.temperature.data = rolling.temperature.data,
#                                                                                                  myform.fin = formula('log.seed ~ mean.temperature')))
#   
#   output_fit_summary.basic = rbind(output_fit_summary.basic, output_fit_summary.temp.basic)
#   
#   
# }
