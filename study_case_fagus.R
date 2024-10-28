library(rcompendium)
#add_description(given = 'Valentin', family = 'Journé',
#                email = 'journe.valentin@gmail.com')
#add_license('CC BY 4.0')
#add_readme_rmd(given = 'Valentin', family = 'Journé')
# add_code_of_conduct(email = 'journe.valentin@gmail.com')
# add_citation(
#   given = 'Valentin',
#   family = 'Journé',
#   organisation = 'Forest Biology Center, Adam Mickiewicz University',
#   open = TRUE,
#   overwrite = FALSE,
#   quiet = FALSE
# )

add_dependencies('.')
devtools::load_all(here::here())
devtools::document()
devtools::check()
#let's make my R compendium 
#remotes::install_deps()
#https://frbcesab.github.io/rcompendium/articles/rcompendium.html
library(weatheRcues)
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
library(patchwork)

functions <- list.files(here("fun"), full.names = T) %>%
  purrr::map(source)




new.folder = FALSE 
if(new.folder == TRUE){
  dir.create(here::here('climate_dailyEOBS'))
  dir.create(here::here('figures'))
}

#access climate zip file 
#https://osf.io/287ms/ and add it in climate_dailyEOBS folder

nfile = check_folder_and_contents(here::here('climate_dailyEOBS'), file_pattern = '.qs')
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
#url.mastreev2 = 'https://raw.githubusercontent.com/JJFoest/MASTREEplus/main/Data/MASTREEplus_2024-06-26_V2.csv'

initial.data.mastree <- data.table::fread(url.mastreev2)

#now use the filtering 
#keep cont, remove index, flower and pollen 
Fagus.seed = initial.data.mastree %>%
  filter(!Variable == "flower" & VarType == "C" & !Unit=="index" & !Variable == "pollen") %>% 
  group_by(Species, VarType) %>% 
  mutate(nMax = n()) %>% 
  filter(Species == "Fagus sylvatica") %>% 
  filter(Year > 1952 & Year < 2023) %>% 
  mutate(sitenewname= paste0(Alpha_Number, "_",Site_number, "_", Species_code)) %>% 
  group_by(sitenewname) %>% 
  mutate(log.seed = log(1+Value),
         scale.seed = scale(Value),
         n = n()) %>% 
  filter(n > 19) %>% #initially 14
  mutate() %>% 
  mutate(Date = paste0( "15/06/",Year)) %>% 
  mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
         plotname.lon.lat = paste0("longitude=",Longitude, "_","latitude=", Latitude)) %>% 
  group_by(plotname.lon.lat) %>%   
  ungroup() %>% 
  as.data.frame()

#get the mothod collection
methods.collection.mv2 = Fagus.seed %>% dplyr::select(sitenewname, Country, Collection_method, Length) %>% distinct()


seed.production.plot = Fagus.seed %>% 
  group_by(plotname.lon.lat) %>% 
  complete(Year = full_seq(Year, 1)) %>% 
  drop_na(Collection_method) %>% 
  ggplot(aes(x=Year, y = log.seed, group = plotname.lon.lat, col = Collection_method))+
  geom_line(na.rm = FALSE, alpha = .1)+
  ggpubr::theme_pubr()+
  ylab('Seed production (log)')+
  facet_grid(Collection_method~., scales = 'free')+
  scale_color_futurama()+
  scale_fill_futurama()+
  theme(legend.position = 'none')


#for methods 
forsummary = Fagus.seed %>% dplyr::select(n, sitenewname, Collection_method) %>% distinct()
summary(forsummary$n)
table(forsummary$Collection_method)



plot_locations<- st_as_sf(Fagus.seed %>% 
                            dplyr::select(Longitude, Latitude, plotname.lon.lat, Collection_method) %>% distinct(), coords = c("Longitude", "Latitude")) %>% 
  st_set_crs(4326)

mapBase <- getMap(resolution = "high") %>%   
  st_as_sf() %>% 
  st_make_valid()%>%  st_crop(mapBase, xmin = -12, xmax = 35, ymin = 35, ymax = 70)

MapsMastree <- ggplot(data=mapBase$geometry)+ 
  geom_sf(fill = "grey90", colour = "black", alpha = .8, linewidth = .25)+ #plot map of France
  xlab(" ")+ ylab(" ")+ #white or aliceblue
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
  scale_color_futurama()+
  scale_fill_futurama()+
  guides(col = "none", fill=guide_legend(title=NULL, override.aes = list(size=8)))



cowplot::save_plot(here("figures/Figure1.png"),MapsMastree+seed.production.plot+plot_annotation(tag_levels = 'a') & 
                     theme(plot.tag = element_text(size = 12)), 
                   ncol = 2, nrow = 2, dpi = 300)

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
run.climwin <-F
#take some time ! more than 1 h 
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
            here('outputs/statistics_absolute_climwin.qs'))}else{
  statistics_absolute_climwin = qs::qread('outputs/statistics_absolute_climwin.qs')
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
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed')+
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
  ylim(0,600)+
  scale_color_futurama()+
  scale_fill_futurama()+
  ggpubr::theme_cleveland()+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  ylab('Days reversed')

cowplot::save_plot(here("figures/FigureS1.png"),climwin.all.sites, ncol = 1.5, nrow = 1.8, dpi = 300)

  
quibble2(statistics_absolute_climwin$WindowOpen, q = c(0.25, 0.5, 0.75))
quibble2(statistics_absolute_climwin$WindowClose, q = c(0.25, 0.5, 0.75))



#it is still useful for me because here I will also work with simple correlation later 
results.moving.site = FULL.moving.climate.analysis(seed.data = Fagus.seed,
                             climate.path = climate.beech.path,
                             refday = 305,
                             lastdays = max(range),
                             rollwin = 1)

#create 50 list - 1 per site 
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
    lastdays = max(range),
    rollwin = 1
  ))
qs::qsave(statistics_csp_method, 
          here('outputs/statistics_csp.qs'))




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
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed')+
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_futurama()+
  scale_fill_futurama()+
  ggpubr::theme_cleveland()+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  ylab('Days reversed')
csp.all.sites

cowplot::save_plot(here("figures/csp_allsites_methods.png"),csp.all.sites, ncol = 1.5, nrow = 1.8, dpi = 300)


quibble2(output_fit_summary.best.csp$window.open, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.best.csp$window.close, q = c(0.25, 0.5, 0.75))

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
    refday = 305,
    knots = NULL
  )
)



ggplot(statistics_psr_method %>% 
         group_by(sitenewname) %>% 
         slice(which.max(r2)))+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
  geom_point(aes(y = window.open, x = sitenewname), size = 2) +
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
  ylab('Days reversed')

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
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed')+
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_futurama()+
  scale_fill_futurama()+
  ggpubr::theme_cleveland()+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  ylab('Days reversed')
psr.all.sites

cowplot::save_plot(here("figures/psr_allsites_methods.png"),psr.all.sites, ncol = 1.5, nrow = 1.8, dpi = 300)


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
    seed.data = Fagus.seed, 
    climate.path = climate.beech.path,
    refday = 305,
    lastdays = max(range),
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
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed')+
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
  ylim(0,600)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_futurama()+
  scale_fill_futurama()+
  ggpubr::theme_cleveland()+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  ylab('Days reversed')
basic.all.sites
cowplot::save_plot(here("figures/basic_allsites_methods.png"),basic.all.sites, ncol = 1.5, nrow = 1.5, dpi = 300)

quibble2(output_fit_summary.basic.best$window.open, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.basic.best$window.close, q = c(0.25, 0.5, 0.75))

#put a long lag 
windows.avg = output_fit_summary.basic.best %>% 
  dplyr::select(sitenewname, window.open,  window.close) %>% 
  mutate(method = 'Signal processing') %>% 
  bind_rows(output_fit_summary.psr.best%>% 
              dplyr::select(sitenewname, window.open, window.close) %>% 
              mutate(method = 'P-spline regression')) %>% 
  bind_rows(statistics_absolute_climwin %>% 
              dplyr::select(sitenewname, WindowOpen, WindowClose) %>% rename(window.open = WindowOpen,
                                                                             window.close = WindowClose) %>% 
              mutate(method = 'Sliding time window')) %>% 
  bind_rows(output_fit_summary.best.csp %>% 
              dplyr::select(sitenewname, window.open , window.close) %>% 
              mutate(method = 'Climate sensitivity profile') )%>%
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

library(RColorBrewer)
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
  scale_colour_manual(values = brewer.pal(6, "Paired")[c(5,6)])+
  scale_colour_manual(values = brewer.pal(6, "Paired")[c(5,6)])+
  ggpubr::theme_pubr()+
  theme(legend.position = c(.8,.8),
        legend.title = element_blank())+
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')

#GENERIC PLOT R2 and AIC
bestr2 = output_fit_summary.basic.best %>% 
  dplyr::select(sitenewname, r2, AIC) %>% 
  mutate(method = 'Signal processing') %>% 
  bind_rows(output_fit_summary.psr.best%>% 
              dplyr::select(sitenewname, r2, AIC) %>% 
              mutate(method = 'P-spline regression')) %>% 
  bind_rows(statistics_absolute_climwin %>% 
              dplyr::select(sitenewname, R2, AIC) %>% rename(r2 = R2) %>% 
              mutate(method = 'Sliding time window')) %>% 
  bind_rows(output_fit_summary.best.csp %>% 
              dplyr::select(sitenewname, r2, AIC) %>% 
              mutate(method = 'Climate sensitivity profile') )

mean_r2 <- bestr2 %>%
  group_by(method) %>%
  summarise(mean_r2 = median(r2, na.rm = TRUE))

averagedensr2method = bestr2 %>% 
  ggplot(aes(x = r2, fill = method, col = method)) +
  geom_density(alpha = 0.2, size = .8) +
  labs(x = "R-squared (r2)", y = "Density") +
  #geom_vline(data = mean_r2, aes(xintercept = mean_r2, color = method),
  #           linetype = "dashed", size = .5) +
  ggpubr::theme_pubclean()+
  geomtextpath::geom_textvline(data = mean_r2, 
                 aes(label = round(mean_r2,3), color = method, xintercept = mean_r2),
                 vjust = -0.3,
                 hjust = 1,
                 fontface = 'bold',
                 linetype = 'dashed',show.legend = FALSE)+
  scale_color_viridis_d(name = 'Method')+
  scale_fill_viridis_d(name = 'Method')+
  guides(col=guide_legend(ncol=2))+
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())

averagedensr2method
cowplot::save_plot(here("figures/averager2.method.png"),averagedensr2method+median.windows.plot+plot_annotation(tag_levels = 'a') & 
                     theme(plot.tag = element_text(size = 12)), ncol = 1.8, nrow = 1.4, dpi = 300)

####################################################################################
##############################################################################
#######################
#Block cross validation 
####################################################################################
##############################################################################
#######################
#CLEANED VERSION
#take 1 day (more 12 hours)
num_blocks <- 5
test_block_nb = 2
num_iterations <- 10  # Specify the number of iterations
output_fit_summary_cv_fin_part2 = NULL
# Loop through each site
save.bc = TRUE
for(i in 1:length(beech.site.all)) {
  site.cv <- Fagus.seed %>% 
    filter(plotname.lon.lat == beech.site.all[i])
  
  # Make 5 blocks of data
  blocks <- create_blocks(site.cv, num_blocks)
  
  # Perform cross-validation multiple times
  for (j in 1:num_iterations) {
    block_indices <- seq_along(blocks)  # Should be the same as block size
    test_indices <- sample(block_indices, size = test_block_nb, replace = FALSE)  # Randomly select 2 validation blocks
    train_indices <- setdiff(block_indices, test_indices)  # Remaining blocks for training 
    
    # Combine the block of training
    train_blocks <- do.call(rbind, blocks[train_indices])  # Combine training blocks
    validation_blocks <- do.call(rbind, blocks[test_indices])  # Combine validation blocks
    
    climate_data <- format_climate_data(site = beech.site.all[i],
                                        path = climate.beech.path, 
                                        scale.climate = T)  
    
    # Run climwin per site
    climwin1 <- climwin_site_days(
      climate_data = climate_data,
      data = train_blocks,  
      site.name = beech.site.all[i],  
      range = range,
      cinterval = 'day',
      refday = c(01, 11),  
      optionwindows = 'absolute',  
      climate_var = 'TMEAN'  
    ) %>% 
      dplyr::mutate(method = 'climwin')
    
    psr1 <- runing_psr_site(bio_data = train_blocks,
                            site = beech.site.all[i],
                            climate_csv = climate_data,
                            tot_days = 600,
                            refday = refday,
                            rollwin = 1,
                            covariates.of.interest = 'TMEAN',
                            matrice = c(3, 1),
                            knots = NULL,
                            tolerancedays = 7,
                            plot = TRUE) %>% 
      dplyr::slice(which.max(r2)) %>% 
      dplyr::mutate(method = 'psr')
    
    moving.site1 <- site.moving.climate.analysis(
      bio_data = train_blocks, 
      climate.data = climate_data, 
      lastdays = lastdays,
      refday = 305,
      myform = formula('log.seed ~ rolling_avg_tmean')
    )
    
    csp1 <- runing_csp_site(Results_CSPsub = moving.site1,
                            data = train_blocks,
                            siteneame.forsub = beech.site.all[i],
                            option.check.name = TRUE,
                            climate_csv = climate_data,
                            refday = 305,
                            lastdays = max(range),
                            rollwin = 1,
                            optim.k = F) %>% 
      dplyr::slice(which.max(r2)) %>% 
      dplyr::mutate(method = 'csp')
    
    basic1 <- runing_basic_cues(lag = 100,
                                threshold = 3,
                                influence = 0,
                                tolerancedays = 7,
                                refday = 305,
                                lastdays = max(range),
                                rollwin = 1,
                                siteforsub = beech.site.all[i],
                                climate_csv = climate_data,
                                Results_CSPsub = moving.site1,
                                data = train_blocks) %>% 
      dplyr::slice(which.max(r2)) %>% 
      dplyr::mutate(method = 'signal')
    
    temp.win <- bind_rows(climwin1, psr1, csp1, basic1) %>% 
      dplyr::select(window.open, window.close, 
                    intercept.estimate, slope.estimate, 
                    method)
    
    yearneed <- 2
    yearperiod <- (min(climate_data$year) + yearneed):max(climate_data$year)
    
    # Apply the function across all years in yearperiod and combine results
    rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                        climate = climate_data, 
                                        yearneed = yearneed, 
                                        refday = refday, 
                                        lastdays = lastdays, 
                                        rollwin = 1, 
                                        variablemoving = 'TMEAN')
    
    output_fit_summary.temp <- purrr::map_dfr(1:nrow(temp.win), ~cross_validation_outputs_windows_modelling(., tible.sitelevel = validation_blocks, 
                                                                                            window_ranges_df = temp.win,
                                                                                            rolling.temperature.data = rolling.temperature.data,
                                                                                            myform.fin = formula('log.seed ~ mean.temperature'))) 
    
    output_fit_summary_cv_fin_part1 <- bind_rows(output_fit_summary_cv_fin_part1, output_fit_summary.temp)
  }
  output_fit_summary_cv_fin_part2 = bind_rows(output_fit_summary_cv_fin_part2, output_fit_summary_cv_fin_part1)
}
if(save.bc == T){
  qs::qsave(output_fit_summary_cv_fin_part2, here(paste0('outputs/outputs_blockcrosstotal_',num_blocks,'_train_',test_block_nb,'.qs' )))
}
#qs::qsave(output_fit_summary_cv_fin_part2, here('outputs/output_fit_summary_cv_fin_part2.qs' ))

output_fit_summary_cv_fin_part2 %>% 
  as_tibble() %>% 
  drop_na(mae) %>% 
  left_join(methods.collection.mv2) %>%
  dplyr::mutate(mae = as.numeric(mae),
                sitenewname = fct_reorder(sitenewname, Collection_method)) %>%
  ggplot()+
  geom_point(aes(x=mae,y=sitenewname, col = Collection_method), alpha = .8, shape = 21)+
  facet_grid(.~method.cues)+
  geom_vline(xintercept = 0)+
  scale_color_futurama()+
  scale_fill_futurama()+
  ggpubr::theme_cleveland()+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  ylab('Days reversed')

block.cross.figure = output_fit_summary_cv_fin_part2 %>% 
  left_join(methods.collection.mv2) %>%
  mutate(method.cues = recode(method.cues,
                           "climwin" = "Sliding time window",
                           "csp" = "Climate sensitivity profile",
                           "psr" = "P-spline regression",
                           'signal' = 'Signal processing')) %>% 
  drop_na(mae) %>% 
  ggplot()+
  gghalves::geom_half_boxplot(aes(x= Collection_method, y = mae, fill = Collection_method), center = T,
                    errorbar.draw=FALSE, 
                    width=0.8, nudge = 0.02, alpha = .9)+
  gghalves::geom_half_violin(aes(x= Collection_method, y = mae, fill = Collection_method),side="r", nudge=0.02, alpha = .65, col= "black")+
  facet_grid(.~method.cues)+
  geom_vline(xintercept = 0)+
  scale_color_futurama()+
  scale_fill_futurama()+
  ggpubr::theme_pubclean()+
  geom_hline(yintercept = 0)+
  theme(legend.position = 'none', legend.title = element_blank(), axis.text.x = element_text(angle = 90, size = 12))+
  ylab('Mean Absolute Error')+
  xlab('')
block.cross.figure
cowplot::save_plot(here("figures/Figure.block.png"),block.cross.figure, ncol = 1.6, nrow = 1.4, dpi = 300)

####################################################################################
##############################################################################
#######################
#data sample size effect on window identification
####################################################################################
##############################################################################
#######################
#now test the frac dataset 
#take more than 10 hours if 10 simulations (climwin takes some times)
# Initialize a list to store all results
result_list_all <- list()

# Define the sample fractions you want to iterate over
fractions <- c(0.3, 0.5, 0.7, 0.9)

# Define the number of iterations for each fraction
n_iter <- 10 #take > 12 and < 24 hours on mac os 

# Outer loop to iterate over the fractions
for (fraction in fractions) {
  
  # Initialize a list to store results for each fraction
  result_list <- list()
  
  # Inner loop to perform the operation multiple times for the given fraction
  for (i in 1:n_iter) {
    # Get the sample list for the given fraction
    sample_result <- run_sampling_fraction_all_methods(fraction = fraction,
                                                       range = range, 
                                                       Fagus.seed = Fagus.seed, 
                                                       beech.site.all = beech.site.all, 
                                                       climate.path = climate.beech.path)
    
    # Append this sample's results to the overall list for this fraction
    result_list[[i]] <- sample_result
  }
  
  # Store the results for this fraction in the overall list, using the fraction as the key
  result_list_all[[paste0("fraction_", fraction * 100)]] <- result_list
}


sampling_data_effect <- map_dfr(names(result_list_all), function(sublist) {
  # extract my sublist 
  new_list <- result_list_all[[sublist]]
  new_list[[1]] <- lapply(new_list[[1]], rename_columns_if_needed)  # Apply renaming to the first list
  
  
  # rbind dplyr 
  datatible.sim <- bind_rows(new_list)
  return(datatible.sim)
})

#qs::qsave(sampling_data_effect, here('outputs/sampling_data_effect.qs'))
final_combined_df = qs::qread(here('final_combined_df.qs'))

summary.test.sampling = final_combined_df %>% 
  mutate(method = recode(method,
                              "climwin" = "Sliding time window",
                              "csp" = "Climate sensitivity profile",
                              "psr" = "P-spline regression",
                              'signal' = 'Signal processing')) %>% 
  group_by(method,sample_fraction) %>%
  summarise(wind.open_median = median(window.open, na.rm = TRUE),
            wind.close_median = median(window.close, na.rm = TRUE),
            wind.close_q25 = quantile(window.close, 0.3, na.rm = TRUE),
            wind.close_q75 = quantile(window.close, 0.7, na.rm = TRUE),
            wind.open_q25 = quantile(window.open, 0.3, na.rm = TRUE),
            wind.open_q75 = quantile(window.open, 0.7, na.rm = TRUE),
            wind.open_mean = mean(window.open),
            wind.close_mean = mean(window.close),
            wind.open_se = se(window.open),
            wind.close_se = se(window.close))



sensi.data.plot = summary.test.sampling %>% 
  dplyr::select(method, sample_fraction, wind.open_median:wind.open_q75) %>% 
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>% 
  separate_wider_delim(typewind, delim = '_', names=c('windows.type', 'metric')) %>% 
  pivot_wider(names_from = 'metric', values_from = value) %>% 
  mutate(sample_fraction = as.numeric(sample_fraction),
         method_sample = paste(method, sample_fraction, sep = " (n=") %>% paste0(")"),
         sample_fraction_factor = as_factor(paste0(sample_fraction*100, '%'))) %>% 
  ggplot(aes(x = sample_fraction_factor, group = windows.type, col = windows.type)) +
  geom_point(aes(y = median), size = 2,
             position = position_dodge(width = 0.5)) +
  geom_pointrange(mapping = aes(y = median, ymin = q25, ymax = q75), 
                  position = position_dodge(width = 0.5)) +
  coord_flip() +
  facet_wrap(.~method, scales = 'free_y')+
  xlab('Method (Sample Size)') + 
  ylab('Days reversed') +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5,6)])+
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5,6)])+
  ggpubr::theme_pubr()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())

sensi.data.plot

cowplot::save_plot(here("figures/sensitivity.method.png"),sensi.data.plot, ncol = 1.4, nrow = 1.4, dpi = 300)




########################################################################
####################################################################################
########################################################################
############################################################
######MAKE simulation study case with climwin 

# Load necessary libraries
library(tidyverse)
library(lubridate)


# Simulate climate data: assume daily data for 20 years
set.seed(123)  
startyear = 1940
years <- startyear:2020
days_per_year <- 365
climate_data_simulated <- data.frame(
  date = seq.Date(from = as.Date(paste0(startyear,"-01-01")), by = "day", length.out = length(years) * days_per_year),
  temp = runif(length(years) * days_per_year, min = 10, max = 30)  # random temperatures
)
  
climate_data_simulated <- climate_data_simulated %>%
  mutate(year = format(date, "%Y"),
         day_of_year = yday(date),
         year = as.numeric(year)) %>%
  mutate(across(c(temp), scale)) %>%
  mutate(across(c(temp), as.vector))

# Select ~ one week in June 
june_week <- climate_data_simulated %>%
  filter(day_of_year >= 150 & day_of_year <= 160)


raw.data.param.alpha = statistics_absolute_climwin$intercept.estimate
raw.data.param.beta = statistics_absolute_climwin$slope.estimate
raw.data.param.sigma = statistics_absolute_climwin$sigma


#generate parameter range values for the simulation
setup.param = parameter.range(raw.data.param.alpha,
                           raw.data.param.beta,
                           raw.data.param.sigma,
                           option = 'min.max')
  

alpha.random.value = runif(1, min = setup.param$alpha[1], max = setup.param$alpha[2])
beta.random.value = runif(1, min = setup.param$beta[1], max = setup.param$beta[2])
sigma.random.value = runif(1, min = setup.param$sigma[1], max = setup.param$sigma[2])



seed_production_simulated <- june_week %>%
  group_by(year) %>%
  summarise(mean_temp_june_week = mean(temp)) %>%
  # Standardize the mean temperature for proper correlation, as in the main analysis
  mutate(mean_temp_scaled = scale(mean_temp_june_week, center = TRUE, scale = TRUE)) %>%
  # Generate log.seed with controlled correlation
  mutate(log.seed = alpha.random.value + beta.random.value * mean_temp_scaled + rnorm(n(), mean = 0, sd = sigma.random.value))%>%
  filter(year > startyear)%>% #to select one year less 
  mutate(year = as.numeric(as.character(year))) %>% 
  mutate(Date = paste0( "15/06/",year)) %>% #need a date, even if absolute window
  mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y")) %>% #make it as in the climate date
  sample_frac(sample.fraction)  

reg.summary.simulation = summary(lm(log.seed~mean_temp_june_week, data = seed_production_simulated))
r2sim.before.climwin = reg.summary.simulation$adj.r.squared

simulation.cimwin = climwin::slidingwin(
  xvar = list(temperature.degree = climate_data_simulated$temp),
  cdate = climate_data_simulated$date,
  bdate = seed_production_simulated$Date2,
  baseline = lm(log.seed~1, data = seed_production_simulated),#i Needed to specify the formula here, if not it is not working properly
  cinterval = 'day',
  range = c(600, 0),
  refday = c(01, 11),
  type = 'absolute',
  stat = "mean",
  cmissing = 'method2',
  func = "lin"
)

generate 1000 datasets
then run simulation for the different methods 

object.result = bind_cols(simulation.cimwin$combos,
          r2.before.climwin = r2sim.before.climwin,
          alpha.random.value.setup = alpha.random.value,
          beta.random.value.setup = beta.random.value,
          sigma.random.value.setup = sigma.random.value)


simulate_seed_production <- function(setup.param, 
                                     june_week, 
                                     climate_data_simulated, 
                                     startyear, 
                                     sample.fraction = 1, 
                                     num_simulations = 10) {
  
  # Container to store results
  results_list <- vector("list", num_simulations)
  
  for (i in seq_len(num_simulations)) {
    
    # Randomly sample alpha, beta, and sigma from the provided ranges to the uniform
    alpha.random.value <- runif(1, min = setup.param$alpha[1], max = setup.param$alpha[2])
    beta.random.value <- runif(1, min = setup.param$beta[1], max = setup.param$beta[2])
    sigma.random.value <- runif(1, min = setup.param$sigma[1], max = setup.param$sigma[2])
    
    # Simulate seed production, need same june week data, coming from the same climate data simulated
    seed_production_simulated <- june_week %>%
      group_by(year) %>%
      summarise(mean_temp_june_week = mean(temp)) %>%
      # Standardize the mean temperature for proper correlation, as in the main analysis
      mutate(mean_temp_scaled = scale(mean_temp_june_week, center = TRUE, scale = TRUE)) %>%
      # Generate log.seed with controlled correlation
      mutate(log.seed = alpha.random.value + beta.random.value * mean_temp_scaled + 
               rnorm(n(), mean = 0, sd = sigma.random.value)) %>%
      filter(year > startyear) %>%
      mutate(year = as.numeric(as.character(year))) %>%
      mutate(Date = paste0("15/06/", year)) %>%
      mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y")) %>%
      sample_frac(sample.fraction)  
    
    # Run regression to get R^2 before climwin
    reg.summary.simulation <- summary(lm(log.seed ~ mean_temp_june_week, data = seed_production_simulated))
    r2sim.before.climwin <- reg.summary.simulation$adj.r.squared
    
    # Perform climwin analysis
    simulation.climwin <- climwin::slidingwin(
      xvar = list(temperature.degree = climate_data_simulated$temp),
      cdate = climate_data_simulated$date,
      bdate = seed_production_simulated$Date2,
      baseline = lm(log.seed ~ 1, data = seed_production_simulated),
      cinterval = 'day',
      range = c(600, 0), #here I used the same as we used in the ms 
      refday = c(01, 11), #main text 
      type = 'absolute',
      stat = "mean",
      cmissing = 'method2',
      func = "lin"
    )
    
    # Bind results and store them
    result <- bind_cols(simulation.climwin$combos,
                        r2.before.climwin = r2sim.before.climwin,
                        alpha.random.value.setup = alpha.random.value,
                        beta.random.value.setup = beta.random.value,
                        sigma.random.value.setup = sigma.random.value)
    
    results_list[[i]] <- result
  }
  
  # Combine all results into a single dataframe
  final_results <- bind_rows(results_list)
  
  return(final_results)
}


simulation1 = simulate_seed_production(setup.param, 
                                     june_week, 
                                     climate_data_simulated, 
                                     startyear, 
                                     sample.fraction = 1, 
                                     num_simulations = 10) 
simulation2 = simulate_seed_production(setup.param, 
                                       june_week, 
                                       climate_data_simulated, 
                                       startyear, 
                                       sample.fraction = 1, 
                                       num_simulations = 100) 
simulation3 = simulate_seed_production(setup.param, 
                                       june_week, 
                                       climate_data_simulated, 
                                       startyear, 
                                       sample.fraction = 1, 
                                       num_simulations = 100) 
simulation4 = simulate_seed_production(setup.param, 
                                       june_week, 
                                       climate_data_simulated, 
                                       startyear, 
                                       sample.fraction = 1, 
                                       num_simulations = 100) 

simulation5 = simulate_seed_production(setup.param, 
                                       june_week, 
                                       climate_data_simulated, 
                                       startyear, 
                                       sample.fraction = 1, 
                                       num_simulations = 100) 





simulation6 = simulate_seed_production(setup.param, 
                                       june_week, 
                                       climate_data_simulated, 
                                       startyear, 
                                       sample.fraction = 1, 
                                       num_simulations = 100) 
simulation7 = simulate_seed_production(setup.param, 
                                       june_week, 
                                       climate_data_simulated, 
                                       startyear, 
                                       sample.fraction = 1, 
                                       num_simulations = 100) 

simulation.all = bind_rows(simulation1, simulation2, simulation3, 
          simulation4, simulation5, simulation6,
          simulation7)

simulation.all %>% 
summarise(wind.open_median = median(WindowOpen, na.rm = TRUE),
          wind.close_median = median(WindowClose, na.rm = TRUE),
          wind.close_q25 = quantile(WindowClose, 0.25, na.rm = TRUE),
          wind.close_q75 = quantile(WindowClose, 0.75, na.rm = TRUE),
          wind.open_q25 = quantile(WindowOpen, 0.25, na.rm = TRUE),
          wind.open_q75 = quantile(WindowOpen, 0.75, na.rm = TRUE)) %>% 
  dplyr::select(wind.open_median:wind.open_q75) %>% 
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>% 
  separate_wider_delim(typewind, delim = '_', names=c('windows.type', 'metric')) %>% 
  pivot_wider(names_from = 'metric', values_from = value) %>% 
  left_join(calendar.climwin %>% mutate(median = days.reversed))

hist(simulation.all$r2.before.climwin)
quantile(simulation.all$r2.before.climwin, 0.75, na.rm = TRUE)
quantile(simulation.all$r2.before.climwin, 0.25, na.rm = TRUE)

#find when correct value 
perfectwin = simulation.all %>% filter(between( WindowOpen, 150, 160),
                                       between( WindowClose, 140, 150))
summary(perfectwin$r2.before.climwin)


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
calendar.climwin = goingbackpastdayscalendar(refday = 304,
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
