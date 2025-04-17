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
library(here)
library(tidyverse)
#https://community.rstudio.com/t/unable-to-install-packages-from-github/124372
#Sys.unsetenv("GITHUB_PAT")

library(climwin)
library(broom)
library(mgcv)
library(here)
#for maps
library(sf)
library("rworldmap")
library("rworldxtra")
#require(rgdal)
library(ggspatial)
library(ggsci)
library(patchwork)

#functions <- list.files(here("fun"), full.names = T) %>%
#  purrr::map(source)

new.folder = FALSE
if (new.folder == TRUE) {
  dir.create(here::here('climate_dailyEOBS'))
  dir.create(here::here('figures'))
}

#access climate zip file
#https://osf.io/287ms/ and add it in climate_dailyEOBS folder

nfile = check_folder_and_contents(
  here::here('climate_dailyEOBS'),
  file_pattern = '.qs'
)
if (nfile$folder_exists == FALSE) {
  stop(print('missing folder climate_dailyEOBS'))
}
if (length(nfile$files_present) == 0) {
  stop(print(
    'missing file to add in the folder of interest name climate_dailyEOBS'
  ))
}
#SHOULD HAVE NO WARNING HERE

# I am generating calendar
#double check values
calendar = goingbackpastdayscalendar(
  refday = 305, #double check later, because climwin start at 1 or 0
  lastdays = 600,
  yearback = 2
) %>%
  filter(YEAR == 1950) %>%
  dplyr::select(MONTHab, DOY, days.reversed) %>%
  mutate(
    datefake = as.Date(DOY, origin = "1948-01-01"),
    day.month = format(datefake, "%m-%d"),
    year = case_when(
      days.reversed < 366 ~ 1,
      days.reversed >= 366 & days.reversed < 730 ~ 2,
      days.reversed >= 730 & days.reversed < 1095 ~ 3,
      TRUE ~ 4
    )
  )

initial.data.mastree <- read.csv(
  '/Users/valentinjourne/Dropbox/Mastree/MASTREEplus_2024-06-26_V2.csv',
  stringsAsFactors = F
)

#load data from JJ Foest, last V2 version
#url.mastreev2 = 'https://raw.githubusercontent.com/JJFoest/MASTREEplus/main/Data/MASTREEplus_2024-06-26_V2.csv'

#initial.data.mastree <- data.table::fread(url.mastreev2)

#now use the filtering
#keep cont, remove index, flower and pollen
Fagus.seed = initial.data.mastree %>%
  filter(
    !Variable == "flower" &
      VarType == "C" &
      !Unit == "index" &
      !Variable == "pollen"
  ) %>%
  group_by(Species, VarType) %>%
  mutate(nMax = n()) %>%
  filter(Species == "Fagus sylvatica") %>%
  filter(Year > 1952 & Year < 2023) %>%
  mutate(
    sitenewname = paste0(Alpha_Number, "_", Site_number, "_", Species_code)
  ) %>%
  group_by(sitenewname) %>%
  mutate(
    log.seed = log(1 + Value),
    scale.seed = scale(Value),
    n = n(),
    scaling01 = (Value - min(Value, na.rm = T)) /
      (max(Value) - min(Value, na.rm = T)),
    scalingbeta = y.transf.betareg(scaling01)
  ) %>%
  filter(n > 19) %>% #initially 14
  mutate() %>%
  mutate(Date = paste0("15/06/", Year)) %>%
  mutate(
    Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
    plotname.lon.lat = paste0(
      "longitude=",
      Longitude,
      "_",
      "latitude=",
      Latitude
    )
  ) %>%
  group_by(plotname.lon.lat) %>%
  ungroup() %>%
  as.data.frame()

#get the mothod collection
methods.collection.mv2 = Fagus.seed %>%
  dplyr::select(sitenewname, Country, Collection_method, Length) %>%
  distinct()


seed.production.plot = Fagus.seed %>%
  group_by(plotname.lon.lat) %>%
  complete(Year = full_seq(Year, 1)) %>%
  drop_na(Collection_method) %>%
  ggplot(aes(
    x = Year,
    y = log.seed,
    group = plotname.lon.lat,
    col = Collection_method
  )) +
  geom_line(na.rm = FALSE, alpha = .1) +
  ggpubr::theme_pubr() +
  ylab('Seed production (log)') +
  facet_grid(Collection_method ~ ., scales = 'free') +
  scale_color_futurama() +
  scale_fill_futurama() +
  theme(legend.position = 'none')


#for methods
forsummary = Fagus.seed %>%
  dplyr::select(n, sitenewname, Collection_method) %>%
  distinct()
summary(forsummary$n)
table(forsummary$Collection_method)


plot_locations <- st_as_sf(
  Fagus.seed %>%
    dplyr::select(Longitude, Latitude, plotname.lon.lat, Collection_method) %>%
    distinct(),
  coords = c("Longitude", "Latitude")
) %>%
  st_set_crs(4326)

mapBase <- getMap(resolution = "high") %>%
  st_as_sf() %>%
  st_make_valid() %>%
  st_crop(mapBase, xmin = -12, xmax = 35, ymin = 35, ymax = 70)

MapsMastree <- ggplot(data = mapBase$geometry) +
  geom_sf(fill = "grey90", colour = "black", alpha = .8, linewidth = .25) + #plot map of France
  xlab(" ") +
  ylab(" ") + #white or aliceblue
  geom_point(
    data = plot_locations %>%
      mutate(x = unlist(map(geometry, 1)), y = unlist(map(geometry, 2))) %>%
      as.data.frame(),
    aes(x = x, y = y, col = Collection_method, fill = Collection_method),
    stroke = 1,
    alpha = .8,
    shape = 21,
    size = 4.3
  ) + #size = rmse,
  scale_size_continuous(
    breaks = c(0.1, 0.3, 0.8),
    range = c(0, 7)
  ) +
  coord_sf(xlim = c(-10, 25), ylim = c(40, 60)) +
  scale_y_continuous(breaks = c(40, 50, 60, 70)) +
  scale_x_continuous(breaks = c(-10, 0, 10, 20)) +
  ylab("Latitude") +
  xlab("Longitude") +
  ggpubr::theme_pubr() +
  scale_color_futurama() +
  scale_fill_futurama() +
  guides(
    col = "none",
    fill = guide_legend(title = NULL, override.aes = list(size = 8))
  )


cowplot::save_plot(
  here("figures/Figure1.png"),
  MapsMastree +
    seed.production.plot +
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(size = 12)),
  ncol = 2,
  nrow = 2,
  dpi = 300
)

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
run.climwin <- T
#take some time ! more than 1 h
if (run.climwin == T) {
  library(furrr) #make it in parallele, https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
  #plan(sequential)#change from multisession because of issues , and it came from the data vector.. so better using sequential and got issue with qs argument
  plan(multisession)

  statistics_absolute_climwin <- future_map_dfr(
    beech.site.all,
    ~ {
      library(weatheRcues)
      # format the climate, .x is the current site
      climate_data <- weatheRcues:::format_climate_data(
        site = .x,
        path = climate.beech.path,
        scale.climate = TRUE
      )
      # Run climwin per site
      weatheRcues::climwin_site_days(
        climate_data = climate_data,
        data = Fagus.seed %>% filter(plotname.lon.lat == .x),
        site.name = .x,
        range = range,
        cinterval = 'day',
        refday = c(01, 11),
        optionwindows = "absolute",
        climate_var = "TMEAN",
        give.clean = T
      )
    },
    .options = furrr_options(seed = TRUE) #setting the seed across the parallel workers.
  )

  qs::qsave(
    statistics_absolute_climwin,
    here('outputs/statistics_absolute_climwin.qs')
  )
} else {
  statistics_absolute_climwin = qs::qread(
    'outputs/statistics_absolute_climwin.qs'
  )
}


mean(statistics_absolute_climwin$window.open)
mean(statistics_absolute_climwin$window.close)
#159-427

climwin.dd = statistics_absolute_climwin %>%
  dplyr::select(sitenewname, window.open, window.close) %>%
  pivot_longer(window.open:window.close, values_to = 'days.reversed') %>%
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed))


climwin.all.sites = climwin.dd %>%
  dplyr::select(sitenewname, name, days.reversed) %>%
  pivot_wider(names_from = "name", values_from = "days.reversed") %>%
  left_join(methods.collection.mv2) %>%
  #arrange(Collection_method) %>%
  mutate(sitenewname = fct_reorder(sitenewname, Collection_method)) %>%
  ggplot() +
  geom_segment(
    aes(
      y = window.open,
      yend = window.close,
      x = sitenewname,
      col = Collection_method
    ),
    size = 2
  ) +
  coord_flip() +
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  ylim(0, 600) +
  scale_color_futurama() +
  scale_fill_futurama() +
  ggpubr::theme_cleveland() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  ylab('Days reversed')

cowplot::save_plot(
  here("figures/FigureS1.png"),
  climwin.all.sites,
  ncol = 1.5,
  nrow = 1.8,
  dpi = 300
)


quibble2(statistics_absolute_climwin$window.open, q = c(0.25, 0.5, 0.75))
quibble2(statistics_absolute_climwin$window.close, q = c(0.25, 0.5, 0.75))


#it is still useful for me because here I will also work with simple correlation later
results.moving.site = FULL.moving.climate.analysis(
  seed.data.all = Fagus.seed,
  climate.path = climate.beech.path,
  refday = 305,
  lastdays = max(range),
  rollwin = 1,
  myform = formula('log.seed ~ TMEAN'),
  model_type = 'lm'
)


#create 50 list - 1 per site
Results_daily = results.moving.site %>%
  arrange(plotname.lon.lat) %>%
  group_by(plotname.lon.lat) %>%
  group_split()

name_daily = results.moving.site %>%
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
  ~ CSP_function_site(
    Results_daily[[.]],
    unique(Results_daily[[.]]$plotname.lon.lat),
    seed.data = Fagus.seed %>% arrange(plotname.lon.lat),
    climate.path = climate.beech.path,
    refday = 305,
    lastdays = max(range),
    rollwin = 1
  )
)
qs::qsave(statistics_csp_method, here('outputs/statistics_csp.qs'))


output_fit_summary.best.csp = statistics_csp_method %>%
  group_by(sitenewname) %>%
  dplyr::slice(which.max(r2))


csp.all.sites = output_fit_summary.best.csp %>%
  left_join(methods.collection.mv2) %>%
  ungroup() %>%
  arrange(Collection_method) %>%
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot() +
  geom_segment(
    aes(
      y = window.open,
      yend = window.close,
      x = sitenewname,
      col = Collection_method
    ),
    size = 2
  ) +
  coord_flip() +
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  ylim(0, 600) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_futurama() +
  scale_fill_futurama() +
  ggpubr::theme_cleveland() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  ylab('Days reversed')
csp.all.sites

cowplot::save_plot(
  here("figures/csp_allsites_methods.png"),
  csp.all.sites,
  ncol = 1.5,
  nrow = 1.8,
  dpi = 300
)


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

site = unique(Fagus.seed$plotname.lon.lat) #[30]
#refday = 305
statistics_psr_method = map_dfr(
  site,
  ~ PSR_function_site(
    site = .x,
    climate.path = climate.beech.path,
    seed.data = Fagus.seed, # Use .x correctly here
    lastdays = max(range),
    matrice = c(3, 1),
    refday = 305,
    knots = NULL,
    tolerancedays = 7
  )
)

qs::qsave(statistics_psr_method, here('outputs/statistics_psr.qs'))


output_fit_summary.psr.best = statistics_psr_method %>%
  group_by(sitenewname) %>%
  slice(which.max(r2))
psr.all.sites = output_fit_summary.psr.best %>%
  left_join(methods.collection.mv2) %>%
  ungroup() %>%
  arrange(Collection_method) %>%
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot() +
  geom_segment(
    aes(
      y = window.open,
      yend = window.close,
      x = sitenewname,
      col = Collection_method
    ),
    size = 2
  ) +
  geom_point(
    aes(y = window.open, x = sitenewname, col = Collection_method),
    size = 2,
    shape = 15
  ) +
  coord_flip() +
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  ylim(0, 600) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_futurama() +
  scale_fill_futurama() +
  ggpubr::theme_cleveland() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  ylab('Days reversed')
psr.all.sites

cowplot::save_plot(
  here("figures/psr_allsites_methods.png"),
  psr.all.sites,
  ncol = 1.5,
  nrow = 1.8,
  dpi = 300
)


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
  ~ basiccues_function_site(
    Results_daily[[.]],
    unique(Results_daily[[.]]$plotname.lon.lat),
    seed.data = Fagus.seed,
    climate.path = climate.beech.path,
    refday = 305,
    lag = 100,
    threshold = 3,
    lastdays = max(range),
    rollwin = 1
  )
)


qs::qsave(statistics_basic_method, here('outputs/statistics_basic.qs'))

output_fit_summary.basic.best = statistics_basic_method %>%
  group_by(sitenewname) %>%
  slice(which.max(r2)) %>%
  ungroup()

basic.all.sites = output_fit_summary.basic.best %>%
  left_join(methods.collection.mv2) %>%
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot() +
  geom_segment(
    aes(
      y = window.open,
      yend = window.close,
      x = sitenewname,
      col = Collection_method
    ),
    size = 2
  ) +
  coord_flip() +
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  ylim(0, 600) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_futurama() +
  scale_fill_futurama() +
  ggpubr::theme_cleveland() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  ylab('Days reversed')
basic.all.sites
cowplot::save_plot(
  here("figures/basic_allsites_methods.png"),
  basic.all.sites,
  ncol = 1.5,
  nrow = 1.5,
  dpi = 300
)

quibble2(output_fit_summary.basic.best$window.open, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.basic.best$window.close, q = c(0.25, 0.5, 0.75))

#put a long lag
windows.avg = output_fit_summary.basic.best %>%
  dplyr::select(sitenewname, window.open, window.close) %>%
  mutate(method = 'Signal processing') %>%
  bind_rows(
    output_fit_summary.psr.best %>%
      dplyr::select(sitenewname, window.open, window.close) %>%
      mutate(method = 'P-spline regression')
  ) %>%
  bind_rows(
    statistics_absolute_climwin %>%
      dplyr::select(sitenewname, window.open, window.close) %>%
      mutate(method = 'Sliding time window')
  ) %>%
  bind_rows(
    output_fit_summary.best.csp %>%
      dplyr::select(sitenewname, window.open, window.close) %>%
      mutate(method = 'Climate sensitivity profile')
  ) %>%
  group_by(method) %>%
  summarise(
    wind.open_median = median(window.open, na.rm = TRUE),
    wind.close_median = median(window.close, na.rm = TRUE),
    wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
    wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
    wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
    wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE),
    wind.open_mean = mean(window.open),
    wind.close_mean = mean(window.close),
    wind.open_se = se(window.open),
    wind.close_se = se(window.close)
  )

library(RColorBrewer)
median.windows.plot = windows.avg %>%
  dplyr::select(method, wind.open_median:wind.open_q75) %>%
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
  separate_wider_delim(
    typewind,
    delim = '_',
    names = c('windows.type', 'metric')
  ) %>%
  pivot_wider(names_from = 'metric', values_from = value) %>%
  ggplot(aes(x = method, group = windows.type, col = windows.type)) +
  geom_point(
    aes(y = median),
    size = 2,
    position = position_dodge(width = 0.2)
  ) +
  geom_pointrange(
    mapping = aes(y = median, ymin = q25, ymax = q75),
    position = position_dodge(width = 0.2)
  ) +
  coord_flip() +
  xlab('') +
  ylab('Days reversed') +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_colour_manual(values = brewer.pal(6, "Paired")[c(5, 6)]) +
  scale_colour_manual(values = brewer.pal(6, "Paired")[c(5, 6)]) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(.8, .8), legend.title = element_blank()) +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')

#GENERIC PLOT R2 and AIC
bestr2 = output_fit_summary.basic.best %>%
  dplyr::select(sitenewname, r2, AIC) %>%
  mutate(method = 'Signal processing') %>%
  bind_rows(
    output_fit_summary.psr.best %>%
      dplyr::select(sitenewname, r2, AIC) %>%
      mutate(method = 'P-spline regression')
  ) %>%
  bind_rows(
    statistics_absolute_climwin %>%
      dplyr::select(sitenewname, R2, AIC) %>%
      rename(r2 = R2) %>%
      mutate(method = 'Sliding time window')
  ) %>%
  bind_rows(
    output_fit_summary.best.csp %>%
      dplyr::select(sitenewname, r2, AIC) %>%
      mutate(method = 'Climate sensitivity profile')
  )

mean_r2 <- bestr2 %>%
  group_by(method) %>%
  summarise(mean_r2 = median(r2, na.rm = TRUE))

averagedensr2method = bestr2 %>%
  ggplot(aes(x = r2, fill = method, col = method)) +
  geom_density(alpha = 0.2, size = .8) +
  labs(x = "R-squared (r2)", y = "Density") +
  #geom_vline(data = mean_r2, aes(xintercept = mean_r2, color = method),
  #           linetype = "dashed", size = .5) +
  ggpubr::theme_pubclean() +
  geomtextpath::geom_textvline(
    data = mean_r2,
    aes(label = round(mean_r2, 3), color = method, xintercept = mean_r2),
    vjust = -0.3,
    hjust = 1,
    fontface = 'bold',
    linetype = 'dashed',
    show.legend = FALSE
  ) +
  scale_color_viridis_d(name = 'Method') +
  scale_fill_viridis_d(name = 'Method') +
  guides(col = guide_legend(ncol = 2)) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

averagedensr2method
cowplot::save_plot(
  here("figures/averager2.method.png"),
  averagedensr2method +
    median.windows.plot +
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(size = 12)),
  ncol = 1.8,
  nrow = 1.4,
  dpi = 300
)

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
num_iterations <- 10 # Specify the number of iterations
#I saved them both (part 1 and part2) just to check if it is working properly
output_fit_summary_cv_fin_part1 = NULL
output_fit_summary_cv_fin_part2 = NULL
# Loop through each site
save.bc = TRUE
for (i in 1:length(beech.site.all)) {
  site.cv <- Fagus.seed %>%
    filter(plotname.lon.lat == beech.site.all[i])

  # Make 5 blocks of data
  blocks <- create_blocks(site.cv, num_blocks)

  # Perform cross-validation multiple times
  for (j in 1:num_iterations) {
    block_indices <- seq_along(blocks) # Should be the same as block size
    test_indices <- sample(block_indices, size = test_block_nb, replace = FALSE) # Randomly select 2 validation blocks
    train_indices <- setdiff(block_indices, test_indices) # Remaining blocks for training

    # Combine the block of training
    train_blocks <- do.call(rbind, blocks[train_indices]) # Combine training blocks
    validation_blocks <- do.call(rbind, blocks[test_indices]) # Combine validation blocks

    climate_data <- format_climate_data(
      site = beech.site.all[i],
      path = climate.beech.path,
      scale.climate = T
    )

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

    psr1 <- runing_psr_site(
      bio_data = train_blocks,
      site = beech.site.all[i],
      climate_csv = climate_data,
      lastdays = max(range),
      refday = refday,
      rollwin = 1,
      covariates.of.interest = 'TMEAN',
      matrice = c(3, 1),
      knots = NULL,
      tolerancedays = 7
    ) %>%
      dplyr::slice(which.max(r2)) %>%
      dplyr::mutate(method = 'psr')

    moving.site1 <- site.moving.climate.analysis(
      bio_data = train_blocks,
      climate_data = climate_data,
      lastdays = lastdays,
      refday = 305,
      myform = formula('log.seed ~ TMEAN'),
      model_type = 'lm'
    )

    csp1 <- runing_csp_site(
      Results_CSPsub = moving.site1,
      data = train_blocks,
      siteneame.forsub = beech.site.all[i],
      option.check.name = TRUE,
      climate_csv = climate_data,
      refday = 305,
      lastdays = max(range),
      rollwin = 1,
      variablemoving = 'TMEAN',
      myform.fin = formula('log.seed ~ TMEAN'),
      model_type = 'lm',
      optim.k = F
    ) %>%
      dplyr::slice(which.max(r2)) %>%
      dplyr::mutate(method = 'csp')

    basic1 <- runing_basic_cues(
      lag = 100,
      threshold = 3,
      influence = 0,
      tolerancedays = 7,
      refday = 305,
      lastdays = max(range),
      rollwin = 1,
      siteforsub = beech.site.all[i],
      climate_csv = climate_data,
      Results_CSPsub = moving.site1,
      myform.fin = formula('log.seed ~ TMEAN'),
      model_type = 'lm',
      data = train_blocks
    ) %>%
      dplyr::slice(which.max(r2)) %>%
      dplyr::mutate(method = 'signal')

    temp.win <- bind_rows(climwin1, psr1, csp1, basic1) %>%
      dplyr::select(
        window.open,
        window.close,
        intercept.estimate,
        slope.estimate,
        method
      )

    yearneed <- 2
    yearperiod <- (min(climate_data$year) + yearneed):max(climate_data$year)

    # Apply the function across all years in yearperiod and combine results
    rolling.data <- map_dfr(
      yearperiod,
      reformat.climate.backtothepast,
      climate = climate_data,
      yearneed = yearneed,
      refday = refday,
      lastdays = lastdays,
      rollwin = 1,
      variablemoving = 'TMEAN'
    )

    output_fit_summary.temp <- purrr::map_dfr(
      1:nrow(temp.win),
      ~ cross_validation_outputs_windows_modelling(
        .,
        tible.sitelevel = validation_blocks,
        window_ranges_df = temp.win,
        rolling.data = rolling.data,
        myform.fin = formula('log.seed ~ TMEAN')
      )
    )

    output_fit_summary_cv_fin_part1 <- bind_rows(
      output_fit_summary_cv_fin_part1,
      output_fit_summary.temp
    )
  }
  output_fit_summary_cv_fin_part2 = bind_rows(
    output_fit_summary_cv_fin_part2,
    output_fit_summary_cv_fin_part1
  )
}
if (save.bc == T) {
  qs::qsave(
    output_fit_summary_cv_fin_part2,
    here(paste0(
      'outputs/outputs_blockcrosstotalv2_',
      num_blocks,
      '_train_',
      test_block_nb,
      '.qs'
    ))
  )
}
#qs::qsave(output_fit_summary_cv_fin_part2, here('outputs/output_fit_summary_cv_fin_part2.qs' ))

output_fit_summary_cv_fin_part2 %>%
  as_tibble() %>%
  drop_na(mae) %>%
  left_join(methods.collection.mv2) %>%
  dplyr::mutate(
    mae = as.numeric(mae),
    sitenewname = fct_reorder(sitenewname, Collection_method)
  ) %>%
  ggplot() +
  geom_point(
    aes(x = mae, y = sitenewname, col = Collection_method),
    alpha = .8,
    shape = 21
  ) +
  facet_grid(. ~ method.cues) +
  geom_vline(xintercept = 0) +
  scale_color_futurama() +
  scale_fill_futurama() +
  ggpubr::theme_cleveland() +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  ylab('Days reversed')

block.cross.figure = output_fit_summary_cv_fin_part2 %>%
  left_join(methods.collection.mv2) %>%
  mutate(
    method.cues = recode(
      method.cues,
      "climwin" = "Sliding time window",
      "csp" = "Climate sensitivity profile",
      "psr" = "P-spline regression",
      'signal' = 'Signal processing'
    )
  ) %>%
  drop_na(mae) %>%
  ggplot() +
  gghalves::geom_half_boxplot(
    aes(x = Collection_method, y = mae, fill = Collection_method),
    center = T,
    errorbar.draw = FALSE,
    width = 0.8,
    nudge = 0.02,
    alpha = .9
  ) +
  gghalves::geom_half_violin(
    aes(x = Collection_method, y = mae, fill = Collection_method),
    side = "r",
    nudge = 0.02,
    alpha = .65,
    col = "black"
  ) +
  facet_grid(. ~ method.cues) +
  geom_vline(xintercept = 0) +
  scale_color_futurama() +
  scale_fill_futurama() +
  ggpubr::theme_pubclean() +
  geom_hline(yintercept = 0) +
  theme(
    legend.position = 'none',
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, size = 12)
  ) +
  ylab('Mean Absolute Error') +
  xlab('')
block.cross.figure
cowplot::save_plot(
  here("figures/Figure.block.png"),
  block.cross.figure,
  ncol = 1.6,
  nrow = 1.4,
  dpi = 300
)

####################################################################################
##############################################################################
#######################
#data sample size effect on window identification
#now test the frac dataset
#take more than 12 hours if 10 iteration (climwin takes some times)
# Initialize a list to store all results
result_list_all <- list()

# Define the sample fractions you want to iterate over
#fractions <- c(0.3, 0.5, 0.7, 0.9)
years_to_extract <- c(5, 10, 15, 20)
# Define the number of iterations for each fraction
n_iter <- 10 #take > 12 and < 24 hours on mac os

# Outer loop to iterate over the fractions
for (num_years in years_to_extract) {
  # Initialize a list to store results for each fraction
  result_list <- list()

  # Inner loop to perform the operation multiple times for the given fraction
  for (i in 1:n_iter) {
    # Get the sample list for the given fraction
    sample_result <- run_sampling_fraction_all_methods(
      years_to_extract = num_years,
      range = range,
      Fagus.seed = Fagus.seed,
      beech.site.all = beech.site.all,
      climate.path = climate.beech.path,
    )

    # Append this sample's results to the overall list for this fraction
    result_list[[i]] <- sample_result
  }

  # Store the results for this fraction in the overall list, using the fraction as the key
  result_list_all[[paste0("years_", num_years)]] <- result_list
}


#updated here 2025 feb
#take few hours for one iteration
#do this only to the longest time series, more than 50 years of observations
subset.long.term = Fagus.seed %>% filter(n > 50)
beech.site.subset.longterm = Fagus.seed %>%
  filter(n > 50) %>%
  dplyr::select(plotname.lon.lat) %>%
  distinct() %>%
  as.matrix()

sample_result = run_sampling_fraction_all_methods(
  years_to_extract = c(5),
  range,
  Fagus.seed = subset.long.term,
  beech.site.all = beech.site.subset.longterm,
  climate.path = climate.beech.path
)

sample_result = results_all_years
#qs::qsave(sample_result, "sample_result.qs")

merged_df <- imap_dfr(sample_result, function(sublist, name) {
  bind_rows(lapply(sublist, function(df) {
    df %>% mutate(source = name) # Add column with list name
  }))
})


summary.sample.size = merged_df %>%
  group_by(method, source, sitenewname) %>%
  mutate(
    r2_highest = ifelse(method %in% c('csp', 'psr', 'signal'), r2 == max(r2), T)
  ) %>%
  dplyr::filter(r2_highest) %>%
  mutate(
    method = recode(
      method,
      "climwin" = "Sliding time window",
      "csp" = "Climate sensitivity profile",
      "psr" = "P-spline regression",
      'signal' = 'Peak signal identification'
    ),
    source = recode(
      source,
      "years_1" = 5,
      "years_2" = 10,
      "years_3" = 15,
      'years_4' = 20
    )
  ) %>%
  group_by(method, source) %>%
  dplyr::select(method, source, window.open, window.close) %>%
  summarise(
    wind.open_median = median(window.open, na.rm = TRUE),
    wind.close_median = median(window.close, na.rm = TRUE),
    wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
    wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
    wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
    wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE),
    wind.open_mean = mean(window.open, na.rm = TRUE),
    wind.close_mean = mean(window.close, na.rm = T),
    wind.open_se = se(window.open),
    wind.close_se = se(window.close)
  )

sensi.data.plot = summary.sample.size %>%
  dplyr::select(method, source, wind.open_median:wind.open_q75) %>%
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
  separate_wider_delim(
    typewind,
    delim = '_',
    names = c('windows.type', 'metric')
  ) %>%
  pivot_wider(names_from = 'metric', values_from = value) %>%
  mutate(
    source = as.numeric(source),
    method_sample = paste(method, source, sep = " (n=") %>% paste0(")"),
    sample_fraction_factor = as_factor(paste0(source, ' years'))
  ) %>%
  ggplot(aes(
    x = sample_fraction_factor,
    group = windows.type,
    col = windows.type
  )) +
  geom_point(
    aes(y = median),
    size = 2,
    position = position_dodge(width = 0.5)
  ) +
  geom_pointrange(
    mapping = aes(y = median, ymin = q25, ymax = q75),
    position = position_dodge(width = 0.5)
  ) +
  coord_flip() +
  facet_wrap(. ~ method, scales = 'free_y') +
  xlab('') +
  ylab('Days reversed') +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
  ggpubr::theme_pubr() +
  theme(legend.position = 'bottom', legend.title = element_blank())
sensi.data.plot
cowplot::save_plot(
  here("figures/sensitivity.method.png"),
  sensi.data.plot,
  ncol = 1.4,
  nrow = 1.4,
  dpi = 300
)

#for text
summary.sample.size %>% filter(source == 20)
summary.sample.size %>%
  filter(method == "Climate sensitivity profile") %>%
  View()
summary.sample.size %>% filter(method == "Sliding time window") %>% View()
summary.sample.size %>%
  filter(method == "Peak signal identification") %>%
  View()

#alternative figure for michal

datatype20.windows = merged_df %>%
  group_by(method, source, sitenewname) %>%
  mutate(
    r2_highest = ifelse(method %in% c('csp', 'psr', 'signal'), r2 == max(r2), T)
  ) %>%
  dplyr::filter(r2_highest) %>%
  mutate(
    method = recode(
      method,
      "climwin" = "Sliding time window",
      "csp" = "Climate sensitivity profile",
      "psr" = "P-spline regression",
      'signal' = 'Peak signal identification'
    ),
    source = recode(
      source,
      "years_1" = 5,
      "years_2" = 10,
      "years_3" = 15,
      'years_4' = 20
    )
  ) %>%
  filter(source == 20) %>% #do it for the 20 years
  left_join(methods.collection.mv2) %>%
  dplyr::select(method, Collection_method, window.open, window.close) %>%
  group_by(method, Collection_method) %>%
  dplyr::select(method, source, window.open, window.close) %>%
  summarise(
    wind.open_median = median(window.open, na.rm = TRUE),
    wind.close_median = median(window.close, na.rm = TRUE),
    wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
    wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
    wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
    wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE),
    wind.open_mean = mean(window.open, na.rm = TRUE),
    wind.close_mean = mean(window.close, na.rm = T),
    wind.open_se = se(window.open),
    wind.close_se = se(window.close)
  )

View(datatype20.windows)


datatype20.windows.figure = datatype20.windows %>%
  dplyr::select(method, Collection_method, wind.open_median:wind.open_q75) %>%
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
  separate_wider_delim(
    typewind,
    delim = '_',
    names = c('windows.type', 'metric')
  ) %>%
  pivot_wider(names_from = 'metric', values_from = value) %>%
  ggplot(aes(x = Collection_method, group = windows.type, col = windows.type)) +
  geom_point(
    aes(y = median),
    size = 2,
    position = position_dodge(width = 0.5)
  ) +
  geom_pointrange(
    mapping = aes(y = median, ymin = q25, ymax = q75),
    position = position_dodge(width = 0.5)
  ) +
  coord_flip() +
  facet_grid(. ~ method, scales = 'free_y') +
  xlab('') +
  ylab('Days reversed') +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
  ggpubr::theme_pubr() +
  theme(legend.position = 'bottom', legend.title = element_blank())
datatype20.windows.figure
#
#
#
#   pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
#   ggplot(aes(x = Collection_method      , col = typewind)) +
#   geom_boxplot(aes(y = value)) +
#   coord_flip() +
#   facet_grid(.~method)+
#   xlab('') +
#   ylab('Days reversed') +
#   theme(legend.position = 'bottom',
#         legend.title = element_blank()) +
#   geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5,6)])+
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5,6)])+
#   ggpubr::theme_pubr()+
#   theme(legend.position = 'bottom',
#         legend.title = element_blank())
cowplot::save_plot(
  here("figures/Figure.20yrsdata.subset.windows.png"),
  datatype20.windows.figure,
  ncol = 2,
  nrow = 1.4,
  dpi = 300
)


####################################################################################
##############################################################################
#######################

#now test the frac dataset
#take more than 12 hours if 10 iteration (climwin takes some times)
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
    sample_result <- run_sampling_fraction_all_methods(
      fraction = fraction,
      range = range,
      Fagus.seed = Fagus.seed,
      beech.site.all = beech.site.all,
      climate.path = climate.beech.path
    )

    # Append this sample's results to the overall list for this fraction
    result_list[[i]] <- sample_result
  }

  # Store the results for this fraction in the overall list, using the fraction as the key
  result_list_all[[paste0("fraction_", fraction * 100)]] <- result_list
}


sampling_data_effect <- map_dfr(names(result_list_all), function(sublist) {
  # extract my sublist
  new_list <- result_list_all[[sublist]]
  new_list[[1]] <- lapply(new_list[[1]], rename_columns_if_needed) # Apply renaming to the first list

  # rbind dplyr
  datatible.sim <- bind_rows(new_list)
  return(datatible.sim)
})

#qs::qsave(sampling_data_effect, here('outputs/sampling_data_effect.qs'))
final_combined_df = qs::qread(here('final_combined_df.qs'))

summary.test.sampling = final_combined_df %>%
  mutate(
    method = recode(
      method,
      "climwin" = "Sliding time window",
      "csp" = "Climate sensitivity profile",
      "psr" = "P-spline regression",
      'signal' = 'Signal processing'
    )
  ) %>%
  group_by(method, sample_fraction) %>%
  summarise(
    wind.open_median = median(window.open, na.rm = TRUE),
    wind.close_median = median(window.close, na.rm = TRUE),
    wind.close_q25 = quantile(window.close, 0.3, na.rm = TRUE),
    wind.close_q75 = quantile(window.close, 0.7, na.rm = TRUE),
    wind.open_q25 = quantile(window.open, 0.3, na.rm = TRUE),
    wind.open_q75 = quantile(window.open, 0.7, na.rm = TRUE),
    wind.open_mean = mean(window.open),
    wind.close_mean = mean(window.close),
    wind.open_se = se(window.open),
    wind.close_se = se(window.close)
  )


sensi.data.plot = summary.test.sampling %>%
  dplyr::select(method, sample_fraction, wind.open_median:wind.open_q75) %>%
  pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
  separate_wider_delim(
    typewind,
    delim = '_',
    names = c('windows.type', 'metric')
  ) %>%
  pivot_wider(names_from = 'metric', values_from = value) %>%
  mutate(
    sample_fraction = as.numeric(sample_fraction),
    method_sample = paste(method, sample_fraction, sep = " (n=") %>%
      paste0(")"),
    sample_fraction_factor = as_factor(paste0(sample_fraction * 100, '%'))
  ) %>%
  ggplot(aes(
    x = sample_fraction_factor,
    group = windows.type,
    col = windows.type
  )) +
  geom_point(
    aes(y = median),
    size = 2,
    position = position_dodge(width = 0.5)
  ) +
  geom_pointrange(
    mapping = aes(y = median, ymin = q25, ymax = q75),
    position = position_dodge(width = 0.5)
  ) +
  coord_flip() +
  facet_wrap(. ~ method, scales = 'free_y') +
  xlab('Method (Sample Size)') +
  ylab('Days reversed') +
  theme(legend.position = 'bottom', legend.title = element_blank()) +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
  ggpubr::theme_pubr() +
  theme(legend.position = 'bottom', legend.title = element_blank())

sensi.data.plot

cowplot::save_plot(
  here("figures/sensitivity.method.png"),
  sensi.data.plot,
  ncol = 1.4,
  nrow = 1.4,
  dpi = 300
)


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
startyear = 1980
years <- startyear:2020
yday <- 365
climate_data_simulated <- data.frame(
  date = seq.Date(
    from = as.Date(paste0(startyear, "-01-01")),
    by = "day",
    length.out = length(years) * yday
  ),
  temp = runif(length(years) * yday, min = 10, max = 30) # random temperatures
)

climate_data_simulated <- climate_data_simulated %>%
  mutate(
    year = format(date, "%Y"),
    yday = yday(date),
    year = as.numeric(year),
    LONGITUDE = NA,
    LATITUDE = NA
  ) %>%
  mutate(across(c(temp), scale)) %>%
  mutate(across(c(temp), as.vector)) %>%
  rename(TMEAN = temp)

#qs::qsave(climate_data_simulated, here::here('outputs/climate_data_simulated.qs'))

# Select ~ one week in June
june_week <- climate_data_simulated %>%
  filter(yday >= 150 & yday <= 160)

raw.data.param.alpha = statistics_absolute_climwin$intercept.estimate
raw.data.param.beta = statistics_absolute_climwin$slope.estimate
raw.data.param.sigma = statistics_absolute_climwin$sigma


#generate parameter range values for the simulation
setup.param = parameter.range(
  raw.data.param.alpha,
  raw.data.param.beta,
  raw.data.param.sigma,
  option = 'min.max'
)

# Container to store results
results_list_fakedata = simulated_fake_bio_data(
  fakeclimatewindow = june_week,
  setup.param = setup.param,
  num_simulations = 1000,
  save_fake_data = TRUE,
  overwrite = FALSE
)


#simulation40years1000sim = simulation_runing_window(climate_data_simulated,
#                                    results_list_fakedata)

simulation_runing_window = function(
  climate_data_simulated,
  results_list_fakedata
) {
  print(paste0("number of simulation:", length(results_list_fakedata)))
  if (
    length(results_list_fakedata) > 10 & length(results_list_fakedata) < 100
  ) {
    print(paste0("It might take few hours"))
  }
  if (length(results_list_fakedata) > 100) {
    print(paste0(
      "It might take some days (because I've made a basic for loop)"
    ))
  }

  fin.sim = NULL
  for (k in 1:length(results_list_fakedata)) {
    #climwin
    # statistics_climwin_method.simulated = climwin_site_days(
    #   climate_data = climate_data_simulated,
    #   data = results_list_fakedata[[k]],
    #   site.name = as.character(k),
    #   range = c(365, 0),
    #   cinterval = 'day',
    #   refday = c(01, 11),
    #   optionwindows = 'absolute',
    #   climate_var = 'TMEAN'
    # ) %>%
    #   dplyr::mutate(method = 'climwin',
    #                 r2 = R2)

    #create object for csp and basic methods
    run.sim.day.res = site.moving.climate.analysis(
      bio_data = as.data.frame(results_list_fakedata[[k]]),
      climate.data = climate_data_simulated,
      lastdays = 365,
      formula('log.seed ~ rolling_avg_tmean'),
      refday = 305,
      yearneed = 1
    )

    statistics_csp_method.simulated = runing_csp_site(
      Results_CSPsub = run.sim.day.res,
      data = as.data.frame(results_list_fakedata[[k]]),
      siteneame.forsub = as.character(k),
      option.check.name = TRUE,
      climate_csv = climate_data_simulated,
      refday = 305,
      lastdays = 365,
      rollwin = 1,
      optim.k = F,
      variablemoving = 'TMEAN',
      yearneed = 1
    ) %>%
      slice(which.max(r2)) %>%
      dplyr::mutate(method = 'csp')

    statistics_basic_method.simulated <- runing_basic_cues(
      lag = 100,
      threshold = 3,
      influence = 0,
      tolerancedays = 7,
      refday = 305,
      lastdays = 365,
      rollwin = 1,
      siteforsub = as.character(k),
      climate_csv = climate_data_simulated,
      Results_CSPsub = run.sim.day.res,
      data = results_list_fakedata[[k]],
      yearneed = 1
    ) %>%
      dplyr::slice(which.max(r2)) %>%
      dplyr::mutate(method = 'signal')

    statistics_psr_method.simulated <- runing_psr_site(
      bio_data = results_list_fakedata[[k]],
      site = as.character(k),
      climate_csv = climate_data_simulated,
      tot_days = 365,
      refday = 305,
      rollwin = 1,
      covariates.of.interest = 'TMEAN',
      matrice = c(3, 1),
      knots = NULL,
      tolerancedays = 7,
      plot = TRUE,
      yearneed = 1
    ) %>%
      dplyr::slice(which.max(r2)) %>%
      dplyr::mutate(method = 'psr')

    #extract just column of interest and previous parameters
    fin.sim.temp = bind_rows(
      #statistics_climwin_method.simulated,
      statistics_csp_method.simulated,
      statistics_basic_method.simulated,
      statistics_psr_method.simulated
    ) %>%
      dplyr::select(
        sitenewname,
        reference.day,
        method,
        window.open,
        window.close,
        slope.estimate,
        intercept.estimate,
        AIC,
        r2
      ) %>%
      dplyr::mutate(
        r2.before.simulate = base::unique(
          results_list_fakedata[[k]]$r2.before.climwin
        ),
        alpha.random.value.setup = base::unique(
          results_list_fakedata[[k]]$alpha.random.value.setup
        ),
        beta.random.value.setup = base::unique(
          results_list_fakedata[[k]]$beta.random.value.setup
        ),
        sigma.random.value.setup = base::unique(
          results_list_fakedata[[k]]$sigma.random.value.setup
        ),
        nb.year.biosimulate = nrow(results_list_fakedata[[k]])
      ) %>%
      dplyr::rename(simulation = sitenewname)

    fin.sim = rbind(fin.sim, fin.sim.temp)
  }
  return(fin.sim)
}


# Load required libraries
simulation40years1000sim = simulation_runing_window(
  climate_data_simulated,
  results_list_fakedata = results_list_fakedata[c(1:2)]
)

# Define the function with parallel processing
simulation_runing_window_parallel <- function(
  climate_data_simulated = climate_data_simulated,
  results_list_fakedata = results_list_fakedata
) {
  library(doParallel)
  library(foreach)
  print(paste0("Number of simulations: ", length(results_list_fakedata)))
  if (
    length(results_list_fakedata) > 10 & length(results_list_fakedata) < 100
  ) {
    print("It might take a few hours.")
  }
  if (length(results_list_fakedata) > 100) {
    print("It might take some days (considering this is parallelized).")
  }

  # Set up parallel backend
  num_cores <- parallel::detectCores() - 3 #just to avoid too much memory :(
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # Use foreach for parallel processing
  #needed to add all function and everything ...
  fin.sim <- foreach::foreach(
    k = 1:length(results_list_fakedata),
    .combine = rbind,
    .packages = c("weatheRcues", 'tidyverse'), #export all the new functions :') and the obect
    .export = c(
      "site.moving.climate.analysis",
      "runing_csp_site",
      "runing_basic_cues",
      "runing_psr_site",
      'reformat.climate.backtothepast',
      'runing.movingwin.analysis',
      'correlation.spearman.se',
      "optimize_and_fit_gam",
      'get_predictions_windows',
      'extract_consecutive_sequences',
      'save_window_ranges',
      "find_concurrent_period",
      "ThresholdingAlgo",
      "replace_0",
      ls(globalenv()) #just use environment .. because some parameters included in env
    )
  ) %dopar%
    {
      library(tidyverse)
      library(weatheRcues)
      # climwin method
      statistics_climwin_method.simulated <- climwin_site_days(
        climate_data = climate_data_simulated,
        data = results_list_fakedata[[k]],
        site.name = as.character(k),
        range = c(365, 0),
        cinterval = 'day',
        refday = c(1, 11),
        optionwindows = 'absolute',
        climate_var = 'TMEAN'
      ) %>%
        dplyr::mutate(method = 'climwin', r2 = R2)

      # csp and basic methods
      run.sim.day.res <- site.moving.climate.analysis(
        bio_data = results_list_fakedata[[k]],
        climate.data = climate_data_simulated,
        lastdays = 365,
        formula('log.seed ~ rolling_avg_tmean'),
        refday = 305,
        yearneed = 1
      )

      statistics_csp_method.simulated <- runing_csp_site(
        Results_CSPsub = run.sim.day.res,
        data = results_list_fakedata[[k]],
        siteneame.forsub = as.character(k),
        option.check.name = TRUE,
        climate_csv = climate_data_simulated,
        refday = 305,
        lastdays = 365,
        rollwin = 1,
        optim.k = F,
        variablemoving = 'TMEAN',
        yearneed = 1
      ) %>%
        dplyr::slice(which.max(r2)) %>%
        dplyr::mutate(method = 'csp')

      statistics_basic_method.simulated <- runing_basic_cues(
        lag = 100,
        threshold = 3,
        influence = 0,
        tolerancedays = 7,
        refday = 305,
        lastdays = 365,
        rollwin = 1,
        siteforsub = as.character(k),
        climate_csv = climate_data_simulated,
        Results_CSPsub = run.sim.day.res,
        data = results_list_fakedata[[k]],
        yearneed = 1
      ) %>%
        dplyr::slice(which.max(r2)) %>%
        dplyr::mutate(method = 'signal')

      #psr method
      statistics_psr_method.simulated <- runing_psr_site(
        bio_data = results_list_fakedata[[k]],
        site = as.character(k),
        climate_csv = climate_data_simulated,
        tot_days = 365,
        refday = 305,
        rollwin = 1,
        covariates.of.interest = 'TMEAN',
        matrice = c(3, 1),
        knots = NULL,
        tolerancedays = 7,
        plot = TRUE,
        yearneed = 1
      ) %>%
        dplyr::slice(which.max(r2)) %>%
        dplyr::mutate(method = 'psr')

      # Combine results into a temporary dataframe
      fin.sim.temp <- dplyr::bind_rows(
        statistics_climwin_method.simulated,
        statistics_csp_method.simulated,
        statistics_basic_method.simulated,
        statistics_psr_method.simulated
      ) %>%
        dplyr::select(
          sitenewname,
          reference.day,
          method,
          window.open,
          window.close,
          slope.estimate,
          intercept.estimate,
          AIC,
          r2
        ) %>%
        dplyr::mutate(
          r2.before.simulate = unique(
            results_list_fakedata[[k]]$r2.before.climwin
          ),
          alpha.random.value.setup = unique(
            results_list_fakedata[[k]]$alpha.random.value.setup
          ),
          beta.random.value.setup = unique(
            results_list_fakedata[[k]]$beta.random.value.setup
          ),
          sigma.random.value.setup = unique(
            results_list_fakedata[[k]]$sigma.random.value.setup
          ),
          nb.year.biosimulate = nrow(results_list_fakedata[[k]]),
          simulation = as.character(k)
        )

      fin.sim.temp
    }

  # Stop parallel processing
  parallel::stopCluster(cl)
  return(fin.sim)
}

#you might have warning on your computer, if you are runing > 100 sim
#started at 1.30pm
#finished 3.45pm ! YOUHOU
#you might have several warning in your session
result_simulations = simulation_runing_window_parallel(
  climate_data_simulated,
  results_list_fakedata
)

#qs::qsave(result_simulations, here('outputs/result_simulations.qs'))

#calendar for sim
calendar.sim = goingbackpastdayscalendar(
  refday = 304,
  lastdays = 400,
  yearback = 2
) %>%
  filter(YEAR == 1950) %>%
  dplyr::select(MONTHab, DOY, days.reversed) %>%
  mutate(
    datefake = as.Date(DOY, origin = "1948-01-01"),
    day.month = format(datefake, "%m-%d"),
    year = case_when(
      days.reversed < 366 ~ 1,
      days.reversed >= 366 & days.reversed < 730 ~ 2,
      days.reversed >= 730 & days.reversed < 1095 ~ 3,
      TRUE ~ 4
    )
  )

#do the round if not issues by merging :")
data.plot.sim = result_simulations %>%
  group_by(method) %>%
  summarise(
    wind.open_median = round(median(window.open, na.rm = TRUE)),
    wind.close_median = round(median(window.close, na.rm = TRUE)),
    wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
    wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
    wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
    wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = starts_with('wind'), names_to = 'typewind') %>%
  separate_wider_delim(
    typewind,
    delim = '_',
    names = c('windows.type', 'metric')
  ) %>%
  pivot_wider(names_from = 'metric', values_from = value) %>%
  ungroup() %>%
  dplyr::left_join(
    calendar.sim %>%
      mutate(median = days.reversed) %>%
      dplyr::select(DOY, median) %>%
      distinct() %>%
      as_tibble(),
    join_by(median)
  ) %>%
  mutate(
    method = recode(
      method,
      "climwin" = "Sliding time window",
      "csp" = "Climate sensitivity profile",
      "psr" = "P-spline regression",
      'signal' = 'Signal processing'
    )
  )


#here with the mean and sd
result_simulations %>%
  group_by(method) %>%
  summarise(
    wind.open_mean = round(mean(window.open, na.rm = TRUE)),
    wind.close_mean = round(mean(window.close, na.rm = TRUE)),
    wind.close_sd = sd(window.close, na.rm = TRUE),
    wind.open_sd = sd(window.open, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = starts_with('wind'), names_to = 'typewind') %>%
  separate_wider_delim(
    typewind,
    delim = '_',
    names = c('windows.type', 'metric')
  ) %>%
  pivot_wider(names_from = 'metric', values_from = value) %>%
  dplyr::left_join(
    calendar.sim %>%
      mutate(mean = days.reversed) %>%
      dplyr::select(DOY, mean) %>%
      distinct() %>%
      as_tibble(),
    join_by(mean)
  )

#now the figures
sim.wind = ggplot(
  data.plot.sim,
  aes(x = method, group = windows.type, col = windows.type)
) +
  geom_point(
    aes(y = median),
    size = 2,
    position = position_dodge(width = 0.2)
  ) +
  geom_pointrange(
    mapping = aes(y = median, ymin = q25, ymax = q75),
    position = position_dodge(width = 0.2)
  ) +
  coord_flip() +
  xlab('') +
  ylab('Days reversed') +
  scale_color_brewer(palette = "Paired") +
  scale_fill_brewer(palette = "Paired") +
  scale_colour_manual(values = brewer.pal(12, "Paired")[c(9, 10)]) +
  scale_colour_manual(values = brewer.pal(12, "Paired")[c(9, 10)]) +
  ggpubr::theme_pubr() +
  theme(legend.position = c(.8, .8), legend.title = element_blank()) +
  geom_hline(yintercept = 156, color = 'black', linetype = 'dotted') +
  geom_hline(yintercept = 146, color = 'black', linetype = 'dotted')

r2sim = ggplot(result_simulations, aes(x = r2.before.simulate)) +
  geom_histogram(
    aes(y = ..density..),
    binwidth = 0.05,
    fill = "#0073C2FF",
    color = "white",
    alpha = 0.8
  ) +
  geom_density(color = "#FC4E07", size = 1) +
  geom_vline(
    aes(xintercept = mean(r2.before.simulate, na.rm = TRUE)),
    color = "blue",
    linetype = "dashed",
    size = 0.8
  ) +
  labs(x = "R2 from simulated relationship", y = "Density") +
  theme_minimal(base_size = 15) +
  ggpubr::theme_cleveland()


cowplot::save_plot(
  here("figures/FigureSim.png"),
  r2sim +
    sim.wind +
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(size = 12)) +
      plot_layout(widths = c(.3, 2)),
  ncol = 1.3,
  nrow = 1,
  dpi = 300
)

#find when correct value
perfectwin = simulation.all %>%
  filter(between(WindowOpen, 150, 160), between(WindowClose, 140, 150))
summary(perfectwin$r2.before.climwin)
