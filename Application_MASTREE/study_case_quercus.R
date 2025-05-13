devtools::load_all(here::here())
library(weatheRcues)
library(tidyverse)
library(sf)
library(ggsci)

################################
#example with all quercus robur petraea and spp (hybrid)
#subset data of interest here
Quercus.seed = initial.data.mastree %>%
  filter(
    !Variable == "flower" &
      VarType == "C" &
      !Unit == "index" &
      !Variable == "pollen"
  ) %>%
  group_by(Species, VarType) %>%
  mutate(nMax = n()) %>%
  filter(str_detect(Species, 'Quercus robur|Quercus petraea|Quercus spp.')) %>%
  filter(
    !str_detect(Country, 'United States of America|Russian Federation (the)')
  ) %>%
  filter(Year > 1952 & Year < 2023) %>%
  mutate(
    sitenewname = paste0(Alpha_Number, "_", Site_number, "_", Species_code)
  ) %>%
  group_by(sitenewname) %>%
  mutate(log.seed = log(1 + Value), scale.seed = scale(Value), n = n()) %>%
  filter(n > 19) %>%
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
  ungroup() %>%
  as.data.frame()

#df for collection method
methods.collection.quercus = Quercus.seed %>%
  dplyr::select(
    sitenewname,
    Country,
    Collection_method,
    Length,
    Longitude,
    Latitude
  ) %>%
  distinct()

#for after
quercus.site.all = unique(Quercus.seed$plotname.lon.lat)
climate.quercus.path <- here::here('climate_dailyEOBS/')


########################################################################################################################
########################################################################################################################
########################################################################################################################
#TEST with climwin
########################################################################################################################
run.climwin <- T
#take some time ! alternative it will load the outputs
if (run.climwin == T) {
  #for parallel to take less time usually ...
  library(furrr)
  plan(multisession)
  #run loop for each sites
  statistics_absolute_climwin_quercus <- future_map_dfr(
    quercus.site.all,
    ~ {
      #for package in progress
      devtools::load_all(here::here())
      library(weatheRcues)
      #format climate each site
      climate_data <- format_climate_data(
        site = .x,
        path = climate.quercus.path,
        scale.climate = T
      )
      #run climwin each sites
      climwin_site_days(
        climate_data = climate_data,
        data = Quercus.seed %>% filter(plotname.lon.lat == .x),
        site.name = .x,
        range = c(600, 0),
        cinterval = 'day',
        refday = c(01, 11),
        optionwindows = 'absolute',
        climate_var = 'TMEAN'
      )
    }
  )

  qs::qsave(
    statistics_absolute_climwin_quercus,
    here::here('outputs/statistics_absolute_climwin_QUERCUS.qs')
  )
} else {
  statistics_absolute_climwin_quercus = qs::qread(
    'outputs/statistics_absolute_climwin_QUERCUS.qs'
  )
}


########################################################################################################################
########################################################################################################################
########################################################################################################################
#TEST with CSP
########################################################################################################################

results.moving.site.quercus = FULL.moving.climate.analysis(
  seed.data.all = Quercus.seed,
  climate.path = climate.quercus.path,
  refday = 305,
  lastdays = 600,
  rollwin = 1
)
Results_daily.quercus = results.moving.site.quercus %>%
  arrange(sitenewname) %>%
  group_by(sitenewname) %>%
  group_split()

name_daily.quercus = results.moving.site.quercus %>%
  arrange(sitenewname) %>%
  group_by(sitenewname) %>%
  group_keys()

# Set names of the list based on group keys
names(Results_daily.quercus) <- apply(name_daily.quercus, 1, paste)

statistics_csp_method.quercus = map_dfr(
  1:length(Results_daily.quercus),
  ~ CSP_function_site(
    Results_daily.quercus[[.]],
    unique(Results_daily.quercus[[.]]$plotname.lon.lat),
    seed.data = Quercus.seed %>% arrange(sitenewname),
    climate.path = climate.beech.path,
    refday = 305,
    lastdays = max(range),
    rollwin = 1
  )
)


########################################################################################################################
########################################################################################################################
########################################################################################################################
#Final figures
########################################################################################################################

#make map
plot_locations.quercus <- st_as_sf(
  Quercus.seed %>%
    dplyr::select(Longitude, Latitude, plotname.lon.lat, Collection_method) %>%
    distinct(),
  coords = c("Longitude", "Latitude")
) %>%
  st_set_crs(4326)

Maps.quercus <- ggplot(data = mapBase$geometry) +
  geom_sf(fill = "grey90", colour = "black", alpha = .8, linewidth = .25) + #plot map of France
  xlab(" ") +
  ylab(" ") + #white or aliceblue
  geom_point(
    data = plot_locations.quercus %>%
      mutate(x = unlist(map(geometry, 1)), y = unlist(map(geometry, 2))) %>%
      as.data.frame(),
    aes(x = x, y = y), #, col = Collection_method, fill = Collection_method
    stroke = 1,
    alpha = .8,
    shape = 21,
    size = 3,
    col = "#3D3B25FF",
    fill = "#3D3B25FF"
  ) +
  scale_size_continuous(
    breaks = c(0.1, 0.3, 0.8),
    range = c(0, 7)
  ) +
  coord_sf(xlim = c(-10, 25), ylim = c(40, 60)) +
  scale_y_continuous(breaks = c(40, 50, 60, 70)) +
  scale_x_continuous(breaks = c(-10, 0, 10, 20)) +
  ylab("Lat") +
  xlab("Lon") +
  ggpubr::theme_pubr() +
  scale_color_futurama() +
  scale_fill_futurama() +
  guides(
    col = "none",
    fill = guide_legend(title = NULL, override.aes = list(size = 8))
  ) +
  theme(legend.position = 'none')
Maps.quercus

#make moving window plot comparison methods csp vs climwin
quercus2method = statistics_csp_method.quercus %>%
  group_by(sitenewname) %>%
  dplyr::slice(which.max(r2)) %>%
  mutate(method = 'climate sensitivity profile') %>%
  ungroup() %>%
  bind_rows(
    statistics_absolute_climwin_quercus %>%
      mutate(method = 'sliding moving window')
  ) %>%
  left_join(methods.collection.quercus) %>%
  mutate(sitenewname = forcats::fct_reorder(sitenewname, Collection_method)) %>%
  ggplot() +
  geom_rect(
    aes(ymin = 62, ymax = 153, xmin = -Inf, xmax = Inf),
    fill = 'yellow',
    col = 'white',
    alpha = 0.01
  ) +
  geom_rect(
    aes(ymin = 428, ymax = 519, xmin = -Inf, xmax = Inf),
    fill = 'yellow',
    col = 'white',
    alpha = 0.01
  ) +
  geom_rect(
    aes(ymin = 246, ymax = 155, xmin = -Inf, xmax = Inf),
    fill = 'darkgreen',
    col = 'white',
    alpha = 0.01
  ) +
  geom_rect(
    aes(ymin = Inf, ymax = 520, xmin = -Inf, xmax = Inf),
    fill = 'darkgreen',
    col = 'white',
    alpha = 0.01
  ) +
  geom_rect(
    aes(ymin = -Inf, ymax = 62, xmin = -Inf, xmax = Inf),
    fill = 'darkgoldenrod',
    col = 'white',
    alpha = 0.01
  ) +
  geom_rect(
    aes(ymin = 427, ymax = 337, xmin = -Inf, xmax = Inf),
    fill = 'darkgoldenrod',
    col = 'white',
    alpha = 0.01
  ) +
  geom_segment(
    aes(
      y = window.open,
      yend = window.close,
      x = sitenewname,
      col = method,
      fill = method
    ),
    size = 2
  ) +
  geom_point(
    aes(
      y = window.open,
      x = sitenewname,
      col = method,
      fill = method
    ),
    size = 1.5,
    shape = 22
  ) +
  coord_flip() +
  geom_hline(yintercept = 133, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
  ylim(0, 600) +
  scale_colour_manual(values = viridis_pal(option = "magma")(10)[c(2, 6)]) +
  scale_fill_manual(values = viridis_pal(option = "magma")(10)[c(2, 6)]) +
  ggpubr::theme_cleveland() +
  theme(legend.position = 'none', legend.title = element_blank()) +
  ylab('Days reversed') +
  facet_grid(method ~ .) +
  theme(
    axis.text.y = element_blank(),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12)
  )
#extract R2 from models
R2quercus = statistics_csp_method.quercus %>%
  group_by(sitenewname) %>%
  dplyr::slice(which.max(r2)) %>%
  dplyr::select(sitenewname, r2) %>%
  mutate(method = 'climate sensivitiy profil') %>%
  bind_rows(
    statistics_absolute_climwin_quercus %>%
      dplyr::select(sitenewname, r2) %>%
      mutate(method = 'sliding moving window')
  ) %>%
  ggplot(aes(x = r2, fill = method)) +
  geom_density(alpha = 0.42) +
  labs(x = expression(R^2), y = "Density") +
  ggpubr::theme_pubclean() +
  theme(legend.position = c(1, 1)) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = viridis_pal(option = "magma")(10)[c(2, 6)]) +
  scale_fill_manual(values = viridis_pal(option = "magma")(10)[c(2, 6)])

#combine all the map
Figure.quercus = cowplot::plot_grid(
  (Maps.quercus / R2quercus),
  quercus2method,
  labels = 'auto'
)
cowplot::save_plot(
  here::here("figures/Figure.quercus.png"),
  Figure.quercus,
  ncol = 1.8,
  nrow = 1.5,
  dpi = 300,
  bg = "white"
)
