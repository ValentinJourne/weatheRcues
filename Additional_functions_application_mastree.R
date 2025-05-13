formatting_mastree_fagus = function(mastreedata) {
  mastreedata %>%
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
      scalingbeta = y_transformation_betareg(scaling01)
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
}

obtained_cleaned_method_collection = function(data) {
  data %>%
    group_by(sitenewname, Country, Collection_method, Length) %>%
    summarise(
      average.log.seed = mean(log.seed, na.rm = T),
      sd.log.seed = sd(log.seed, na.rm = T)
    ) %>%
    dplyr::select(
      sitenewname,
      Country,
      Collection_method,
      Length,
      average.log.seed,
      sd.log.seed
    ) %>%
    mutate(
      Collection_method = factor(
        dplyr::recode(
          as_factor(Collection_method),
          "seed count" = "Seed count",
          "seed trap" = "Seed trap",
          "visual crop assessment" = "Visual crop",
        ),
        levels = c("Seed count", "Seed trap", "Visual crop")
      )
    ) %>%
    ungroup() %>%
    distinct()
}
