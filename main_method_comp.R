library(tidyverse)
library(climwin)
initial.data.mastree <- read.csv('/Users/vjourne/Library/CloudStorage/Dropbox/Mastree/MASTREEplus_2022-02-03_V1.csv',stringsAsFactors = F)


y.transf.betareg <- function(y){
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}
#now use the filtering 
#keep cont, remove index, flower and pollen 
Fagus.seed = initial.data.mastree %>%
  filter(!Variable == "flower" & VarType == "C" & !Unit=="index" & !Variable == "pollen") %>% 
  group_by(Species, VarType) %>% 
  mutate(nMax = n()) %>% 
  filter(Species == "Fagus sylvatica") %>% 
  #filter(Unit == "seeds/m2" | Unit == "seeds/individual") %>% 
  filter(Year > 1951 & Year < 2020) %>% 
  mutate(sitenewname= paste0(Alpha_Number, "_",Site_number, "_", Species_code)) %>% 
  group_by(sitenewname) %>% 
  mutate(ScaledSeedProduction = (Value - min(Value, na.rm = T))/(max(Value)-min(Value, na.rm = T)),
         ScaledSeedProductionBR = y.transf.betareg(ScaledSeedProduction),
         log.seed = log(1+Value),
         scale.seed = scale(Value),
         n = n()) %>% 
  filter(n > 14) %>% 
  mutate() %>% 
  mutate(Date = paste0( "15/06/",Year)) %>% 
  as.data.frame() %>% 
  mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y"))

plot(Fagus.seed$log.seed, Fagus.seed$scale.seed)
hist(Fagus.seed$log.seed)

#some option 
#here for climwin 
optionwindows = 'absolute'
referenceclimwin = c(01, 10)
#going back to one year before + 6 month more or less
range = c((round(365/2)+365), 0)

#climate path 
#temperature.path = here::here('~/Documents/GITprojects/breakdownmast/climatedailyERA5')
climate.beech.path <- "/Users/vjourne/Library/CloudStorage/Dropbox/MastreeForemast/climateSiteDaily"
#beech site unique name 
beech.site.all = unique(Fagus.seed$sitenewname)

cinterval = "day"

# Define a function to process each site
climwin_site_days <- function(site.name) {
  
  # Filter climate data for the site
  climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = site.name)
  
  # Load the biological data for the site
  bio_data <- Fagus.seed %>%
    filter(sitenewname == site.name) %>%
    as.data.frame()
  
  # Reformat the climate data
  climate_csv <- read_csv(climate_beech_unique, col_types = cols(...1 = col_skip())) %>%
    mutate(Tmean = base::scale(Tmean),
           Prp = base::scale(Prp),
           across(c(Tmean, Prp), as.vector),
           date = as.Date(paste0(day,'/',month,'/', year), format = "%d/%m/%Y"),
           yday = lubridate::yday(date))
  
  # Set the range for the sliding window
  range <- c((round(365/2) + 365), 0)
  
  # Run the climwin analysis
  climwin_output <- slidingwin(
    xvar = list(temperature.degree = climate_csv$Tmean),
    cdate = climate_csv$date,
    bdate = bio_data$Date2,
    baseline = lm(log.seed ~ 1, data = bio_data),
    cinterval = cinterval,
    range = range,
    type = optionwindows,
    refday = c(01, 11), # Absolute based on 01-11 (1st November)
    stat = "mean",
    cmissing = 'method2',
    func = "lin"
  )
  
  # Extract the summary for the best model
  broom_summary <- broom::tidy(climwin_output[[1]]$BestModel) %>%
    filter(term == 'climate') %>%
    select(-term)
  
  # Extract performance statistics and combine with site information
  statistics <- bind_cols(
    sitenewname = site.name,
    climwin_output$combos, # Extract statistics for all variants
    performance::model_performance(climwin_output[[1]]$BestModel) %>% as.data.frame(),
    broom_summary
  )
  
  return(statistics)
}


#option to run in parallel to make it faster
run.climwin = F
if(run.climwin==T){
  library(furrr)#make it in parallele, https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
  plan(multisession)#background R sessions (on current machine) , the computer is flying (alt cluster)
  statistics_absolute_climwin <- future_map_dfr(beech.site.all, climwin_site_days)
  qs::qsave(statistics_absolute_climwin, 
            here('statistics_absolute_climwin.weeks.qs'))}else{
  statistics_absolute_climwin = qs::qread('statistics_absolute_climwin.qs')
}
mean(statistics_absolute_climwin$WindowOpen)
mean(statistics_absolute_climwin$WindowClose)
#159-427
plotall(dataset = climwin_output.relative[[1]]$Dataset, 
               bestmodel = climwin_output.relative[[1]]$BestModel, 
               bestmodeldata = climwin_output.relative[[1]]$BestModelData,
               title=climwin_output.relative$combos$climate)


goingbackpastdayscalendar <-function(refday = 274, 
                                     lastdays = 1095, 
                                     yearback = 3){
  monthstart = c('-01-01')
  DATE = seq(as.Date(paste0("1946", monthstart)), as.Date(paste0("1953", monthstart)), by="days")
  MONTH =  format(as.Date(DATE, format="%Y-%m-%d"),"%m") %>% as.numeric()
  MONTHab = month.abb[MONTH]
  YEAR =  format(as.Date(DATE, format="%Y-%m-%d"),"%Y") %>% as.numeric()
  
  DOY = yday(DATE)
  dfata = data.frame(DATE,YEAR,MONTHab, DOY)
  yearperiod = 1946:1953
  sizevec = length(unique(YEAR))-yearback
  refday = refday
  vectotemp = NULL
  
  for(k in 1:sizevec){
    #year +3, becusse of the 36 month analysis 
    yearsref = yearperiod[k]
    yearrefminusOne <- yearsref-yearback #PREVIOUS = 3
    
    tt <- dfata %>% 
      filter(YEAR <= yearsref & YEAR >= yearrefminusOne) %>% 
      mutate(referenceFin = ifelse(YEAR == yearsref & DOY == refday, 1,
                                   ifelse(YEAR == yearsref & DOY > refday, NA, 0))) %>% 
      filter(!is.na(referenceFin)) %>% 
      as.data.frame()
    #create sequence going back 365 month before 
    seqDays <- seq(1,nrow(tt),1)
    newsequance <- rep(seqDays)
    
    ttup <- tt %>% 
      mutate(days.reversed = rev(newsequance))%>% 
      filter(days.reversed< lastdays )
    ttupfin = ttup %>%
      arrange(days.reversed)  %>% 
      mutate(YEAR = max(YEAR))  
    
    vectotemp <- rbind(vectotemp, ttupfin) 
  }
  
  return(vectotemp)
}


calendar = goingbackpastdayscalendar(refday = 305,
                          lastdays = 547, 
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

climwin.dd = statistics_absolute_climwin %>% 
  dplyr::select(sitenewname, WindowOpen, WindowClose) %>% 
  pivot_longer(WindowOpen:WindowClose, values_to = 'days.reversed') %>% 
  left_join(calendar %>% dplyr::select(MONTHab, DOY, days.reversed)) 

ggplot(climwin.dd %>% dplyr::select(sitenewname, name, days.reversed) %>% 
         pivot_wider(names_from = "name", values_from = "days.reversed"))+
  geom_segment(aes(y = WindowOpen, yend = WindowClose, x = sitenewname), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,547)

quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}
  
quibble2(statistics_absolute_climwin$WindowOpen, q = c(0.25, 0.5, 0.75))
quibble2(statistics_absolute_climwin$WindowClose, q = c(0.25, 0.5, 0.75))


# Function to process each year to get past moving climate 
reformat.climate.backtothepast <- function(yearsref = 2000, 
                                           climate = climate, 
                                           yearneed = 2, 
                                           refday = 274, 
                                           lastdays = 1095, 
                                           rollwin = 1, 
                                           variablemoving = 'temperature.degree') {
  # Print parameter values - to see if no shit 
  print(paste("rollwin - size window:", rollwin))
  print(paste("refday - reference day:", refday))
  print(paste("lastdays - last day:", lastdays))
  print(paste("variablemoving - variable to move:", variablemoving))
  print(paste("yearneed - adjust to maturation fruit time :", yearneed, 'years'))
  
  if (!variablemoving %in% names(climate)) {
    warning(paste("Warning: Column", variablemoving, "not found in the dataset."))
    return(NULL)
  }
  
  yearrefminusOne <- yearsref - yearneed
  tt <- climate %>%
    filter(year <= yearsref & year >= yearrefminusOne) %>%
    mutate(referenceFin = ifelse(year == yearsref & yday == refday, 1, ifelse(year == yearsref & yday > refday, NA, 0))) %>%
    filter(!is.na(referenceFin)) %>%
    as.data.frame()
  
  # Create sequence going back lastdays days before the reference day
  seqDays <- seq(1, nrow(tt), 1)
  newsequance <- rep(seqDays)
  
  ttup <- tt %>%
    mutate(days.reversed = rev(newsequance)) %>%
    filter(days.reversed < lastdays)
  
  #use !!sym; convert a string, here my variable name, to a symbol
  ttupfin <- ttup %>%
    arrange(days.reversed) %>%
    mutate(rolling_avg_tmean = zoo::rollmeanr(!!sym(variablemoving), k = rollwin, fill = NA, align = 'right')) %>%
    mutate(year = max(year)) %>%
    dplyr::select(sitenewname, year, date, yday, days.reversed, rolling_avg_tmean)
  
  return(ttupfin)
}

#se correlation
correlation.spearman.se <- function (cor, n) 
{
  se.cor <- sqrt((1 - cor^2)^2 * (1 + cor^2/2)/(n - 3))
  return(se.cor)
}

library(broom)
runing.movingwin.analysis = function(Fagus.seed = Fagus.seed,
                                     rolling.temperature.data = rolling.temperature.data,
                                     method = 'spearman',
                                     covariates.of.interest = 'rolling_avg_tmean',
                                     myform = formula('log.seed~rolling_avg_tmean')){
  
  #merge data seed to moving climate
  tible.sitelevel = Fagus.seed %>% #site = bio_data 
    rename(year = Year) %>% 
    left_join(rolling.temperature.data) %>% 
    drop_na(!!sym(covariates.of.interest))
  
  #define correlation - calculate correlation and se, extract also p value 
  n = tible.sitelevel %>% dplyr::select(year, sitenewname) %>% distinct() %>% nrow()
  
  correlation.all <- tible.sitelevel %>% 
    nest(data = -days.reversed) %>%
    mutate(correlation = map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$estimate)) %>% 
    mutate(pvalue.cor = map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$p.value))
  
  
  cortemp = correlation.all %>% 
    unnest(c(correlation, pvalue.cor)) %>% 
    dplyr::select(days.reversed, correlation, pvalue.cor) %>% 
    dplyr::mutate(correlation.se = correlation.spearman.se(.$correlation, n)) %>% 
    dplyr::mutate(sitenewname = unique(tible.sitelevel$sitenewname))
  
  #use purr for iteration 
  #here broom package used, because it is simple lm (and not glmmTMB, need to adjust then )
  fitted_models = tible.sitelevel %>%
    nest(data = -days.reversed) %>%
    mutate(model = purrr::map(data, ~lm(myform, data = ., na.action = na.omit)),
           tidied = purrr::map(model, tidy),
           glanced = purrr::map(model, glance),
           augmented = purrr::map(model, augment))
  
  slope = fitted_models %>%
    unnest(tidied) %>% 
    filter(str_detect(term, as.character(myform)[3])) %>% 
    dplyr::select(days.reversed,term, estimate, std.error, p.value) %>% 
    mutate(sitenewname  = unique(tible.sitelevel$sitenewname)) %>% 
    left_join(fitted_models %>%
                unnest(glanced) %>% 
                dplyr::select(days.reversed,r.squared, logLik) %>% 
                mutate(sitenewname  = unique(tible.sitelevel$sitenewname))) %>% 
    dplyr::select(sitenewname, days.reversed, everything()) %>% 
    left_join(cortemp)
  
  return(slope)
}


#now lets try the other methods I updated - based on thakery method nature 2016 
site.moving.climate.analysis <- function(site.name, 
                                         Fagus.seed, 
                                         climate.beech.path, 
                                         lastdays, 
                                         myform) {
  # Load the biological data for the site
  bio_data <- Fagus.seed %>%
    filter(sitenewname == site.name) %>%
    as.data.frame() 
  
  # Climate load and format
  climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = site.name)
  # Reformat the climate data
  climate_csv <- read_csv(climate_beech_unique, col_types = cols(...1 = col_skip())) %>%
    mutate(Tmean = base::scale(Tmean),
           Prp = base::scale(Prp),
           across(c(Tmean, Prp), as.vector),
           date = as.Date(paste0(day,'/',month,'/', year), format = "%d/%m/%Y"),
           yday = lubridate::yday(date))
  
  # Define the year period
  yearneed <- 2
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  
  # Apply the function across all years in yearperiod and combine results
  rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                      climate = climate_csv, 
                                      yearneed = yearneed, 
                                      refday = 305, 
                                      lastdays = lastdays, 
                                      rollwin = 1, 
                                      variablemoving = 'Tmean')
  
  results.moving.site = runing.movingwin.analysis(Fagus.seed = Fagus.seed,
                                                  rolling.temperature.data = rolling.temperature.data,
                                                  method = 'spearman',
                                                  covariates.of.interest = 'rolling_avg_tmean',
                                                  myform = myform)
  
  return(results.moving.site)
}


# Main function to process all sites
FULL.moving.climate.analysis <- function(Fagus.seed = Fagus.seed,
                                         climate.beech.path = climate.beech.path,
                                         refday = 305,
                                         lastdays = max(range),
                                         rollwin = 1) {
  al.sites <- unique(Fagus.seed$sitenewname)
  
  results.moving <- map_dfr(al.sites, 
                            ~site.moving.climate.analysis(site.name = .x, 
                                                          Fagus.seed = Fagus.seed, 
                                                          climate.beech.path = climate.beech.path,                            
                                                          lastdays = lastdays,
                                                          myform = formula('log.seed~rolling_avg_tmean')))
  
  return(results.moving)
}

#it is still useful for me because here I will also work with simple correlation
results.moving.site = FULL.moving.climate.analysis(Fagus.seed = Fagus.seed,
                             climate.beech.path = climate.beech.path,
                             refday = 305,
                             lastdays = max(range),
                             rollwin = 1)


ggplot(results.moving.site)+
  geom_bar(aes(y = r.squared, x = days.reversed), size = 2, stat="identity" ) +
  geom_vline(xintercept = 133, color = 'red')+
  geom_vline(xintercept = 498, color = 'red')

#ok now use the other function for identificatuon cues 
#function from Nature Tackery paper 
#function to find concurrent window eg sequences of days
find_concurrent_period=function(temp_window,pred_C){
  
  #create dummy index
  idx=1
  #find differences between submitted dates and note whic are greater than 1
  diff_id = c(1,which(diff(temp_window)>1))
  
  if(length(diff_id)>1){if((which(diff(temp_window)>1))[1]==1){diff_id=diff_id[-1]}}
  if(length(temp_window)>1){
    #find all the  series of sequentioal dates and give unique id
    for(rv in 1:length(diff_id)){
      
      if(rv==length(diff_id)){
        idx=c(idx,rep(rv,length((diff_id[rv]+1):length(temp_window))))
      }else{
        idx=c(idx,rep(rv,length((diff_id[rv]+1):(diff_id[rv+1]))))
      }
      
    }
  }
  
  #from the estimated coefficients and the concurrent window indices (idx) find which period has the most extreme average. call that the period to use
  mx_coef = which.max(tapply(abs(pred_C$fit[temp_window]),idx,mean))
  
  return(temp_window[idx==mx_coef])
  
}

get_predictions_windows <- function(slope_gam, rs_gam, temporary) {
  pred_S <- predict(slope_gam, se.fit = TRUE, type = "response")
  pred_R <- predict(rs_gam, se.fit = TRUE, type = "response")
  
  coef_range_l <- which(((pred_S$fit - 1.96 * pred_S$se.fit) - 100) <= quantile(pred_S$fit - 100, 0.025))
  coef_range_u <- which(((pred_S$fit + 1.96 * pred_S$se.fit) - 100) >= quantile(pred_S$fit - 100, 0.975))
  r_sq_range <- which(((pred_R$fit + 1.96 * pred_R$se.fit) - 100) >= quantile(pred_R$fit - 100, 0.975))
  
  temp_window_l <- coef_range_l[is.element(coef_range_l, r_sq_range)]
  temp_window_u <- coef_range_u[is.element(coef_range_u, r_sq_range)]
  
  if (length(temp_window_l) == 0) temp_window_l <- coef_range_l
  if (length(temp_window_u) == 0) temp_window_u <- coef_range_u
  
  days <- find_concurrent_period(temp_window = c(temp_window_l, temp_window_u), pred_C = pred_S)
  list(pred_S = pred_S, pred_R = pred_R, days = days)
}

optimize_and_fit_gam <- function(temporary, optim.k = TRUE, plots = F, k = 20) {
  if (optim.k) {
    # Function to check significance in k.check
    is_significant <- function(check) {
      p_values <- check[,'p-value']
      all(p_values < 0.05) # Returns TRUE if all p-values are significant
    }
    
    # Range of k values to try
    k_values <- seq(10, 365)  # Adjust to one per day 
    optimal_k <- NULL
    
    for (k in 1:length(k_values)) {
      # Fit the model
      slope_gam <- gam(slope ~ s(day, k = k_values[k], bs = "cr"), data = temporary)
      
      # Perform k.check
      check <- k.check(slope_gam)
      pvalue <- check[,'p-value']
      kfin <- check[,"k'"]
      
      # Check if all p-values are significant
      if (!is_significant(check)) {
        optimal_k <- k_values[k]
        break
      }
    }
    if (!is.null(optimal_k)) {
      cat("Optimal k found:", optimal_k, "\n")
      slope_gam <- gam(slope ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      rs_gam <- gam(r_s ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      k <- optimal_k
    } else {
      cat("No optimal k found within the specified range. Using default k = -1.\n")
      slope_gam <- gam(slope ~ s(day, k = -1, bs = "cr"), data = temporary)
      rs_gam <- gam(r_s ~ s(day, k = -1, bs = "cr"), data = temporary)
      k <- -1
    }
  } else {
    slope_gam <- gam(slope ~ s(day, k = k, bs = "cr"), data = temporary)
    rs_gam <- gam(r_s ~ s(day, k = k, bs = "cr"), data = temporary)
    k <- k
  }
  
  if (plots) {
    results <- cowplot::plot_grid(
      draw(slope_gam, residuals = TRUE) + ylab('Slope (partial effect)'), 
      draw(rs_gam, residuals = TRUE) + ylab('R2 (partial effect)')
    )
    print(results)
  }
  
  list(slope_gam = slope_gam, rs_gam = rs_gam, k = k)
  
}

#extract sequence 
extract_consecutive_sequences <- function(values, keep_all = FALSE) {
  # Initialize variables
  sequences <- list()
  current_sequence <- integer(0)
  
  # Iterate through the values
  #if does not match afte +1 then it will create new sequence length 
  for (i in seq_along(values)) {
    if (i == 1 || values[i] == values[i - 1] + 1) {
      current_sequence <- c(current_sequence, values[i])
    } else {
      sequences <- c(sequences, list(current_sequence))
      current_sequence <- values[i]  # Start new sequence
    }
  }
  
  # Add the last sequence
  sequences <- c(sequences, list(current_sequence))
  
  if (!keep_all) {
    # Find the longest sequence
    longest_sequence <- sequences[[which.max(lengths(sequences))]]
    return(longest_sequence)
  } else {
    return(sequences)
  }
}

# save windows ranges 
save_window_ranges <- function(sequences_days) {
  # Initialize a list to store the results
  window_ranges <- list()
  
  # Check the length of sequences_days
  if (length(sequences_days) == 1) {
    # If only one sequence, save it directly
    windowsclose <- sequences_days[[1]][1]
    windowsopen <- tail(sequences_days[[1]], n = 1)
    
    # Store the result in the list
    window_ranges[[1]] <- c(windowsclose, windowsopen)
  } else {
    # Loop through each sequence and save the windows
    for (p in 1:length(sequences_days)) {
      windowsclose <- sequences_days[[p]][1]
      windowsopen <- tail(sequences_days[[p]], n = 1)
      
      # Store the result in the list
      window_ranges[[p]] <- c(windowsclose, windowsopen)
    }
  }
  
  # Convert the list to a data frame for easier handling
  window_ranges_df <- do.call(rbind, lapply(window_ranges, function(x) {
    data.frame(windowsclose = x[2], windowsopen = x[1])
  }))
  
  return(window_ranges_df)
}

# Define the function to extract sequences with tolerance 
#before I made simple function to extract full sequence, but here it is more tricky, because they have gaps in between 
extract_sequences_auto <- function(vec, tolerance) {
  # Sort the vector to handle sequences in ascending order
  vec <- sort(vec)
  
  # Initialize variables
  sequences <- list()
  current_seq <- c(vec[1])
  
  # Iterate over the vector to find sequences
  for (i in 2:length(vec)) {
    # Check if the current element is within tolerance of the last element in current_seq
    if (vec[i] - current_seq[length(current_seq)] <= tolerance) {
      current_seq <- c(current_seq, vec[i])
    } else {
      # Save the current sequence if it's not empty
      if (length(current_seq) > 1) {
        sequences[[length(sequences) + 1]] <- current_seq
      }
      # Start a new sequence
      current_seq <- c(vec[i])
    }
  }
  
  # Add the last sequence if it's not empty
  if (length(current_seq) > 1) {
    sequences[[length(sequences) + 1]] <- current_seq
  }
  
  return(sequences)
}


#create 61 list - 1 per site 
Results_CSP = results.moving.site %>%
  group_by(sitenewname) %>%
  group_split()

nameCSP =   results.moving.site %>% 
  group_by(sitenewname) %>%
  group_keys()

# Set names of the list based on group keys
names(Results_CSP) <- apply(nameCSP, 1, paste)


output_fit_summary = NULL
rollwin = 1
lastdays = max(range)
myform.fin = formula('log.seed ~ mean.temperature')
refday = 305
library(mgcv)


for(i in 1:length(Results_CSP)){
  #now make subset and run regression
  list_slope <- as.list(Results_CSP[[i]]$estimate)
  list_rs <- as.list(Results_CSP[[i]]$r.squared)
  
  day = seq(rollwin,(lastdays-1),1)
  slope = list_slope
  r_s = list_rs
  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), 
                          day = day, 
                          r_s = unlist(r_s))
  #use the function to optimize k 
  results <- optimize_and_fit_gam(temporary, optim.k = F, plots = F, k = (12+6)) #I will specify just number of month 
  #one year and 6 month for 1.5 years 
  #get days windows 
  #need to add rollwin because it start at 1... 
  days = get_predictions_windows(slope_gam = results[[1]], 
                                 rs_gam = results[[2]], 
                                 temporary)$days
  
  sequences_days = extract_consecutive_sequences(days, keep_all = TRUE)
    
  #subset the data, 
  #scale seed production and do logit transformation (will use this logit for climwin later)
  tible.sitelevel = Fagus.seed %>% #site = bio_data 
    dplyr::filter(sitenewname == unique(Results_CSP[[i]]$sitenewname))
  
  #get the climate id
  # Climate load and format
  climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = unique(Results_CSP[[i]]$sitenewname))
  # Reformat the climate data
  climate_csv <- read_csv(climate_beech_unique, col_types = cols(...1 = col_skip())) %>%
    mutate(Tmean = base::scale(Tmean),
           Prp = base::scale(Prp),
           across(c(Tmean, Prp), as.vector),
           date = as.Date(paste0(day,'/',month,'/', year), format = "%d/%m/%Y"),
           yday = lubridate::yday(date))
  
  # Define the year period
  yearneed <- 2
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  
  
  # Apply the function across all years in yearperiod and combine results
  #here the map dfr will basically do same as aplly by runing the function over all the time period we want 
  rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                      climate = climate_csv, 
                                      yearneed = yearneed, 
                                      refday = 305, 
                                      lastdays = lastdays, 
                                      rollwin = 1, 
                                      variablemoving = 'Tmean')
  
  window_ranges_df <- save_window_ranges(sequences_days) %>% 
    mutate(windows.sequences.number = 1:nrow(.))
  

  #now for me makes an uggly loop
  # Initialize a list to store climate data for each window
  #climate_windows_best_list <- list()
  
  # Loop through each window range in window_ranges_df
  #i need to do this because sometimes the algo identify more sequences of days 
  for (z in 1:nrow(window_ranges_df)) {
    windowsopen <- window_ranges_df$windowsopen[z]
    windowsclose <- window_ranges_df$windowsclose[z]
    window_number <- window_ranges_df$windows.sequences.number[z]
    
    # Filter the rolling temperature data according to the current window range
    climate_windows_best <- rolling.temperature.data %>%
      dplyr::filter(days.reversed >= windowsopen & days.reversed <= windowsclose) %>%
      group_by(sitenewname, year) %>%
      summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(window_number = window_number)  # Add the window number for identification
    
    #need to change year colnames :( 
    fit.best.model = tible.sitelevel %>% 
      rename(year = Year) %>% 
      left_join(climate_windows_best) %>%
      lm(myform.fin, data = .)
    
    intercept = tidy(fit.best.model)$estimate[1]
    intercept.se = tidy(fit.best.model)$std.error[1]
    estimate.model = tidy(fit.best.model)$estimate[2]
    estimate.se.model = tidy(fit.best.model)$std.error[2]
    pvalue.model = tidy(fit.best.model)$p.value[2]
    r2 = glance(fit.best.model)$r.squared
    AIC = glance(fit.best.model)$AIC
    nobs = glance(fit.best.model)$nobs
    nsequence.id = z
    
    output_fit_summary.temp <- data.frame(sitenewname = unique(Results_CSP[[i]]$sitenewname),
                                          reference.day = refday,
                                          windows.size = rollwin,
                                          knots.number = results$k,
                                          window.open = windowsopen,
                                          window.close = windowsclose,
                                          intercept = intercept,
                                          intercept.se = intercept.se,
                                          estimate = estimate.model,
                                          estimate.se = estimate.se.model,
                                          pvalue = pvalue.model,
                                          r2 = r2,
                                          AIC = AIC,
                                          nobs = nobs,
                                          nsequence.id = nsequence.id)
    output_fit_summary = rbind(output_fit_summary, output_fit_summary.temp)
    
  }
}

output_fit_summary.best = output_fit_summary %>% 
  group_by(sitenewname) %>% 
  slice(which.max(r2))

ggplot(output_fit_summary.best)+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,547)

#no subset 
ggplot(output_fit_summary)+
  geom_segment(aes(y = window.open, yend = window.close, x = sitenewname), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,547)

#change name here between closing and opning I mixed them 
quibble2(output_fit_summary.best$window.open, q = c(0.25, 0.5, 0.75))
quibble2(output_fit_summary.best$window.close, q = c(0.25, 0.5, 0.75))


#now simple correlaation and extract window size based on the correlation 
sub.correlation.10days <- results.moving.site %>%
  group_by(sitenewname) %>%
  arrange(days.reversed) %>%
  mutate(
    max_corr_date = days.reversed[which.max(abs(correlation))]
  ) %>%
  filter(days.reversed >= max_corr_date - 10 & days.reversed <= max_corr_date + 10) %>%
  ungroup() %>%
  group_by(sitenewname) %>%
  filter(days.reversed == min(days.reversed) | days.reversed == max(days.reversed)) %>%
  ungroup()
  
#plot with simple correlation, extract top highest correlation , and then 10+- days before highest correlation
ggplot(sub.correlation.10days %>% dplyr::select(sitenewname, days.reversed) %>% 
         group_by(sitenewname) %>% 
         mutate(windows = 1:2) %>% 
         mutate(windows = ifelse(windows == 1, 'close', 'open')) %>% 
         pivot_wider(names_from = "windows", values_from = "days.reversed"))+
  geom_segment(aes(y = open, yend = close, x = sitenewname), size = 2) +
  coord_flip()+
  geom_hline(yintercept = 133, color = 'red')+
  geom_hline(yintercept = 498, color = 'red')+
  ylim(0,547)

stat.correlation = sub.correlation.10days %>% dplyr::select(sitenewname, days.reversed) %>% 
  group_by(sitenewname) %>% 
  mutate(windows = 1:2) %>% 
  mutate(windows = ifelse(windows == 1, 'close', 'open')) %>% 
  pivot_wider(names_from = "windows", values_from = "days.reversed")
quibble2(stat.correlation$open, q = c(0.25, 0.5, 0.75))
quibble2(stat.correlation$close, q = c(0.25, 0.5, 0.75))

#now test last method P spline 
site = '2114_1_FAGSYL'
tot_days = max(range)
covariates.of.interest = 'Tmean'
knots = NULL
output_fit_summary.psr = NULL


bio_data = Fagus.seed %>% 
  dplyr::filter(sitenewname ==site ) 
selectrows <- bio_data
bio_data$seeds <- round(bio_data$log.seed, 3)

climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = site)
# Reformat the climate data
climate_csv <- read_csv(climate_beech_unique, col_types = cols(...1 = col_skip())) %>%
  mutate(Tmean = base::scale(Tmean),
         Prp = base::scale(Prp),
         across(c(Tmean, Prp), as.vector),
         date = as.Date(paste0(day,'/',month,'/', year), format = "%d/%m/%Y"),
         yday = lubridate::yday(date))

# need climate data to be arranged with year as row
# need to reduce climate dataframe to only year, yday and temp
climate2 <- data.frame(year = climate_csv$year, yday = climate_csv$yday, temp = climate_csv$Tmean)
tempmat <- climate2 %>% spread(, key = yday, value = temp)
tempmat <- tempmat[,-1]
#number years monitoring seeds 
ny<-length(bio_data$Year)
nt<-tot_days-1

## Formatting data
#based on Simmonds et al code
index.matrix=matrix(1:nt,ny,nt,byrow=TRUE)


# Define the year period
yearneed <- 2
#will fiter the year period here to the year needeed 
yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)

# Apply the function across all years in yearperiod and combine results
rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                    climate = climate_csv, 
                                    yearneed = yearneed, 
                                    refday = 305, 
                                    lastdays = lastdays, 
                                    rollwin = 1, 
                                    variablemoving = covariates.of.interest)

#merge data seed to moving climate
tible.sitelevel = bio_data %>% #site = bio_data 
  rename(year = Year) %>% 
  left_join(rolling.temperature.data) %>% 
  drop_na(!!sym('rolling_avg_tmean'))

climate2 <- data.frame(year = tible.sitelevel$year, 
                       yday = tible.sitelevel$days.reversed, 
                       temp = tible.sitelevel$rolling_avg_tmean)

covariate.matrix = climate2 %>%
  spread(key = yday, value = temp) %>%
  dplyr::select(-year) %>%                    
  as.matrix()    
covariate.matrix <- unname(as.matrix(covariate.matrix))

#second order, but from https://link.springer.com/article/10.1007/s00484-007-0141-4
#it is possible to adjust between 1-3 (instead of 2) like c(1,1) instead of c(2,1)
#here I specify to 0 instead of 1 as in SImmonds et al, meaning that I have no penaties on slope
#if using 1 instead of 0, the model is not able to identify a cues 
#and I think it is normal and ok to be more flexible because we might have multiple peaks with time 
#So if I want to specify K = 30, for example it is working 
#not if knots are too low 
#which make sense when looking https://link.springer.com/article/10.1007/s00484-011-0472-z 
if(is.null(knots)){
  ifelse(ny < 21, K <- 8, K <- 10)
  if(ny < 11){K <- 6}
  model<-gam(seeds~s(index.matrix,k=K,m=c(2,0),bs="ps",by=covariate.matrix), 
             data = bio_data, method="GCV.Cp")
  summary(model)
}else{
K=knots
model<-gam(seeds~s(index.matrix,k=K,m=c(2,1),bs="ps",by=covariate.matrix), 
           data = bio_data, method="GCV.Cp")
summary(model)}

plotted <- plot(model)
coefs <- data.frame(fit = plotted[[1]]$fit)
#"The most important days of temperature influence on phenology in the year may be identified by as those 
#with partial coefficients greater than or less than zero by more than two times their standard error."
marker <- c(which(coefs$fit > round((1.96*sd(coefs$fit))+mean(coefs$fit),2)), 
            which(coefs$fit < round((1.96*sd(coefs$fit))-mean(coefs$fit)),2))
window <- round(plotted[[1]]$x[find_concurrent_period(marker, coefs)])

#here does not work because not consecutive 
#extract_consecutive_sequences(window, keep_all = T) 
sequences_days = extract_sequences_auto(window, tolerance = 20) 
#i guess now will do the same shit as the other methods 
window_ranges_df <- save_window_ranges(sequences_days) %>% 
  mutate(windows.sequences.number = 1:nrow(.))

for (z in 1:nrow(window_ranges_df)) {
  windowsopen <- window_ranges_df$windowsopen[z]
  windowsclose <- window_ranges_df$windowsclose[z]
  window_number <- window_ranges_df$windows.sequences.number[z]
  
  # Filter the rolling temperature data according to the current window range
  climate_windows_best <- rolling.temperature.data %>%
    dplyr::filter(days.reversed >= windowsopen & days.reversed <= windowsclose) %>%
    group_by(sitenewname, year) %>%
    summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(window_number = window_number)  # Add the window number for identification
  
  #need to change year colnames :( 
  fit.best.model = bio_data %>% 
    rename(year = Year) %>% 
    left_join(climate_windows_best) %>%
    lm(myform.fin, data = .)
  
  intercept = tidy(fit.best.model)$estimate[1]
  intercept.se = tidy(fit.best.model)$std.error[1]
  estimate.model = tidy(fit.best.model)$estimate[2]
  estimate.se.model = tidy(fit.best.model)$std.error[2]
  pvalue.model = tidy(fit.best.model)$p.value[2]
  r2 = glance(fit.best.model)$r.squared
  AIC = glance(fit.best.model)$AIC
  nobs = glance(fit.best.model)$nobs
  nsequence.id = z
  
  output_fit_summary.temp.psr <- data.frame(sitenewname = unique(bio_data$sitenewname),
                                        reference.day = refday,
                                        windows.size = rollwin,
                                        knots.number = K,
                                        window.open = windowsopen,
                                        window.close = windowsclose,
                                        intercept = intercept,
                                        intercept.se = intercept.se,
                                        estimate = estimate.model,
                                        estimate.se = estimate.se.model,
                                        pvalue = pvalue.model,
                                        r2 = r2,
                                        AIC = AIC,
                                        nobs = nobs,
                                        nsequence.id = nsequence.id)
  output_fit_summary.psr = rbind(output_fit_summary.psr, output_fit_summary.temp.psr)
  
}

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
#copy from https://stackoverflow.com/questions/19170130/combining-random-forests-built-with-different-training-sets-in-r
my_combine <- function (...) 
{
  pad0 <- function(x, len) c(x, rep(0, len - length(x)))
  padm0 <- function(x, len) rbind(x, matrix(0, nrow = len - 
                                              nrow(x), ncol = ncol(x)))
  rflist <- list(...)
  areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
  if (any(!areForest)) 
    stop("Argument must be a list of randomForest objects")
  rf <- rflist[[1]]
  classRF <- rf$type == "classification"
  trees <- sapply(rflist, function(x) x$ntree)
  ntree <- sum(trees)
  rf$ntree <- ntree
  nforest <- length(rflist)
  haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
  vlist <- lapply(rflist, function(x) rownames(importance(x)))
  numvars <- sapply(vlist, length)
  if (!all(numvars[1] == numvars[-1])) 
    stop("Unequal number of predictor variables in the randomForest objects.")
  for (i in seq_along(vlist)) {
    if (!all(vlist[[i]] == vlist[[1]])) 
      stop("Predictor variables are different in the randomForest objects.")
  }
  haveForest <- sapply(rflist, function(x) !is.null(x$forest))
  if (all(haveForest)) {
    nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
    rf$forest$nrnodes <- nrnodes
    rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
    rf$forest$nodestatus <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$nodestatus, nrnodes)))
    rf$forest$bestvar <- do.call("cbind", lapply(rflist, 
                                                 function(x) padm0(x$forest$bestvar, nrnodes)))
    rf$forest$xbestsplit <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$xbestsplit, nrnodes)))
    rf$forest$nodepred <- do.call("cbind", lapply(rflist, 
                                                  function(x) padm0(x$forest$nodepred, nrnodes)))
    tree.dim <- dim(rf$forest$treemap)
    if (classRF) {
      rf$forest$treemap <- array(unlist(lapply(rflist, 
                                               function(x) apply(x$forest$treemap, 2:3, pad0, 
                                                                 nrnodes))), c(nrnodes, 2, ntree))
    }
    else {
      rf$forest$leftDaughter <- do.call("cbind", lapply(rflist, 
                                                        function(x) padm0(x$forest$leftDaughter, nrnodes)))
      rf$forest$rightDaughter <- do.call("cbind", lapply(rflist, 
                                                         function(x) padm0(x$forest$rightDaughter, nrnodes)))
    }
    rf$forest$ntree <- ntree
    if (classRF) 
      rf$forest$cutoff <- rflist[[1]]$forest$cutoff
  }
  else {
    rf$forest <- NULL
  }
  #
  #Tons of stuff removed here...
  #
  if (classRF) {
    rf$confusion <- NULL
    rf$err.rate <- NULL
    if (haveTest) {
      rf$test$confusion <- NULL
      rf$err.rate <- NULL
    }
  }
  else {
    rf$mse <- rf$rsq <- NULL
    if (haveTest) 
      rf$test$mse <- rf$test$rsq <- NULL
  }
  rf
}

#try smooth with spline 
# spline_fit <- smooth.spline(data_beforesign$day, data_beforesign$combo_r2_slope, spar = 0.4)  # Adjust spar for more or less smoothness
# 
# df <- data_beforesign %>%
#   mutate(smoothed_combined_score = predict(spline_fit, data_beforesign$day)$y)
# ggplot(df, aes(x = day)) +
#   geom_line(aes(y = combo_r2_slope), size = 1) +
#   geom_line(aes(y=smoothed_combined_score), col = 'red')


data_beforesign = temporary %>% 
   mutate(combo_r2_slope = slope*(r_s)) 



ggplot(data_beforesign, aes(x = day)) +
  geom_line(aes(y = combo_r2_slope), size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')



ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

y <- c(659, 613, 586, 642, 685, 695, 691, 733, 638, 708, 706, 691, 703, 712, 715, 700, 693, 682, 717, 711, 722, 700, 704, 704, 715, 691, 670, 684, 689, 711, 680, 692, 686, 710, 702, 699, 702, 715, 691, 670, 684, 689, 711, 673, 688, 699, 677, 701, 680, 692, 686, 710, 702, 699, 717, 691, 703, 712, 715, 700, 693, 682, 717, 711, 722, 700, 704, 704, 715, 691, 670, 684, 689, 711, 680, 692, 686, 710, 702, 699, 702, 706, 651, 712, 722, 734, 705, 714, 691, 704, 704, 669, 712, 715, 689, 715, 691, 670, 684, 689, 711, 651, 712, 722, 734, 705, 714, 691, 686, 676, 690, 693, 702, 694, 682, 693, 724, 693, 707, 684, 687, 705, 680, 680, 705, 680, 693, 700, 704, 704, 669, 712, 715, 689, 715, 691, 670, 684, 689, 711, 673, 678, 699, 677, 680, 682, 676, 690, 658, 675, 663, 667, 682, 673, 675, 656, 652, 563, 544, 542, 540, 532, 538, 505, 526, 565, 629, 720, 713, 720, 720, 773, 732, 740, 695, 689, 723, 685, 726, 710, 684, 693, 715, 692, 683, 712, 707, 693, 699, 717, 703, 687, 682, 690, 716, 708, 713, 700, 676, 708, 691, 717, 711, 722, 688, 695, 641, 666, 638, 639, 600, 635, 609, 653, 671, 649, 716, 708, 713, 700, 676, 708, 691, 717, 711, 722, 700, 704, 704, 669, 712, 715, 689, 715, 691, 670, 684, 689, 711, 673, 688, 699, 677, 701, 680, 692, 686, 710, 702, 699, 717, 691, 703, 712, 715, 700, 693, 682, 717, 711, 722, 700, 704, 704, 715, 691, 670, 684, 689, 711, 680, 692, 686, 710, 702, 699, 702, 715, 691, 670, 684, 689, 711, 673, 688, 699, 677, 701, 680, 692, 686, 710, 702, 699, 717, 691, 703, 712, 715, 700, 693, 682, 717, 711, 722, 700, 704, 704, 715, 691, 670, 684, 769, 767, 740, 752, 686, 710, 702, 699, 702, 706, 651, 712, 722, 734, 705, 714, 691, 704, 704, 669, 712, 715, 689, 715, 691, 670, 684, 689, 711, 704, 669, 712, 715, 689, 715, 691, 670, 684, 689, 711, 673, 688, 699, 677, 701, 680, 692, 686, 710, 702, 699, 717, 691, 703, 712, 715, 700, 693, 682, 717, 711, 722, 700, 704, 704, 715, 691, 670, 684, 689, 711, 680, 692, 686, 710, 702, 699, 702, 715, 691, 670, 684, 689, 711, 673, 688, 699, 677, 701, 680, 692, 686, 710, 702, 699, 665, 630, 808, 686, 787, 781, 796, 815, 786, 793, 664, 717, 691, 703, 712, 715, 700, 693, 682, 717, 711, 722, 700, 704, 704, 715, 691, 670, 684, 689, 711, 680, 692, 686, 710, 702, 699, 702, 706, 651, 712, 722, 734, 705, 714, 691, 704, 704, 669, 712, 715, 689, 715, 691, 670, 684, 689, 711)
y = data_beforesign$combo_r2_slope

#put a long lag 
lag       <- 100
threshold <- 3
influence <- 0
result <- ThresholdingAlgo(y,lag,threshold,influence)

# Plot result
par(mfrow = c(2,1),oma = c(2,2,0,0) + 0.1,mar = c(0,0,2,1) + 0.2)
plot(1:length(y),y,type="l",ylab="",xlab="") 
lines(1:length(y),result$avgFilter,type="l",col="cyan",lwd=2)
lines(1:length(y),result$avgFilter+threshold*result$stdFilter,type="l",col="green",lwd=2)
lines(1:length(y),result$avgFilter-threshold*result$stdFilter,type="l",col="green",lwd=2)
plot(result$signals,type="S",col="red",ylab="",xlab="",ylim=c(-1.5,1.5),lwd=2)


replace_0 <- function(vec, threshold) {
  # Convert any isolated 1 to 0
  vec_rle <- rle(vec)
  vec_rle$values[vec_rle$lengths == 1 & vec_rle$values == 1] <- 0
  vec <- inverse.rle(vec_rle)
  
  # Replace 0 sequences based on the threshold
  vec_rle <- rle(vec)
  for (i in seq_along(vec_rle$values)) {
    if (vec_rle$values[i] == 0 && vec_rle$lengths[i] <= threshold) {
      if (i > 1 && i < length(vec_rle$values) && vec_rle$values[i-1] == vec_rle$values[i+1]) {
        vec_rle$values[i] <- vec_rle$values[i-1]
      }
    }
  }
  
  # Reconstruct the vector
  result <- inverse.rle(vec_rle)
  
  return(result)
}


plot(replace_0(result$signals, threshold = 10),type="S",col="red",ylab="",xlab="",ylim=c(-1.5,1.5),lwd=2)

#my short ml 
list_slope <- as.list(Results_CSP[[i]]$estimate)
list_rs <- as.list(Results_CSP[[i]]$r.squared)

day = seq(rollwin,(lastdays-1),1)
slope = list_slope
r_s = list_rs
#k = nrow(site)

temporary <- data.frame(sitenewname = unique(Results_CSP[[1]]$sitenewname),
                        slope = unlist(slope), 
                        day = day, 
                        r_s = unlist(r_s))


optimization.ml <- function(data,
                            num_partitions = 10,
                            partition_size = 0.8) {

  # Train Random Forest models
  #https://stackoverflow.com/questions/19170130/combining-random-forests-built-with-different-training-sets-in-r
  #https://stackoverflow.com/questions/70842819/is-there-an-r-function-in-the-ranger-package-that-can-combine-two-random-forest
  library(caret) #for parititon 
  library(ranger) #for rf (it takes weights in it eventually )
  
  set.seed(9999)
  #create x partition 
  partitions <- list()
  #do not specify list T in partitioon 
  for (i in 1:num_partitions) {
    trainIndex <- createDataPartition(data$day, p = partition_size, times = 1, list = F)
    partitions[[i]] <- trainIndex
  }
  # Train the model on each partition and make predictions
  trainData = list()
  #testData = list()
  #rf_model_slope = list()
  rf_model = list()
  predictions.ranger.data=NULL
  
  for (i in 1:num_partitions) {
    # Subset data into training and testing sets
    trainData[[i]] <- data[partitions[[i]], ]
    rf_model[[i]] <- ranger::ranger(
      #formula = slope ~ day, #this one max slope and take into account r2 as weights 
      formula = r_s ~ day, #i ends up with this, just because it maximise r2 
      data = trainData[[i]],
      #case.weights = trainData[[i]]$r_s,    # Include R2 as weights
      importance = "impurity"    
    )
    #do not forget to make prediction for the full dataset ! 
    predictions.ranger.data.temp = data.frame(sitenewname = unique(data$sitenewname),
                                                                   day = data$day, 
               actual_slope = data$slope, 
               actual_r2 = data$r_s,
               predicted_slope = predict(rf_model[[i]], data = data)$predictions,
               iter = i)
    predictions.ranger.data = rbind(predictions.ranger.data, predictions.ranger.data.temp)
    
  }
  
  #use ranger to make prediction on the full sampel   
  #have same actual r2 and day and slope 
  test.pred = predictions.ranger.data %>% 
    group_by(sitenewname, day, actual_r2, actual_slope) %>% 
    summarise(pred.avg = mean(predicted_slope))

  return(test.pred)
}

test.pred = optimization.ml(temporary,
                            num_partitions = 20,
                            partition_size = 0.8)

ggplot(test.pred, aes(x = day)) +
  geom_line(aes(y = actual_r2, color = "Actual Slope"), size = 1, linetype = "dashed") +
  geom_line(aes(y = pred.avg), size = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')
data.prediction = test.pred


find.windows.optimal.valleys = function(data.prediction, #the one use for windows identification 
                                        rolling.temperature.data,#the data with rolling 
                                        quantilehigh = 0.9,
                                        minpeakdistance = 15,
                                        plot = FALSE){
  
  if(unique(data.prediction$sitenewname) != unique(rolling.temperature.data$sitenewname)){
    stop("Error: Not the same name between the prediction and the formatted data :(")
  }
  
  library(pracma)
  #use r pracma
  data_for_peaks <- data.prediction %>%
    ungroup() %>% 
    dplyr::select(day, pred.avg)
  #find peaks with pracma fuction 
  peaks <- findpeaks(abs(data_for_peaks$pred.avg), 
                     minpeakheight = quantile(abs(data_for_peaks$pred.avg), quantilehigh),
                     minpeakdistance = minpeakdistance) 
  peaks_df <- data_for_peaks[peaks[,2],] %>%
    mutate(peak_id = row_number())%>% 
    mutate(type = 'peak') 
  WinCloValley = data_for_peaks[peaks[,3],] %>%
    mutate(id = row_number())%>% 
    rename(windowsclose = 1) %>% 
    dplyr::select(windowsclose, id)
  WinOpValley = data_for_peaks[peaks[,4],] %>%
    mutate(id = row_number()) %>% 
    rename(windowsopen = 1,) %>% 
    dplyr::select(windowsopen, id)
  
  window_ranges_df = left_join(WinCloValley, WinOpValley) %>% 
    rename(windows.sequences.number = id)
  
  output_fit_summary.new = NULL
  for (z in 1:nrow(window_ranges_df)) {
    windowsopen <- window_ranges_df$windowsopen[z]
    windowsclose <- window_ranges_df$windowsclose[z]
    window_number <- window_ranges_df$windows.sequences.number[z]
    
    # Filter the rolling temperature data according to the current window range
    #does not change anything if sign different, just important the date smaller and bigger than
    climate_windows_best <- rolling.temperature.data %>%
      dplyr::filter(days.reversed <= windowsopen & days.reversed >= windowsclose) %>%
      group_by(sitenewname, year) %>%
      summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(window_number = window_number)  # Add the window number for identification
    
    #need to change year colnames :( 
    fit.best.model = bio_data %>% 
      rename(year = Year) %>% 
      left_join(climate_windows_best) %>%
      lm(myform.fin, data = .)
    
    intercept = tidy(fit.best.model)$estimate[1]
    intercept.se = tidy(fit.best.model)$std.error[1]
    estimate.model = tidy(fit.best.model)$estimate[2]
    estimate.se.model = tidy(fit.best.model)$std.error[2]
    pvalue.model = tidy(fit.best.model)$p.value[2]
    r2 = glance(fit.best.model)$r.squared
    AIC = glance(fit.best.model)$AIC
    nobs = glance(fit.best.model)$nobs
    nsequence.id = z
    
    output_fit_summary.temp.new <- data.frame(sitenewname = unique(bio_data$sitenewname),
                                              reference.day = refday,
                                              windows.size = rollwin,
                                              window.open = windowsopen,
                                              window.close = windowsclose,
                                              intercept = intercept,
                                              intercept.se = intercept.se,
                                              estimate = estimate.model,
                                              estimate.se = estimate.se.model,
                                              pvalue = pvalue.model,
                                              r2 = r2,
                                              AIC = AIC,
                                              nobs = nobs,
                                              nsequence.id = nsequence.id)
    output_fit_summary.new = rbind(output_fit_summary.new, output_fit_summary.temp.new)
    
  }
  if(plot == T){
    plott = ggplot(data_for_peaks, aes(x = day)) +
      geom_line(aes(y = pred.avg), size = 1) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = 'none')+
      geom_vline(xintercept = c(peaks_df$day), col = 'red')+
      geom_vline(xintercept = c(WinCloValley$windowsclose), col = 'green')+
      geom_vline(xintercept = c(WinOpValley$windowsopen), col = 'blue')
    print(plott)
  }
  return(output_fit_summary.new)
}



ttt = find.windows.optimal.valleys(data.prediction, #the one use for windows identification 
                                   rolling.temperature.data,#the data with rolling 
                                   quantilehigh = 0.9,
                                   minpeakdistance = 15, 
                                   plot = TRUE)

ttt







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