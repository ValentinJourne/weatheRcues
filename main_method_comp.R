initial.data.mastree <- read.csv('/Users/vjourne/Library/CloudStorage/Dropbox/Mastree/MASTREEplus_2022-02-03_V1.csv',stringsAsFactors = F)

library(tidyverse)
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
