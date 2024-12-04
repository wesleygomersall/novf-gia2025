library(data.table)
library(tidyverse)

LEN_CUTOFF <- 180 # suggested cutoff: minimum length of beta barrel
FILE <- ''

mydata <- fread(FILE) %>%  
  mutate(length = nchar(protein)) %>% 
  filter(length > LEN_CUTOFF) %>% 
  filter(str_detect(protein, "^M")) # start with Met

mydata %>%
  ggplot(aes(x=length)) + 
  geom_histogram() + 
  ggtitle("Histogram of protein lengths after filtering") + 
  theme_bw()

mydata %>% 
  ggplot(aes(x=MEF)) + 
  geom_histogram() + 
  ggtitle("Histogram of protein MEF after filtering") + 
  theme_bw()
