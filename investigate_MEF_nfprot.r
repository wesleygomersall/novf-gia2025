rm(list=ls())

library(data.table)
library(tidyverse)

setwd('~/bioinfo/NFluor/novf-gia2025/')

LEN_CUTOFF <- 180 # suggested cutoff: minimum length of beta barrel
BLUE_MEF_FILE <- 'blue/blue_MEF_FINAL.csv'
RED_MEF_FILE <- 'red/red_MEF_FINAL.csv'

bluemefdata <- fread(BLUE_MEF_FILE) %>%  
  mutate(length = nchar(protein)) %>% 
  filter(length > LEN_CUTOFF) %>% 
  filter(str_detect(protein, "^M")) # start with Met
redmefdata <- fread(RED_MEF_FILE) %>%  
  mutate(length = nchar(protein)) %>% 
  filter(length > LEN_CUTOFF) %>% 
  filter(str_detect(protein, "^M")) # start with Met

bluemefdata %>%
  ggplot(aes(x=length)) + 
  geom_histogram() + 
  ggtitle("Blue_MEF protein lengths after filtering") + 
  theme_bw()
redmefdata %>%
  ggplot(aes(x=length)) + 
  geom_histogram() + 
  ggtitle("Red_MEF protein lengths after filtering") + 
  theme_bw()

bluemefdata %>% 
  ggplot(aes(x=MEF)) + 
  geom_histogram() + 
  ggtitle("Blue_MEF brightness histogram after filtering") + 
  theme_bw()
redmefdata %>% 
  ggplot(aes(x=MEF)) + 
  geom_histogram() + 
  ggtitle("Red_MEF brightness histogram after filtering") + 
  theme_bw()
