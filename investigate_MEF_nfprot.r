rm(list=ls())

library(data.table)
library(tidyverse)
library(patchwork)
library(ggExtra)
library(zoo)

setwd('~/bioinfo/NFluor/novf-gia2025/')

# for normalizing graphs
normalize <-function(m){
  (m - min(m))/(max(m)-min(m))
}

LEN_CUTOFF <- 0 # suggested cutoff: originally 180 minimum length of beta barrel
BLUE_MEF_FILE <- 'blue/blue_MEF_FINAL.csv'
BLUE_COUNTS <- 'blue/blue_protcountsraw.tsv'
RED_MEF_FILE <- 'red/red_MEF_FINAL.csv'
RED_COUNTS <- 'red/red_protcountsraw.tsv'
FPBASE_DATA <- 'all_proteins.csv'

# create a list of the protein sequences from FPbase
fpbase_proteins <- fread(FPBASE_DATA) %>% 
  select(c("name", "seq", "slug", "states.0.brightness", "doi", 
           "genbank", "ipg_id", "states.0.brightness", "states.0.em_max", 
           "states.0.ex_max", "states.0.ext_coeff", "states.0.lifetime", 
           "states.0.maturation", "states.0.name", "states.0.pka", 
           "states.0.qy", "states.0.slug" ))
  
bluemefdata <- fread(BLUE_MEF_FILE) %>%  
  mutate(length = nchar(protein)) %>% 
  filter(length > LEN_CUTOFF) %>% 
  filter(str_detect(protein, "^M")) # start with Met
redmefdata <- fread(RED_MEF_FILE) %>%  
  mutate(length = nchar(protein)) %>% 
  filter(length > LEN_CUTOFF) %>% 
  filter(str_detect(protein, "^M")) # start with Met
# for checking for duplicates
# n_blue_occur <- data.frame(table(bluemefdata$protein))
# n_red_occur <- data.frame(table(redmefdata$protein))

allmefs <- full_join(bluemefdata, redmefdata, by = 'protein') %>% # this is an outer join
  mutate(MEF.x = na.aggregate(na_if(MEF.x, NA), FUN = min)) %>% 
  mutate(MEF.y = na.aggregate(na_if(MEF.y, NA), FUN = min)) # make the NAs the minimum value

# generate list of the proteins from FPbase which exist in this data
blue_found_prots <- data.frame()
red_found_prots <- data.frame()
all_found_prots <- data.frame()

for (variant in bluemefdata$protein)
{
  if (variant %in% fpbase_proteins$seq){
    found <- fpbase_proteins[match(variant, fpbase_proteins$seq)]
    blue_found_prots <- rbind(blue_found_prots, found)
  }
}
blue_found_prots <- merge(x = blue_found_prots, y = bluemefdata[ , c("protein", "MEF")], 
      by.x = "seq", by.y = "protein", all.x=TRUE)

for (variant in redmefdata$protein)
{
  if (variant %in% fpbase_proteins$seq){
    found <- fpbase_proteins[match(variant, fpbase_proteins$seq)]
    red_found_prots <- rbind(red_found_prots, found)
  }
}
red_found_prots <- merge(x = red_found_prots, y = redmefdata[ , c("protein", "MEF")], 
      by.x = "seq", by.y = "protein", all.x=TRUE)

length(intersect(blue_found_prots$seq, red_found_prots$seq))

for (variant in allmefs$protein)
{
  if (variant %in% fpbase_proteins$seq){
    found <- fpbase_proteins[match(variant, fpbase_proteins$seq)]
    all_found_prots <- rbind(all_found_prots, found)
  }
}
all_found_prots <- merge(x = all_found_prots, y = bluemefdata[ , c("protein", "MEF")], 
      by.x = "seq", by.y = "protein", all.x=TRUE) %>% 
  rename("MEF.blue" = "MEF")
all_found_prots <- merge(x = all_found_prots, y = redmefdata[ , c("protein", "MEF")], 
      by.x = "seq", by.y = "protein", all.x=TRUE) %>% 
  rename("MEF.red" = "MEF")

normalize <-function(myvec){
  (myvec - min(myvec))/(max(myvec)-min(myvec))
}

all_found_prots <- all_found_prots %>% 
  mutate(normMEF.blue =  (MEF.blue - min(MEF.blue, na.rm = TRUE)) / (max(MEF.blue, na.rm = TRUE) - min(MEF.blue, na.rm = TRUE))) %>% 
  mutate(normMEF.red =  (MEF.red - min(MEF.red, na.rm = TRUE)) / (max(MEF.red, na.rm = TRUE) - min(MEF.red, na.rm = TRUE))) %>% 
  mutate(predicted.color = ifelse(is.na(normMEF.blue), "Red", 
                                  ifelse(is.na(normMEF.red), "Blue", 
                                         ifelse(normMEF.blue > normMEF.red, "Blue", 
                                                ifelse(normMEF.red > normMEF.blue, "Red", "Inconclusive")))))

bluecounts <- fread(BLUE_COUNTS)
redcounts <- fread(RED_COUNTS)
totalnum_blues <- sum(bluecounts$totalreads)
totalnum_reds <- sum(redcounts$totalreads)
sum(bluecounts$totalfrac)
sum(redcounts$totalfrac)

# bluecounts <- bluecounts %>% 
  # mutate(bin1raw = bin1 * totalnum_blues)

redvsblue <- allmefs %>% 
  select(protein, MEF.x, MEF.y) %>% 
  rename(RED_MEFL = MEF.y) %>% 
  rename(BLUE_MEFL = MEF.x)

redvsblue <- left_join(redvsblue, bluecounts, by= "protein", relationship = "many-to-many") %>% 
  rename(blue_bin1 = bin1) %>% 
  rename(blue_bin2 = bin2) %>% 
  rename(blue_bin3 = bin3) %>% 
  rename(blue_bin4 = bin4) %>% 
  rename(blue_bin5 = bin5) %>% 
  rename(blue_bin6 = bin6) %>% 
  rename(blue_bin7 = bin7) %>% 
  rename(blue_bin8 = bin8) %>% 
  rename(blue_bin9 = bin9) %>% 
  rename(blue_tot = totalreads)

redvsblue <- left_join(redvsblue, redcounts, by= "protein", relationship = "many-to-many") %>% 
  rename(red_bin1 = bin1) %>% 
  rename(red_bin2 = bin2) %>% 
  rename(red_bin3 = bin3) %>% 
  rename(red_bin4 = bin4) %>% 
  rename(red_bin5 = bin5) %>% 
  rename(red_bin6 = bin6) %>% 
  rename(red_bin7 = bin7) %>% 
  rename(red_bin8 = bin8) %>% 
  rename(red_bin9 = bin9) %>% 
  rename(red_tot = totalreads)

redvsblue <- redvsblue %>% group_by(protein) %>%
  summarise(blue_bin1 = sum(blue_bin1, na.rm = T),
            blue_bin2 = sum(blue_bin2, na.rm = T), 
            blue_bin3 = sum(blue_bin3, na.rm = T), 
            blue_bin4 = sum(blue_bin4, na.rm = T), 
            blue_bin5 = sum(blue_bin5, na.rm = T), 
            blue_bin6 = sum(blue_bin6, na.rm = T), 
            blue_bin7 = sum(blue_bin7, na.rm = T), 
            blue_bin8 = sum(blue_bin8, na.rm = T), 
            blue_bin9 = sum(blue_bin9, na.rm = T),
            blue_tot = sum(blue_tot, na.rm = T), 
            blue_MEFL = mean(BLUE_MEFL, na.rm = T),
            red_bin1 = sum(red_bin1, na.rm = T),
            red_bin2 = sum(red_bin2, na.rm = T), 
            red_bin3 = sum(red_bin3, na.rm = T), 
            red_bin4 = sum(red_bin4, na.rm = T), 
            red_bin5 = sum(red_bin5, na.rm = T), 
            red_bin6 = sum(red_bin6, na.rm = T), 
            red_bin7 = sum(red_bin7, na.rm = T), 
            red_bin8 = sum(red_bin8, na.rm = T), 
            red_bin9 = sum(red_bin9, na.rm = T),
            red_tot = sum(red_tot, na.rm = T), 
            red_MEFL = mean(RED_MEFL, na.rm = T))

write.csv(redvsblue,"redvsblue.csv", row.names = FALSE)

# PLOTS ########################################################

blue_len_hist <- bluemefdata %>%
  ggplot(aes(x=length)) + 
  geom_histogram(color = 'black', fill = 'blue') + 
  ylim(range(0,5000)) + 
  # scale_y_continuous(name = "Count", breaks = c(0, 4, 32, 256, 2048), labels = c("0", "4", "32", "256", "2048"), trans = 'log2') + 
  xlab("Protein Length") + 
  ylab("Count") + 
  theme_bw()
red_len_hist <- redmefdata %>%
  ggplot(aes(x=length)) + 
  geom_histogram(color = 'black', fill = 'red') + 
  ylim(range(0,5000)) + 
  # scale_y_continuous(name = "Count", breaks = c(0, 4, 32, 256, 2048), labels = c("0", "4", "32", "256", "2048"), trans = 'log2') + 
  xlab("Protein Length") + 
  ylab("") + 
  theme_bw()

bluebrightnesses <- bluemefdata %>% 
  mutate(MEFnorm = normalize(MEF)) %>% 
  ggplot(aes(x=MEFnorm)) + 
  geom_histogram(color = 'black', fill = 'blue') + 
  ylim(range(0,2050)) + 
  scale_y_continuous(name = "Count (log)", breaks = c(0, 4, 32, 256, 2048), labels = c("0", "4", "32", "256", "2048"), trans = 'log2') + 
  # ggtitle("Blue brightnesses") + 
  xlab("Normalized Blue Intensity") + 
  ylab("Count (log)") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))  
redbrightnesses <- redmefdata %>% 
  mutate(MEFnorm = normalize(MEF)) %>% 
  ggplot(aes(x=MEFnorm)) + 
  geom_histogram(color = 'black', fill = 'red') + 
  ylim(range(0,2050)) + 
  scale_y_continuous(name = "", breaks = c(0, 4, 32, 256, 2048), labels = c("0", "4", "32", "256", "2048"), trans = 'log2') + 
  # ggtitle("Red brightnesses") + 
  xlab("Normalized Red Intensity") + 
  ylab("") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))  

# normalize this plot for display.

mefplotnorm <- allmefs %>% 
  mutate(infpbase = protein %in% fpbase_proteins$seq) %>% 
  mutate(MEF.x = normalize(MEF.x)) %>% 
  mutate(MEF.y = normalize(MEF.y)) %>% 
  arrange(infpbase) %>% 
  ggplot(aes(x=MEF.x, y=MEF.y, colour = infpbase)) + 
  annotate("rect",xmin=0,xmax=0.25,ymin=0.25,ymax=1,alpha=0.3,fill="red") +
  annotate("rect",xmin=0.25,xmax=1,ymin=0,ymax=0.25,alpha=0.3,fill="blue") +
  annotate("rect",xmin=0.25,xmax=1,ymin=0.25,ymax=1,alpha=0.3,fill="purple") +
  scale_color_manual(values=c("black", "black")) + #made these the same color, but change later to show fpbase proteins
  geom_point(size = 1, alpha = 0.5) +
  theme_bw() + 
  xlab("Normalized Blue Intensity") + 
  ylab("Normalized Red Intensity") +
  theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.margin = margin(1, 0, 0, 0, "cm")) # trying to make patchwork look nicer

# original version
# mefplot <- allmefs %>% 
  # mutate(infpbase = protein %in% fpbase_proteins$seq) %>% 
  # arrange(infpbase) %>% 
  # ggplot(aes(x=MEF.x, y=MEF.y, colour = infpbase)) + 
  # scale_color_manual(values=c("black", "red")) +
  # geom_point() +
  # theme_bw() + 
  # xlab("Blue MEF") + 
  # ylab("Red MEF") +
  # theme(legend.position = "none") 
# ggMarginal(mefplot, type = "density", size = 5) # plot mefplot with density

# now time to separate the proteins found only in the red or blue sorts and see if there is interesting stuff to correlate there.. 

# blue_found_prots %>% 
  # ggplot(aes(x = states.0.em_max, y = states.0.brightness, colour = MEF)) + 
  # geom_point()
# red_found_prots %>% 
  # ggplot(aes(x = states.0.em_max, y = states.0.brightness, colour = MEF)) + 
  # geom_point()

blue_found <- all_found_prots %>% 
  mutate(blue = ifelse(is.na(normMEF.blue), "No", ifelse(normMEF.blue < 0.25, "No", "Yes"))) %>% 
  arrange(blue) %>% 
  ggplot(aes(x = states.0.ex_max, y = states.0.em_max, colour = blue, alpha = blue)) + 
  scale_color_manual(values=c("black", "blue")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  geom_point(size = 1) +
  theme_bw() + 
  ylab("Emission wavelength (nm)") +
  xlab("Excitation wavelength (nm)") + 
  theme(legend.position = "none") 

red_found <- all_found_prots %>% 
  mutate(red = ifelse(is.na(normMEF.red), "No", ifelse(normMEF.red < 0.25, "No", "Yes"))) %>% 
  arrange(red) %>% 
  ggplot(aes(x = states.0.ex_max, y = states.0.em_max, colour = red, alpha = red)) + 
  scale_color_manual(values=c("black", "red")) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  geom_point(size = 1) +
  theme_bw() + 
  # ylab("Emission wavelength (nm)") +
  ylab("") + 
  xlab("Excitation wavelength (nm)") + 
  theme(legend.position = "none")

# bluebrightnesses + redbrightnesses - ggMarginal(mefplotnorm, type = "density", size = 5) + plot_layout(ncol = 1) + plot_annotation(tag_levels = "A")

# Generate Plots
pdf(file = "length_hist.pdf", width = 8, height = 4) # width and height in inches
blue_len_hist + red_len_hist + plot_annotation(tag_levels = "A")
dev.off()

pdf(file = "brightness_hist.pdf", width = 8, height = 4) 
bluebrightnesses + redbrightnesses + plot_annotation(tag_levels = "A")
dev.off()

# pdf(file = "red_vs_blue.pdf", width = 8, height = 8) 
# mefplotnorm
png(file = "red_vs_blue.png", width = 5.5, height = 5.5, res = 1200, units = "in") 
ggMarginal(mefplotnorm, type = "density", size = 5) # plot mefplot with density
dev.off()

pdf(file = "FPbase_proteins.pdf", width = 8, height = 4) 
blue_found + red_found + plot_annotation(tag_levels = "A")
dev.off()