library(ggplot2)
library(dplyr)
library(rstudioapi)
library(stringr)
library(cowplot)
library(DescTools)
library(ggrepel)

####### LOAD TRIALS #######

tracks <- choose.files(default = "", caption = "Select files",
                      multi = TRUE, filters = Filters,
                      index = nrow(Filters))

background <- choose.files(default = "", caption = "Select files",
                      multi = TRUE, filters = Filters,
                      index = nrow(Filters))

list1 <- tracks[1:10]
list2 <- tracks[11:20]
list3 <- tracks[21:29]

blist1 <- background[1:10]
blist2 <- background[11:20]
blist3 <- background[21:29]

tracks_1v <- list()
tracks_2v <- list()
tracks_3v <- list()

background_1v <- list()
background_2v <- list()
background_3v <- list()

for (file_path in list1){
  data <- read.csv(file_path, header = TRUE)
  tracks_1v <- append(tracks_1v, list(data))
}

for (file_path in list2){
  data <- read.csv(file_path, header = TRUE)
  tracks_2v <- append(tracks_2v, list(data))
}

for (file_path in list3){
  data <- read.csv(file_path, header = TRUE)
  tracks_3v <- append(tracks_3v, list(data))
}

for (file_path in blist1){
  data <- read.csv(file_path, header = TRUE)
  background_1v <- append(background_1v, list(data))
}

for (file_path in blist2){
  data <- read.csv(file_path, header = TRUE)
  background_2v <- append(background_2v, list(data))
}

for (file_path in blist3){
  data <- read.csv(file_path, header = TRUE)
  background_3v <- append(background_3v, list(data))
}

setwd(dirname(getActiveDocumentContext()$path))
getwd()

############  ANALYSIS  ###########################

############# 1V ###############################

EFS_1V <- function(t1, b1){
  #TRACKMATE DATA
  csv_frame <- t1
  CELLS <- csv_frame[order(csv_frame$TRACK_ID, csv_frame$FRAME),]
  CELLS <- CELLS[order(CELLS$FRAME),]
  #BACKGROUND SPOTS
  spot_data <- b1
  SPOTS <- spot_data[-c(1),]
  SPOTS <- as.data.frame(apply(SPOTS, 2, as.numeric))
  SPOTS <- SPOTS %>% mutate(BG_Avg = rowMeans(SPOTS))
  SPOTS$FRAME <- (1:nrow(SPOTS))-1
  SPOTS <- SPOTS[, c(101,102)]
  #COMBINED CELL AND SPOT DATA
  DATA <- left_join(CELLS, SPOTS)
  DATA <- DATA[-c(6:9)]
  #BACKGROUND SUBTRACTION
  DATA$Avg_after_BGsub <- NA
  DATA <- DATA %>%
    mutate(Avg_after_BGsub = MEAN_INTENSITY_CH1 - BG_Avg)
  #CTR & EXP DATA SPLIT
  CTR <- subset(DATA, FRAME<25)
  EXP <- subset(DATA, FRAME >= 25)
  
  ##################### 1V PRE/CONTROL #################################################
  CTR <- CTR[order(CTR$TRACK_ID, CTR$FRAME),]
  
  #INTENSITY IN 80TH PERCENTILE
  quantile(CTR$Avg_after_BGsub, probs=.80)
  Q <- quantile(CTR$Avg_after_BGsub, probs=.80)
  
  
  
  CTR$Top20 <- NA
  CTR$Top20 <- ifelse(CTR$Avg_after_BGsub >= Q, "Top20val",
                      ifelse(CTR$Avg_after_BGsub, "No"))
  
  #FIRST POINT
  ctr_first <- CTR %>% group_by(TRACK_ID) %>%
    slice_head(n=1)
  
  ctr_first <- ctr_first[c("TRACK_ID", "Avg_after_BGsub")]
  ctr_first <- ctr_first %>% rename(point1 = Avg_after_BGsub)
  ctr <- left_join(CTR, ctr_first)
  
  #PERCENT CHANGE FROM FIRST POINT
  ctr <- ctr %>% group_by(TRACK_ID) %>% mutate(Change = Avg_after_BGsub-point1,
                                               PctChg = (Change/point1)*100)
  
  ctr_FindBase <- ctr
  ctr_FindBase <- ctr_FindBase[order(ctr_FindBase$TRACK_ID, ctr_FindBase$FRAME),]
  
  #https://stackoverflow.com/questions/61770151/counting-rows-until-a-condition-is-met-in-r-nas-before-the-condition-is-met
  ctr_FindBase <- ctr_FindBase %>%
    group_by(TRACK_ID) %>%
    group_by(temp = cumsum(PctChg >= 10), .add = TRUE) %>%
    mutate(out = row_number()-1) %>%
    group_by(TRACK_ID) %>%
    mutate(out = replace(out, row_number()<which.max(PctChg >= 10), NA))
  
  #https://stackoverflow.com/questions/24477748/r-count-na-by-group
  ctr_FindBase_count <- aggregate(out ~ TRACK_ID, ctr_FindBase, function(x) {sum(is.na(x))}, na.action=NULL)
  
  ctr_Lth <- ctr %>% group_by(TRACK_ID) %>% mutate(lgth = length(TRACK_ID))
  ctr_Lth <- ctr_Lth[c("TRACK_ID", "lgth")]
  ctr_Lth <- ctr_Lth %>% group_by(TRACK_ID) %>% slice_head(n=1)
  
  ctr_FindBase_Bn <- left_join(ctr_FindBase_count, ctr_Lth)
  
  for(i in 1:nrow(ctr_FindBase_Bn)){
    if(ctr_FindBase_Bn[i, "out"] == 0){
      ctr_FindBase_Bn[i, "out"] = ctr_FindBase_Bn[i, "lgth"]
    }
  }
  
  ctr_FindBase_Bn <- ctr_FindBase_Bn[-c(3)]
  ctr_FindBase_Bn <- ctr_FindBase_Bn %>% rename(meancount = out)
  
  #MEAN OF FIRST N COLUMNS - BASELINE AVG
  #https://stackoverflow.com/questions/64203431/subset-number-of-rows-per-group-based-on-a-value-in-another-dataframe
  ctr_base <- ctr_FindBase %>%
    left_join(ctr_FindBase_Bn) %>%
    group_by(TRACK_ID) %>%
    arrange(TRACK_ID, FRAME) %>%
    slice(1:meancount[1])
  
  ctr_base2 <- aggregate(Avg_after_BGsub ~ TRACK_ID, ctr_base, FUN=mean)
  ctr_base2 <- ctr_base2 %>% rename(base_avg = Avg_after_BGsub)
  
  ctr_base3 <- full_join(ctr_FindBase, ctr_base2)
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  
  #NORMALIZATION
  ctr_base3 <- ctr_base3 %>% mutate(baseline = Avg_after_BGsub/base_avg)
  ctr_base3 <- ctr_base3 %>% mutate(time = (((FRAME-24)*5)/60))
  
  
  #LAG PCT CHANGE
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Prev_int = lag(baseline),
                                                           Change_lag = baseline-Prev_int,
                                                           PctChg_lag = (Change_lag/Prev_int)*100)
  
  #LEAD INTENSITY GROUP
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Lead_int = lead(baseline))
  
  ctr_base3$Lead_int[is.na(ctr_base3$Lead_int)] <- 0
  
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Change_lead = (baseline - Lead_int))
  
  #DETERMINE SPIKE
  ctr_base3$Peak <- NA
  ctr_base3$minchange <- NA
  ctr_base3$PctChg_lag[is.na(ctr_base3$PctChg_lag)] <- 0
  ctr_base3$Change_lag[is.na(ctr_base3$Change_lag)] <- 0
  
  for (i in 1:nrow(ctr_base3)){
    if (ctr_base3[i, "Change_lag"] >0 && ctr_base3[i, "Change_lead"] >0){
      ctr_base3[i, "Peak"] = "Peak"
    } else{
      ctr_base3[i, "Peak"] = "No"
    }
  }
  
  counter = 0
  for (i in 1:nrow(ctr_base3)){
    if(ctr_base3[i, "Peak"] == "No"){
      ctr_base3[i, "minchange"] = counter
    } else{
      counter = counter+1
    }
  }
  
  ctr_spikedata <- ctr_base3 %>% mutate(leadingPCT = lead(PctChg_lag))
  ctr_spikedata1 <- aggregate(leadingPCT ~ minchange, ctr_spikedata, FUN=sum)
  
  ctr_spikedata2 <- full_join(ctr_spikedata1, ctr_base3)
  ctr_spikedata2 <- ctr_spikedata2[order(ctr_spikedata2$TRACK_ID, ctr_spikedata2$FRAME),]
  ctr_spikedata2<- ctr_spikedata2[c(3:24, 1, 2)]
  
  ctr_spikedata2 <- ctr_spikedata2 %>% mutate(PKpct = lag(leadingPCT))
  ctr_spikedata2$itsaspike <- NA
  ctr_spikedata2$PKpct[is.na(ctr_spikedata2$PKpct)] <- 0
  
  for (i in 1:nrow(ctr_spikedata2)){
    if(ctr_spikedata2[i, "Peak"] == "Peak"){
      if(ctr_spikedata2[i, "PKpct"] >= 10 || ctr_spikedata2[i, "PctChg_lag"] >= 10){
        ctr_spikedata2[i, "itsaspike"] = "SPIKE"
      } else{
        ctr_spikedata2[i, "itsaspike"] = "Not"
      }
    }else{
      ctr_spikedata2[i, "itsaspike"] = "Not"
    }
  }
  
  ###OUTLIER TEST####
  
  summary(ctr_spikedata2$baseline)
  ctr_mean <- mean(ctr_spikedata2$baseline)
  ctr_sd <- sd(ctr_spikedata2$baseline)
  
  ctr_spikedata2 <- ctr_spikedata2 %>% mutate(Zscore = (baseline-ctr_mean)/(ctr_sd))
  
  ctr_spikedata2 <- ctr_spikedata2[!(ctr_spikedata2$TRACK_ID %in% ctr_spikedata2$TRACK_ID[ctr_spikedata2$Zscore > 3]),]
  ctr_spikedata2 <- ctr_spikedata2[!(ctr_spikedata2$TRACK_ID %in% ctr_spikedata2$TRACK_ID[ctr_spikedata2$Zscore < -3]),]
  
  
  print(length(unique(ctr_base3$TRACK_ID)))
  print(length(unique(ctr_spikedata2$TRACK_ID)))
  
  #COUNT SPIKES/CELL (From Lag)
  ctr_numspikes <- sum(str_count(ctr_spikedata2$itsaspike, "SPIKE"))
  ctr_NUM <- count(ctr_spikedata2, itsaspike)
  
  ctr_celltotal <- length(unique(ctr_spikedata2$TRACK_ID))
  
  ctr_Activity <- (ctr_numspikes)/(ctr_celltotal*2)
  
  ##################### 1V POST/EXPERIMENTAL #################################################
  EXP <- EXP[order(EXP$TRACK_ID, EXP$FRAME),]
  
  #INTENSITY IN 80TH PERCENTILE
  quantile(EXP$Avg_after_BGsub, probs=.80)
  Q <- quantile(EXP$Avg_after_BGsub, probs=.80)
  
  
  
  EXP$Top20 <- NA
  EXP$Top20 <- ifelse(EXP$Avg_after_BGsub >= Q, "Top20val",
                      ifelse(EXP$Avg_after_BGsub, "No"))
  
  #FIRST POINT
  exp_first <- EXP %>% group_by(TRACK_ID) %>%
    slice_head(n=1)
  
  exp_first <- exp_first[c("TRACK_ID", "Avg_after_BGsub")]
  exp_first <- exp_first %>% rename(point1 = Avg_after_BGsub)
  EXP <- left_join(EXP, exp_first)
  
  #PERCENT CHANGE FROM FIRST POINT
  EXP <- EXP %>% group_by(TRACK_ID) %>% mutate(Change = Avg_after_BGsub-point1,
                                               PctChg = (Change/point1)*100)
  
  exp_FindBase <- EXP
  exp_FindBase <- exp_FindBase[order(exp_FindBase$TRACK_ID, exp_FindBase$FRAME),]
  
  #https://stackoverflow.com/questions/61770151/counting-rows-until-a-condition-is-met-in-r-nas-before-the-condition-is-met
  exp_FindBase <- exp_FindBase %>%
    group_by(TRACK_ID) %>%
    group_by(temp = cumsum(PctChg >= 10), .add = TRUE) %>%
    mutate(out = row_number()-1) %>%
    group_by(TRACK_ID) %>%
    mutate(out = replace(out, row_number()<which.max(PctChg >= 10), NA))
  
  #https://stackoverflow.com/questions/24477748/r-count-na-by-group
  exp_FindBase_count <- aggregate(out ~ TRACK_ID, exp_FindBase, function(x) {sum(is.na(x))}, na.action=NULL)
  
  exp_Lth <- EXP %>% group_by(TRACK_ID) %>% mutate(lgth = length(TRACK_ID))
  exp_Lth <- exp_Lth[c("TRACK_ID", "lgth")]
  exp_Lth <- exp_Lth %>% group_by(TRACK_ID) %>% slice_head(n=1)
  
  exp_FindBase_Bn <- left_join(exp_FindBase_count, exp_Lth)
  
  for(i in 1:nrow(exp_FindBase_Bn)){
    if(exp_FindBase_Bn[i, "out"] == 0){
      exp_FindBase_Bn[i, "out"] = exp_FindBase_Bn[i, "lgth"]
    }
  }
  
  exp_FindBase_Bn <- exp_FindBase_Bn[-c(3)]
  exp_FindBase_Bn <- exp_FindBase_Bn %>% rename(meancount = out)
  
  #MEAN OF FIRST N COLUMNS - BASELINE AVG
  #https://stackoverflow.com/questions/64203431/subset-number-of-rows-per-group-based-on-a-value-in-another-dataframe
  exp_base <- exp_FindBase %>%
    left_join(exp_FindBase_Bn) %>%
    group_by(TRACK_ID) %>%
    arrange(TRACK_ID, FRAME) %>%
    slice(1:meancount[1])
  
  exp_base2 <- aggregate(Avg_after_BGsub ~ TRACK_ID, exp_base, FUN=mean)
  exp_base2 <- exp_base2 %>% rename(base_avg = Avg_after_BGsub)
  
  exp_base3 <- full_join(exp_FindBase, exp_base2)
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  
  #NORMALIZATION
  exp_base3 <- exp_base3 %>% mutate(baseline = Avg_after_BGsub/base_avg)
  exp_base3 <- exp_base3 %>% mutate(time = (((FRAME-25)*5)/60))
  

  #LAG PCT CHANGE
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Prev_int = lag(baseline),
                                                           Change_lag = baseline-Prev_int,
                                                           PctChg_lag = (Change_lag/Prev_int)*100)
  
  #LEAD INTENSITY GROUP
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Lead_int = lead(baseline))
  
  exp_base3$Lead_int[is.na(exp_base3$Lead_int)] <- 0
  
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Change_lead = (baseline - Lead_int))
  
  #DETERMINE SPIKE
  exp_base3$Peak <- NA
  exp_base3$minchange <- NA
  exp_base3$PctChg_lag[is.na(exp_base3$PctChg_lag)] <- 0
  exp_base3$Change_lag[is.na(exp_base3$Change_lag)] <- 0
  
  for (i in 1:nrow(exp_base3)){
    if (exp_base3[i, "Change_lag"] >0 && exp_base3[i, "Change_lead"] >0){
      exp_base3[i, "Peak"] = "Peak"
    } else{
      exp_base3[i, "Peak"] = "No"
    }
  }
  
  counter = 0
  for (i in 1:nrow(exp_base3)){
    if(exp_base3[i, "Peak"] == "No"){
      exp_base3[i, "minchange"] = counter
    } else{
      counter = counter+1
    }
  }
  
  exp_spikedata <- exp_base3 %>% mutate(leadingPCT = lead(PctChg_lag))
  exp_spikedata1 <- aggregate(leadingPCT ~ minchange, exp_spikedata, FUN=sum)
  
  exp_spikedata2 <- full_join(exp_spikedata1, exp_base3)
  exp_spikedata2 <- exp_spikedata2[order(exp_spikedata2$TRACK_ID, exp_spikedata2$FRAME),]
  exp_spikedata2<- exp_spikedata2[c(3:24, 1, 2)]
  
  exp_spikedata2 <- exp_spikedata2 %>% mutate(PKpct = lag(leadingPCT))
  exp_spikedata2$itsaspike <- NA
  exp_spikedata2$PKpct[is.na(exp_spikedata2$PKpct)] <- 0
  
  for (i in 1:nrow(exp_spikedata2)){
    if(exp_spikedata2[i, "Peak"] == "Peak"){
      if(exp_spikedata2[i, "PKpct"] >= 10 || exp_spikedata2[i, "PctChg_lag"] >= 10){
        exp_spikedata2[i, "itsaspike"] = "SPIKE"
      } else{
        exp_spikedata2[i, "itsaspike"] = "Not"
      }
    }else{
      exp_spikedata2[i, "itsaspike"] = "Not"
    }
  }
  
  ###OUTLIER TEST####
  
  summary(exp_spikedata2$baseline)
  exp_mean <- mean(exp_spikedata2$baseline)
  exp_sd <- sd(exp_spikedata2$baseline)
  
  exp_spikedata2 <- exp_spikedata2 %>% mutate(Zscore = (baseline-exp_mean)/(exp_sd))
  
  exp_spikedata2 <- exp_spikedata2[!(exp_spikedata2$TRACK_ID %in% exp_spikedata2$TRACK_ID[exp_spikedata2$Zscore > 3]),]
  exp_spikedata2 <- exp_spikedata2[!(exp_spikedata2$TRACK_ID %in% exp_spikedata2$TRACK_ID[exp_spikedata2$Zscore < -3]),]
  
  print(length(unique(exp_base3$TRACK_ID)))
  print(length(unique(exp_spikedata2$TRACK_ID)))
  
  #COUNT SPIKES/CELL (From Lag)
  exp_numspikes <- sum(str_count(exp_spikedata2$itsaspike, "SPIKE"))
  exp_NUM <- count(exp_spikedata2, itsaspike)
  
  exp_celltotal <- length(unique(exp_spikedata2$TRACK_ID))
  
  exp_Activity <- (exp_numspikes)/(exp_celltotal*2)
  
  
  ######## INDIVIDUAL SPIKE DATA ###################
  
  Cspike <- subset(ctr_spikedata2, itsaspike == "SPIKE")
  Cspike2 <- Cspike[c('TRACK_ID', 'itsaspike')]
  Cspike3 <- Cspike2 %>% group_by(TRACK_ID) %>%
    mutate(counts = str_count(itsaspike, "SPIKE"))
  Cspike4 <- aggregate(counts ~ TRACK_ID, Cspike3, FUN=sum)
  
  Espike <- subset(exp_spikedata2, itsaspike == "SPIKE")
  Espike2 <- Espike[c('TRACK_ID', 'itsaspike')]
  Espike3 <- Espike2 %>% group_by(TRACK_ID) %>%
    mutate(counts = str_count(itsaspike, "SPIKE"))
  Espike4 <- aggregate(counts ~ TRACK_ID, Espike3, FUN=sum)
  
  Cspikes <<- Cspike4
  Espikes <<- Espike4
  
  
  ############ MAKE DATA FRAME #################
  STATS <- c("1V_CTR_SPM", "1V_EXP_SPM", "1V_CTR_SPIKE_TOTAL", "1V_CTR_CELL_TOTAL", "1V_EXP_SPIKE_TOTAL", "1V_EXP_CELL_TOTAL")
  COUNT <- c(ctr_Activity, exp_Activity, ctr_numspikes, ctr_celltotal, exp_numspikes, exp_celltotal)
  
  Stats_1V_funcdata <<- data.frame(STATS, COUNT)
  
}

activity_list_1v <- list()
ctr_spikelist_1v <- list()
exp_spikelist_1v <- list()

for (i in 1:length(tracks_1v)) {
  results <- EFS_1V(tracks_1v[[i]], background_1v[[i]])
  activity_list_1v[[i]] <- assign(paste0("v1_ACTIVITY_", i, sep = ''), Stats_1V_funcdata, envir = .GlobalEnv)
  ctr_spikelist_1v[[i]] <- assign(paste0("v1_CTR_Spikes_", i, sep = ''), Cspikes, envir = .GlobalEnv)
  exp_spikelist_1v[[i]] <- assign(paste0("v1_EXP_Spikes_", i, sep = ''), Espikes, envir = .GlobalEnv)
  
}

STATS_1V <- NULL

for (i in 1:length(activity_list_1v)) {
  common <- c("STATS", "COUNT")
  stats1 <- activity_list_1v[[i]]
  
  if (is.null(STATS_1V)) {
    STATS_1V <- stats1
  } else {
    STATS_1V <- full_join(STATS_1V, stats1, by = common)
  }
}

####ORDER & SPLIT
STATS_1V <- STATS_1V[order(STATS_1V$STATS),]

STATS_1V_ctr <-STATS_1V[STATS_1V$STATS == "1V_CTR_SPM",]
STATS_1V_exp <-STATS_1V[STATS_1V$STATS == "1V_EXP_SPM",]

s1V_ctr_mean <- mean(STATS_1V_ctr$COUNT)
s1V_exp_mean <- mean(STATS_1V_exp$COUNT)

s1V_ctr_sd <- sd(STATS_1V_ctr$COUNT)
s1V_exp_sd <- sd(STATS_1V_exp$COUNT)

### T-TEST ###

t_test_1V <- t.test(STATS_1V_ctr$COUNT, STATS_1V_exp$COUNT, paired = TRUE)
Pval_1V <- t_test_1V$p.value

############ MAKE DATA FRAME #################

MEANS_1V <- c("1V_CTR_SPM", "1V_EXP_SPM", "1V_CTR_StDv", "1V_EXP_StDv", "1V_Pval")
COUNT_1V <- c(s1V_ctr_mean, s1V_exp_mean, s1V_ctr_sd, s1V_exp_sd, Pval_1V)

Stats_1VA <- data.frame(MEANS_1V, COUNT_1V)

################################################################
###############################################################
###############################################################
###############################################################

########          2v        #######################################

###################################################################

##################################################################

EFS_2V <- function(t2, b2){
  #TRACKMATE DATA
  csv_frame <- t2
  CELLS <- csv_frame[order(csv_frame$TRACK_ID, csv_frame$FRAME),]
  CELLS <- CELLS[order(CELLS$FRAME),]
  #BACKGROUND SPOTS
  spot_data <- b2
  SPOTS <- spot_data[-c(1),]
  SPOTS <- as.data.frame(apply(SPOTS, 2, as.numeric))
  SPOTS <- SPOTS %>% mutate(BG_Avg = rowMeans(SPOTS))
  SPOTS$FRAME <- (1:nrow(SPOTS))-1
  SPOTS <- SPOTS[, c(101,102)]
  #COMBINED CELL AND SPOT DATA
  DATA <- left_join(CELLS, SPOTS)
  DATA <- DATA[-c(6:9)]
  #BACKGROUND SUBTRACTION
  DATA$Avg_after_BGsub <- NA
  DATA <- DATA %>%
    mutate(Avg_after_BGsub = MEAN_INTENSITY_CH1 - BG_Avg)
  #CTR & EXP DATA SPLIT
  CTR <- subset(DATA, FRAME<25)
  EXP <- subset(DATA, FRAME >= 25)
  
  ##################### 2V PRE/CONTROL #################################################
  CTR <- CTR[order(CTR$TRACK_ID, CTR$FRAME),]
  
  #INTENSITY IN 80TH PERCENTILE
  quantile(CTR$Avg_after_BGsub, probs=.80)
  Q <- quantile(CTR$Avg_after_BGsub, probs=.80)
  
  
  
  CTR$Top20 <- NA
  CTR$Top20 <- ifelse(CTR$Avg_after_BGsub >= Q, "Top20val",
                      ifelse(CTR$Avg_after_BGsub, "No"))
  
  #FIRST POINT
  ctr_first <- CTR %>% group_by(TRACK_ID) %>%
    slice_head(n=1)
  
  ctr_first <- ctr_first[c("TRACK_ID", "Avg_after_BGsub")]
  ctr_first <- ctr_first %>% rename(point1 = Avg_after_BGsub)
  ctr <- left_join(CTR, ctr_first)
  
  #PERCENT CHANGE FROM FIRST POINT
  ctr <- ctr %>% group_by(TRACK_ID) %>% mutate(Change = Avg_after_BGsub-point1,
                                               PctChg = (Change/point1)*100)
  
  ctr_FindBase <- ctr
  ctr_FindBase <- ctr_FindBase[order(ctr_FindBase$TRACK_ID, ctr_FindBase$FRAME),]
  
  #https://stackoverflow.com/questions/61770151/counting-rows-until-a-condition-is-met-in-r-nas-before-the-condition-is-met
  ctr_FindBase <- ctr_FindBase %>%
    group_by(TRACK_ID) %>%
    group_by(temp = cumsum(PctChg >= 10), .add = TRUE) %>%
    mutate(out = row_number()-1) %>%
    group_by(TRACK_ID) %>%
    mutate(out = replace(out, row_number()<which.max(PctChg >= 10), NA))
  
  #https://stackoverflow.com/questions/24477748/r-count-na-by-group
  ctr_FindBase_count <- aggregate(out ~ TRACK_ID, ctr_FindBase, function(x) {sum(is.na(x))}, na.action=NULL)
  
  ctr_Lth <- ctr %>% group_by(TRACK_ID) %>% mutate(lgth = length(TRACK_ID))
  ctr_Lth <- ctr_Lth[c("TRACK_ID", "lgth")]
  ctr_Lth <- ctr_Lth %>% group_by(TRACK_ID) %>% slice_head(n=1)
  
  ctr_FindBase_Bn <- left_join(ctr_FindBase_count, ctr_Lth)
  
  for(i in 1:nrow(ctr_FindBase_Bn)){
    if(ctr_FindBase_Bn[i, "out"] == 0){
      ctr_FindBase_Bn[i, "out"] = ctr_FindBase_Bn[i, "lgth"]
    }
  }
  
  ctr_FindBase_Bn <- ctr_FindBase_Bn[-c(3)]
  ctr_FindBase_Bn <- ctr_FindBase_Bn %>% rename(meancount = out)
  
  #MEAN OF FIRST N COLUMNS - BASELINE AVG
  #https://stackoverflow.com/questions/64203431/subset-number-of-rows-per-group-based-on-a-value-in-another-dataframe
  ctr_base <- ctr_FindBase %>%
    left_join(ctr_FindBase_Bn) %>%
    group_by(TRACK_ID) %>%
    arrange(TRACK_ID, FRAME) %>%
    slice(1:meancount[1])
  
  ctr_base2 <- aggregate(Avg_after_BGsub ~ TRACK_ID, ctr_base, FUN=mean)
  ctr_base2 <- ctr_base2 %>% rename(base_avg = Avg_after_BGsub)
  
  ctr_base3 <- full_join(ctr_FindBase, ctr_base2)
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  
  #NORMALIZATION
  ctr_base3 <- ctr_base3 %>% mutate(baseline = Avg_after_BGsub/base_avg)
  ctr_base3 <- ctr_base3 %>% mutate(time = (((FRAME-24)*5)/60))
  
  
  #LAG PCT CHANGE
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Prev_int = lag(baseline),
                                                           Change_lag = baseline-Prev_int,
                                                           PctChg_lag = (Change_lag/Prev_int)*100)
  
  #LEAD INTENSITY GROUP
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Lead_int = lead(baseline))
  
  ctr_base3$Lead_int[is.na(ctr_base3$Lead_int)] <- 0
  
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Change_lead = (baseline - Lead_int))
  
  #DETERMINE SPIKE
  ctr_base3$Peak <- NA
  ctr_base3$minchange <- NA
  ctr_base3$PctChg_lag[is.na(ctr_base3$PctChg_lag)] <- 0
  ctr_base3$Change_lag[is.na(ctr_base3$Change_lag)] <- 0
  
  for (i in 1:nrow(ctr_base3)){
    if (ctr_base3[i, "Change_lag"] >0 && ctr_base3[i, "Change_lead"] >0){
      ctr_base3[i, "Peak"] = "Peak"
    } else{
      ctr_base3[i, "Peak"] = "No"
    }
  }
  
  counter = 0
  for (i in 1:nrow(ctr_base3)){
    if(ctr_base3[i, "Peak"] == "No"){
      ctr_base3[i, "minchange"] = counter
    } else{
      counter = counter+1
    }
  }
  
  ctr_spikedata <- ctr_base3 %>% mutate(leadingPCT = lead(PctChg_lag))
  ctr_spikedata1 <- aggregate(leadingPCT ~ minchange, ctr_spikedata, FUN=sum)
  
  ctr_spikedata2 <- full_join(ctr_spikedata1, ctr_base3)
  ctr_spikedata2 <- ctr_spikedata2[order(ctr_spikedata2$TRACK_ID, ctr_spikedata2$FRAME),]
  ctr_spikedata2<- ctr_spikedata2[c(3:24, 1, 2)]
  
  ctr_spikedata2 <- ctr_spikedata2 %>% mutate(PKpct = lag(leadingPCT))
  ctr_spikedata2$itsaspike <- NA
  ctr_spikedata2$PKpct[is.na(ctr_spikedata2$PKpct)] <- 0
  
  for (i in 1:nrow(ctr_spikedata2)){
    if(ctr_spikedata2[i, "Peak"] == "Peak"){
      if(ctr_spikedata2[i, "PKpct"] >= 10 || ctr_spikedata2[i, "PctChg_lag"] >= 10){
        ctr_spikedata2[i, "itsaspike"] = "SPIKE"
      } else{
        ctr_spikedata2[i, "itsaspike"] = "Not"
      }
    }else{
      ctr_spikedata2[i, "itsaspike"] = "Not"
    }
  }
  
  ###OUTLIER TEST####
  
  summary(ctr_spikedata2$baseline)
  ctr_mean <- mean(ctr_spikedata2$baseline)
  ctr_sd <- sd(ctr_spikedata2$baseline)
  
  ctr_spikedata2 <- ctr_spikedata2 %>% mutate(Zscore = (baseline-ctr_mean)/(ctr_sd))
  
  ctr_spikedata2 <- ctr_spikedata2[!(ctr_spikedata2$TRACK_ID %in% ctr_spikedata2$TRACK_ID[ctr_spikedata2$Zscore > 3]),]
  ctr_spikedata2 <- ctr_spikedata2[!(ctr_spikedata2$TRACK_ID %in% ctr_spikedata2$TRACK_ID[ctr_spikedata2$Zscore < -3]),]
  
  print(length(unique(ctr_base3$TRACK_ID)))
  print(length(unique(ctr_spikedata2$TRACK_ID)))
  
  #COUNT SPIKES/CELL (From Lag)
  ctr_numspikes <- sum(str_count(ctr_spikedata2$itsaspike, "SPIKE"))
  ctr_NUM <- count(ctr_spikedata2, itsaspike)
  
  ctr_celltotal <- length(unique(ctr_spikedata2$TRACK_ID))
  
  ctr_Activity <- (ctr_numspikes)/(ctr_celltotal*2)
  
  ##################### 2V POST/EXPERIMENTAL #################################################
  EXP <- EXP[order(EXP$TRACK_ID, EXP$FRAME),]
  
  #INTENSITY IN 80TH PERCENTILE
  quantile(EXP$Avg_after_BGsub, probs=.80)
  Q <- quantile(EXP$Avg_after_BGsub, probs=.80)
  
  
  
  EXP$Top20 <- NA
  EXP$Top20 <- ifelse(EXP$Avg_after_BGsub >= Q, "Top20val",
                      ifelse(EXP$Avg_after_BGsub, "No"))
  
  #FIRST POINT
  exp_first <- EXP %>% group_by(TRACK_ID) %>%
    slice_head(n=1)
  
  exp_first <- exp_first[c("TRACK_ID", "Avg_after_BGsub")]
  exp_first <- exp_first %>% rename(point1 = Avg_after_BGsub)
  EXP <- left_join(EXP, exp_first)
  
  #PERCENT CHANGE FROM FIRST POINT
  EXP <- EXP %>% group_by(TRACK_ID) %>% mutate(Change = Avg_after_BGsub-point1,
                                               PctChg = (Change/point1)*100)
  
  exp_FindBase <- EXP
  exp_FindBase <- exp_FindBase[order(exp_FindBase$TRACK_ID, exp_FindBase$FRAME),]
  
  #https://stackoverflow.com/questions/61770151/counting-rows-until-a-condition-is-met-in-r-nas-before-the-condition-is-met
  exp_FindBase <- exp_FindBase %>%
    group_by(TRACK_ID) %>%
    group_by(temp = cumsum(PctChg >= 10), .add = TRUE) %>%
    mutate(out = row_number()-1) %>%
    group_by(TRACK_ID) %>%
    mutate(out = replace(out, row_number()<which.max(PctChg >= 10), NA))
  
  #https://stackoverflow.com/questions/24477748/r-count-na-by-group
  exp_FindBase_count <- aggregate(out ~ TRACK_ID, exp_FindBase, function(x) {sum(is.na(x))}, na.action=NULL)
  
  exp_Lth <- EXP %>% group_by(TRACK_ID) %>% mutate(lgth = length(TRACK_ID))
  exp_Lth <- exp_Lth[c("TRACK_ID", "lgth")]
  exp_Lth <- exp_Lth %>% group_by(TRACK_ID) %>% slice_head(n=1)
  
  exp_FindBase_Bn <- left_join(exp_FindBase_count, exp_Lth)
  
  for(i in 1:nrow(exp_FindBase_Bn)){
    if(exp_FindBase_Bn[i, "out"] == 0){
      exp_FindBase_Bn[i, "out"] = exp_FindBase_Bn[i, "lgth"]
    }
  }
  
  exp_FindBase_Bn <- exp_FindBase_Bn[-c(3)]
  exp_FindBase_Bn <- exp_FindBase_Bn %>% rename(meancount = out)
  
  #MEAN OF FIRST N COLUMNS - BASELINE AVG
  #https://stackoverflow.com/questions/64203431/subset-number-of-rows-per-group-based-on-a-value-in-another-dataframe
  exp_base <- exp_FindBase %>%
    left_join(exp_FindBase_Bn) %>%
    group_by(TRACK_ID) %>%
    arrange(TRACK_ID, FRAME) %>%
    slice(1:meancount[1])
  
  exp_base2 <- aggregate(Avg_after_BGsub ~ TRACK_ID, exp_base, FUN=mean)
  exp_base2 <- exp_base2 %>% rename(base_avg = Avg_after_BGsub)
  
  exp_base3 <- full_join(exp_FindBase, exp_base2)
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  
  #NORMALIZATION
  exp_base3 <- exp_base3 %>% mutate(baseline = Avg_after_BGsub/base_avg)
  exp_base3 <- exp_base3 %>% mutate(time = (((FRAME-25)*5)/60))
  
  #LAG PCT CHANGE
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Prev_int = lag(baseline),
                                                           Change_lag = baseline-Prev_int,
                                                           PctChg_lag = (Change_lag/Prev_int)*100)
  
  #LEAD INTENSITY GROUP
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Lead_int = lead(baseline))
  
  exp_base3$Lead_int[is.na(exp_base3$Lead_int)] <- 0
  
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Change_lead = (baseline - Lead_int))
  
  #DETERMINE SPIKE
  exp_base3$Peak <- NA
  exp_base3$minchange <- NA
  exp_base3$PctChg_lag[is.na(exp_base3$PctChg_lag)] <- 0
  exp_base3$Change_lag[is.na(exp_base3$Change_lag)] <- 0
  
  for (i in 1:nrow(exp_base3)){
    if (exp_base3[i, "Change_lag"] >0 && exp_base3[i, "Change_lead"] >0){
      exp_base3[i, "Peak"] = "Peak"
    } else{
      exp_base3[i, "Peak"] = "No"
    }
  }
  
  counter = 0
  for (i in 1:nrow(exp_base3)){
    if(exp_base3[i, "Peak"] == "No"){
      exp_base3[i, "minchange"] = counter
    } else{
      counter = counter+1
    }
  }
  
  exp_spikedata <- exp_base3 %>% mutate(leadingPCT = lead(PctChg_lag))
  exp_spikedata1 <- aggregate(leadingPCT ~ minchange, exp_spikedata, FUN=sum)
  
  exp_spikedata2 <- full_join(exp_spikedata1, exp_base3)
  exp_spikedata2 <- exp_spikedata2[order(exp_spikedata2$TRACK_ID, exp_spikedata2$FRAME),]
  exp_spikedata2<- exp_spikedata2[c(3:24, 1, 2)]
  
  exp_spikedata2 <- exp_spikedata2 %>% mutate(PKpct = lag(leadingPCT))
  exp_spikedata2$itsaspike <- NA
  exp_spikedata2$PKpct[is.na(exp_spikedata2$PKpct)] <- 0
  
  for (i in 1:nrow(exp_spikedata2)){
    if(exp_spikedata2[i, "Peak"] == "Peak"){
      if(exp_spikedata2[i, "PKpct"] >= 10 || exp_spikedata2[i, "PctChg_lag"] >= 10){
        exp_spikedata2[i, "itsaspike"] = "SPIKE"
      } else{
        exp_spikedata2[i, "itsaspike"] = "Not"
      }
    }else{
      exp_spikedata2[i, "itsaspike"] = "Not"
    }
  }
  
  ###OUTLIER TEST####
  
  summary(exp_spikedata2$baseline)
  exp_mean <- mean(exp_spikedata2$baseline)
  exp_sd <- sd(exp_spikedata2$baseline)
  
  exp_spikedata2 <- exp_spikedata2 %>% mutate(Zscore = (baseline-exp_mean)/(exp_sd))
  
  exp_spikedata2 <- exp_spikedata2[!(exp_spikedata2$TRACK_ID %in% exp_spikedata2$TRACK_ID[exp_spikedata2$Zscore > 3]),]
  exp_spikedata2 <- exp_spikedata2[!(exp_spikedata2$TRACK_ID %in% exp_spikedata2$TRACK_ID[exp_spikedata2$Zscore < -3]),]
  
  print(length(unique(exp_base3$TRACK_ID)))
  print(length(unique(exp_spikedata2$TRACK_ID)))
  
  #COUNT SPIKES/CELL (From Lag)
  exp_numspikes <- sum(str_count(exp_spikedata2$itsaspike, "SPIKE"))
  exp_NUM <- count(exp_spikedata2, itsaspike)
  
  exp_celltotal <- length(unique(exp_spikedata2$TRACK_ID))
  
  exp_Activity <- (exp_numspikes)/(exp_celltotal*2)
  
 
  ######## INDIVIDUAL SPIKE DATA ###################
  
  Cspike <- subset(ctr_spikedata2, itsaspike == "SPIKE")
  Cspike2 <- Cspike[c('TRACK_ID', 'itsaspike')]
  Cspike3 <- Cspike2 %>% group_by(TRACK_ID) %>%
    mutate(counts = str_count(itsaspike, "SPIKE"))
  Cspike4 <- aggregate(counts ~ TRACK_ID, Cspike3, FUN=sum)
  
  Espike <- subset(exp_spikedata2, itsaspike == "SPIKE")
  Espike2 <- Espike[c('TRACK_ID', 'itsaspike')]
  Espike3 <- Espike2 %>% group_by(TRACK_ID) %>%
    mutate(counts = str_count(itsaspike, "SPIKE"))
  Espike4 <- aggregate(counts ~ TRACK_ID, Espike3, FUN=sum)
  
  Cspikes <<- Cspike4
  Espikes <<- Espike4
  
  ############ MAKE DATA FRAME #################
  STATS <- c("2V_CTR_SPM", "2V_EXP_SPM", "2V_CTR_SPIKE_TOTAL", "2V_CTR_CELL_TOTAL", "2V_EXP_SPIKE_TOTAL", "2V_EXP_CELL_TOTAL")
  COUNT <- c(ctr_Activity, exp_Activity, ctr_numspikes, ctr_celltotal, exp_numspikes, exp_celltotal)
  
  Stats_2V_funcdata <<- data.frame(STATS, COUNT)
}


activity_list_2v <- list()
ctr_spikelist_2v <- list()
exp_spikelist_2v <- list()

for (i in 1:length(tracks_2v)) {
  results <- EFS_2V(tracks_2v[[i]], background_2v[[i]])
  activity_list_2v[[i]] <- assign(paste0("v2_ACTIVITY_", i, sep = ''), Stats_2V_funcdata, envir = .GlobalEnv)
  ctr_spikelist_2v[[i]] <- assign(paste0("v2_CTR_Spikes_", i, sep = ''), Cspikes, envir = .GlobalEnv)
  exp_spikelist_2v[[i]] <- assign(paste0("v2_EXP_Spikes_", i, sep = ''), Espikes, envir = .GlobalEnv)
  
}

STATS_2V <- NULL

for (i in 1:length(activity_list_2v)) {
  common <- c("STATS", "COUNT")
  stats2 <- activity_list_2v[[i]]
  
  if (is.null(STATS_2V)) {
    STATS_2V <- stats2
  } else {
    STATS_2V <- full_join(STATS_2V, stats2, by = common)
  }
}

####ORDER & SPLIT
STATS_2V <- STATS_2V[order(STATS_2V$STATS),]

STATS_2V_ctr <-STATS_2V[STATS_2V$STATS == "2V_CTR_SPM",]
STATS_2V_exp <-STATS_2V[STATS_2V$STATS == "2V_EXP_SPM",]

s2V_ctr_mean <- mean(STATS_2V_ctr$COUNT)
s2V_exp_mean <- mean(STATS_2V_exp$COUNT)

s2V_ctr_sd <- sd(STATS_2V_ctr$COUNT)
s2V_exp_sd <- sd(STATS_2V_exp$COUNT)

### T-TEST ###

t_test_2V <- t.test(STATS_2V_ctr$COUNT, STATS_2V_exp$COUNT, paired = TRUE)
Pval_2V <- t_test_2V$p.value

############ MAKE DATA FRAME #################

MEANS_2V <- c("2V_CTR_SPM", "2V_EXP_SPM", "2V_CTR_StDv", "2V_EXP_StDv", "2V_Pval")
COUNT_2V <- c(s2V_ctr_mean, s2V_exp_mean, s2V_ctr_sd, s2V_exp_sd, Pval_2V)

Stats_2VA <- data.frame(MEANS_2V, COUNT_2V)



################################################################
###############################################################
###############################################################
###############################################################

########          3v        #######################################

###################################################################

##################################################################


EFS_3V <- function(t3, b3){
  #TRACKMATE DATA
  csv_frame <- t3
  CELLS <- csv_frame[order(csv_frame$TRACK_ID, csv_frame$FRAME),]
  CELLS <- CELLS[order(CELLS$FRAME),]
  #BACKGROUND SPOTS
  spot_data <- b3
  SPOTS <- spot_data[-c(1),]
  SPOTS <- as.data.frame(apply(SPOTS, 2, as.numeric))
  SPOTS <- SPOTS %>% mutate(BG_Avg = rowMeans(SPOTS))
  SPOTS$FRAME <- (1:nrow(SPOTS))-1
  SPOTS <- SPOTS[, c(101,102)]
  #COMBINED CELL AND SPOT DATA
  DATA <- left_join(CELLS, SPOTS)
  DATA <- DATA[-c(6:9)]
  #BACKGROUND SUBTRACTION
  DATA$Avg_after_BGsub <- NA
  DATA <- DATA %>%
    mutate(Avg_after_BGsub = MEAN_INTENSITY_CH1 - BG_Avg)
  #CTR & EXP DATA SPLIT
  CTR <- subset(DATA, FRAME<25)
  EXP <- subset(DATA, FRAME >= 25)
  
  ##################### 3V PRE/CONTROL #################################################
  CTR <- CTR[order(CTR$TRACK_ID, CTR$FRAME),]
  
  #INTENSITY IN 80TH PERCENTILE
  quantile(CTR$Avg_after_BGsub, probs=.80)
  Q <- quantile(CTR$Avg_after_BGsub, probs=.80)
  
  
  
  CTR$Top20 <- NA
  CTR$Top20 <- ifelse(CTR$Avg_after_BGsub >= Q, "Top20val",
                      ifelse(CTR$Avg_after_BGsub, "No"))
  
  #FIRST POINT
  ctr_first <- CTR %>% group_by(TRACK_ID) %>%
    slice_head(n=1)
  
  ctr_first <- ctr_first[c("TRACK_ID", "Avg_after_BGsub")]
  ctr_first <- ctr_first %>% rename(point1 = Avg_after_BGsub)
  ctr <- left_join(CTR, ctr_first)
  
  #PERCENT CHANGE FROM FIRST POINT
  ctr <- ctr %>% group_by(TRACK_ID) %>% mutate(Change = Avg_after_BGsub-point1,
                                               PctChg = (Change/point1)*100)
  
  ctr_FindBase <- ctr
  ctr_FindBase <- ctr_FindBase[order(ctr_FindBase$TRACK_ID, ctr_FindBase$FRAME),]
  
  #https://stackoverflow.com/questions/61770151/counting-rows-until-a-condition-is-met-in-r-nas-before-the-condition-is-met
  ctr_FindBase <- ctr_FindBase %>%
    group_by(TRACK_ID) %>%
    group_by(temp = cumsum(PctChg >= 10), .add = TRUE) %>%
    mutate(out = row_number()-1) %>%
    group_by(TRACK_ID) %>%
    mutate(out = replace(out, row_number()<which.max(PctChg >= 10), NA))
  
  #https://stackoverflow.com/questions/24477748/r-count-na-by-group
  ctr_FindBase_count <- aggregate(out ~ TRACK_ID, ctr_FindBase, function(x) {sum(is.na(x))}, na.action=NULL)
  
  ctr_Lth <- ctr %>% group_by(TRACK_ID) %>% mutate(lgth = length(TRACK_ID))
  ctr_Lth <- ctr_Lth[c("TRACK_ID", "lgth")]
  ctr_Lth <- ctr_Lth %>% group_by(TRACK_ID) %>% slice_head(n=1)
  
  ctr_FindBase_Bn <- left_join(ctr_FindBase_count, ctr_Lth)
  
  for(i in 1:nrow(ctr_FindBase_Bn)){
    if(ctr_FindBase_Bn[i, "out"] == 0){
      ctr_FindBase_Bn[i, "out"] = ctr_FindBase_Bn[i, "lgth"]
    }
  }
  
  ctr_FindBase_Bn <- ctr_FindBase_Bn[-c(3)]
  ctr_FindBase_Bn <- ctr_FindBase_Bn %>% rename(meancount = out)
  
  #MEAN OF FIRST N COLUMNS - BASELINE AVG
  #https://stackoverflow.com/questions/64203431/subset-number-of-rows-per-group-based-on-a-value-in-another-dataframe
  ctr_base <- ctr_FindBase %>%
    left_join(ctr_FindBase_Bn) %>%
    group_by(TRACK_ID) %>%
    arrange(TRACK_ID, FRAME) %>%
    slice(1:meancount[1])
  
  ctr_base2 <- aggregate(Avg_after_BGsub ~ TRACK_ID, ctr_base, FUN=mean)
  ctr_base2 <- ctr_base2 %>% rename(base_avg = Avg_after_BGsub)
  
  ctr_base3 <- full_join(ctr_FindBase, ctr_base2)
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  
  #NORMALIZATION
  ctr_base3 <- ctr_base3 %>% mutate(baseline = Avg_after_BGsub/base_avg)
  ctr_base3 <- ctr_base3 %>% mutate(time = (((FRAME-24)*5)/60))
  
  #LAG PCT CHANGE
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Prev_int = lag(baseline),
                                                           Change_lag = baseline-Prev_int,
                                                           PctChg_lag = (Change_lag/Prev_int)*100)
  
  #LEAD INTENSITY GROUP
  ctr_base3 <- ctr_base3[order(ctr_base3$TRACK_ID, ctr_base3$FRAME),]
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Lead_int = lead(baseline))
  
  ctr_base3$Lead_int[is.na(ctr_base3$Lead_int)] <- 0
  
  ctr_base3 <- ctr_base3 %>% group_by(TRACK_ID) %>% mutate(Change_lead = (baseline - Lead_int))
  
  #DETERMINE SPIKE
  ctr_base3$Peak <- NA
  ctr_base3$minchange <- NA
  ctr_base3$PctChg_lag[is.na(ctr_base3$PctChg_lag)] <- 0
  ctr_base3$Change_lag[is.na(ctr_base3$Change_lag)] <- 0
  
  for (i in 1:nrow(ctr_base3)){
    if (ctr_base3[i, "Change_lag"] >0 && ctr_base3[i, "Change_lead"] >0){
      ctr_base3[i, "Peak"] = "Peak"
    } else{
      ctr_base3[i, "Peak"] = "No"
    }
  }
  
  counter = 0
  for (i in 1:nrow(ctr_base3)){
    if(ctr_base3[i, "Peak"] == "No"){
      ctr_base3[i, "minchange"] = counter
    } else{
      counter = counter+1
    }
  }
  
  ctr_spikedata <- ctr_base3 %>% mutate(leadingPCT = lead(PctChg_lag))
  ctr_spikedata1 <- aggregate(leadingPCT ~ minchange, ctr_spikedata, FUN=sum)
  
  ctr_spikedata2 <- full_join(ctr_spikedata1, ctr_base3)
  ctr_spikedata2 <- ctr_spikedata2[order(ctr_spikedata2$TRACK_ID, ctr_spikedata2$FRAME),]
  ctr_spikedata2<- ctr_spikedata2[c(3:24, 1, 2)]
  
  ctr_spikedata2 <- ctr_spikedata2 %>% mutate(PKpct = lag(leadingPCT))
  ctr_spikedata2$itsaspike <- NA
  ctr_spikedata2$PKpct[is.na(ctr_spikedata2$PKpct)] <- 0
  
  for (i in 1:nrow(ctr_spikedata2)){
    if(ctr_spikedata2[i, "Peak"] == "Peak"){
      if(ctr_spikedata2[i, "PKpct"] >= 10 || ctr_spikedata2[i, "PctChg_lag"] >= 10){
        ctr_spikedata2[i, "itsaspike"] = "SPIKE"
      } else{
        ctr_spikedata2[i, "itsaspike"] = "Not"
      }
    }else{
      ctr_spikedata2[i, "itsaspike"] = "Not"
    }
  }
  
  ###OUTLIER TEST####
  
  summary(ctr_spikedata2$baseline)
  ctr_mean <- mean(ctr_spikedata2$baseline)
  ctr_sd <- sd(ctr_spikedata2$baseline)
  
  ctr_spikedata2 <- ctr_spikedata2 %>% mutate(Zscore = (baseline-ctr_mean)/(ctr_sd))
  
  ctr_spikedata2 <- ctr_spikedata2[!(ctr_spikedata2$TRACK_ID %in% ctr_spikedata2$TRACK_ID[ctr_spikedata2$Zscore > 3]),]
  ctr_spikedata2 <- ctr_spikedata2[!(ctr_spikedata2$TRACK_ID %in% ctr_spikedata2$TRACK_ID[ctr_spikedata2$Zscore < -3]),]
  
  print(length(unique(ctr_base3$TRACK_ID)))
  print(length(unique(ctr_spikedata2$TRACK_ID)))
  
  #COUNT SPIKES/CELL (From Lag)
  ctr_numspikes <- sum(str_count(ctr_spikedata2$itsaspike, "SPIKE"))
  ctr_NUM <- count(ctr_spikedata2, itsaspike)
  
  ctr_celltotal <- length(unique(ctr_spikedata2$TRACK_ID))
  
  ctr_Activity <- (ctr_numspikes)/(ctr_celltotal*2)
  
  ##################### 3V POST/EXPERIMENTAL #################################################
  EXP <- EXP[order(EXP$TRACK_ID, EXP$FRAME),]
  
  #INTENSITY IN 80TH PERCENTILE
  quantile(EXP$Avg_after_BGsub, probs=.80)
  Q <- quantile(EXP$Avg_after_BGsub, probs=.80)
  
  
  
  EXP$Top20 <- NA
  EXP$Top20 <- ifelse(EXP$Avg_after_BGsub >= Q, "Top20val",
                      ifelse(EXP$Avg_after_BGsub, "No"))
  
  #FIRST POINT
  exp_first <- EXP %>% group_by(TRACK_ID) %>%
    slice_head(n=1)
  
  exp_first <- exp_first[c("TRACK_ID", "Avg_after_BGsub")]
  exp_first <- exp_first %>% rename(point1 = Avg_after_BGsub)
  EXP <- left_join(EXP, exp_first)
  
  #PERCENT CHANGE FROM FIRST POINT
  EXP <- EXP %>% group_by(TRACK_ID) %>% mutate(Change = Avg_after_BGsub-point1,
                                               PctChg = (Change/point1)*100)
  
  exp_FindBase <- EXP
  exp_FindBase <- exp_FindBase[order(exp_FindBase$TRACK_ID, exp_FindBase$FRAME),]
  
  #https://stackoverflow.com/questions/61770151/counting-rows-until-a-condition-is-met-in-r-nas-before-the-condition-is-met
  exp_FindBase <- exp_FindBase %>%
    group_by(TRACK_ID) %>%
    group_by(temp = cumsum(PctChg >= 10), .add = TRUE) %>%
    mutate(out = row_number()-1) %>%
    group_by(TRACK_ID) %>%
    mutate(out = replace(out, row_number()<which.max(PctChg >= 10), NA))
  
  #https://stackoverflow.com/questions/24477748/r-count-na-by-group
  exp_FindBase_count <- aggregate(out ~ TRACK_ID, exp_FindBase, function(x) {sum(is.na(x))}, na.action=NULL)
  
  exp_Lth <- EXP %>% group_by(TRACK_ID) %>% mutate(lgth = length(TRACK_ID))
  exp_Lth <- exp_Lth[c("TRACK_ID", "lgth")]
  exp_Lth <- exp_Lth %>% group_by(TRACK_ID) %>% slice_head(n=1)
  
  exp_FindBase_Bn <- left_join(exp_FindBase_count, exp_Lth)
  
  for(i in 1:nrow(exp_FindBase_Bn)){
    if(exp_FindBase_Bn[i, "out"] == 0){
      exp_FindBase_Bn[i, "out"] = exp_FindBase_Bn[i, "lgth"]
    }
  }
  
  exp_FindBase_Bn <- exp_FindBase_Bn[-c(3)]
  exp_FindBase_Bn <- exp_FindBase_Bn %>% rename(meancount = out)
  
  #MEAN OF FIRST N COLUMNS - BASELINE AVG
  #https://stackoverflow.com/questions/64203431/subset-number-of-rows-per-group-based-on-a-value-in-another-dataframe
  exp_base <- exp_FindBase %>%
    left_join(exp_FindBase_Bn) %>%
    group_by(TRACK_ID) %>%
    arrange(TRACK_ID, FRAME) %>%
    slice(1:meancount[1])
  
  exp_base2 <- aggregate(Avg_after_BGsub ~ TRACK_ID, exp_base, FUN=mean)
  exp_base2 <- exp_base2 %>% rename(base_avg = Avg_after_BGsub)
  
  exp_base3 <- full_join(exp_FindBase, exp_base2)
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  
  #NORMALIZATION
  exp_base3 <- exp_base3 %>% mutate(baseline = Avg_after_BGsub/base_avg)
  exp_base3 <- exp_base3 %>% mutate(time = (((FRAME-25)*5)/60))
  
  #LAG PCT CHANGE
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Prev_int = lag(baseline),
                                                           Change_lag = baseline-Prev_int,
                                                           PctChg_lag = (Change_lag/Prev_int)*100)
  
  #LEAD INTENSITY GROUP
  exp_base3 <- exp_base3[order(exp_base3$TRACK_ID, exp_base3$FRAME),]
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Lead_int = lead(baseline))
  
  exp_base3$Lead_int[is.na(exp_base3$Lead_int)] <- 0
  
  exp_base3 <- exp_base3 %>% group_by(TRACK_ID) %>% mutate(Change_lead = (baseline - Lead_int))
  
  #DETERMINE SPIKE
  exp_base3$Peak <- NA
  exp_base3$minchange <- NA
  exp_base3$PctChg_lag[is.na(exp_base3$PctChg_lag)] <- 0
  exp_base3$Change_lag[is.na(exp_base3$Change_lag)] <- 0
  
  for (i in 1:nrow(exp_base3)){
    if (exp_base3[i, "Change_lag"] >0 && exp_base3[i, "Change_lead"] >0){
      exp_base3[i, "Peak"] = "Peak"
    } else{
      exp_base3[i, "Peak"] = "No"
    }
  }
  
  counter = 0
  for (i in 1:nrow(exp_base3)){
    if(exp_base3[i, "Peak"] == "No"){
      exp_base3[i, "minchange"] = counter
    } else{
      counter = counter+1
    }
  }
  
  exp_spikedata <- exp_base3 %>% mutate(leadingPCT = lead(PctChg_lag))
  exp_spikedata1 <- aggregate(leadingPCT ~ minchange, exp_spikedata, FUN=sum)
  
  exp_spikedata2 <- full_join(exp_spikedata1, exp_base3)
  exp_spikedata2 <- exp_spikedata2[order(exp_spikedata2$TRACK_ID, exp_spikedata2$FRAME),]
  exp_spikedata2<- exp_spikedata2[c(3:24, 1, 2)]
  
  exp_spikedata2 <- exp_spikedata2 %>% mutate(PKpct = lag(leadingPCT))
  exp_spikedata2$itsaspike <- NA
  exp_spikedata2$PKpct[is.na(exp_spikedata2$PKpct)] <- 0
  
  for (i in 1:nrow(exp_spikedata2)){
    if(exp_spikedata2[i, "Peak"] == "Peak"){
      if(exp_spikedata2[i, "PKpct"] >= 10 || exp_spikedata2[i, "PctChg_lag"] >= 10){
        exp_spikedata2[i, "itsaspike"] = "SPIKE"
      } else{
        exp_spikedata2[i, "itsaspike"] = "Not"
      }
    }else{
      exp_spikedata2[i, "itsaspike"] = "Not"
    }
  }
  
  ###OUTLIER TEST####
  
  summary(exp_spikedata2$baseline)
  exp_mean <- mean(exp_spikedata2$baseline)
  exp_sd <- sd(exp_spikedata2$baseline)
  
  exp_spikedata2 <- exp_spikedata2 %>% mutate(Zscore = (baseline-exp_mean)/(exp_sd))
  
  exp_spikedata2 <- exp_spikedata2[!(exp_spikedata2$TRACK_ID %in% exp_spikedata2$TRACK_ID[exp_spikedata2$Zscore > 3]),]
  exp_spikedata2 <- exp_spikedata2[!(exp_spikedata2$TRACK_ID %in% exp_spikedata2$TRACK_ID[exp_spikedata2$Zscore < -3]),]
  
  print(length(unique(exp_base3$TRACK_ID)))
  print(length(unique(exp_spikedata2$TRACK_ID)))
  
  #COUNT SPIKES/CELL (From Lag)
  exp_numspikes <- sum(str_count(exp_spikedata2$itsaspike, "SPIKE"))
  exp_NUM <- count(exp_spikedata2, itsaspike)
  
  exp_celltotal <- length(unique(exp_spikedata2$TRACK_ID))
  
  exp_Activity <- (exp_numspikes)/(exp_celltotal*2)
  
  
  ######## INDIVIDUAL SPIKE DATA ###################
  
  Cspike <- subset(ctr_spikedata2, itsaspike == "SPIKE")
  Cspike2 <- Cspike[c('TRACK_ID', 'itsaspike')]
  Cspike3 <- Cspike2 %>% group_by(TRACK_ID) %>%
    mutate(counts = str_count(itsaspike, "SPIKE"))
  Cspike4 <- aggregate(counts ~ TRACK_ID, Cspike3, FUN=sum)
  
  Espike <- subset(exp_spikedata2, itsaspike == "SPIKE")
  Espike2 <- Espike[c('TRACK_ID', 'itsaspike')]
  Espike3 <- Espike2 %>% group_by(TRACK_ID) %>%
    mutate(counts = str_count(itsaspike, "SPIKE"))
  Espike4 <- aggregate(counts ~ TRACK_ID, Espike3, FUN=sum)
  
  Cspikes <<- Cspike4
  Espikes <<- Espike4
  
  ############ MAKE DATA FRAME #################
  STATS <- c("3V_CTR_SPM", "3V_EXP_SPM", "3V_CTR_SPIKE_TOTAL", "3V_CTR_CELL_TOTAL", "3V_EXP_SPIKE_TOTAL", "3V_EXP_CELL_TOTAL")
  COUNT <- c(ctr_Activity, exp_Activity, ctr_numspikes, ctr_celltotal, exp_numspikes, exp_celltotal)
  
  Stats_3V_funcdata <<- data.frame(STATS, COUNT)
}


activity_list_3v <- list()
ctr_spikelist_3v <- list()
exp_spikelist_3v <- list()

for (i in 1:length(tracks_3v)) {
  results <- EFS_3V(tracks_3v[[i]], background_3v[[i]])
  activity_list_3v[[i]] <- assign(paste0("v3_ACTIVITY_", i, sep = ''), Stats_3V_funcdata, envir = .GlobalEnv)
  ctr_spikelist_3v[[i]] <- assign(paste0("v3_CTR_Spikes_", i, sep = ''), Cspikes, envir = .GlobalEnv)
  exp_spikelist_3v[[i]] <- assign(paste0("v3_EXP_Spikes_", i, sep = ''), Espikes, envir = .GlobalEnv)
  
}

STATS_3V <- NULL

for (i in 1:length(activity_list_3v)) {
  common <- c("STATS", "COUNT")
  stats3 <- activity_list_3v[[i]]
  
  if (is.null(STATS_3V)) {
    STATS_3V <- stats3
  } else {
    STATS_3V <- full_join(STATS_3V, stats3, by = common)
  }
}

####ORDER & SPLIT
STATS_3V <- STATS_3V[order(STATS_3V$STATS),]

STATS_3V_ctr <-STATS_3V[STATS_3V$STATS == "3V_CTR_SPM",]
STATS_3V_exp <-STATS_3V[STATS_3V$STATS == "3V_EXP_SPM",]

s3V_ctr_mean <- mean(STATS_3V_ctr$COUNT)
s3V_exp_mean <- mean(STATS_3V_exp$COUNT)

s3V_ctr_sd <- sd(STATS_3V_ctr$COUNT)
s3V_exp_sd <- sd(STATS_3V_exp$COUNT)

### T-TEST ###

t_test_3V <- t.test(STATS_3V_ctr$COUNT, STATS_3V_exp$COUNT, paired = TRUE)
Pval_3V <- t_test_3V$p.value

############ MAKE DATA FRAME #################

MEANS_3V <- c("3V_CTR_SPM", "3V_EXP_SPM", "3V_CTR_StDv", "3V_EXP_StDv", "3V_Pval")
COUNT_3V <- c(s3V_ctr_mean, s3V_exp_mean, s3V_ctr_sd, s3V_exp_sd, Pval_3V)

Stats_3VA <- data.frame(MEANS_3V, COUNT_3V)




################# 1,2,3 COMBINED ######################

###############################
##############################


stats_1v <- Stats_1VA
stats_1v <- stats_1v %>% rename("MEANS" = "MEANS_1V")
stats_1v <- stats_1v %>% rename("COUNT" = "COUNT_1V")

stats_2v <- Stats_2VA
stats_2v <- stats_2v %>% rename("MEANS" = "MEANS_2V")
stats_2v <- stats_2v %>% rename("COUNT" = "COUNT_2V")

stats_3v <- Stats_3VA
stats_3v <- stats_3v %>% rename("MEANS" = "MEANS_3V")
stats_3v <- stats_3v %>% rename("COUNT" = "COUNT_3V")

STATS_all <- full_join(stats_1v, stats_2v)
STATS_all <- full_join(STATS_all, stats_3v)

### SUBSET and STDV ##
stats_sd <- STATS_all[grep("StDv", STATS_all$MEANS),]
stats_pval <- STATS_all[grep("Pval", STATS_all$MEANS),]
STATS_base <- STATS_all[grep("SPM", STATS_all$MEANS),]

stats_sd <- stats_sd %>% rename(SD = COUNT)

STATS_base[STATS_base == "1V_CTR_SPM"] <- "1V_CTR"
STATS_base[STATS_base == "2V_CTR_SPM"] <- "2V_CTR"
STATS_base[STATS_base == "3V_CTR_SPM"] <- "3V_CTR"
STATS_base[STATS_base == "1V_EXP_SPM"] <- "1V_EXP"
STATS_base[STATS_base == "2V_EXP_SPM"] <- "2V_EXP"
STATS_base[STATS_base == "3V_EXP_SPM"] <- "3V_EXP"

stats_sd[stats_sd == "1V_CTR_StDv"] <- "1V_CTR"
stats_sd[stats_sd == "2V_CTR_StDv"] <- "2V_CTR"
stats_sd[stats_sd == "3V_CTR_StDv"] <- "3V_CTR"
stats_sd[stats_sd == "1V_EXP_StDv"] <- "1V_EXP"
stats_sd[stats_sd == "2V_EXP_StDv"] <- "2V_EXP"
stats_sd[stats_sd == "3V_EXP_StDv"] <- "3V_EXP"

STATS_base <- full_join(STATS_base, stats_sd)


### GRAPH & CSV OUTPUT ###
STATS_base <- STATS_base %>% mutate(Run = c(1,1,2,2,3,3))

STATS_base <- STATS_base %>% rename(RUN = MEANS,
                                    EFS = Run)
STATS_base[STATS_base == "1V_CTR"] <- "c_PRE"
STATS_base[STATS_base == "2V_CTR"] <- "c_PRE"
STATS_base[STATS_base == "3V_CTR"] <- "c_PRE"
STATS_base[STATS_base == "1V_EXP"] <- "POST"
STATS_base[STATS_base == "2V_EXP"] <- "POST"
STATS_base[STATS_base == "3V_EXP"] <- "POST"

sn <- c(10,10,10,10,9,9)

STATS_base2 <- STATS_base %>% mutate(SamNum = sn)
STATS_base2 <- STATS_base2 %>% mutate(SEM = SD/(sqrt(SamNum)))

ggplot(STATS_base2, aes(x=EFS, y=COUNT, fill=RUN))+
  geom_bar(stat='identity', position=position_dodge(), color="black", width=0.5)+
  geom_errorbar(aes(ymin=COUNT-SEM, ymax=COUNT+SEM), width =.2,
                position = position_dodge(.5), size=1)+
  scale_fill_manual(values = c("azure3", "darkslategray"),
                    labels = c("Pre", "Post"))+
  scale_x_discrete(limits = c("1 V/cm", "2 V/cm", "3 V/cm"))+
  scale_y_continuous(limits=c(0,1.2), breaks = seq(0, 1.2, 0.2), expand = expansion(mult = c(0, 0)))+
  labs(x = "",
       y = "Spiking Activity (Spikes/min)",
       title = "Spiking Activity",
       subtitle = "",
       fill = "")+
  theme_cowplot(12)+
  theme(plot.title = element_text(size= 32, hjust=0.5),
        plot.subtitle = element_text(size=18, hjust=0.5),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=26, vjust=2),
        axis.text.x = element_text(size=32, vjust=-2),
        axis.text.y = element_text(size=26),
        axis.line = element_line(size=1.25),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=18),
        legend.position = "top",
        legend.justification = "center")

write.csv(STATS_base2,'stats_base.csv', row.names = FALSE)
write.csv(Stats_1VA,'Stats_1Vcm.csv', row.names = FALSE)
write.csv(Stats_2VA,'Stats_2Vcm.csv', row.names = FALSE)
write.csv(Stats_3VA,'Stats_3Vcm.csv', row.names = FALSE)



# ACTIVE POPULATION 

lengths_1v_ctr <- data.frame()
lengths_1v_exp <- data.frame()

for (i in 1:length(ctr_spikelist_1v)) {
  ID_length <- length(ctr_spikelist_1v[[i]][["TRACK_ID"]])
  lengths_1v_ctr <- rbind(lengths_1v_ctr, data.frame(Length = ID_length))
}

for (i in 1:length(exp_spikelist_1v)) {
  ID_length <- length(exp_spikelist_1v[[i]][["TRACK_ID"]])
  lengths_1v_exp <- rbind(lengths_1v_exp, data.frame(Length = ID_length))
}

lengths_2v_ctr <- data.frame()
lengths_2v_exp <- data.frame()

for (i in 1:length(ctr_spikelist_2v)) {
  ID_length <- length(ctr_spikelist_2v[[i]][["TRACK_ID"]])
  lengths_2v_ctr <- rbind(lengths_2v_ctr, data.frame(Length = ID_length))
}

for (i in 1:length(exp_spikelist_2v)) {
  ID_length <- length(exp_spikelist_2v[[i]][["TRACK_ID"]])
  lengths_2v_exp <- rbind(lengths_2v_exp, data.frame(Length = ID_length))
}

lengths_3v_ctr <- data.frame()
lengths_3v_exp <- data.frame()

for (i in 1:length(ctr_spikelist_3v)) {
  ID_length <- length(ctr_spikelist_3v[[i]][["TRACK_ID"]])
  lengths_3v_ctr <- rbind(lengths_3v_ctr, data.frame(Length = ID_length))
}

for (i in 1:length(exp_spikelist_3v)) {
  ID_length <- length(exp_spikelist_3v[[i]][["TRACK_ID"]])
  lengths_3v_exp <- rbind(lengths_3v_exp, data.frame(Length = ID_length))
}

lengths_1v_ctr$Trial <- 1:10
lengths_1v_exp$Trial <- 1:10
lengths_2v_ctr$Trial <- 11:20
lengths_2v_exp$Trial <- 11:20
lengths_3v_ctr$Trial <- 21:29
lengths_3v_exp$Trial <- 21:29

active_1v <- merge(lengths_1v_ctr, lengths_1v_exp, by = 'Trial', all=TRUE)
active_1v <- active_1v %>% rename(CTR = Length.x, EXP = Length.y)

active_2v <- merge(lengths_2v_ctr, lengths_2v_exp, by = 'Trial', all=TRUE)
active_2v <- active_2v %>% rename(CTR = Length.x, EXP = Length.y)

active_3v <- merge(lengths_3v_ctr, lengths_3v_exp, by = 'Trial', all=TRUE)
active_3v <- active_3v %>% rename(CTR = Length.x, EXP = Length.y)

active_all <- full_join(active_1v, active_2v)
active_all <- full_join(active_all, active_3v)



active_all <- active_all %>% mutate(PctDiff = ((EXP-CTR)/CTR)*100)
active_all <- active_all %>% mutate(PctRd = round(PctDiff))

change_mean <- mean(active_all$PctDiff)


conditions <- active_all[c(1:3)]
condC <- conditions[-c(3)]
condC <- condC %>% rename(AC = CTR)
condC <- condC %>% mutate(COND = "CTR")

condE <- conditions[-c(2)]
condE <- condE %>% rename(AC = EXP)
condE <- condE %>% mutate(COND = "EXP")

conditions <- rbind(condC, condE)


ggplot(conditions, aes(x=COND, y=AC))+
  geom_boxplot(aes(fill=COND), size=1)+
  stat_boxplot(geom = 'errorbar', width = 0.3, size=1.5)+
  scale_fill_manual(values = c("dodgerblue2", "chartreuse3"))+
  geom_point(color="black", fill="deeppink2", size=3, shape=21)+
  geom_line(aes(group=Trial),linewidth = 0.6)+
  scale_x_discrete(labels = c("Control", "Post-EFS"))+
  labs(x = "Trial Condition",
       y = "Active Cells (>= 1 Spike)",
       title = "Effects of EFS on Active Cell Population",
       subtitle = "",
       fill = "")+
  theme_cowplot(12)+
  theme(plot.title = element_text(size= 26, hjust=0.5, face="bold"),
        plot.subtitle = element_text(size=18, hjust=0.5),
        axis.title.x = element_text(size=26, face="bold"),
        axis.title.y = element_text(size=26, face="bold", vjust=2),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=26),
        axis.line = element_line(size=1.25),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=18),
        legend.position = "none")+
  scale_y_continuous(breaks=seq(0,400,50))

