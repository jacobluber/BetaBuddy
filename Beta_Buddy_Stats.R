library(ggplot2)
library(dplyr)
library(cowplot)
library(rstudioapi)
library(tidyr)


#### EXPERIMENTAL DATA ######
#TRACKMATE DATA
csv_frame <- read.csv('experimentjup.csv',header=TRUE)
CELLS <- csv_frame[order(csv_frame$TRACK_ID, csv_frame$FRAME),]
CELLS <- CELLS[order(CELLS$FRAME),]

#BACKGROUND SPOTS
spot_data <- read.csv('backgroundsub.csv',header=TRUE)
SPOTS <- spot_data[-c(1),]
SPOTS <- as.data.frame(apply(SPOTS, 2, as.numeric))
SPOTS <- SPOTS %>% mutate(BG_Avg = rowMeans(SPOTS))
SPOTS$FRAME <- (1:nrow(SPOTS))-1
SPOTS <- SPOTS[, c(101,102)]

#COMBINED CELL AND SPOT DATA
DATA <- left_join(CELLS, SPOTS)

#BACKGROUND SUBTRACTION
DATA$Avg_after_BGsub <- NA
DATA <- DATA %>%
  mutate(Avg_after_BGsub = MEAN_INTENSITY_CH1 - BG_Avg)

#AVERAGE FIRST 5
avg_5 <- DATA %>% 
  group_by(TRACK_ID) %>%
  slice_head(n=5)
avg_5 <- aggregate(avg_5$Avg_after_BGsub, list(avg_5$TRACK_ID), FUN=mean)
DATA$F5avg <- NA
z <- length(unique(DATA$TRACK_ID))
for(i in 1:z){DATA$F5avg[which(DATA$TRACK_ID==i)] <- avg_5$x[which(avg_5$Group.1==i)]}
DATA %>%
  mutate(
    baseline = as.numeric(DATA$Avg_after_BGsub)/as.numeric(DATA$F5avg),
  ) -> DATA
print(length(unique(DATA$TRACK_ID)))

#KEEP ONLY TRACK IDS WITH MAX NUMBER OF FRAMES
f = length(unique(DATA$FRAME))

DATA %>%
   group_by(TRACK_ID) %>%
   filter(n() == f) -> DATA
print(length(unique(DATA$TRACK_ID)))

#VARIANCE
var_data <- DATA[order(DATA$TRACK_ID),]
var_data$variance <- NA
var_group <- aggregate(var_data$baseline, list(var_data$TRACK_ID), FUN=var)
for(i in 1:z){var_data$variance[which(var_data$TRACK_ID==i)] <- var_group$x[which(var_group$Group.1==i)]}
VAR <- var_data

#TOP 20% OF VARIANCE
n <- 20
VgrpPCT <- var_group[var_group$x > quantile(var_group$x, prob=1-n/100),]
VgrpPCT <- VgrpPCT %>% rename(TRACK_ID = Group.1)
VAR20 <- left_join(VAR, VgrpPCT)
VAR20 <- VAR20 %>% rename(Top20Var = x)
EXP <- VAR20

#FILTER BASELINE COLUMNS
EXPt <- EXP %>% rename(EXP_base = baseline,
                      EXP_track = TRACK_ID)
EXPt <- EXPt[-c(6:12)]


#### GRAPH PLOTS ####

#Raw Experimental Data Line Plot
RAW_exp <- ggplot(CELLS, aes(x=FRAME, y=MEAN_INTENSITY_CH1, group=TRACK_ID, col=factor(TRACK_ID)))+
  geom_line(size=0.7)+
  labs(title = 'RAW EXPERIMENTAL MEAN INTENSITIES',
       x = 'Frame',
       y = 'Mean Intensity')+
  theme_cowplot(12)+
  theme(legend.position='none',
        plot.title = element_text(size= 18, hjust=0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        plot.background = element_rect(fill="white"))+
  scale_color_hue(l = 30)
ggsave("Raw_plot.jpg",  width=12, height=8, dpi=300)
RAW_exp

#Normalized Experimental Line Plot
exp20 <- ggplot(EXPt, aes(x=FRAME, y=EXP_base, group=EXP_track, col=factor(Top20Var)))+
  geom_line(size=0.7)+
  labs(title = 'Normalized Experimental Intensities',
       x = 'Frame',
       y = 'Mean Intensity',
  caption='*Top 20% of Cells with the Greatest Variance are Highlighted')+
  theme_cowplot(12)+
  theme(legend.position='none',
        plot.title = element_text(size= 18, hjust=0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        plot.background = element_rect(fill="white"))+
  scale_color_hue(l = 30)
ggsave("line_plot.jpg",  width=12, height=8, dpi=300)
exp20


#EXP Dot Plot
Dot <- ggplot(EXPt, aes(x=POSITION_X, y=POSITION_Y, col=factor(Top20Var)))+
  geom_point(size=2)+
  labs(title = "Position of Experimental Cells Over Time",
       x="X Position",
       y="Y Position",
       color="Mean Intensity",
  caption='*Top 20% of Cells with the Greatest Variance are Highlighted')+
  theme_cowplot(12)+
  theme(plot.title = element_text(size= 18, hjust=0.5),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position='none',
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill="gray99", color="black"),
        strip.background = element_rect(fill="navy"),
        strip.text.x = element_text(color="white", face="bold"),
        plot.background = element_rect(fill="white"))+
  scale_color_hue(l = 30)+
  scale_y_reverse()+
  facet_wrap(~ EXPt$FRAME,ncol=5,scales='fixed')
ggsave("Dot_Plot.jpg",  width=10, height=8, dpi=300)
Dot


###WRITE CSV FILE - FINAL STATS###
EXP_csv <- EXP
EXP_csv <- EXP_csv %>% rename(Normalized_Mean_Intensity = baseline,
                              Track_ID = TRACK_ID,
                              Frame = FRAME,
                              Mean_Intensity_Raw = MEAN_INTENSITY_CH1,
                              Median_Intensity_Raw = MEDIAN_INTENSITY_CH1,
                              Max_Intensity_Raw = MAX_INTENSITY_CH1,
                              Min_Intensity_Raw = MIN_INTENSITY_CH1,
                              Area = AREA,
                              Background_Sample_Average = BG_Avg,
                              Mean_Intensity_Post_Bckgd_Subtraction = Avg_after_BGsub,
                              First_5_Pts_Avg = F5avg,
                              Normalized_Mean_Intensity = baseline,
                              Variance = variance,
                              Top_20_Percent_Variance = Top20Var)

EXP_csv$Top_20_Percent_Variance <- ifelse(EXP_csv$Top_20_Percent_Variance>0,"Yes")
EXP_csv$Top_20_Percent_Variance <- EXP_csv$Top_20_Percent_Variance %>% replace_na("No")

write.csv(EXP_csv,'Experimental_Data.csv', row.names = FALSE)


print("Stats Run Complete!")
