library(CausalImpact)
library(tidyverse)
library(lubridate)
library(data.table)
library(reshape2)

### BEFORE RUNNING CODE, CHANGE THIS PATH TO THE NEW LOCAL PATH
setwd("/Users/andrewzhang/Desktop/gerstein/wearableSensorData/software/")

rm(list=ls())

### DEFINING FUNCTIONS TO REPLACE NA VALUES WITH NEIGHBORING MEANS

replacena1 <- function(l) {
  stopifnot(is.numeric(l))
  indx <- is.na(l)
  l[indx] <- vapply(which(indx), function(x) {mean(c(l[x - 1], l[x + 1])) }, FUN.VALUE = double(1))
  l
}

replacena2 <- function(l) {
  stopifnot(is.numeric(l))
  indx <- is.na(l)
  l[indx] <- vapply(which(indx), function(x) {if( is.na(l[x-1]) ){ l[x+1] } else if( is.na(l[x+1])) { l[x-1] } else mean(c(l[x - 1], l[x + 1])) }, FUN.VALUE = double(1))
  l
}

replacena.while <- function(x) {
  x <- replacena1(x)
  success <- FALSE
  while (!success) {
    x <- replacena2(x)
    success <- sum(is.na(x))==0
    #print(sum(is.na(x)))
  }
  return(x)
}


### READING AND PROCESSING APPLE WATCH DATA

watch <- read.csv("parsed_full_export.csv", sep=",", header=FALSE)
watch <- watch %>% separate(V1, c("Date", "Time", "Time.Zone"), sep = " ")

# set interval you want to use
interval <- 5 #minutes

date1 <- as.POSIXct("2019-12-27 00:00:00")
date2 <- as.POSIXct("2020-03-21 00:00:00") # same original dates for combined.data
int <- interval(date1, date2) 

colnames(watch)[4:6] <- c("Variable", "Units", "Value")
watch$Variable <- as.character(watch$Variable)
watch$Variable <- substr(watch$Variable, 25, nchar(watch$Variable))


### READING AND PROCESSING GLUCOSE / INSULIN SENSOR DATA
med07.full.glucose <- read.table("MED07_Full_Glucose.csv", as.is=T, header=T, sep="\t")
med07.full.insulin <- read.table("MED07_Full_Insulin.csv", as.is=T, header=T, sep="\t")

# change the dates to POSIXct formatting. column name could vary depending on your file
med07.full.glucose$DateTime <- as.POSIXct(med07.full.glucose$Timestamp..YYYY.MM.DDThh.mm.ss.,format = "%Y-%m-%dT%H:%M:%S")
med07.full.insulin$DateTime <- as.POSIXct(med07.full.insulin$EventDateTime,format = "%Y-%m-%dT%H:%M:%S")

# align the time to the time interval
med07.full.glucose$AlignTime <- align.time(med07.full.glucose$DateTime, interval*60)
med07.full.insulin$AlignTime <- align.time(med07.full.insulin$DateTime, interval*60)
med07.full.insulin <- med07.full.insulin[-15354:-15359, ] # remove skipped values corresponding to daylight savings

# create a vector of time intervals that spans the whole study period
total.period <- format(seq.POSIXt(as.POSIXct("2019-12-20 00:00:00"), as.POSIXct("2020-03-21 00:00:00"), by = "5 min"),
                       "%Y-%m-%d %H:%M:%S")


### CREATING A DATA FRAME TO EVENTUALLY HOLD COMBINED DATA
combined.data <- data.frame(Time=total.period, Glucose=NA, IOB=NA)

# match the aligned time in the studies, to the vector that spans the whole study
combined.data[match(as.character(med07.full.glucose$AlignTime), as.character(combined.data$Time)),]$Glucose <- med07.full.glucose$Glucose.Value..mg.dL.
combined.data[match(as.character(med07.full.insulin$AlignTime), as.character(combined.data$Time)),]$IOB <- med07.full.insulin$IOB

# change the data to be numeric
combined.data$Glucose <- as.numeric(combined.data$Glucose)
combined.data$IOB <- as.numeric(combined.data$IOB)

# remove rows corresponding to missing glucose data
combined.data <- combined.data[-1:-2017, ]

# using the the replace NA functions to do the neighboring means
combined.data$Glucose <- replacena.while(combined.data$Glucose)
combined.data$IOB <- replacena.while(combined.data$IOB)

# grab only the date
combined.data$date <- as.Date(as.character(combined.data$Time))

# normalize glucose and IOB values so they can be displayed on the same plot
combined.data$Glucose.scaled <- scale(combined.data$Glucose)
combined.data$IOB.scaled <- scale(combined.data$IOB)

# extract hour for variance analysis
combined.data$HOUR <- strftime(combined.data$Time, format="%H")
combined.data$date.HOUR <- paste(combined.data$date, combined.data$HOUR, sep=" ")

# aggregate and calculate variance of metrics per hour
ag_hourly_variance <- aggregate(combined.data[, 2:3], list(combined.data$date.HOUR), FUN=var)

# match back to original data frame
combined.data$Glucose.variance <- ag_hourly_variance$Glucose[match(combined.data$date.HOUR, ag_hourly_variance$Group.1)]
combined.data$IOB.variance <- ag_hourly_variance$IOB[match(combined.data$date.HOUR, ag_hourly_variance$Group.1)]


### PROCESSING AND ADDING APPLE WATCH VARIABLES TO COMBINED DATA

# Heart Rate Variability
watch.HRV <- watch[which(watch$Variable=="HeartRateVariabilitySDNN"),]
watch.HRV$AlignTime <- align.time(as.POSIXct(paste(watch.HRV$Date,watch.HRV$Time, sep=" ")),interval*60)
watch.HRV.new <- watch.HRV[watch.HRV$AlignTime %within% int,]

combined.data$HRV <- NA
combined.data[match(as.character(watch.HRV.new$AlignTime), as.character(combined.data$Time)),]$HRV <- as.numeric(as.character(watch.HRV.new$Value))
if(is.na(combined.data$HRV[1])){
  combined.data$HRV[1] <- na.omit(combined.data$HRV)[1]
}
if(is.na(combined.data$HRV[length(combined.data$HRV)])){
  combined.data$HRV[length(combined.data$HRV)] <- na.omit(combined.data$HRV)[length(na.omit(combined.data$HRV))]
}
combined.data$HRV <- replacena.while(combined.data$HRV)

# StepCount
watch.StepCount <- watch[which(watch$Variable=="StepCount"),]
watch.StepCount$AlignTime <- align.time(as.POSIXct(paste(watch.StepCount$Date,watch.StepCount$Time,sep=" ")),interval*60)
watch.StepCount.new <- watch.StepCount[watch.StepCount$AlignTime %within% int,]
combined.data$StepCount <- NA
combined.data[match(as.character(watch.StepCount.new$AlignTime), as.character(combined.data$Time)),]$StepCount <- as.numeric(as.character(watch.StepCount.new$Value))
if(is.na(combined.data$StepCount[1])){
  combined.data$StepCount[1] <- na.omit(combined.data$StepCount)[1]
}
if(is.na(combined.data$StepCount[length(combined.data$StepCount)])){
  combined.data$StepCount[length(combined.data$StepCount)] <- na.omit(combined.data$StepCount)[length(na.omit(combined.data$StepCount))]
}
combined.data$StepCount <- replacena.while(combined.data$StepCount)

# ActiveEnergyBurned 
watch.ActiveEnergyBurned <- watch[which(watch$Variable=="ActiveEnergyBurned"),]
watch.ActiveEnergyBurned$AlignTime <- align.time(as.POSIXct(paste(watch.ActiveEnergyBurned$Date,watch.ActiveEnergyBurned$Time,sep=" ")),interval*60)
watch.ActiveEnergyBurned.new <- watch.ActiveEnergyBurned[watch.ActiveEnergyBurned$AlignTime %within% int,]
combined.data$ActiveEnergyBurned <- NA
combined.data[match(as.character(watch.ActiveEnergyBurned.new$AlignTime), as.character(combined.data$Time)),]$ActiveEnergyBurned <- as.numeric(as.character(watch.ActiveEnergyBurned.new$Value))
if(is.na(combined.data$ActiveEnergyBurned[1])){
  combined.data$ActiveEnergyBurned[1] <- na.omit(combined.data$ActiveEnergyBurned)[1]
}
if(is.na(combined.data$ActiveEnergyBurned[length(combined.data$ActiveEnergyBurned)])){
  combined.data$ActiveEnergyBurned[length(combined.data$ActiveEnergyBurned)] <- na.omit(combined.data$ActiveEnergyBurned)[length(na.omit(combined.data$ActiveEnergyBurned))]
}
combined.data$ActiveEnergyBurned <- replacena.while(combined.data$ActiveEnergyBurned)

# BasalEnergyBurned
watch.BasalEnergyBurned <- watch[which(watch$Variable=="BasalEnergyBurned"),]
watch.BasalEnergyBurned$AlignTime <- align.time(as.POSIXct(paste(watch.BasalEnergyBurned$Date,watch.BasalEnergyBurned$Time,sep=" ")),interval*60)
watch.BasalEnergyBurned.new <- watch.BasalEnergyBurned[watch.BasalEnergyBurned$AlignTime %within% int,]
combined.data$BasalEnergyBurned <- NA
combined.data[match(as.character(watch.BasalEnergyBurned.new$AlignTime), as.character(combined.data$Time)),]$BasalEnergyBurned <- as.numeric(as.character(watch.BasalEnergyBurned.new$Value))
if(is.na(combined.data$BasalEnergyBurned[1])){
  combined.data$BasalEnergyBurned[1] <- na.omit(combined.data$BasalEnergyBurned)[1]
}
if(is.na(combined.data$BasalEnergyBurned[length(combined.data$BasalEnergyBurned)])){
  combined.data$BasalEnergyBurned[length(combined.data$BasalEnergyBurned)] <- na.omit(combined.data$BasalEnergyBurned)[length(na.omit(combined.data$BasalEnergyBurned))]
}
combined.data$BasalEnergyBurned <- replacena.while(combined.data$BasalEnergyBurned)

# DistanceWalkingRunning
watch.DistanceWalkingRunning <- watch[which(watch$Variable=="DistanceWalkingRunning"),]
watch.DistanceWalkingRunning$AlignTime <- align.time(as.POSIXct(paste(watch.DistanceWalkingRunning$Date,watch.DistanceWalkingRunning$Time,sep=" ")),interval*60)
watch.DistanceWalkingRunning.new <- watch.DistanceWalkingRunning[watch.DistanceWalkingRunning$AlignTime %within% int,]
combined.data$DistanceWalkingRunning <- NA
combined.data[match(as.character(watch.DistanceWalkingRunning.new$AlignTime), as.character(combined.data$Time)),]$DistanceWalkingRunning <- as.numeric(as.character(watch.DistanceWalkingRunning.new$Value))
if(is.na(combined.data$DistanceWalkingRunning[1])){
  combined.data$DistanceWalkingRunning[1] <- na.omit(combined.data$DistanceWalkingRunning)[1]
}
if(is.na(combined.data$DistanceWalkingRunning[length(combined.data$DistanceWalkingRunning)])){
  combined.data$DistanceWalkingRunning[length(combined.data$DistanceWalkingRunning)] <- na.omit(combined.data$DistanceWalkingRunning)[length(na.omit(combined.data$DistanceWalkingRunning))]
}
combined.data$DistanceWalkingRunning <- replacena.while(combined.data$DistanceWalkingRunning)

# HeartRate
watch.HeartRate <- watch[which(watch$Variable=="HeartRate"),]
watch.HeartRate$AlignTime <- align.time(as.POSIXct(paste(watch.HeartRate$Date,watch.HeartRate$Time,sep=" ")),interval*60)
watch.HeartRate.new <- watch.HeartRate[watch.HeartRate$AlignTime %within% int,]
combined.data$HeartRate <- NA
combined.data[match(as.character(watch.HeartRate.new$AlignTime), as.character(combined.data$Time)),]$HeartRate <- as.numeric(as.character(watch.HeartRate.new$Value))
if(is.na(combined.data$HeartRate[1])){
  combined.data$HeartRate[1] <- na.omit(combined.data$HeartRate)[1]
}
if(is.na(combined.data$HeartRate[length(combined.data$HeartRate)])){
  combined.data$HeartRate[length(combined.data$HeartRate)] <- na.omit(combined.data$HeartRate)[length(na.omit(combined.data$HeartRate))]
}
combined.data$HeartRate <- replacena.while(combined.data$HeartRate)

# FlightsClimbed
watch.FlightsClimbed <- watch[which(watch$Variable=="FlightsClimbed"),]
watch.FlightsClimbed$AlignTime <- align.time(as.POSIXct(paste(watch.FlightsClimbed$Date,watch.FlightsClimbed$Time,sep=" ")),interval*60)
watch.FlightsClimbed.new <- watch.FlightsClimbed[watch.FlightsClimbed$AlignTime %within% int,]
combined.data$FlightsClimbed <- NA
combined.data[match(as.character(watch.FlightsClimbed.new$AlignTime), as.character(combined.data$Time)),]$FlightsClimbed <- as.numeric(as.character(watch.FlightsClimbed.new$Value))
if(is.na(combined.data$FlightsClimbed[1])){
  combined.data$FlightsClimbed[1] <- na.omit(combined.data$FlightsClimbed)[1]
}
if(is.na(combined.data$FlightsClimbed[length(combined.data$FlightsClimbed)])){
  combined.data$FlightsClimbed[length(combined.data$FlightsClimbed)] <- na.omit(combined.data$FlightsClimbed)[length(na.omit(combined.data$FlightsClimbed))]
}
combined.data$FlightsClimbed <- replacena.while(combined.data$FlightsClimbed)

# AppleExerciseTime
watch.AppleExerciseTime <- watch[which(watch$Variable=="AppleExerciseTime"),]
watch.AppleExerciseTime$AlignTime <- align.time(as.POSIXct(paste(watch.AppleExerciseTime$Date,watch.AppleExerciseTime$Time,sep=" ")),interval*60)
watch.AppleExerciseTime.new <- watch.AppleExerciseTime[watch.AppleExerciseTime$AlignTime %within% int,]
combined.data$AppleExerciseTime <- NA
combined.data[match(as.character(watch.AppleExerciseTime.new$AlignTime), as.character(combined.data$Time)),]$AppleExerciseTime <- as.numeric(as.character(watch.AppleExerciseTime.new$Value))
if(is.na(combined.data$AppleExerciseTime[1])){
  combined.data$AppleExerciseTime[1] <- na.omit(combined.data$AppleExerciseTime)[1]
}
if(is.na(combined.data$AppleExerciseTime[length(combined.data$AppleExerciseTime)])){
  combined.data$AppleExerciseTime[length(combined.data$AppleExerciseTime)] <- na.omit(combined.data$AppleExerciseTime)[length(na.omit(combined.data$AppleExerciseTime))]
}
combined.data$AppleExerciseTime <- replacena.while(combined.data$AppleExerciseTime)


### EXPLORATORY DATA ANALYSIS AND VISUALIZATION

# plot of glucose and insulin over the entire 12 week period
tmp <- combined.data[, c("Time", "Glucose.scaled", "IOB.scaled")]
tmp <- tmp %>% pivot_longer(c('Glucose.scaled', 'IOB.scaled'), names_to = "AppleWatchVariable", values_to="Value")
tmp %>% ggplot() + 
  geom_line(mapping = aes(x = Time, y = Value, group = AppleWatchVariable, color = AppleWatchVariable)) + 
  ggtitle("Glucose and Insulin vs. Time, Full 12 Weeks") + ylab("Z-score") + xlab("Time") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# plot of hourly variance against time for glucose and IOB, full 12 weeks
tmp <- combined.data[, c("Time", "Glucose.variance", "IOB.variance")]
tmp <- tmp %>% pivot_longer(c("Glucose.variance", "IOB.variance"), names_to = "AppleWatchVariable", values_to="Value")
ggplot(tmp) + geom_line(mapping = aes(x = Time, y = Value, group = AppleWatchVariable, color = AppleWatchVariable)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  ggtitle("Glucose and Insulin Hourly Variance vs. Time, Full 12 Weeks") + 
  ylab("Variance") + xlab("Time")


### CAUSAL IMPACT ANALYSES

# causal impact of exercise intervention, using the raw, untransformed values of glucose and IOB
period1 <- combined.data
intervention.exercise <- 4032
y <- (period1$Glucose)
x1 <- period1$IOB
post.period <- c(intervention.exercise, dim(period1)[1])
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA
y <- as.numeric(y)
x1 <- as.numeric(x1)
df.tmp <- data.frame(y,x1)
ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1, ss, niter = 1000)

impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response),
                       model.args = list(nseasons = 7, season.duration = 1,dynamic.regression=T))
plot(impact) + ggtitle("12 Week Causal Impact Analysis Of Exercise Using 5-Minute Glucose and IOB Readings")
summary(impact,"report")

# causal impact of exercise on hourly variance
intervention.exercise <- 361  # exercise regime begins on the day corresopnding to row 361 of ag_hourly_variance
ag_hourly_variance <- ag_hourly_variance[-2040, ] # remove the last row which has NAs
y <- ag_hourly_variance$Glucose
x1 <- ag_hourly_variance$IOB
post.period <- c(intervention.exercise, nrow(ag_hourly_variance))
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA
y <- as.numeric(y)
x1 <- as.numeric(x1)

ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1, ss, niter = 1000)

impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response),
                       model.args = list(nseasons = 7, season.duration = 1,dynamic.regression=T))
plot(impact) + ggtitle("12 Week Causal Impact Analysis of Exercise Using Hourly Glucose and IOB Variance")
summary(impact,"report")

# now let's try to aggregate by day, using mean value
ag_daily_mean <- aggregate(combined.data[, 2:3], list(combined.data$date), FUN=mean)
colnames(ag_daily_mean)[-1] <- paste(colnames(ag_daily_mean)[-1], "Mean", sep=".")
intervention.exercise <- 16    # exercise regime begins on day 16
y <- ag_daily_mean$Glucose.Mean
x1 <- ag_daily_mean$IOB.Mean
post.period <- c(intervention.exercise, nrow(ag_daily_mean))
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA
y <- as.numeric(y)
x1 <- as.numeric(x1)

ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1, ss, niter = 1000)

impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response),
                       model.args = list(nseasons = 7, season.duration = 1,dynamic.regression=T))
plot(impact) + ggtitle("12 Week Causal Impact Analysis of Exercise Using 24 Hour Glucose and IOB Mean")
ggsave(filename="Causal_Impact_Daily_Mean.pdf", device="pdf", width=16, height=10)
summary(impact,"report")

# Look at significant predictors of daily glucose variance
ag_daily_variance <- aggregate(combined.data[, c(2:3, 11:18)], list(combined.data$date), FUN=var)
colnames(ag_daily_variance)[c(-1, -11)] <- paste(colnames(ag_daily_variance)[c(-1, -11)], "Variance", sep=".")
y <- (ag_daily_variance$Glucose.Variance)
x1 <- as.numeric(ag_daily_variance$IOB.Variance)
x2 <- as.numeric(ag_daily_variance$HRV.Variance)
x3 <- as.numeric(ag_daily_variance$StepCount.Variance)
x4 <- as.numeric(ag_daily_variance$ActiveEnergyBurned.Variance)
x5 <- as.numeric(ag_daily_variance$BasalEnergyBurned.Variance)
x6 <- as.numeric(ag_daily_variance$DistanceWalkingRunning.Variance)
x7 <- as.numeric(ag_daily_variance$HeartRate.Variance)
x8 <- as.numeric(ag_daily_variance$FlightsClimbed.Variance)
x9 <- as.numeric(ag_daily_variance$AppleExerciseTime)
lm1 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9)

# significant variables include IOB and FlightsClimbed (p-value < 0.001)
# read in some new total insulin data we have
new_insulin <- read.table("MED07_Newest_Insulin.csv", sep=",", header=TRUE)
new_insulin <- new_insulin[which(!is.na(new_insulin$Total.Insulin)), ]
new_insulin$Date <- as.character(mdy(new_insulin$Date))

colnames(ag_daily_variance)[1] <- "Date"
ag_daily_variance$Date <- as.character(ag_daily_variance$Date)
ag_daily_variance$new_daily_insulin <- new_insulin$Total.Insulin[match(ag_daily_variance$Date, new_insulin$Date)]
x1 <- as.numeric(ag_daily_variance$new_daily_insulin)
lm2 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9) 
# our original IOB is better at predicting glucose than the new Total.Insulin by a factor of 100!
# thus, we'll just use IOB instead of this new insulin for the following analysis


### TARGET RANGE ANALYSIS. Clinically, the target glucose level is between 70-180.
### Here, we are interested in seeing what fraction of glucose readings fall within the target range per day

glucose_in_target <- aggregate(combined.data[, 2], list(combined.data$date), FUN=function(x){y <- sum(x > 70 & x < 180)/length(x); return(y)})
colnames(glucose_in_target) <- c("Date", "PercentGlucoseInTargetRange")
ggplot(glucose_in_target) + geom_line(mapping = aes(x = Date, y = PercentGlucoseInTargetRange, group = 1)) +
  geom_vline(xintercept = as.numeric(glucose_in_target$Date[15]), color="red", linetype="dashed") + 
  ggtitle("Percent Daily Glucose Readings Above Target Range, Full 12 Weeks") +
  ylab("Percent Glucose Readings") + xlab("24 Hour Bins")

# redo causal impact analysis, but consider this "percent in target range" as response variable

intervention.exercise <- 16
period1 <- data.frame(glucose_in_target$PercentGlucoseInTargetRange, ag_daily_mean$IOB.Mean)[-86, ]
colnames(period1) <- c("PercentGlucoseInTargetRange", "IOB.Mean")
y <- period1$PercentGlucoseInTargetRange
x1 <- period1$IOB.Mean

post.period <- c(intervention.exercise, dim(period1)[1])
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA
y <- as.numeric(y)
x1 <- as.numeric(x1)
df.tmp <- data.frame(y,x1)

ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1, ss, niter = 1000)

impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response))
plot(impact) + ggtitle("Causal Impact Analysis: Percent Daily Glucose Readings In Target Range, Full 12 Weeks")
summary(impact,"report")


### causal impact analysis, with percent of glucose readings ABOVE target range as response variable

glucose_above_target <- aggregate(combined.data[, 2], list(combined.data$date), FUN=function(x){y <- sum(x > 180)/length(x); return(y)})
colnames(glucose_above_target) <- c("Date", "PercentGlucoseAboveTargetRange")
ag_daily_variance$glucose.above.target <- glucose_above_target$PercentGlucoseAboveTargetRange

period1 <- ag_daily_variance[-86, ]
y <- (period1$glucose.above.target)
x1 <- period1$IOB
post.period <- c(intervention.exercise, dim(period1)[1])
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA
y <- as.numeric(y)
x1 <- as.numeric(x1)
df.tmp <- data.frame(y,x1)

ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1, ss, niter = 1000)

impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response))
plot(impact) + ggtitle("Causal Impact Analysis: Percent Daily Glucose Readings Above Target Range, Full 12 Weeks")
summary(impact,"report")


### SOME CORRELATION ANALYSES FOR COVARIATES

# getting correlation matrix between apple watch variables, aggregated by variance over 24 hours
fwrite(cor(period1[, c(-1, -11)]), file="cor.csv", sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

# plotting the relationship between glucose.in.target and IOB.mean
ag_daily_mean$Group.1 <- as.character(ag_daily_mean$Group.1)
tmp <- data.frame(glucose_in_target$PercentGlucoseInTargetRange, ag_daily_mean$IOB.Mean)
ggplot(tmp, aes(x=ag_daily_mean.IOB.Mean, y = glucose_in_target.PercentGlucoseInTargetRange)) + geom_point() + geom_smooth(method=lm) + 
  ggtitle("Daily Fraction of Glucose In Target vs. Daily Mean IOB, Full 12 Weeks")

# plotting the relationship between glucose.in.target and IOB.variance
tmp <- data.frame(glucose_in_target$PercentGlucoseInTargetRange, ag_daily_variance$IOB.Variance)
ggplot(tmp, aes(x=ag_daily_variance.IOB.Variance, y = glucose_in_target.PercentGlucoseInTargetRange)) + geom_point() + geom_smooth(method=lm) + 
  ggtitle("Daily Fraction of Glucose In Target vs. Daily Variance of IOB, Full 12 Weeks")

# trying some new covariates
ag_daily_mean <- aggregate(combined.data[, c(2:3, 11:18)], list(combined.data$date), FUN=mean)
ag_daily_sum <- aggregate(combined.data[, c(2:3, 11:18)], list(combined.data$date), FUN=sum)
ag_daily_min <- aggregate(combined.data[, c(2:3, 11:18)], list(combined.data$date), FUN=min)
ag_daily_max <- aggregate(combined.data[, c(2:3, 11:18)], list(combined.data$date), FUN=max)
ag_daily_range <- ag_daily_max - ag_daily_min 
ag_daily_range$Group.1 <- ag_daily_max$Group.1

y <- as.numeric(glucose_in_target$PercentGlucoseInTargetRange)
x1 <- as.numeric(ag_daily_variance$IOB.Variance)
x2 <- as.numeric(ag_daily_variance$HRV.Variance)
x3 <- as.numeric(ag_daily_variance$StepCount.Variance)
x4 <- as.numeric(ag_daily_variance$ActiveEnergyBurned.Variance)
x5 <- as.numeric(ag_daily_variance$BasalEnergyBurned.Variance)
x6 <- as.numeric(ag_daily_variance$DistanceWalkingRunning.Variance)
x7 <- as.numeric(ag_daily_variance$HeartRate.Variance)
x8 <- as.numeric(ag_daily_variance$FlightsClimbed.Variance)
cor(cbind(y, x1, x2, x3, x4, x5, x6, x7, x8)[-86, ]) 
lm1 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8)

y <- as.numeric(glucose_in_target$PercentGlucoseInTargetRange)
x1 <- as.numeric(ag_daily_mean$IOB)
x2 <- as.numeric(ag_daily_mean$HRV)
x3 <- as.numeric(ag_daily_mean$StepCount)
x4 <- as.numeric(ag_daily_mean$ActiveEnergyBurned)
x5 <- as.numeric(ag_daily_mean$BasalEnergyBurned)
x6 <- as.numeric(ag_daily_mean$DistanceWalkingRunning)
x7 <- as.numeric(ag_daily_mean$HeartRate)
x8 <- as.numeric(ag_daily_mean$FlightsClimbed)
cor(cbind(y, x1, x2, x3, x4, x5, x6, x7, x8)[-86, ])
lm2 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8)

y <- as.numeric(glucose_in_target$PercentGlucoseInTargetRange)
x1 <- as.numeric(ag_daily_min$IOB)
x2 <- as.numeric(ag_daily_min$HRV)
x3 <- as.numeric(ag_daily_min$StepCount)
x4 <- as.numeric(ag_daily_min$ActiveEnergyBurned)
x5 <- as.numeric(ag_daily_min$BasalEnergyBurned)
x6 <- as.numeric(ag_daily_min$DistanceWalkingRunning)
x7 <- as.numeric(ag_daily_min$HeartRate)
x8 <- as.numeric(ag_daily_min$FlightsClimbed)
cor(cbind(y, x1, x2, x3, x4, x5, x6, x7, x8)[-86, ]) 
lm3 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8)

y <- as.numeric(glucose_in_target$PercentGlucoseInTargetRange)
x1 <- as.numeric(ag_daily_max$IOB) 
x2 <- as.numeric(ag_daily_max$HRV) # significant
x3 <- as.numeric(ag_daily_max$StepCount)
x4 <- as.numeric(ag_daily_max$ActiveEnergyBurned)
x5 <- as.numeric(ag_daily_max$BasalEnergyBurned)
x6 <- as.numeric(ag_daily_max$DistanceWalkingRunning)
x7 <- as.numeric(ag_daily_max$HeartRate)
x8 <- as.numeric(ag_daily_max$FlightsClimbed)
cor(cbind(y, x1, x2, x3, x4, x5, x6, x7, x8)[-86, ])
lm4 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8)

y <- as.numeric(glucose_in_target$PercentGlucoseInTargetRange)
x1 <- as.numeric(ag_daily_range$IOB)
x2 <- as.numeric(ag_daily_range$HRV)
x3 <- as.numeric(ag_daily_range$StepCount)
x4 <- as.numeric(ag_daily_range$ActiveEnergyBurned)
x5 <- as.numeric(ag_daily_range$BasalEnergyBurned)
x6 <- as.numeric(ag_daily_range$DistanceWalkingRunning)
x7 <- as.numeric(ag_daily_range$HeartRate)
x8 <- as.numeric(ag_daily_range$FlightsClimbed)
cor(cbind(y, x1, x2, x3, x4, x5, x6, x7, x8)[-86, ])
lm5 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8)

y <- as.numeric(glucose_in_target$PercentGlucoseInTargetRange)
x1 <- as.numeric(ag_daily_sum$IOB) # significant
x2 <- as.numeric(ag_daily_sum$HRV)
x3 <- as.numeric(ag_daily_sum$StepCount)
x4 <- as.numeric(ag_daily_sum$ActiveEnergyBurned)
x5 <- as.numeric(ag_daily_sum$BasalEnergyBurned)
x6 <- as.numeric(ag_daily_sum$DistanceWalkingRunning)
x7 <- as.numeric(ag_daily_sum$HeartRate)
x8 <- as.numeric(ag_daily_sum$FlightsClimbed)
cor(cbind(y, x1, x2, x3, x4, x5, x6, x7, x8)[-86, ])
lm6 <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8)

### Obtaining inclusion probabilities for model

GetInclusionProbabilities <- function(bsts.object) {
  burn <- SuggestBurn(0.1, bsts.object)
  beta <- bsts.object$coefficients
  beta <- beta[-(1:burn), , drop = FALSE]
  inclusion.prob <- colMeans(beta != 0)
  index <- order(inclusion.prob)
  inclusion.prob <- inclusion.prob[index]
  return(data.frame(predictor = names(inclusion.prob),
                    inclusion.prob = inclusion.prob))
}

### CAUSAL IMPACT WITH NEW COVARIATES, "PERCENT IN TARGET RANGE" AS RESPONSE VAR

intervention.exercise <- 16
period1 <- data.frame(glucose_in_target$PercentGlucoseInTargetRange, ag_daily_mean$IOB,
                      ag_daily_range$HRV, ag_daily_variance$IOB.Variance, ag_daily_min$DistanceWalkingRunning)
colnames(period1) <- c("PercentGlucoseInTargetRange", "IOB.Mean", "HRV.Range", "IOB.Variance", "DistanceWalkingRunning.Min")

period1$PercentGlucoseAboveTargetRange <- aggregate(combined.data[, 2], list(combined.data$date), FUN=function(x){y <- sum(x > 180)/length(x); return(y)})$x
period1 <- period1[-86, ] # last row has some NAs

y <- as.numeric(period1$PercentGlucoseInTargetRange)
x1 <- as.numeric(period1$IOB.Mean)
x2 <- as.numeric(period1$HRV.Range)
x3 <- as.numeric(period1$IOB.Variance)
x4 <- as.numeric(period1$DistanceWalkingRunning.Min)

post.period <- c(intervention.exercise, dim(period1)[1])
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA

ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1+x2+x3+x4, ss, niter = 1000)
GetInclusionProbabilities(bsts.model)
impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response))
plot(impact) + ggtitle("Causal Impact Analysis: Percent Daily Glucose Readings In Target Range, Full 12 Weeks")
ggsave(filename="Figure_5b_Causal_Impact_Glucose_In_Target.pdf", device="pdf", width=16, height=10)
summary(impact,"report")


# causal impact analysis, with percent of glucose readings ABOVE target range as response variable
y <- (period1$PercentGlucoseAboveTargetRange)

post.period <- c(intervention.exercise, dim(period1)[1])
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA
y <- as.numeric(y)

ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1+x2+x3+x4, ss, niter = 1000)
GetInclusionProbabilities(bsts.model)
impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response))
plot(impact) + ggtitle("Causal Impact Analysis: Percent Daily Glucose Readings Above Target Range, Full 12 Weeks")
ggsave(filename="Figure_5c_Causal_Impact_Glucose_Above_Target.pdf", device="pdf", width=16, height=10)
summary(impact,"report")


### RECREATE PLOT OF GLUCOSE, INSULIN OVER 12 WEEKS, NOW WITH SIGNIFICANT COVARIATES TOO 
ag_daily_variance$Date <- as.Date(ag_daily_variance$Date)
tmp <- data.frame(combined.data$Time, combined.data$Glucose, combined.data$IOB)
tmp$mean_daily_IOB <- ag_daily_mean$IOB[match(combined.data$date, ag_daily_mean$Group.1)]
tmp$range_daily_HRV <- ag_daily_range$HRV[match(combined.data$date, ag_daily_range$Group.1)]
tmp$variance_daily_IOB <- ag_daily_variance$IOB.Variance[match(combined.data$date, ag_daily_variance$Date)]
tmp$min_daily_DistanceWalkingRunning <- ag_daily_min$DistanceWalkingRunning[match(combined.data$date, ag_daily_min$Group.1)]
tmp[, -1] <- scale(tmp[, -1])
colnames(tmp) <- c("Time", "Glucose", "IOB", "Daily IOB Mean", "Daily HRV Range",
                   "Daily IOB Variance", "Daily Distance Walking and Running Min")

tmp <- tmp %>% pivot_longer(c('Glucose', 'IOB', 'Daily IOB Mean', 'Daily HRV Range',
                              'Daily IOB Variance', 'Daily Distance Walking and Running Min'), names_to = "Variable", values_to="Value")
tmp %>% ggplot() + 
  geom_line(mapping = aes(x = Time, y = Value, group = Variable, color = Variable)) + 
  ggtitle("Glucose, Insulin, and Significant Covariates vs. Time, Full 12 Weeks") + ylab("Value") + xlab("Time") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(filename="Figure_5a_Covariates_12_Weeks.pdf", device="pdf", width=16, height=8)



