rm(list=ls())
library(CausalImpact)
library(tidyverse)
library(lubridate)
library(data.table)
library(reshape2)

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


# command line arguments: 
# 1) filename of catted 4-column file (date, time, variable, value)
# YYYY-MM-DD HH:MM:SS
# 2) length of interval to align data to in minutes
# 3) start date of analysis in YYYY-MM-DD
# 4) start time of analysis in HH:MM:SS
# 5) end date of analysis in YYYY-MM-DD
# 6) end time of analysis in HH:MM:SS
# 7) date of intervention in YYYY-MM-DD
# 8) time of intervention in HH:MM:SS
# 9) dependent variable

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9)
{
  print("Incorrect number of arguments supplied.")
  quit()
}

if (file.exists(args[1]))
{
  combined.data <- fread(args[1])
  if (ncol(combined.data) == 4)
  {
    colnames(combined.data) <- c("Date", "Time", "Variable", "Value")
  } else {
    print("Incorrect number of columns detected in file. Make sure the file has Date, Time, Variable, and Value columns.")
    quit()
  }
} else {
  print("File not found. Did you forget to include the path with the file?")
  quit()
}

interval <- as.integer(args[2])
if (is.na(interval) | interval < 1)
{
  print("Invalid interval found. Please enter a positive integer value for interval.")
  quit()
}

date1 <- as.POSIXct(paste(args[3], args[4]))
date2 <- as.POSIXct(paste(args[5], args[6]))
intervention.date <- paste(args[7], args[8])
int <- interval(date1, date2) 

if (!identical(which(duplicated(combined.data)), integer(0))) {
  combined.data <- combined.data[-which(duplicated(combined.data)), ] 
}

vars <- unique(combined.data$Variable)

total.period <- format(seq.POSIXt(date1, date2, by = paste(interval, "min")),
                       "%Y-%m-%d %H:%M:%S")
aligned.data <- data.frame(Time=total.period)

print("ALIGNING AND CLEANING DATA...")
for (i in 1:length(vars))
{
  tmp <- combined.data[which(combined.data$Variable == vars[i]), ]
  tmp$AlignTime <- align.time(as.POSIXct(paste(tmp$Date,tmp$Time, sep=" ")), interval*60)
  tmp.new <- tmp[tmp$AlignTime %within% int,]
  if (nrow(tmp.new) > 0) {
    aligned.data[, vars[i]] <- NA
    aligned.data[match(as.character(tmp.new$AlignTime), as.character(aligned.data$Time)), vars[i]] <- as.numeric(as.character(tmp.new$Value))
    if (is.na(aligned.data[1, vars[i]])){
      aligned.data[1, vars[i]] <- na.omit(aligned.data[, vars[i]])[1]
    }
    if(is.na(aligned.data[nrow(aligned.data), vars[i]])){
      aligned.data[nrow(aligned.data), vars[i]] <- na.omit(aligned.data[, vars[i]])[length(na.omit(aligned.data[, vars[i]]))]
    }
    aligned.data[, vars[i]] <- replacena.while(aligned.data[, vars[i]])
  }
  print(vars[i])
}

print("ANALYZING CAUSAL IMPACT...")
## CAUSAL IMPACT ANALYSIS
intervention.row <- which(as.character(aligned.data$Time) == as.character(intervention.date))
y <- aligned.data[, args[9]] # dependent variable
x1 <- aligned.data[, -which(names(aligned.data) %in% c(args[9], "Time"))] # covariates
post.period <- c(intervention.row, dim(aligned.data)[1])
post.period.response <- y[post.period[1] : post.period[2]]
y[post.period[1] : post.period[2]] <- NA
y <- as.numeric(y)
x1 <- as.matrix(x1)
df.tmp <- data.frame(y,x1)
ss <- AddLocalLinearTrend(list(), y)
ss <- AddLocalLevel(list(), y)
bsts.model <- bsts(y~x1, state.specification=ss, data=df.tmp, niter = 1000)

impact <- CausalImpact(bsts.model = bsts.model,
                       post.period.response = as.numeric(post.period.response),
                       model.args = list(nseasons = 7, season.duration = 1,dynamic.regression=T))

print("SAVING ANALYSIS RESULTS TO DIRECTORY")
title <- paste("Causal Impact of Intervention On", args[9])
p <- plot(impact) + ggtitle(title)
ggsave(filename = "CausalImpact.pdf", plot = p, device = "pdf", width = 10, height = 7)

title <- "CausalImpactSummary.txt"
fileConn <- file(title)
writeLines(capture.output(summary(impact, "report")), fileConn)
close(fileConn)

print("DONE!")
