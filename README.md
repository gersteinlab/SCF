# SCF

This repository contains two separate scripts written in R, along with the data sets used for our analysis and testing. The first script, titled "Andrew_SCF_Analysis.R," contains the back end code for the diabetes-related clinical sensor data analyses in the Results section. With a quick adjustment to the working directory set, this code can be run interactively via the R console. The second script, titled "software.R," functions as a command line utility which automates some of the key analyses conducted in the first script. 

Both scripts require that the user have installed the following packages: CausalImpact, tidyverse, lubridate, data.table, reshape2

The R script "software.R" requires a single, complete data file with four columns corresponding to Date, Time, Variable, and Value, in that order. Prior to running, the script and its corresponding data file must be placed in the same directory. Running the script requires the following 9 user-inputted command-line arguments: 

1) filename of 4-column file (date, time, variable, value)
2) length of interval to align data to in minutes
3) start date of analysis in YYYY-MM-DD
4) start time of analysis in HH:MM:SS
5) end date of analysis in YYYY-MM-DD
6) end time of analysis in HH:MM:SS
7) date of intervention in YYYY-MM-DD
8) time of intervention in HH:MM:SS
9) dependent variable name
