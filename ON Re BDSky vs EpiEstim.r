#### Producing a Ontario Re figure
#### Load required libraries #####
# If not installed uncomment.
# %%R
# install.packages("devtools")
# library(devtools)

# devtools::install_github("laduplessis/bdskytools")
library(coda)
library(bdskytools)
library(beastio)
library(lubridate)
library(plotly) 
library(grid)
library(zoo)
library(tidyverse)
library(EpiEstim)

ON_youngest_tip = 2022.7615296803654
all_youngest_tip = ON_youngest_tip


#### Load ON-2022 logfiles ####
fname <- "beast_xmls/ON-2022/ON_combined_run_1_3_4.log" 
ON_trace  <- beastio::readLog(fname, burnin=0.1)

#### Get origin esitmates and enter tip dates ####
ON_origin_HPD <-  beastio::getHPDMedian(ON_trace[, "origin_BDSKY_Serial"])
ON_origin_med <- ON_origin_HPD[2]

#### Get esitmates of Re for relevent timeframes ####
ON_Re_sky    <- beastio::getLogFileSubset(ON_trace, "reproductiveNumber_BDSKY_Serial")
ON_Re_hpd    <- t(beastio::getHPDMedian(ON_Re_sky))
ON_gridTimes  <- seq(0, ON_origin_med, length.out=100)  
ON_Re_gridded <- mcmc(bdskytools::gridSkyline(ON_Re_sky, ON_trace[, "origin_BDSKY_Serial"], ON_gridTimes))
ON_Re_gridded_hpd <- t(getHPDMedian(ON_Re_gridded))
ON_times <- ON_youngest_tip - ON_gridTimes


#### Figure
# Convert to dates Month Year
ON_dates = lubridate::date_decimal(ON_times)

# Data for origin KDE plots
ON_origin_all = ON_youngest_tip - ON_trace[, "origin_BDSKY_Serial"]

#### Plotting figure Epidimiological estimates
x_min_date = floor_date(min(ON_dates),'month') 
x_max_date = max(ON_dates)
x_max_date = ceiling_date(x_max_date ,'month')
date_ticks = seq(as.POSIXct(x_min_date),as.POSIXct(x_max_date),by='months')
y_min = 0
y_max = max(ON_Re_gridded_hpd)
x_min_decimal = lubridate::decimal_date(x_min_date)
x_max_decimal = lubridate::decimal_date(x_max_date)



plotSkyline(ON_dates, ON_Re_gridded_hpd, type='smooth', xlab="", xaxt='n',   ylab=expression("Ontario R"["t"]), lwd=2,
            col='blue', fill=adjustcolor("blue", alpha=0.25), xlim = c(x_min_date, x_max_date), ylims=c(y_min, y_max))
axis.POSIXct(side=1,at=date_ticks,
             format = '%y-%b', las=2)


#### From Case inicdence ####
#Read data
df0 <- read.csv("data/data_case_counts.csv")
str(df0)


incidence_df<-(df0
     %>% na.omit()
     %>% rename("dates"="sampled_date", "I"="case_counts")
     %>% mutate(dates=as.Date(dates, format = "%d-%b-%y"))
     %>% complete(dates = seq.Date(min(dates), max(dates), by="day"))
     %>% replace_na(list(I=0))
) 

#### Plot of WGSs and case incidence ####
WGS_metadata <- read.csv("data/PHO_metadata.csv")
WGS_date_counts = as.data.frame(table(WGS_metadata$Collection_Date))


WGS_date_counts <-(WGS_date_counts
                   %>% rename('dates'="Var1", "Sequenced"="Freq")
                   %>% mutate(dates=as.Date(dates, format = "%Y-%m-%d"))
                   %>% complete(dates = seq.Date(min(incidence_df$dates), max(incidence_df$dates), by="day"))
                   %>% replace_na(list(WGS=0))
)

case_and_WGS_dates = merge(incidence_df,WGS_date_counts, by='dates')
min_date = as.Date(min(ON_dates))
case_and_WGS_dates = (case_and_WGS_dates
                      %>% complete(dates = seq.Date(min(incidence_df$dates), max(incidence_df$dates), by="day"))
                      %>% replace_na(list(I=0,Sequenced=0))
)
case_and_WGS_dates$Not_Sequenced <- case_and_WGS_dates$I - case_and_WGS_dates$Sequenced
dates <- case_and_WGS_dates$dates
case_and_WGS_dates_t = t(case_and_WGS_dates[,-1:-2])
colnames(case_and_WGS_dates_t) <- dates

# As barplot
par(mar = c(6, 4, 1, 1))
ticks = barplot(case_and_WGS_dates_t,  col=c("green","lightblue"),xlab="",
                ylab='Count', xaxt='n')
ind <- seq(2, nrow(case_and_WGS_dates), by = 7)
axis(1, ticks[ind], labels = format(dates[ind], "%b-%d"),las=2)
legend(x = "topright", text.font = 4, lty=c(1,1), lwd = 2,
       col= c("green","lightblue"),
       legend=c("Sequenced","Not Sequenced"))
title(xlab="Month-Day", mgp=c(4.5,0.75,0), family="Calibri Light",cex.lab=1.2)

## We need i) the time window(s) over which to estimate R and ii) information on the distribution of the serial interval.
## specify the mean and standard deviation of the serial interval. In that case an offset gamma distribution is used for the serial interval.

## Resources::
## Anne Cori: https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
## plot.estimate: https://rdrr.io/cran/EpiEstim/man/plot.estimate_R.html

## plot the incidence dist:
# library(incidence)
# plot(as.incidence(df$I, dates = (df$dates)))

## we only specify the mean and standard deviation of the serial interval. In that case an offset gamma distribution is used for the serial interval. In the following example, we use the mean and standard deviation of the serial interval for flu

## Info of SI is based on the following paper: https://jamanetwork.com/journals/jama/article-abstract/197988
## transmission study of Hand, Foot and mouth disease (EV-71), link below. They use an odd phrasing but they have a serial interval estimate “The median transmission interval was 3 days and the mean (SD) interval was 3.7 (2.6) days.”
config_low <- make_config(list(mean_si = 2, 
                               std_si = 2.6))
res_parametric_si_low <- estimate_R(incidence_df, 
                                    method="parametric_si",
                                    config = config_low
)
config_mid <- make_config(list(mean_si = 3.7, 
                               std_si = 2.6))
res_parametric_si_mid <- estimate_R(incidence_df, 
                                    method="parametric_si",
                                    config = config_mid
)
config_high <- make_config(list(mean_si = 7, 
                                std_si = 2.6))
res_parametric_si_high <- estimate_R(incidence_df, 
                                    method="parametric_si",
                                    config = config_high
)


## Index of time vector for R(t)
t_endvec<-res_parametric_si_mid[["R"]][["t_end"]]


## output df of R(t) with CI:
Re_from_incidence_df <- data.frame(dates=res_parametric_si_mid[["dates"]][t_endvec],
                                   low_R=res_parametric_si_low[["R"]][["Mean(R)"]],
                                   low_R_lower=res_parametric_si_low[["R"]][["Quantile.0.05(R)"]],
                                   low_R_upper=res_parametric_si_low[["R"]][["Quantile.0.975(R)"]],
                                   mid_R=res_parametric_si_mid[["R"]][["Mean(R)"]],
                                   mid_R_lower=res_parametric_si_mid[["R"]][["Quantile.0.05(R)"]],
                                   mid_R_upper=res_parametric_si_mid[["R"]][["Quantile.0.975(R)"]],
                                   high_R=res_parametric_si_high[["R"]][["Mean(R)"]],
                                   high_R_lower=res_parametric_si_high[["R"]][["Quantile.0.05(R)"]],
                                   high_R_upper=res_parametric_si_high[["R"]][["Quantile.0.975(R)"]]
)


## plot R(t) with the confidence interval
plot(x=Re_from_incidence_df$dates, y=Re_from_incidence_df$low_R, type='l',col='orange',
     xlim = as.Date(c(x_min_date, x_max_date)), ylim=c(0,6))
polygon(c(Re_from_incidence_df$dates, rev(Re_from_incidence_df$dates)),
        c(Re_from_incidence_df$low_R_upper ,rev(Re_from_incidence_df$low_R_lower)),
        col=adjustcolor("orange", alpha=0.5) , border=NA)
lines(x=Re_from_incidence_df$dates, y=Re_from_incidence_df$mid_R, type='l',col='red',
      xlim = as.Date(c(x_min_date, x_max_date)))
polygon(c(Re_from_incidence_df$dates, rev(Re_from_incidence_df$dates)),
        c(Re_from_incidence_df$mid_R_upper ,rev(Re_from_incidence_df$mid_R_lower)),
        col=adjustcolor('red', alpha=0.5) , border=NA)
lines(x=Re_from_incidence_df$dates, y=Re_from_incidence_df$high_R, type='l',col='brown',
      xlim = as.Date(c(x_min_date, x_max_date)))
polygon(c(Re_from_incidence_df$dates, rev(Re_from_incidence_df$dates)),
        c(Re_from_incidence_df$high_R_upper ,rev(Re_from_incidence_df$high_R_lower)),
        col=adjustcolor('brown', alpha=0.5) , border=NA)



#### Comparing Epiestim Re  with BDSky's ####
x_min_date <- as.POSIXct('2022-08-01 UTC')
x_max_date <- as.POSIXct('2022-11-01 UTC')
date_ticks = seq(x_min_date,x_max_date,by='months')
nf <- layout( matrix(c(1,2,3,4), ncol=1, nrow=4), heights = c(1,1,1,1.5))
y_max = max(c(ON_Re_gridded_hpd,Re_from_incidence_df$high_R_upper))
par(mar = c(0.75, 4, 0.25, 1))
Re_from_incidence_df$dates = as.POSIXct(Re_from_incidence_df$dates)
plotSkyline(ON_dates, ON_Re_gridded_hpd, type='step', xlab="", xaxt='n',
            ylab=expression("Ontario R"["t"]), lwd=2,
            col='blue', fill=adjustcolor("blue", alpha=0.25), xlim = c(x_min_date, x_max_date), ylims=c(y_min, y_max))
lines(x=Re_from_incidence_df$dates, y=Re_from_incidence_df$high_R, type='l',col='brown')
polygon(c(Re_from_incidence_df$dates, rev(Re_from_incidence_df$dates)),
        c(Re_from_incidence_df$high_R_upper ,rev(Re_from_incidence_df$high_R_lower)),
        col=adjustcolor('brown', alpha=0.25) , border=NA)
axis.POSIXct(side=1,at=date_ticks,
             labels = F, las=2)
abline(h = c(2,3,4,5,6,7), col = "grey", lty = "dotted")
abline(h = 1, col = "black", lty = "dashed")
legend(x = "topright", text.font = 4, lty=c(1,1,1,1), lwd = 2,
       col= c('blue',"brown", "red",'darkorange3'),
       legend=c("BDSky", "EpiEstim (High-SI)","EpiEstim (Mid-SI)", "EpiEstim (Low-SI)"))

par(mar = c(0.75, 4, 0, 1))
plotSkyline(ON_dates, ON_Re_gridded_hpd, type='step', xlab="", xaxt='n',
            ylab=expression("Ontario R"["t"]), lwd=2,
            col='blue', fill=adjustcolor("blue", alpha=0.25), xlim = c(x_min_date, x_max_date), ylims=c(y_min, y_max))
lines(x=Re_from_incidence_df$dates, y=Re_from_incidence_df$mid_R, type='l',col='red')
polygon(c(Re_from_incidence_df$dates, rev(Re_from_incidence_df$dates)),
        c(Re_from_incidence_df$mid_R_upper ,rev(Re_from_incidence_df$mid_R_lower)),
        col=adjustcolor('red', alpha=0.25) , border=NA)
axis.POSIXct(side=1,at=date_ticks,
             labels = F, las=2)
abline(h = c(2,3,4,5,6,7), col = "grey", lty = "dotted")
abline(h = 1, col = "black", lty = "dashed")

par(mar = c(0.75, 4, 0, 1))
plotSkyline(ON_dates, ON_Re_gridded_hpd, type='step', xlab="", xaxt='n',
            ylab=expression("Ontario R"["t"]), lwd=2,
            col='blue', fill=adjustcolor("blue", alpha=0.25), xlim = c(x_min_date, x_max_date), ylims=c(y_min, y_max))
lines(x=Re_from_incidence_df$dates, y=Re_from_incidence_df$low_R, type='l',col='darkorange3')
polygon(c(Re_from_incidence_df$dates, rev(Re_from_incidence_df$dates)),
        c(Re_from_incidence_df$low_R_upper ,rev(Re_from_incidence_df$low_R_lower)),
        col=adjustcolor('darkorange3', alpha=0.25) , border=NA)
abline(h = c(2,3,4,5,6,7), col = "grey", lty = "dotted")
abline(h = 1, col = "black", lty = "dashed")
axis.POSIXct(side=1,at=date_ticks,
             labels = F, las=2)

par(mar = c(5.25, 4, 0, 1))

case_and_WGS_dates2 <-(case_and_WGS_dates
                       %>% complete(dates = seq.Date(as.Date(x_min_date), as.Date(x_max_date), by="day"))
                       %>% replace_na(list(I=0,Sequenced=0,Not_Sequenced=0))
)
case_and_WGS_dates2_t = t(case_and_WGS_dates2[,-1:-2])
colnames(case_and_WGS_dates2_t) <- case_and_WGS_dates2$dates
ticks =barplot(case_and_WGS_dates2_t,  col=c("green","lightblue"),xlab="",
               ylab='Count', xaxt='n')
abline(h = seq(0,max(case_and_WGS_dates2$I),by=2), col = "grey", lty = "dotted")
ind <- match(as.Date(date_ticks),case_and_WGS_dates2$dates)
axis(1, ticks[ind], labels = format(case_and_WGS_dates2$dates[ind], "%b-%d"),las=1)
legend(x = "topright", text.font = 4, lty=c(1,1), lwd = 2,
       col= c("green","lightblue"),
       legend=c("Sequenced","Not Sequenced"))
title(xlab="Month-Day", mgp=c(2.5,0.75,0), family="Calibri Light",cex.lab=1.2)




