#### Producing a Combined Re figure for the JH and PHO dataset
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

ON_youngest_tip = 2022.7615296803654
MD_youngest_tip = 2022.6957762557079
all_youngest_tip = min(c(ON_youngest_tip, MD_youngest_tip))


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

#### Load MD-2022 logfiles ####
fname <- "beast_xmls/MD-2022/MD_combined_run_4to6.log"  
MD_trace  <- beastio::readLog(fname, burnin=0.1)

#### Get origin esitmates and enter tip dates ####
MD_origin_HPD <-  beastio::getHPDMedian(MD_trace[, "origin_BDSKY_Serial"])
MD_origin_med <- MD_origin_HPD[2]

#### Get esitmates of Re for relevent timeframes ####
MD_Re_sky    <- beastio::getLogFileSubset(MD_trace, "reproductiveNumber_BDSKY_Serial")
MD_Re_hpd    <- t(beastio::getHPDMedian(MD_Re_sky))
MD_gridTimes  <- seq(0, MD_origin_med, length.out=100)  
MD_Re_gridded <- mcmc(bdskytools::gridSkyline(MD_Re_sky, MD_trace[, "origin_BDSKY_Serial"], MD_gridTimes))
MD_Re_gridded_hpd <- t(getHPDMedian(MD_Re_gridded))
MD_times <- MD_youngest_tip - MD_gridTimes

#### Figure
# Convert to dates Month Year
ON_dates = lubridate::date_decimal(ON_times)
MD_dates = lubridate::date_decimal(MD_times)

# End of Social Contact restriction Ontario and Mayland
MD_NPI_lift_begin = as.POSIXct("2021-03-01", format="%Y-%m-%d")
MD_NPI_lift_end = as.POSIXct("2021-08-13", format="%Y-%m-%d")
ON_NPI_lift_begin = as.POSIXct("2022-01-31", format="%Y-%m-%d")
ON_NPI_lift_end = as.POSIXct("2022-03-14", format="%Y-%m-%d")

# Data for origin KDE plots
ON_origin_all = ON_youngest_tip - ON_trace[, "origin_BDSKY_Serial"]
MD_origin_all = MD_youngest_tip - MD_trace[, "origin_BDSKY_Serial"]


#### Plotting figure Epidimiological estimates
x_min_date = floor_date(lubridate::date_decimal(min(MD_origin_all)),'month') 
x_max_date = max(c(ON_dates, MD_dates))
x_max_date = ceiling_date(x_max_date ,'month')
date_ticks = seq(as.POSIXct(x_min_date),as.POSIXct(x_max_date),by='months')
y_min = 0
y_max = max(c(ON_Re_gridded_hpd, MD_Re_gridded_hpd))
x_min_decimal = lubridate::decimal_date(x_min_date)
x_max_decimal = lubridate::decimal_date(x_max_date)

nf <- layout( matrix(c(1,2,3), ncol=1, nrow=3), heights = c(1,1,1.5))
#### Density plot of origin dates ###
par(mar = c(0.75, 4, 0, 1))
den_ON_origin =density(ON_origin_all)
plot(den_ON_origin,xlim = c(x_min_decimal, x_max_decimal), col='blue', lwd=2,
     ylab='Epidemic Origin',xlab="",main="", xaxt='n')
polygon(den_ON_origin, col=adjustcolor("blue", alpha=0.25))
den_MD_origin =density(MD_origin_all)
lines(den_MD_origin, col="orange", lwd=2)
polygon(den_MD_origin, col=adjustcolor("orange", alpha=0.25))
axis(side=1,at=lubridate::decimal_date(date_ticks), las=2, labels = F)
rect(xleft =lubridate::decimal_date(MD_NPI_lift_begin),
     xright =lubridate::decimal_date(MD_NPI_lift_end),
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("grey", alpha=0.25))
rect(xleft = lubridate::decimal_date(ON_NPI_lift_begin), 
     xright = lubridate::decimal_date(ON_NPI_lift_end), 
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("grey", alpha=0.25))
abline(h = c(2,4,6,8), col = "grey", lty = "dotted")
abline(v = c(min(ON_times),  min(MD_times)), 
       col= c('blue', "orange"), lwd = 2, lty = "dotdash")
legend(x = "topright", text.font = 4, lty=c(1,1,4), lwd = 2,
       col= c('blue', "orange",'black'),
       legend=c("ON-2022", "MD-2022", 'Median Origin'))

par(mar = c(0.75, 4, 0, 1))
plotSkyline(ON_dates, ON_Re_gridded_hpd, type='smooth', xlab="", xaxt='n',   ylab=expression("Ontario R"["t"]), lwd=2,
            col='blue', fill=adjustcolor("blue", alpha=0.25), xlim = c(x_min_date, x_max_date), ylims=c(y_min, y_max))
abline(h = c(2,3,4,5,6), col = "grey", lty = "dotted")
abline(h = 1, col = "black", lty = "dashed")
abline(v = c(min(ON_dates), min(MD_dates)), 
       col= c('blue',"orange"), lwd = 2, lty = "dotdash")
axis.POSIXct(side=1,at=date_ticks,
             format = '%y-%b', las=2, labels = F)
rect(xleft = ON_NPI_lift_begin, xright = ON_NPI_lift_end, ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("grey", alpha=0.25))
#title(xlab="Year-Month", mgp=c(4,0.75,0), family="Calibri Light",cex.lab=1.2)

par(mar = c(5.25, 4, 0, 1))
plotSkyline(MD_dates, MD_Re_gridded_hpd, type='smooth', xlab="", xaxt='n',   ylab=expression("Maryland R"["t"]), lwd=2,
            col="orange", fill=adjustcolor("orange", alpha=0.25), xlim = c(x_min_date, x_max_date), ylims=c(y_min, y_max))
abline(h = c(2,3,4,5,6), col = "grey", lty = "dotted")
abline(h = 1, col = "black", lty = "dashed")
abline(v = c(min(ON_dates), min(MD_dates)), 
       col= c('blue', "orange"), lwd = 2, lty = "dotdash")
rect(xleft = MD_NPI_lift_begin, xright = MD_NPI_lift_end, ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("grey", alpha=0.25))
axis.POSIXct(side=1,at=date_ticks,
             format = '%y-%b', las=2)
title(xlab="Year-Month", mgp=c(4,0.75,0), family="Calibri Light",cex.lab=1.2)

#### Max Re values

index = which.max(ON_Re_gridded_hpd[2,])
ON_Re_max = ON_Re_gridded_hpd[,index]
index = which.max(MD_Re_gridded_hpd[2,])
MD_Re_max = MD_Re_gridded_hpd[,index]


#### Origin HPD intervals as dates
ON_origin_HPD_date <- lubridate::date_decimal(ON_youngest_tip - ON_origin_HPD)
MD_origin_HPD_date <- lubridate::date_decimal(MD_youngest_tip - MD_origin_HPD)

#### TMRCA HPD intervals as dates
ON_TMRCA_HPD = beastio::getHPDMedian(ON_trace[, 'Tree.height'])
MD_TMRCA_HPD = beastio::getHPDMedian(MD_trace[, 'Tree.height'])

#### Plotting figure Evolutionary estimates
# Data for origin KDE plots
ON_TMRCA_all = ON_youngest_tip - ON_trace[, 'Tree.height']
MD_TMRCA_all = MD_youngest_tip - MD_trace[, 'Tree.height']


nf <- layout( matrix(c(1), ncol=1, nrow=1))
den_ON_TMRCA =density(ON_TMRCA_all)
plot(den_ON_TMRCA,xlim = c(x_min_decimal, x_max_decimal), col='blue', lwd=2,
     ylab='TMRCA',xlab="",main="", xaxt='n')
polygon(den_ON_TMRCA, col=adjustcolor("blue", alpha=0.25))
den_MD_TMRCA =density(MD_TMRCA_all)
lines(den_MD_TMRCA, col="orange", lwd=2)
polygon(den_MD_TMRCA, col=adjustcolor("orange", alpha=0.25))
axis(side=1,at=lubridate::decimal_date(date_ticks), las=2, labels = zoo::as.yearmon(date_ticks))
rect(xleft =lubridate::decimal_date(MD_NPI_lift_begin),
     xright =lubridate::decimal_date(MD_NPI_lift_end),
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("grey", alpha=0.25))
rect(xleft = lubridate::decimal_date(ON_NPI_lift_begin), 
     xright = lubridate::decimal_date(ON_NPI_lift_end), 
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA, col = adjustcolor("grey", alpha=0.25))
abline(h = c(1,2,3,4,5,6), col = "grey", lty = "dotted")
abline(v = c(median(ON_TMRCA_all), median(MD_TMRCA_all)), 
       col= c('blue', "orange"), lwd = 2, lty = "dotdash")
legend(x = "topleft", text.font = 4, lty=c(1,1,4), lwd = 2,
       col= c('blue',"orange",'black'),
       legend=c("ON-2022", "MD-2022", 'Median TMRCA'))



ON_TMRCA_HPD_date <- lubridate::date_decimal(ON_youngest_tip - ON_TMRCA_HPD)
MD_TMRCA_HPD_date <- lubridate::date_decimal(MD_youngest_tip - MD_TMRCA_HPD)



