library(lubridate)
library(copula)
library(QRM)

# load in the data
data = read.csv('testingSet.csv',head=TRUE)

# convert on chock times to the POSIXct format 
on_chock = as.POSIXct(strptime(data$on_chocks_time, "%Y-%m-%d %H:%M:%S"), tz = "UTC")
on_chock_date = date(on_chock)
# find the days in the test set
on_chock_date_unique = unique(on_chock_date) 
on_chock=c()

# set up the pinball loss function
pinball <- function(pred, act, q) {
  I = 1*(act>pred)
  l <- (act - pred)*q*I + (pred - act)*(1-q)*(1-I)
  return(mean(l))
}

# load in the file contains the segment each passenger in the test set belongs to
leaf_no = read.csv('leaf-testing.csv',head=FALSE)

# load in the parameters of the gamma distributions. 
coefs = read.csv('coef.csv',head=TRUE)
gumbel_coef = matrix(NA,dim(data)[1],2)
for (i in 1: dim(data)[1]) {
  gumbel_coef[i,] = as.numeric(coefs[which(coefs$leaf == leaf_no[i,2]),c(1,3)])
}

# simulate with Gaussian copula
# we add correlations to the connection times of passengers arriving on the same day same flight
# we use Gaussian copula with correlation coefficient equals 0.5 to add the correlation
# the copula parameter 0.5 is obtained by cross validation
simulation_result = matrix(NA,dim(data)[1],500)
set.seed(201)
for (i in 1:length(on_chock_date_unique)) { # for each day
  day = on_chock_date_unique[i]
  this_day_pax = which(on_chock_date == day) # find passengers arriving on day i
  flight_no = unique(data$ib_flight_no[this_day_pax]) # find the flights arriving at the airport on day i
  n_flight = length(flight_no)
  for (j in 1:n_flight) { # for each flight
    this_flight_pax = this_day_pax[which(data[this_day_pax,9] == flight_no[j])] # find the passengers arriving on flight j
    if (length(this_flight_pax) > 1) {
      u = rCopula(500, normalCopula(0.5, dim = length(this_flight_pax))) # sample from Gaussian copula
      for (l in 1:500){
        simulation_result[this_flight_pax,l] = qgamma (u[l,], shape = gumbel_coef[this_flight_pax,1], scale = gumbel_coef[this_flight_pax,2]) # sample from the marginal gamma distribution
      }
    } else {
      simulation_result[this_flight_pax,] = rgamma(500, shape = gumbel_coef[this_flight_pax,1], scale = gumbel_coef[this_flight_pax,2])
    }
  }
}

# code for running independent simulation without Copulas
# simulation_result = matrix(NA,dim(data)[1],500)
# set.seed(201)
# for (i in 1:dim(simulation_result)[1]){
#   simulation_result[i,] = rgamma(500, shape = gumbel_coef[i,1], scale = gumbel_coef[i,2])
# }

# calculate passengers arriving time as on chock time + connection time
simulation_result2 = matrix(NA,dim(data),500)
on_chocks_time = strptime(data$on_chocks_time, "%Y-%m-%d %H:%M:%S", tz = "UTC")
for (i in 1:dim(data)[1]){
  a = rep(on_chocks_time[i],500)
  a$min = a$min + round(simulation_result[i,])
  simulation_result2[i,] = as.character(a)
}
actual = data$Delta

# pinball losses of the predictions for passengers' connection times in the test set.
connectionTimes05 = matrix(NA,dim(data),1)
connectionTimes25 = matrix(NA,dim(data),1)
connectionTimes50 = matrix(NA,dim(data),1)
connectionTimes75 = matrix(NA,dim(data),1)
connectionTimes95 = matrix(NA,dim(data),1)

for (i in 1:dim(data)[1]){
  connectionTimes25[i] = qgamma(0.25, shape = gumbel_coef[i,1], scale = gumbel_coef[i,2])
  connectionTimes75[i] = qgamma(0.75, shape = gumbel_coef[i,1], scale = gumbel_coef[i,2])
  connectionTimes05[i] = qgamma(0.05, shape = gumbel_coef[i,1], scale = gumbel_coef[i,2])
  connectionTimes95[i] = qgamma(0.95, shape = gumbel_coef[i,1], scale = gumbel_coef[i,2])
  connectionTimes50[i] = qgamma(0.50, shape = gumbel_coef[i,1], scale = gumbel_coef[i,2])
}
actual = data$Delta

# Calculate the hit rate of the 50% and 90% central prediction intervals
s = 0
for (i in 1:dim(data)[1]) {
  if (connectionTimes25[i] <= actual[i] & connectionTimes75[i] >= actual[i]) s=s+1
}
s/dim(data)[1]
s = 0
for (i in 1:dim(data)[1]) {
  if (connectionTimes05[i] <= actual[i] & connectionTimes95[i] >= actual[i]) s=s+1
}
s/dim(data)[1]

# pinball losses of the 0.05, 0.25, 0.5, 0.75 and 0.95 quantiles 
pinball(connectionTimes05, actual, 0.05)
pinball(connectionTimes25, actual, 0.25)
pinball(connectionTimes50, actual, 0.50)
pinball(connectionTimes75, actual, 0.75)
pinball(connectionTimes95, actual, 0.95)

# average pinball loss of the five quantiles
mean(c(pinball(connectionTimes05, actual, 0.05),
       pinball(connectionTimes25, actual, 0.25),
       pinball(connectionTimes50, actual, 0.50),
       pinball(connectionTimes75, actual, 0.75),
       pinball(connectionTimes95, actual, 0.95) ))

# calculate the MAE of the point forecasts (the median)
mean(abs(connectionTimes50- actual))

# free some variables to save some space
coefs=c()
gumbel_coef=c()
indimedian=c()
indip25=c()
indip75=c()
indip05=c()
indip95=c()
leaf_no=c()
simulation_result=c()
a=c()

# calculate the actual passenger flows
n1 <- length(on_chock_date_unique)
df = rep(ymd_hms(paste(on_chock_date_unique[1], '04:30:00'), tz = "UTC"),79*n1)
for (i in 1:n1) {
  df[((i-1)*79+1) : (i*79)] <- seq(ymd_hms(paste(on_chock_date_unique[i], '04:30:00'), tz = "UTC"), by = '15 min',length.out=(79))
}
df = as.POSIXct(strptime(df, "%Y-%m-%d %H:%M:%S"), tz = "UTC")
confTime = strptime(data$local_conform_time, "%Y-%m-%d %H:%M:%S", tz = "UTC")
confTime1 = with(confTime, as.POSIXct(ceiling(as.numeric(confTime)/(15*60)) * (15*60), origin = "1970-01-01"))
countPax1 = rep(NA,78*length(on_chock_date_unique))
for (i in 1:length(on_chock_date_unique)) {
  this_date <- on_chock_date_unique[i]
  n <- which(on_chock_date==this_date)
  this_confTime <- confTime1[n]
  this_confTime <- c(as.POSIXct(paste(this_date,"04:45:00")), this_confTime, as.POSIXct(paste(this_date+1,"00:00:00")))
  countPax = data.frame(table(cut(this_confTime, breaks = "15 mins")))
  countPax$Freq[1] = countPax$Freq[1] - 1
  countPax$Freq[78] = countPax$Freq[78] - 1
  countPax1[((i-1)*78+1) : (i*78)] = countPax$Freq
}

# calculate the predicted passenger flows
predictCountPax <- matrix(NA,78*length(on_chock_date_unique),500)
countPax <- rep(NA,(length(df)-1))
for (j in 1:500) {
  confTime = strptime(simulation_result2[,j], "%Y-%m-%d %H:%M:%S", tz = "UTC")
  confTime1 = with(confTime, as.POSIXct(ceiling(as.numeric(confTime)/(15*60)) * (15*60), origin = "1970-01-01"))
  for (i in 1:length(on_chock_date_unique)) {
    this_date <- on_chock_date_unique[i]
    n <- which(on_chock_date==this_date)
    this_confTime <- confTime1[n]
    this_confTime <- c(as.POSIXct(paste(this_date,"04:45:00")), this_confTime, as.POSIXct(paste(this_date+1,"00:00:00")))
    countPax = data.frame(table(cut(this_confTime, breaks = "15 mins")))
    countPax_Freq = countPax$Freq[1:78]
    countPax_Freq[1] = countPax_Freq[1] - 1
    countPax_Freq[78] = countPax_Freq[78] - 1
    predictCountPax[((i-1)*78+1) : (i*78),j] = countPax_Freq
  }
}

# Calculate the quantiles of the predicted passenger flows
p50 = matrix(NA, dim(predictCountPax)[1],1)
p25 = matrix(NA, dim(predictCountPax)[1],1)
p75 = matrix(NA, dim(predictCountPax)[1],1)
p05 = matrix(NA, dim(predictCountPax)[1],1)
p95 = matrix(NA, dim(predictCountPax)[1],1)
for (i in 1:dim(predictCountPax)[1]) {
  p50[i] = median(predictCountPax[i,])
  p25[i] = quantile(predictCountPax[i,],0.25)
  p75[i] = quantile(predictCountPax[i,],0.75)
  p05[i] = quantile(predictCountPax[i,],0.05)
  p95[i] = quantile(predictCountPax[i,],0.95)
}

# find out the peak hours: 
# If the sum of the actual number of passengers arriving in a given 15-minute interval and those of the previous three 15 minutes is above 200, this 15 minutes interval can be regarded as part of a peak hour
peakornot <- rep(NA, length(countPax1))
for (i in 3:length(peakornot)) {
  if ( i%%78 == 0) next
  if (sum(countPax1[(i-2):i]) >= 200) peakornot[i] <- 1
}
n1 = which(!is.na(peakornot)) 
p05_peak = p05[n1]
p25_peak = p25[n1]
p75_peak = p75[n1]
p95_peak = p95[n1]
p50_peak = p50[n1]
countPax2 = countPax1[n1]

# hit rate of the central 50% and 90% prediction intervals
s = 0
for (i in 1:length(countPax2)) {
  if (p25_peak[i] <= countPax2[i] & p75_peak[i] >= countPax2[i]) s=s+1
}  
s/length(countPax2)  
s = 0
for (i in 1:length(countPax2)) {
  if (p05_peak[i] <= countPax2[i] & p95_peak[i] >= countPax2[i]) s=s+1
}  
s/length(countPax2)  

# calculate the pinball losses of the 0.05, 0.25, 0.5, 0.75 and 0.95 quantiles
pinball(p05_peak, countPax2, 0.05) 
pinball(p25_peak, countPax2, 0.25) 
pinball(p50_peak, countPax2, 0.50) 
pinball(p75_peak, countPax2, 0.75) 
pinball(p95_peak, countPax2, 0.95) 

# calculate the average pinball loss of the five quantiles
mean(c(pinball(p05_peak, countPax2, 0.05),
       pinball(p25_peak, countPax2, 0.25),
       pinball(p50_peak, countPax2, 0.50),
       pinball(p75_peak, countPax2, 0.75), 
       pinball(p95_peak, countPax2, 0.95)))

# calculate the MAE of the point forecasts (or the median)
mean(abs(p50_peak - countPax2))

