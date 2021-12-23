library(readr)
library(zoo)
library(dlm)
library(TSA)
library(lubridate)
library(tseries)
library(xts)
library(forecast)
library(MLmetrics)
library(gtools)
library(caTools)
library(gplots)
library(fma)
library(expsmooth)
library(fpp)

s <- read_csv("C:/Users/Tilakteja/Downloads/monthly-sunspots.txt")
View(s)

#sunspot plot
plot(s$Sunspots,type="o")

#converting data into time series form
s.ts <- matrix(s$Sunspots, nrow =2820 , ncol = 1)
s.ts <- as.vector(t(s.ts))
s.ts <- ts(s.ts,start=c(1749,01), frequency=12)

su.ts <- matrix(s$Sunspots, nrow =2820 , ncol = 1)
su.ts <- as.vector(t(s.ts))
su.ts <- ts(s.ts,start=c(1,1),end = c(20,12), frequency=12)
su.ts
length(su.ts)

suu.ts <- matrix(s$Sunspots, nrow =2820 , ncol = 1)
suu.ts <- as.vector(t(s.ts))
suu.ts <- ts(s.ts,start=c(1,1),end = c(10,12), frequency=12)
suu.ts

plot(ts(as.vector(su.ts)),type='o',ylab='Sunspots')
plot(decompose(su.ts))   #seasonality                       
plot.ts(s)

#eda
library(forecast)
ggseasonplot(suu.ts)

plot(su.ts, main = "", ylab = "sunspot", xlab = "First 20 data ", 
     bty = "l", ylim = c(0, 150))

#ACF PACF plots
tsdisplay(s.ts,main = "ACF and PACF of Sunspots")

#Aggregate() Function in R Splits the data into subsets,
# computes summary statistics for each subsets and returns the result in a group by form.
aggregate(su.ts, 4, median)
AG <- aggregate(su.ts, 4, mean)
autoplot(AG)
summary(su.ts)
aggregate(su.ts, 4, mean)
AGG <- aggregate(su.ts, 4, mean)
autoplot(AGG)


#timeseries
library(boot)
library(forecast)
#Variance stabilizing transform
sun1 <- 2*(sqrt(sunspot.year+1)-1)
sun1
autoplot(sun1,xlab = "Year",ylab = "Sample",main = "Time Series")

#Auto Arima for original data
auto.arima(a.ts)

#making data into stationary form
library(tseries)
library(forecast)
autoplot(s.ts)

#unit ratio test
adf.test(s.ts)
#Kwiatkowski-Phillips-Schmidt-Shin (KPSS) tests are used for testing a null hypothesis that an observable
# time series is stationary around a deterministic trend against the alternative of a unit root.
kp=kpss.test(s.ts)
kp
s.ts_d <- diff(s.ts, differences = 1)
adf.test(s.ts_d, k=12)  #d=1
kpss=kpss.test(s.ts_d)
kpss

autoplot(s.ts_d)  # now it is stationary


pacf(s.ts_d) #p=2
acf(s.ts_d)  #q=1
#Arima of sunspots
tsmod <- Arima(y = s.ts_d, order = c(2,0,1))
tsmod
auto.arima(s.ts)
predict(tsmod,n.ahead = 120,newxreg = NULL,se.fit = TRUE)
forecast(tsmod,h=120)
s.ts_d
#div data for testing
set.seed(3000)
dt=sort(sample(nrow(s),nrow(s)*.7))
train<-s[dt,]
test<-s[-dt,]
train
test
tsdisplay(a.ts)
#converting train data into timeseries 
a.ts <- matrix(train$Sunspots, nrow =1973 , ncol = 1)
a.ts <- as.vector(t(a.ts))
a.ts <- ts(a.ts,start=c(1819,08), frequency=12)
a.ts
autoplot(a.ts)
length(a.ts)
#Auto Arima for train data
tsmo_train <- auto.arima(a.ts)
tsmo_train
autoplot(a.ts)
adf.test(a.ts)
KPSS=kpss.test(a.ts)
KPSS

pacf(a.ts) #p=5
acf(a.ts)  #q=2

tsmod_train <- Arima(y = a.ts, order = c(5,0,2))
tsmod_train
pre <- predict(tsmo_train,n.ahead = 120,newxreg = NULL,se.fit = TRUE)
pre
pree <- round(pre$pred,1)
pree
forecast(tsmod_train,h=120)

# converting test data into time series form
b.ts <- matrix(test$Sunspots, nrow =847 , ncol = 1)
b.ts <- as.vector(t(b.ts))
b.ts <- ts(b.ts,start=c(1749,01), frequency=12)
autoplot(b.ts)
adf.test(b.ts)

mape = (1/240) * sum(abs(s.ts-pre$pred) /s.ts) * 100

mean(abs((s.ts-pre)/s.ts)) * 100
MAPE(y_pred = exp(pre$pred.values) , y_true = s.ts)

aa <- a.ts[1:120]
aa

length(pree)
length(aa)

MAPE(aa,pree)
tsdisplay(a.ts,main = "ACF,PACF,OBSERVED PLOTS")
