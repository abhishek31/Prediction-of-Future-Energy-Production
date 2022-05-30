library(readxl)
library(astsa)
library(tseries)
library(xts)
library(forecast)

data<-read.csv(file="Electric_Production.csv", header=TRUE,  sep=",")
data_full=ts(data$IPG2211A2N,start=c(1985,1), frequency = 12)

############################################################################
## Data Exploration 
############################################################################
plot(data_full)
acf2(data_full)
plot(data_full, xlab = "Time", ylab = "IP index")
title('IP index vs Time relation')

## ACF proves that model is not stationary
## The data is clearly showing an increasing trend over the years.
## Also the variance of the data is not constant.

## Decompose dataset into various components to see trend, seasonality and the random
## components
plot(decompose(data_full))

## Hence differencing and logarithmic transformation is required to convert it 
## to a stationary time series.

#There is trend and increasing variance so we take diff and log of the data
plot(data_full)
plot(diff(data_full))
plot(diff(log(data_full)))

# we stick to taking the difference of the log of the data.

acf2(diff(log(data_full)))

## Augmented DF test to check the stationarity of data
adf.test(diff(log(data_full)))
## ADF test indicates the data is stationary now

############################################################################
## Splitting the data into train and test sets in 80:20 ratio 
############################################################################

train_data<-data[1:324,]
test_data<-data[325:397,]

## Build Y variable using train dataset

y <- ts((train_data$IPG2211A2N), start=c(1985,1),frequency = 12)
plot(y)

## Build Y1 variable using test dataset

y1 <- ts((test_data$IPG2211A2N), start=c(2012,1),frequency = 12)
plot(y1)

############################################################################
## Sarima models
############################################################################

#one step ahead prediction for auto-arima
library(forecast)
library(forecast)
#one step ahead prediction for auto arima
model <- auto.arima(y,stationary=F, seasonal = TRUE)
model.fit <- Arima(y1, model=model)
model$arma
onestep.for <- fitted(model.fit)
onestep.for
## Prediction visualization
autoplot(forecast::forecast(model, h = 73),xlab = "Time", ylab = "IP index")

#MAPE - Mean absolute percentage error.
sum(abs((y1-onestep.for)/y1))/73
#2.17%

#fitting the best  manually selected Sarima models 

##Model 1
sarima(log(data_full), 2, 1, 3, 2, 1, 1, S=12)

model<-Arima(log(y),order=c(2,1,3),seasonal=list(order=c(2,1,1),period=12),lambda=0)
model.fit <- Arima((y1), model=model)
onestep.for <- fitted(model.fit)

#MAPE - Mean absolute percentage error.
sum(abs((y1-onestep.for)/y1))/73
#1.97%

##Model 2
sarima(log(data_full), 2, 1, 3, 1, 1, 1, S=12)

model<-Arima(log(y),order=c(2,1,3),seasonal=list(order=c(1,1,1),period=12),lambda=0)
model.fit <- Arima((y1), model=model)
onestep.for <- fitted(model.fit)

#MAPE - Mean absolute percentage error.
sum(abs((y1-onestep.for)/y1))/73
#2.05%

##Model 3
sarima(log(data_full), 3, 1, 3, 2, 1, 1, S=12)

model<-Arima(log(y),order=c(3,1,3),seasonal=list(order=c(2,1,1),period=12),lambda=0)
model.fit <- Arima((y1), model=model)
onestep.for <- fitted(model.fit)

#MAPE - Mean absolute percentage error.
sum(abs((y1-onestep.for)/y1))/73
#1.98%

##Model 4
sarima(log(data_full), 3, 1, 1, 4, 1, 0, S=12)

model<-Arima(log(y),order=c(3,1,1),seasonal=list(order=c(4,1,0),period=12),lambda=0)
model.fit <- Arima((y1), model=model)
onestep.for <- fitted(model.fit)

#MAPE - Mean absolute percentage error.
sum(abs((y1-onestep.for)/y1))/73
#2.00%

### Comparison of  the coefficients of the Sarima models
## The coefficients of the Sarima model (2, 1, 3, 2, 1, 1, S=12) 
## with the lowest MAPE have the lowest standard errors among other models and hence, 
##this model is selected as the best performing model.

##############################################################################
## Prediction using Prophet 
##############################################################################

library(zoo)
library(prophet)
library(ggplot2)

## Build Dataframe from training data
df<-data.frame(date=as.Date(as.yearmon(time(y))),Y=as.matrix(y))

## Change name of columns as prophet takes input ds as the name of input columns 
## and y as variable, we are predicting.
names(df) = c("ds", "y")
df
qplot(ds,y,data=df)

## Model building
m = prophet(df)

## Forecasting
future = make_future_dataframe(m, periods = 73, include_history = TRUE)
forecast = predict(m, future)
forecast$yhat
#MAPE - Mean absolute percentage error.
MAPE_prophet<-sum(abs((test_data$IPG2211A2N-forecast$yhat[325:397])/test_data$IPG2211A2N))/73
MAPE_prophet
## 24.10%
## Visualize forecast
plot(m, forecast)
prophet_plot_components(m, forecast)


##############################################################################
## Prediction using BSTS 
##############################################################################

library(bsts)
library(dplyr)
library(lubridate)

### Run the bsts model
ss <- AddLocalLinearTrend(list(), y)
ss <- AddSeasonal(ss, y, nseasons = 12)
bsts.model <- bsts(y, state.specification = ss, niter = 500, ping=0)
### Get a suggested number of burn-ins
burn <- SuggestBurn(0.1, bsts.model)

### Predict
p <- predict.bsts(bsts.model, horizon = 73, burn = burn, quantiles = c(.025, .975))

### Actual versus predicted
d2 <- data.frame(
  # fitted values and predictions
  c(as.numeric(-colMeans(bsts.model$one.step.prediction.errors[-(1:burn),])+y),  
    as.numeric(p$mean)),
  # actual data and dates 
  as.numeric(data_full),
  as.Date(time(data_full)))
names(d2) <- c("Fitted", "Actual", "Date")

### MAPE (mean absolute percentage error)
MAPE_BSTS <- sum(abs((d2$Actual[325:397]-d2$Fitted[325:397])/d2$Actual))/73
MAPE_BSTS
## 19.15%

### 95% forecast credible interval
posterior.interval <- cbind.data.frame(
  as.numeric(p$interval[1,]),
  as.numeric(p$interval[2,]), 
  subset(d2, year(Date)>2011)$Date)
names(posterior.interval) <- c("LL", "UL", "Date")

### Join intervals to the forecast
d3 <- left_join(d2, posterior.interval, by="Date")

### Plot actual versus predicted with credible intervals for the holdout period
ggplot(data=d3, aes(x=Date)) +
  geom_line(aes(y=Actual, colour = "Actual"), size=1.2) +
  geom_line(aes(y=Fitted, colour = "Fitted"), size=1.2, linetype=2) +
  theme_bw() + theme(legend.title = element_blank()) + ylab("") + xlab("") +
  geom_vline(xintercept=as.numeric(as.Date("2011-12-01")), linetype=2) + 
  geom_ribbon(aes(ymin=LL, ymax=UL), fill="grey", alpha=0.5) +
  ggtitle(paste0("BSTS -- Holdout MAPE = ", round(100*MAPE_BSTS,2), "%")) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
############################################################################
##                                END
############################################################################

df_viz <- data.frame(Models=c("Auto - ARIMA", "SARIMA(2,1,3)(2,1,1)", "SARIMA(2,1,3)(1,1,1)","SARIMA(3,1,3)(2,1,1)","SARIMA(3,1,1)(4,1,0)","Prophet","LSTM","BSTS"),
                MAPE=c(2.17, 1.97, 2.05,1.98,2.00,24.10,0.05,20.13))
library(ggplot2)
# Basic barplot
p<-ggplot(data=df_viz, aes(x=Models, y=MAPE)) +
  geom_bar(stat="identity", fill="steelblue",width=0.25)+
  theme_minimal(base_size = 10)
p
