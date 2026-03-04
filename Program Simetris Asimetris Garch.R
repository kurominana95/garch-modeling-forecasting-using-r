# Program ARIMA-GARCH dan ARIMA-GJR GARCH Skripsi

# Identitas:
# Nama : Najwa Khoir Aldawiyah
# NIM : 188221045
# Program Studi : S-1 Statistika
# Angkatan : 2022

# Import Library
library(rugarch)
library(aTSA)
library(FinTS) 
library(lmtest)
library(forecast)
library(TSA)
library(tseries)
library(xts)
library(openxlsx)
library(ggplot2)
library(MLmetrics)

#------ Set Working Directory -----
setwd("C:/Users/NAJWA KHOIR ALDWYAH/OneDrive/Dokumen/Statistika/Skripsi")

#------ Import Data ------
stock <- read.xlsx("Data Skripsi Garch.xlsx", sheet=1)

#----- Preprocessing Data -----
stock$Date <- as.Date(as.numeric(stock$Date), origin = "1899-12-30")
jii <- xts(stock$JII, order.by=stock$Date)
length(jii)

##---- Split Train-Test ----
nTrain <- 0.9*length(jii);round(nTrain)
rJiiTrain <- jii[1:round(nTrain)]
rJiiTest <- jii[(round(nTrain)+1):length(jii)];length(rJiiTest)
View(rJiiTrain)

#---- Visualisasi ----
df_train <- data.frame(Date = as.Date(time(rJiiTrain),frac=1), value=as.numeric(rJiiTrain),set="Training")
df_test <- data.frame(Date = as.Date(time(rJiiTest),frac=1), value=as.numeric(rJiiTest),set="Testing")
df_all <- rbind(df_train,df_test)
ggplot(df_all, aes(x = Date, y = value, color = set)) +
  geom_line(size = 0.6) +
  scale_color_manual(values = c("Training" = "black", "Testing" = "blue")) +
  labs(title = "Weekly Closing Price of \n Jakarta Islamic Index",
       x = "Date", y = "Price") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",        
      size = 16
    ),
    legend.position = "bottom",   
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

#---- Stasioneritas ----
adf.test(returnJii)
adf.test(rJiiTrain)

#----- Box Jenkins ARIMA ----
##---- Identifikasi Model ----
acf(rJiiTrain)
pacf(rJiiTrain)
##---- Differencing ----
rJiiTrain.1 <- diff(rJiiTrain)
adf.test(rJiiTrain.1[-1,])
acf(rJiiTrain.1[-1,],main="ACF Differencing 1",lag.max = 15)
pacf(rJiiTrain.1[-1,],main="PACF Differencing 1",lag.max = 15)
eacf(rJiiTrain.1[-1,])
rJiiTrain.1 <- rJiiTrain.1[-1,]
rJiiTrain.2 <- diff(rJiiTrain.1)
acf(rJiiTrain.2[-1,],main="ACF Differencing 2",lag.max = 15)
pacf(rJiiTrain.2[-1],main="ACF Differencing 2",lag.max = 15)

##---- Modeling -----
modelarima <- arima(rJiiTrain, order=c(4,2,0), method="ML", include.mean=FALSE)
coeftest(modelarima)
AIC(modelarima)
MSE(fitted(modelarima),rJiiTrain)
RMSE(fitted(modelarima),rJiiTrain)
MAE(fitted(modelarima),rJiiTrain)
MAPE(fitted(modelarima),rJiiTrain)

##---- Uji Diagnostik ----
checkresiduals(modelarima$residuals)
print(Box.test(modelarima$residuals, lag=12, type="Ljung-Box"))
print(Box.test(modelarima$residuals, lag=24, type="Ljung-Box"))
print(Box.test(modelarima$residuals, lag=36, type="Ljung-Box"))
ks.test(modelarima$residuals,"pnorm",mean(modelarima$residuals),sd(modelarima$residuals))
arch.test(modelarima)

#----- GARCH Modeling -----
## ----- Model Simetri GARCH -----
save_garch <- function(X, p, q, type = "norm"){
  spc <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(p, q)),
                    mean.model = list(armaOrder = c(4, 0),
                                      include.mean = F),
                    distribution.model = type)
  mod <- ugarchfit(spec = spc, data = X)
  return(list(model=mod,spec=spc))
}
g10 <- save_garch(X = rJiiTrain.2[-1], 
                  p = 1, 
                  q = 1)
(g10$model)
infocriteria(g10$model)
signbias(g10$model)
ArchTest(residuals(g10$model,standardize=TRUE),lags=4)
ArchTest(residuals(g10$model,standardize=TRUE),lags=8)
ArchTest(residuals(g10$model,standardize=TRUE),lags=12)

## ----- Model Asimetri GARCH -----
save_agarch <- function(X, p, q, type = "norm"){
  spc <- ugarchspec(variance.model = list(model = "fGARCH",submodel="GJRGARCH",
                                          garchOrder = c(p, q)),
                    mean.model = list(armaOrder = c(4, 0),
                                      include.mean = F),
                    distribution.model = type)
  mod <- ugarchfit(spec = spc, data = X)
  return(list(model=mod,spec=spc))
}

g01 <- save_agarch(X = rJiiTrain.2[-1], 
                   p = 1, 
                   q = 1)
g01$model
infocriteria(g01$model)
signbias(g01$model)
ArchTest(residuals(g01$model,standardize=TRUE),lags=4)
ArchTest(residuals(g01$model,standardize=TRUE),lags=8)
ArchTest(residuals(g01$model,standardize=TRUE),lags=12)

# --- Fitting ----
## ----- Fit ARIMA AR(4) ------
modelarima.fix <- arima(rJiiTrain, order=c(4,2,0), method="ML", include.mean=FALSE)
fit.arima <- fitted(modelarima.fix)
isTRUE(length(fit.arima)=length(rJiiTrain))

## ----- Fit ARIMA GARCH (1,1) ------
coef(g10$model)
omega <- 43.33
alpha <- 0.544
beta  <- 0.251
residu <- residuals(modelarima.fix)
var_resid <- residu^2
var_garch <- rep(0, length(rJiiTrain))
var_garch[1] <- var(var_resid)  # bisa pakai variance awal
for(i in 2:length(rJiiTrain)){
  var_garch[i] <- omega + alpha * var_resid[i-1] + beta * var_garch[i-1]
}
std_residual <- (residu - mean(residu)) / sd(residu)
resid_model <- std_residual * sqrt(var_garch)
resid_model <- xts(resid_model,order.by=time(rJiiTrain))
arima_garch_fit <- as.numeric(fit.arima) + as.numeric(resid_model)

resi_arimag <- rJiiTrain - arima_garch_fit
mape_arimag <- mean(abs(resi_arimag/rJiiTrain));mape_arimag
MAPE(arima_garch_fit,rJiiTrain)
RMSE(arima_garch_fit,rJiiTrain)
MSE(arima_garch_fit,rJiiTrain)

df_garch <- data.frame(
  waktu = time(rJiiTrain),
  actual = as.numeric(rJiiTrain),
  fitted = as.numeric(arima_garch_fit)
)
write.xlsx(df_garch, "Actual vs Fit ARIMA-GARCH.xlsx")

ggplot(df_garch, aes(x = waktu)) +
  geom_line(aes(y = actual, color = "Actual"), size = 1) +
  geom_line(aes(y = fitted, color = "Fitted"), size = 0.8) +
  labs(title = "Aktual vs Fitted ARIMA-GARCH",
       x = "Waktu", y = "Nilai", color = "Series") +
  theme_minimal()

## ----- Fit ARIMA-GJR GARCH (1,1) -----
## ---- Fit ARIMA-GJR GARCH------
modelarima.fix1 <- arima(rJiiTrain, order=c(4,2,0), method="ML", include.mean=FALSE,fixed=c(-0.93497,-0.78055,-0.53590,-0.23175))
summary(modelarima.fix)
omega.gjr <- 28.853135
alpha.gjr <- 0.096802
beta.gjr <- 0.666569
eta.gjr <- 0.999972
residu <- residuals(modelarima.fix)
var_resid <- residu^2
var_gjr <- rep(NA,length(rJiiTrain))
var_gjr[1] <- 0
for(i in 2:length(rJiiTrain)){
  indicator <- ifelse(residu[i-1] < 0, 1, 0)
  var_gjr[i] <- omega.gjr + (alpha.gjr * var_resid[i-1]) + (eta.gjr * var_resid[i-1] * indicator) + (beta.gjr * var_gjr[i-1])
}
std_residual <- (residu - mean(residu)) / sd(residu)
resid_model_gjr <- std_residual * sqrt(var_gjr)
resid_model_gjr <- xts(resid_model_gjr,order.by=time(rJiiTrain))
arima_gjr_fit <- as.numeric(fit.arima) + as.numeric(resid_model_gjr)

resi_arimagjr <- rJiiTrain - arima_gjr_fit;resi_arimagjr
mape_arimagjr <- mean(abs(resi_arimagjr/rJiiTrain));mape_arimagjr
MAPE(arima_gjr_fit,rJiiTrain)
RMSE(arima_gjr_fit,rJiiTrain)
MAE(arima_gjr_fit,rJiiTrain)

df_gjr <- data.frame(
  waktu = time(rJiiTrain),
  actual = as.numeric(rJiiTrain),
  fitted = as.numeric(arima_gjr_fit)
)
write.xlsx(df_garch, "Actual vs Fit ARIMA-GJR GARCH.xlsx")

ggplot(df_gjr, aes(x = waktu)) +
  geom_line(aes(y = actual, color = "Actual"), size = 1) +
  geom_line(aes(y = fitted, color = "Fitted"), size = 1) +
  labs(title = "Aktual vs Fitted ARIMA-GJR GARCH",
       x = "Waktu", y = "Nilai", color = "Series") +
  theme_minimal()

# ----- Prediksi Out Sample -----
## ---- Prediksi ARIMA -----
fore.arima <- predict(modelarima, n.ahead=length(rJiiTest))$pred
resi.fore.arima <- as.numeric(rJiiTest)-as.numeric(fore.arima)
MAPE(fore.arima,rJiiTest)
MSE(fore.arima,rJiiTest)
MAPE(fore.arima,rJiiTest)

## ---- Prediksi ARIMA-GARCH -----
var_resi_fore <- resi.fore.arima^2
var_garch_fore <- rep(NA, length(rJiiTest))
var_garch_fore[1] <- omega + alpha*tail(var_resid,1) + beta*tail(var_garch,1)
var_garch_fore
for(i in 2:length(rJiiTest)){
  var_garch_fore[i] <- omega + alpha * var_resi_fore[i-1] + beta * var_garch_fore[i-1]
}
std_residual_fore <- (resi.fore.arima - mean(resi.fore.arima)) / sd(resi.fore.arima)
resid_model_garch_fore <- std_residual_fore * sqrt(var_garch_fore)
resid_model_garch_fore <- xts(resid_model_garch_fore,order.by=time(rJiiTest))
arima_garch_fore <- as.numeric(fore.arima) + as.numeric(resid_model_garch_fore)

resi_arima_garch_fore <- rJiiTest - arima_garch_fore
mape_arima_garch_fore <- mean(abs(resi_arima_garch_fore/rJiiTest));mape_arima_garch_fore
MAPE(arima_garch_fore,rJiiTest)
RMSE(arima_garch_fore,rJiiTest)
MSE(arima_garch_fore,rJiiTest)
MAE(arima_garch_fore,rJiiTest)

df_fore_garch <- data.frame (
  waktu = time(rJiiTest),
  actual = as.numeric(rJiiTest),
  predict = as.numeric(arima_garch_fore),
  arima = as.numeric(fore.arima)
)
ggplot(df_fore_garch, aes(x = waktu)) +
  geom_line(aes(y = actual, color = "Actual"),
            size = 1.1, alpha = 0.95) +
  geom_line(aes(y = arima, color = "ARIMA"),
            size = 0.9, alpha = 0.85) +
  geom_line(aes(y = predict, color = "ARIMA-GARCH"),
            size = 0.8, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Actual" = "black",      
      "ARIMA" = "#009E73",        
      "ARIMA-GARCH" = "yellow"  
    )
  ) +
  labs(
    x = "Time",
    y = "Price",
    color = "Series"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.5, "cm")
  )


## ----- Prediksi ARIMA GJR GARCH ------
omega.gjr <- 28.853135
alpha.gjr <- 0.096802
beta.gjr <- 0.666569
eta.gjr <- 0.999972
var_resi_fore <- resi.fore.arima^2
var_gjr_fore <- rep(0,length(rJiiTest))
var_gjr_fore[1] <- omega.gjr + (alpha.gjr * tail(var_resid,1)) + (eta.gjr * tail(var_resid,1) * (ifelse(tail(residu,1) < 0, 1, 0))) + (beta.gjr * tail(var_gjr,1))
for(i in 2:length(rJiiTest)){
  indicator <- ifelse(resi.fore.arima[i-1] < 0, 1, 0)
  var_gjr_fore[i] <- omega.gjr + (alpha.gjr * var_resi_fore[i-1]) + (eta.gjr * var_resi_fore[i-1] * indicator) + (beta.gjr * var_gjr_fore[i-1])
}
std_residual_fore <- (resi.fore.arima - mean(resi.fore.arima)) / sd(resi.fore.arima)
resid_model_gjr_fore <- std_residual_fore * sqrt(var_gjr_fore)
resid_model_gjr_fore <- xts(resid_model_gjr_fore,order.by=time(rJiiTest))
arima_gjr_fore <- as.numeric(fore.arima) + as.numeric(resid_model_gjr_fore)

resi_arima_gjr_fore <- rJiiTest - arima_gjr_fore
mape_arima_gjr_fore <- mean(abs(resi_arima_gjr_fore/rJiiTest));mape_arima_gjr_fore
MAPE(arima_gjr_fore,rJiiTest)
RMSE(arima_gjr_fore,rJiiTest)
MSE(arima_gjr_fore,rJiiTest)
MAE(arima_gjr_fore,rJiiTest)

df_fore_gjr <- data.frame (
  waktu = time(rJiiTest),
  actual = as.numeric(rJiiTest),
  predict = as.numeric(arima_gjr_fore)
)
ggplot(df_fore_gjr,aes(x = waktu))+
  geom_line(aes(y = actual, color = "Actual"), size = 1) +
  geom_line(aes(y = predict, color = "Predict"), size = 1) +
  labs(title = "Aktual vs Predict ARIMA-GJR GARCH",
       x = "Waktu", y = "Nilai", color = "Series") +
  theme_minimal()


# --- compare testing (all) ----
df_fore_test_all <- data.frame (
  waktu = time(rJiiTest),
  actual = as.numeric(rJiiTest),
  predict_gjr = as.numeric(arima_gjr_fore),
  predict_garch = as.numeric(arima_garch_fore)
)

ggplot(df_fore_test_all, aes(x = waktu)) +
  geom_line(aes(y = actual, color = "Actual"), size = 1) +
  geom_line(aes(y = predict_gjr, color = "Predict GJR-GARCH"), size = 1) +
  geom_line(aes(y = predict_garch, color = "Predict GARCH"), size = 1) +
  scale_color_manual(
    values = c(
      "Actual" = "black",
      "Predict GJR-GARCH" = "blue",
      "Predict GARCH" = "red"
    )
  ) +
  labs(
    title = "Actual vs Forecast (GARCH vs GJR-GARCH)",
    x = "Time",
    y = "Price",
    color = "JII"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

