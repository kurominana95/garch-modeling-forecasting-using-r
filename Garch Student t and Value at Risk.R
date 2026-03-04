#--- Import Library -----
library(rugarch)
library(aTSA)
library(FinTS) 
library(lmtest)
library(forecast)
library(TSA)
library(tseries)
library(xts)
library(openxlsx)
library(tidyverse)
library(nortest)
library(PerformanceAnalytics)
library(ggplot2)
library(MLmetrics)
library(fGarch)

#------ Set Working Directory -----
setwd("C:/Users/NAJWA KHOIR ALDWYAH/OneDrive/Dokumen/Statistika/Semester 7/Teori resiko")

#------ Import Data ------
stock <- read.xlsx("Stock.xlsx", sheet=2)

#----- Preprocessing Data -----
stock$Date <- as.Date(as.numeric(stock$Date), origin = "1899-12-30")
ptba <- xts(stock$PTBA, order.by=stock$Date)

#------return------
return_ptba <- CalculateReturns(ptba, method = "log")
return_ptba <- return_ptba[-1]
hist(return_ptba)

#----- Plot histogram -----
hist(return_ptba, 
     breaks = 30,
     freq = FALSE,
     col = "coral", 
     border = "black",
     main = "Histogram Return PTBA",
     cex.main = 1.8,
     xlab = "Return Log PTBA")
lines(density(return_ptba, na.rm = TRUE), 
      col = "black", 
      lwd = 3)

#---- Split Train-Test ----
nTrain <- 0.9*length(return_ptba);round(nTrain)
rptbaTrain <- return_ptba[1:round(nTrain)]
rptbaTest <- return_ptba[(round(nTrain)+1):length(return_ptba)];length(rptbaTest)
ks.test(return_ptba,"pnorm",mean(return_ptba),sd(return_ptba))

#---- Visualisasi ----
plot(return_ptba)
write.xlsx(return_ptba,"return_ptba.xlsx")
plot(rptbaTrain)
df_train <- data.frame(Date = as.Date(time(rptbaTrain),frac=1), value=as.numeric(rptbaTrain),set="Training")
df_test <- data.frame(Date = as.Date(time(rptbaTest),frac=1), value=as.numeric(rptbaTest),set="Testing")
df_all <- rbind(df_train,df_test)
ggplot(df_all, aes(x = Date, y = value, color = set)) +
  geom_line(size = 0.6) +
  scale_color_manual(values = c("Training" = "black", "Testing" = "coral")) +
  labs(title = "Time Series Plot of Return PTBA",
       x = "Date", y = "Return of PTBA") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",        
      size = 18
    ),
    legend.position = "bottom",   
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

#---- Stasioneritas ----
adf.test(rptbaTrain)
adf.test(rptbaTest)

rptbaTrain <- rptbaTrain[-1]

#----- Box Jenkins ARIMA ----
##---- Identifikasi Model ----
acf(rptbaTrain)
pacf(rptbaTrain)

##---- Differencing ----
rptbaTrain.1 <- diff(rptbaTrain)
adf.test(rptbaTrain.1[-1])
acf(rptbaTrain.1[-1,])
pacf(rptbaTrain.1[-1,])
eacf(rptbaTrain.1[-1,])
rptbaTrain.1 <- rptbaTrain.1[-1,]

##---- Modeling -----
modelarima <- arima(rptbaTrain, order=c(0,1,1), method="ML", include.mean=FALSE)
coeftest(modelarima)
summary(modelarima)
AIC(modelarima)
MSE(fitted(modelarima),rptbaTrain)
MAE(fitted(modelarima),rptbaTrain)

##---- Uji Diagnostik ----
checkresiduals(modelarima$residuals)
print(Box.test(modelarima$residuals, lag=1, type="Ljung-Box"))
jarque.bera.test(modelarima$residuals)
ks.test(modelarima$residuals,"pnorm",mean(modelarima$residuals),sd(modelarima$residuals))
arch.test(modelarima)

#----- ARCH GARCH -----
## ----- Model ARCH GARCH -----
save_garch <- function(X, p, q, type = "std"){
  spc <- ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(p, q), variance.targeting=FALSE),
                    mean.model = list(armaOrder = c(0, 1),
                                      include.mean = F),
                    distribution.model = type)
  mod <- ugarchfit(spec = spc, data = X)
  return(list(model=mod,spec=spc))
}
g10 <- save_garch(X = rptbaTrain.1, p = 1, q = 1)
infocriteria(g10$model)
g10$model
signbias(g10$model)
ArchTest(residuals(g10$model,standardize=TRUE),lags=4)
ArchTest(residuals(g10$model,standardize=TRUE),lags=8)
ArchTest(residuals(g10$model,standardize=TRUE),lags=12)


#------Var ECF-------
library(e1071)
b=matrix(nrow = 7,ncol=1)
colnames(b) <- c("Return Harga Aset")
rownames(b)<- c("Mean","Std.dev","Variansi","Maksimum","Minimum","Skewness","Kurtosis")
b[1,]=mean(return_ptba)
b[2,]=sd(return_ptba)
b[3,]=var(return_ptba)
b[4,]=min(return_ptba)
b[5,]=max(return_ptba)
b[6,]= skewness(return_ptba)
b[7,]= kurtosis(return_ptba)
print(b)

index(return_ptba)[which.min(return_ptba)]
index(return_ptba)[which.max(return_ptba)]

Calculate_ECF <- function(alpha=0.05, mean, sd, skew, kurt, hold=1, nilai_investasi=1000000)
{
  par_kuantil = qnorm(alpha, mean = 0, sd = 1)
  par_skew = skew
  par_kurt = kurt
  ECF = par_kuantil + (((par_kuantil^2 - 1)*par_skew)/6) + (((par_kuantil^3 - (3*par_kuantil))*(par_kurt-3))/24) -
    (((2*par_kuantil^3 - (5*par_kuantil))*par_skew^2)/36)
  
  var = nilai_investasi*(mean+ECF*sd)*sqrt(hold)
  risk = var/nilai_investasi
  return(list(var = var, risk = risk, ecf = ECF))
}

ecf_5 = Calculate_ECF(alpha=0.05,mean=b[1,],sd=b[2,],skew=b[6,],kurt=b[7,])
ecf_1 = Calculate_ECF(alpha=0.01,mean=b[1,],sd=b[2,],skew=b[6,],kurt=b[7,])
ecf_5$var
ecf_5$risk
ecf_1$var
ecf_1$risk

return_ptba
-0.984407
ma_full <- arima(return_ptba,order=c(0,1,1),fixed = c(-0.984407), include.mean = FALSE, method = "ML")
summary(ma_full)
fore_ma <- predict(ma_full,h=1)$pred
fore_ma
fit_ma <- fitted(ma_full)

spec_fixed <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,1), include.mean = FALSE),
  distribution.model = "std",
  
  fixed.pars = list(
    ma1    = -0.984406,
    omega  = 0.0001663,
    alpha1 = 0.2306235,
    beta1  = 0.4248633,
    shape  = 2.5996696
  )
)
fit <- ugarchfit(spec_fixed, data = return_ptba, solver = "hybrid")
sigma_t <- sigma(fit);as.numeric(sigma_t)
et <- as.numeric(tail((return_ptba),1) - as.numeric(tail(fit_ma,1)));et
sigma_hat <- 0.0001663 + 0.4248633*et^2 + 0.2306235*as.numeric(tail(sigma_t,1));sigma_hat
coef_fit  <- coef(g10$model)
nu_hat    <- if("shape" %in% names(coef_fit)) coef_fit["shape"] else NA ;nu_hat

mu     <- mu_hat
sigma  <- sigma_hat
nu     <- nu_hat
mu;sigma;nu

# -----------------------------------------
# Fungsi untuk menghitung VaR & ES Student-t
# -----------------------------------------

compute_VaR_ES <- function(alpha, mu, sigma, nu){
  
  # Quantile raw (non-standardized)
  q_raw <- qt(alpha, df = nu)
  # Standardized t-quantile
  q_std <- q_raw * sqrt((nu - 2)/nu)
  # VaR parametric (1-step forecast)
  VaR  <- mu + sigma * q_std
  # pdf at quantile
  f_q  <- dt(q_raw, df = nu)
  # Standardized ES (left tail)
  ES_std <- - ( f_q * (nu + q_raw^2) ) / ( alpha * (nu - 1) ) * sqrt((nu - 2)/nu)
  # ES parametric final
  ES  <- mu + sigma * ES_std
  return(list(VaR = VaR, ES = ES))
}

# -----------------------------------------
# Hitung untuk alpha 5% dan 1%
# -----------------------------------------

res_5  <- compute_VaR_ES(alpha = 0.05, mu = b[1,], sigma = b[2,], nu = 2.59967)
res_1  <- compute_VaR_ES(alpha = 0.01,  mu = b[1,], sigma = b[2,], nu = 2.59967)

# -----------------------------------------
# Tampilkan hasil
# -----------------------------------------

cat("===========================================\n")
cat("   PARAMETRIC STUDENT-t VaR & ES (STATIC)\n")
cat("===========================================\n\n")

cat("Alpha = 5% (VaR 95%)\n")
cat("VaR_5% =", res_5$VaR, "\n")
cat("ES_5%  =", res_5$ES,  "\n\n")

cat("Alpha = 1% (VaR 99%)\n")
cat("VaR_1% =", res_1$VaR, "\n")
cat("ES_1%  =", res_1$ES,  "\n")

#----Evaluasi prediksi----
mean_model_fix <- arima(rptbaTrain,order=c(0,1,1),fixed = c(-0.984407), include.mean = FALSE, method = "ML")
summary(mean_model_fix)
fitted_values <- as.numeric(fitted(mean_model_fix));length(fitted_values)
fore_test <- predict(mean_model_fix, n.ahead = 27)
fore_test$pred
MSE(fore_test$pred,rptbaTest)
MAE(fore_test$pred,rptbaTest)

## ----- Fit ARIMA GARCH (1,1) ------
omega  = 0.000166
alpha = 0.230624
beta  = 0.42486
shape  = 2.5996696
residu <- residuals(mean_model_fix)
var_resid <- residu^2
var_garch <- rep(0, length(rptbaTrain))
var_garch[1] <- var(var_resid)  # bisa pakai variance awal
for(i in 2:length(rptbaTrain)){
  var_garch[i] <- omega + alpha * var_resid[i-1] + beta * var_garch[i-1]
}
std_residual <- (residu - mean(residu)) / sd(residu)
resid_model <- std_residual * sqrt(var_garch)
resid_model <- xts(resid_model,order.by=time(rptbaTrain))
arima_garch_fit <- as.numeric(fitted_values) + as.numeric(resid_model)
resi_arimag <- rptbaTrain - arima_garch_fit
mape_arimag <- mean(abs(resi_arimag)/rptbaTrain);mape_arimag
MSE(arima_garch_fit,rptbaTrain)
MAE(arima_garch_fit,rptbaTrain)
smape(rptbaTrain,arima_garch_fit)

df_garch <- data.frame(
  waktu = time(rptbaTrain),
  actual = as.numeric(rptbaTrain),
  fitted = as.numeric(arima_garch_fit)
)

ggplot(df_garch, aes(x = waktu)) +
  geom_line(aes(y = actual, color = "Actual"), size = 1) +
  geom_line(aes(y = fitted, color = "Fitted"), linetype = "dashed", size = 0.8) +
  labs(title = "Aktual vs Fitted ARIMA-GARCH",
       x = "Waktu", y = "Nilai", color = "Series") +
  theme_minimal()

## ---- Fore ARIMA-GARCH
resi.fore.arima <- as.numeric(rptbaTest)-as.numeric(fore_test$pred)
var_resi_fore <- resi.fore.arima^2
var_resi_fore
var_garch_fore <- rep(0, length(rptbaTest))
var_garch_fore[1] <- tail(var_garch, 1) 
var_garch_fore
for(i in 2:length(rptbaTest)){
  var_garch_fore[i] <- omega + alpha * var_resi_fore[i-1] + beta * var_garch_fore[i-1]
}
std_residual_fore <- (resi.fore.arima - mean(resi.fore.arima)) / sd(resi.fore.arima)
resid_model_garch_fore <- std_residual_fore * sqrt(var_garch_fore)
resid_model_garch_fore <- xts(resid_model_garch_fore,order.by=time(rptbaTest))
arima_garch_fore <- as.numeric(fore_test$pred) + as.numeric(resid_model_garch_fore)
resi_arima_garch_fore <- rptbaTest - arima_garch_fore
mape_arima_garch_fore <- mean(abs(resi_arima_garch_fore)/rptbaTest);mape_arima_garch_fore
MSE(arima_garch_fore,rptbaTest)
MAE(arima_garch_fore,rptbaTest)
RMSE(arima_garch_fore,rptbaTest)
smape(rptbaTest, arima_garch_fore)

df_fore_garch <- data.frame (
  waktu = time(rptbaTest),
  actual = as.numeric(rptbaTest),
  predict = as.numeric(arima_garch_fore)
)
ggplot(df_fore_garch, aes(x = waktu)) +
  geom_line(aes(y = actual, color = "Actual"), size = 1) +
  geom_line(aes(y = predict, color = "Predict"), size = 1) +
  scale_color_manual(values = c("Actual" = "black", 
                                "Predict" = "#CD5B45")) +
  labs(title = "Actual vs Predict IMA(1,1)-GARCH-t (1,1)",
       x = "Date", y = "Return", color = "Series") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",        
      size = 16
    ),
    legend.position = "bottom",   
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

miu_t = -0.001219732
sigma2_t = omega + alpha*tail(var_resi_fore,1) + beta*tail(var_garch_fore,1)
sigma2_t
sigma_t = sqrt(sigma2_t)
sigma_t

var_st_95 = compute_VaR_ES(0.05,mu = miu_t,sigma = sigma_t, nu = 2.59967)$VaR
var_st_99 = compute_VaR_ES(0.01,mu = miu_t,sigma = sigma_t, nu = 2.59967)$VaR
var_st_95
var_st_99
var_ecf_95 = Calculate_ECF(0.05,mean=miu_t,sd=sigma_t,skew=b[6,],kurt = b[7,])$risk
var_ecf_99 = Calculate_ECF(0.01,mean=miu_t,sd=sigma_t,skew=b[6,],kurt = b[7,])$risk
var_ecf_95
var_ecf_99
