#ustalanie sciezki plikow
setwd("D:\\Rynki finansowe\\Quant invest")

#dolaczanie niezbednych biblioteki
library(quantmod)
library(fPortfolio)
library(PerformanceAnalytics)

#pobieranie notowan
x<-read.csv2("Quant_Invest_Fundusze.csv",header=T,stringsAsFactors = F)

#notowania
notowania<-xts(x[,-1],order.by=as.Date(x[,1]))
notowania<-na.locf(notowania)
storage.mode(notowania)<-"numeric"

#dane do konstrukcji wag portfela obejmuja okres od 2000 do 2018 roku
stopyzwrotu<-na.omit(100*diff(log(notowania)))["2000/2018"]

#benchmarki
wzorcowe<-notowania
benchmarki <- na.omit(100*diff(log(wzorcowe)))["2000/2018"]

#stopa wolna od ryzyka
rf <- read.csv2("5-Year-Poland-Bond-Yield.csv", header=T,stringsAsFactors = F)
rf <- xts(rf[,2],order.by=as.Date(rf$Data))
storage.mode(rf)<-"numeric"
rf_stacjonarne <- diff(rf)["2000/2018"]

#Elementy konstrukcji portfela Markowitza
scenarios <- dim(benchmarki)[1] #liczba obserwacji
assets <- dim(benchmarki)[2] #liczba aktywow w portfelu
data_ts <- as.timeSeries(benchmarki) #przeksztalcenie na szereg typu ts
spec <- portfolioSpec() #wywoâ‰¥anie standardowej specyfikacji portfela
setSolver(spec) <- "solveRquadprog" #wybor metody optymalizacji
setNFrontierPoints(spec) <- 20  #liczba punktow na granicy efektywnoÃºci
#Wykluczamy krotka sprzedaz
constraints <- c("LongOnly")
portfolioConstraints(data_ts, spec, constraints)
frontier <- portfolioFrontier(data_ts, spec, constraints)
#print(frontier)
frontierPlot(object=frontier)
tailoredFrontierPlot(object = frontier)
#weightsPlot(frontier, col = rainbow(assets))

#portfel maksymalizujacy wsk. Sharpe
maxsharpe<-tangencyPortfolio(data_ts)
wagi_sharpe<-maxsharpe@portfolio@portfolio$weights
wagi_sharpe #portfel Sharpe

benchmark_Markowitz<-t(wagi_sharpe%*%t(benchmarki["2000/2018"]))
benchmark_Markowitz<-xts(benchmark_Markowitz,order.by=as.Date(index(benchmarki["2000/2018"])))/10 #przyrownanie funduszy
names(benchmark_Markowitz)<- c("Benchmark")

srednia_wazona <- function(stopy, wagi) {
  if (dim(stopy)[2]==length(wagi) ){
    sumy = xts(0, as.Date("1900-01-01","%Y-%m-%d")) 
    for (i in 1:nrow(stopy)){
      suma = 0
      suma = sum(stopy[i,]*wagi)
      sumy<-c(sumy,xts(suma,as.Date(index(stopy[i]),"%Y-%m-%d")))
    }
    return(sumy[-1,])
  } 
}

#### Optymalizacja poprzez Algorytm Genetyczny
portfolio_returns = function(x, p.name="Portfel_optymalny") {
  port.returns = 0
  for (i in 1:length(x)) {
    port.returns = port.returns + stopyzwrotu[,i] * x[i]
  }
  names(port.returns) = p.name
  return (port.returns)
}

sharpe = function(x) {
  port.returns = portfolio_returns(x)
  return (mean(port.returns-rf_stacjonarne)/sqrt(var(port.returns)))
}

constraint = function(x) {
  boundary_constr = (sum(x)-1)**2   # "suma x = 1" ograniczenie
  
  for (i in 1:length(x)) {
    boundary_constr = boundary_constr + 
      max(c(0,x[i]-1))**2 +  # "x <= 1" ograniczenie
      max(c(0,-x[i]))**2     # "x >= 0" ograniczenie
  }
  return (boundary_constr)
}

zmienna <- -1000
zmienna2 <- 1000

obj = function(x) {
  # maksymalizacja wskaŸnika Sharpe oraz przeliczenie SRRI
  odchylenie <- sd(to.period(portfolio_returns(as.vector(x))["2014/2018"], period = 'weeks')*52)^0.5
  if(odchylenie>=2 & odchylenie <10){
   tmp <- (-sharpe(x))+100*constraint(x) 
  } else {
    tmp <- 45000
  }
  return (tmp)
}

library(GA)
ga_res = ga(
  type="real-valued", 
  
  # "ga" funkcja optymalizuje, wiêc mno¿ymy funkcjê celu przez -1
  function(x){-obj(x)},
  lower = rep(0,ncol(stopyzwrotu)), # x_i >= 0
  upper = rep(1,ncol(stopyzwrotu)), # x_i <= 1
  maxiter = 50000, # maksymalna liczba iteracji
  run=50, # Jezeli maksymalne dopasowanie pozostaje takie samo przez 50 iteracji
  # wtedy algorytm siê zatrzymuje
  parallel=TRUE,  # Pozakuje jak algorytm pracuje
  monitor=TRUE,
  seed=1 # u¿yteczne dla replikacji wyników
)

 # Wagi przyporzadkowane do funduszy
solution = as.vector(summary(ga_res)$solution)
cbind(names(stopyzwrotu), solution)

#wykres osi¹gnieæ portfela optymalnego na tle poszczególnych funduszy
par(mfrow=c(1,1))
optimal_returns = portfolio_returns(solution, p.name="Portfel optymalny")
plot(cumsum(optimal_returns),main="",type="l",lwd=3, ylim=c(-80, 200),major.ticks="quarters",grid.ticks.on="quarters")
lines(cumsum(stopyzwrotu[,1]),col="blue")
lines(cumsum(stopyzwrotu[,2]),col="red")
lines(cumsum(stopyzwrotu[,3]),col="green3")
lines(cumsum(stopyzwrotu[,4]),col="violet", lwd=3)
lines(cumsum(stopyzwrotu[,5]),col="peru")
lines(cumsum(stopyzwrotu[,6]),col="orange")
lines(cumsum(stopyzwrotu[,7]),col="yellow")
lines(cumsum(benchmark_Markowitz),col="black",lwd=2)
addLegend(legend.loc="topleft", legend.names=c("Portfel optymalny", colnames(notowania), names(benchmark_Markowitz)), 
          col=c("black","blue","red","green3","violet","peru","orange","yellow","black"),lwd=c(5,5,5,5,5,5,5,5,5))
#odchylenie standardowe portfela optymalnego z ostatnich 5-ciu lat (SRRI)
Std.dev.portfolio <- (sd(to.period(optimal_returns["2014/2018"], period = 'weeks'))*52)^0.5
print(Std.dev.portfolio)

optimal_returns_cumulative <- cumsum(optimal_returns)
benchmark_Markowitz_cumulative <- cumsum(benchmark_Markowitz)

#stopy zwrotu skumulowanych stóp zwrotu (5 i 1 roczna efektywna)
portf.return.5y <- as.numeric(optimal_returns_cumulative["2018-12-31"])-as.numeric(optimal_returns_cumulative[,1]["2013-12-31"])
print(portf.return.5y)
portf.return.1y = ((1+portf.return.5y/100)^(1/5)-1)*100
print(portf.return.1y)

#wyniki portfela optymalnego na tle benchmarku
plot(optimal_returns_cumulative,main="", ylim=c(-50,150), major.ticks="quarters",grid.ticks.on="quarters")
lines(benchmark_Markowitz_cumulative,col="blue", lwd=3)
addLegend(legend.loc="topleft", legend.names=c("Portfel optymalny", "Benchmark"), 
          col=c("black","blue"),lwd=c(5,5))

#daty w 2019
daty <- seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by="day")
daty<- daty[!weekdays(daty) %in% c('sobota','niedziela')]
swieta_polskie <- as.Date(c("2019-01-01", "2019-04-22", "2019-05-01", "2019-05-03", "2019-06-20",
                "2019-08-15", "2019-11-01", "2019-11-11", "2019-12-25", "2019-12-26"))
daty <- daty[which(!daty %in% swieta_polskie)]
dni_pracujace <- length(daty)
dni_calosc <- c(index(optimal_returns_cumulative), daty)

library(forecast)
#roznicowanie -znalezenie optymlnego stopnia roznicowaniania portfela
ndiffs(optimal_returns_cumulative) #optimal_returns

#ustalenie stopnia opóŸnieñ 
Acf(optimal_returns)
#wynika z wykresu, ze opoznienie autoregresji wyniesie 1
Pacf(optimal_returns)
#opoznienie sredniej ruchomej wyniesie 0

#ustalenie hidden layer w sieciach neuronowych
#Pacf(stopyzwrotu$AP)
#Pacf(stopyzwrotu$ARR)
#Pacf(stopyzwrotu$ARW)
#Pacf(stopyzwrotu$G)
#Pacf(stopyzwrotu$OP)
#Pacf(stopyzwrotu$ORR)
#Pacf(stopyzwrotu$ORW)

#for(i in 1:ncol(stopyzwrotu)){
#print(auto.arima(stopyzwrotu[,i]))
#}

#przygotowanie danych do ML
dataML <- xts(cbind(optimal_returns, stopyzwrotu)) #ML - machine learning
index_data <- round(0.75*nrow(dataML))

train <- dataML[1:index_data,]
test <- dataML[-c(1:index_data),]

#pierwszy sposob -Arima

arima.fit <- Arima(optimal_returns_cumulative, order=c(1,1,0))
arima.fit.prediction <- forecast(arima.fit, nrow(test))
plot(arima.fit.prediction)
arima.fit.prediction.xts <- xts(arima.fit.prediction$mean, order.by=index(test))
MSE.arima <- sum((diff(arima.fit.prediction.xts)-test$Portfel.optymalny[-1,]))^2/nrow(test)

#drugi spsosob -uczenie maszynowe Fitting Generalized Linear Models
lm.fit <- glm(Portfel.optymalny~., data=train)
#summary(lm.fit)
pr.lm <- predict(lm.fit,test)
MSE.lm <- sum((pr.lm - test$Portfel.optymalny)^2)/nrow(test)
#print(MSE.lm)

#trzecie -sieci nauronowe

#skalowanie metoda max min
maxs <- apply(dataML, 2, max) 
mins <- apply(dataML, 2, min)
scaled <- as.data.frame(scale(dataML, center = mins, scale = maxs - mins))

trainNN <- scaled[1:index_data,]
testNN <- scaled[-c(1:index_data),]

library(neuralnet)
trainNNdataset <-neuralnet(Portfel.optymalny~AP+ARR+ARW+G+OP+ORR+ORW, data=trainNN, hidden=c(4,5), linear.output = TRUE)
plot(trainNNdataset)
#hidden 4,5 zosta³o wybrane na podstawie dopasowania metod¹ Pacf (badania stopnia autoregresji - trzeba to jakos ³adnie uj¹æ)

predicted.dataNN <- compute(trainNNdataset,testNN)

predicted.dataNN2 <- predicted.dataNN$net.result*(max(dataML$Portfel.optymalny)-min(dataML$Portfel.optymalny))+min(dataML$Portfel.optymalny)
test.r <- (testNN$Portfel.optymalny)*(max(dataML$Portfel.optymalny)-min(dataML$Portfel.optymalny))+min(dataML$Portfel.optymalny)
MSE.nn <- sum((test.r - predicted.dataNN2)^2)/nrow(test)
MSE.nn
MSE.lm
MSE.arima
#najmniejszy b³¹d ma glm


plot(testNN$Portfel.optymalny,predicted.dataNN2,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
#abline(-1,1,lwd=2) #real?????
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

plot(testNN$Portfel.optymalny,pr.lm,col='blue',main='Real vs predicted lm',pch=18, cex=0.7)
#abline(0,1,lwd=2) #real?????
legend('bottomright',legend='LM',pch=18,col='blue', bty='n', cex=.95)

plot(testNN$Portfel.optymalny, diff(arima.fit.prediction.xts), col='green',main='Real vs predicted ARIMA',pch=18, cex=0.7)

#ARIMA prediction
prognoza_ARIMA_ <-forecast(arima.fit, 251)
prognoza_ARIMA <- xts(prognoza_ARIMA_$mean, order.by=daty)
plot(diff(prognoza_ARIMA), col='green')
prognoza_ARIMA_calosc <- c(optimal_returns_cumulative, prognoza_ARIMA)
plot(prognoza_ARIMA_calosc)
rect(xleft="2019-01-01",xright ="2019-12-31", ybottom=-Inf, ytop=Inf, col = "blue" )
#nie dziala rect, jak macie czas to spróbujcie prosze zaznaczyc okres 2019

#GLM pediction
pr_lm_1y <- xts(pr.lm[1:251], order.by = daty)
plot(pr_lm_1y, col='blue')
plot(cumsum(pr_lm_1y))

prognoza_GLM_calosc <- c(optimal_returns, pr_lm_1y)
plot(prognoza_GLM_calosc)
plot(cumsum(prognoza_GLM_calosc))
rect(xleft="2019-01-01",xright ="2019-12-31", ybottom=-Inf, ytop=Inf, col = "blue" )
#nie dziala rect, jak macie czas to spróbujcie prosze zaznaczyc okres 2019


#NN prediction
prediction_dataNN_1y <- xts(predicted.dataNN2[1:251], order.by = daty)
plot(prediction_dataNN_1y, col='red')
plot(cumsum(prediction_dataNN_1y))

prognoza_NN_calosc <- c(optimal_returns, prediction_dataNN_1y)
plot(prognoza_NN_calosc)
plot(cumsum(prognoza_NN_calosc))
rect(xleft="2019-01-01",xright ="2019-12-31", ybottom=-Inf, ytop=Inf, col = "blue" )
#nie dziala rect, jak macie czas to spróbujcie prosze zaznaczyc okres 2019

###CROSS-VALIDATION

library(boot)

#ARIMA
k <- 12
cv.arima.error <-vector()
for(i in 1:k){
  index_data_CV_a <- sample(1:nrow(optimal_returns_cumulative),round(0.9*nrow(optimal_returns_cumulative)))
  train.cv_a <- optimal_returns_cumulative[index_data_CV_a,]
  test.cv_a <- optimal_returns_cumulative[-index_data_CV_a,]
  arima.cv <- Arima(train.cv_a, order=c(1,1,0))
  arima.prediction.cv <- forecast(arima.cv, nrow(test.cv_a))
  arima.prediction.cv <- xts(arima.prediction.cv$mean, order.by = index(test.cv_a))
  cv.arima.error[i] <- sum(na.omit(diff(test.cv_a)-diff(arima.prediction.cv)))^2/(nrow(arima.prediction.cv)-1)
}
cv.arima.error

mean(cv.arima.error)
boxplot(cv.arima.error,xlab='MSE CV',col='cyan',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for ARIMA',horizontal=TRUE)

#linear model
set.seed(200)
lm.fit.CV <- glm(Portfel.optymalny~., data=dataML)
#cv.glm.lm <- cv.glm(dataML,lm.fit.CV,K=12)$delta[1]
cv.glm.lm <- cv.glm(dataML,lm.fit.CV,K=12)

boxplot(cv.glm.lm$delta[1],xlab='MSE CV',col='cyan',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for ARIMA',horizontal=TRUE)

#neural networks
set.seed(450)
cv.error <- NULL
k <- 12

library(plyr) 
pbar <- create_progress_bar('text')
pbar$init(k)

for(i in 1:k){
  index_data_CV <- sample(1:nrow(dataML),round(0.9*nrow(dataML)))
  train.cv <- scaled[index_data_CV,]
  test.cv <- scaled[-index_data_CV,]
  
  nn.cv <- neuralnet(Portfel.optymalny~AP+ARR+ARW+G+OP+ORR+ORW,data=train.cv,hidden=c(4,4),linear.output=T)
  
  pr.nn.cv <- compute(nn.cv,test.cv)
  pr.nn.cv <- pr.nn.cv$net.result*(max(dataML$Portfel.optymalny)-min(dataML$Portfel.optymalny))+min(dataML$Portfel.optymalny)
  
  test.cv.r <- (test.cv$Portfel.optymalny)*(max(dataML$Portfel.optymalny)-min(dataML$Portfel.optymalny))+min(dataML$Portfel.optymalny)
  
  cv.error[i] <- sum((test.cv.r - pr.nn.cv)^2)/nrow(test.cv)
  
  pbar$step()
}

mean(cv.error)
cv.error

boxplot(cv.error,xlab='MSE CV',col='cyan',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for NN',horizontal=TRUE)


  
#Prognoza portfela i benchmarku 1-rocznego metod¹ ARIMA
#portfel_vs_benchmark <- cbind(optimal_returns_cumulative, benchmarki_wazone_cumulative)
#colnames(portfel_vs_benchmark) <- c("Skumulowane stopy z portfela", "Skumulowane stopy z benchmarku")
#predicted_portfel=forecast(auto.arima(portfel_vs_benchmark$`Skumulowane stopy z portfela`),dni_pracujace)
#predicted_benchmark=forecast(auto.arima(portfel_vs_benchmark$`Skumulowane stopy z benchmarku`),dni_pracujace)

#par(mfrow=c(2,1))
#plot(predicted_portfel, main="")
#lines(fitted(pred2),col="red")
#plot(predicted_benchmark, col="red", main="")
#portfel_prognozowany <- data.frame(fitted(predicted_portfel))
#benchmark_prognozowany <- data.frame(fitted(predicted_benchmark))
#portfel_vs_benchmark_prognozowany <- cbind(fitted(portfel_prognozowany), fitted(benchmark_prognozowany))
#names(portfel_vs_benchmark_prognozowany) <- c("Portfel prognozowany", "Benchmark prognozowany")
#as.xts(portfel_vs_benchmark_prognozowany, order.by=as.Date())
#predicted_portfel$mean

#statystyki bazowe
library("xtable")
statystyki_stopy <- basicStats(portfel_vs_benchmark) #portfel optymalny
#xtable(statystyki_stopy)
portfel_prognozowany_stat <- basicStats(portfel_prognozowany[4800:4956,]) #portfel prognozowany
#xtable(portfel_prognozowany_stat)
benchmark_prognozowany_stat <- basicStats(benchmark_prognozowany[4800:4956,]) #benchmark prognozowany
#xtable(benchmark_prognozowany_stat)

#statystyki portfela optymalnego

table.AnnualizedReturns(diff(portfel_vs_benchmark), geometric=FALSE)
#xtable(table.AnnualizedReturns(diff(portfel_vs_benchmark, geometric=FALSE)))

dane.full<-na.omit(na.locf(merge(optimal_returns,benchmarki_wazone)))
table.CAPM(dane.full[,1]/100,dane.full[,2]/100,Rf=0)
#xtable(table.CAPM(optimal_returns,benchmarki_wazone,scale=252,Rf=0))

#statystyki portfela prognozowanego
dane.prog<-diff(portfel_vs_benchmark_prognozowany)

table.AnnualizedReturns(dane.prog/100, geometric=FALSE)
#xtable(table.AnnualizedReturns(dane.prog, geometric=FALSE))
table.CAPM(dane.prog[,1]/100,dane.prog[,2]/100,Rf=0)
#xtable(table.CAPM(dane.prog[,1]/100,dane.prog[,2]/100,Rf=0))

### u¿yæ analize casuality
#i okresle przyczynowosc
#causality(model.var,cause="GSPC.Adjusted") #stopa zwrotu wplywa na wolumen
#causality(model.var,cause="GSPC.Volume") #wolumen wplywa na stope zwrotu

#na podst danych rynkowych wygenerowac dane scenariusze ekonomiczne
