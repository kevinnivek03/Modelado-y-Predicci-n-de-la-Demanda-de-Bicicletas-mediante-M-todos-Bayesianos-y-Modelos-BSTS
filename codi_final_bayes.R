# ============================
# 1. Cargar librerías
# ============================

library(ggplot2)
library(skimr)
library(dplyr)

# ============================
# 2. Cargar datos
# ============================
day <- read.csv("day.csv")

# ============================
# 3. Estructura y resumen
# ============================
str(day)
summary(day)
skim(day)

# ============================
# 4. Histograma general de cnt
# ============================
windows()
hist(day$cnt, col="lightblue", border="white",
     main="Histograma de la demanda diaria (cnt)",
     xlab="Número de bicicletas por día")

# ============================
# 5. Serie de tiempo
# ============================
day$fecha <- as.Date(day$dteday)
ggplot(day, aes(x = fecha, y = cnt)) +
  geom_line() +
  labs(title = "Demanda diaria de bicicletas en el tiempo",
       x = "Fecha", y = "Número de bicicletas (cnt)")

# ============================
# 6. STL (descomposición)
# ============================
ts_cnt <- ts(day$cnt, frequency = 365)
windows()
plot(stl(ts_cnt, s.window = "periodic"),
     main = "Descomposición STL de la demanda diaria de bicicletas")

# ============================
# 7. Scatter cnt vs registered
# ============================
ggplot(day, aes(x = registered, y = cnt)) +
  geom_point(alpha=0.4, color="skyblue3") +
  geom_smooth(method="loess", color="red") +
  labs(title = "Relación entre usuarios y demanda",
       x = "registered", y = "cnt")

# ============================
# 8. Boxplot por día laboral
# ============================
day$workingday <- factor(day$workingday,
                         labels=c("No laboral","Laboral"))

ggplot(day, aes(x = workingday, y = cnt)) +
  geom_boxplot(fill="lightblue", alpha=0.6) +
  labs(title="Demandas cnt por workingday")

# ============================
# 9. Matriz de correlación
# ============================
vars <- day %>% select(temp, atemp, hum, windspeed, casual, registered, cnt)
cor(vars)

# ==========================================
# 1. Matriz de diseño X y variable respuesta y
# ==========================================

# Selección de predictores relevantes
vars_modelo <- c("registered", "temp", "atemp", "hum", "windspeed", "workingday")

# Matriz diseño con intercepto
X <- model.matrix(~ registered + temp + atemp + hum + windspeed + workingday, data = day)

# Variable respuesta
y <- day$cnt

n <- nrow(X)
p <- ncol(X)

# Vista rápida
head(X)

# ==========================================
# 2. Definir parámetros de las prioris
# ==========================================

b0 <- rep(0, p)              # vector de ceros
B0 <- 0.001 * diag(p)        # matriz casi no informativa
c0 <- 0.001                  # IG shape
d0 <- 0.001                  # IG rate

# ==========================================
# 3. Implementar el Gibbs para β y σ²
# ==========================================

library(mvtnorm)

Gibbs_Regresion <- function(cant_aleat, beta_ini, sigma2_ini){
  
  beta_aleat <- vector("list", cant_aleat + 1)
  sigma2_aleat <- rep(NA, cant_aleat + 1)
  
  # iniciales
  beta_aleat[[1]] <- beta_ini
  sigma2_aleat[1] <- sigma2_ini
  
  for(i in 1:cant_aleat){
    
    # ----- Condicional de beta | sigma2 -----
    Bn <- solve(B0 + (1/sigma2_aleat[i]) * t(X) %*% X)
    bn <- Bn %*% (B0 %*% b0 + (1/sigma2_aleat[i]) * t(X) %*% y)
    
    beta_aleat[[i+1]] <- as.numeric(rmvnorm(1, mean = bn, sigma = Bn))
    
    # ----- Condicional de sigma2 | beta -----
    resid <- y - X %*% beta_aleat[[i+1]]
    cn <- c0 + n/2
    dn <- d0 + 0.5 * t(resid) %*% resid
    sigma2_aleat[i+1] <- 1 / rgamma(1, shape = cn, rate = dn)
  }
  
  return(list(beta = beta_aleat, sigma2 = sigma2_aleat))
}

# ==========================================
# 4. Ejecutar el muestreador
# ==========================================

set.seed(123)
resultados <- Gibbs_Regresion(
  cant_aleat = 10000,
  beta_ini = rep(0, p),
  sigma2_ini = 1
)

# Extraer cadenas
beta_chain <- do.call(rbind, resultados$beta)
sigma2_chain <- resultados$sigma2

# Resumen inicial
colMeans(beta_chain)
mean(sigma2_chain)

# ==========================================
# 5. Graficar histogramas de las betas
# ==========================================

windows()
par(mfrow=c(3,3))
for(j in 1:p){
  hist(beta_chain[,j], col="lightblue", border="white",
       main=paste("Posterior β", j-1),
       xlab="")
}
# ==========================================
# 6. Diagnósticos MCMC
# ==========================================

library(coda)

cadenas_par = cbind(beta_chain, sigma2_chain)

# Cada columna es un parámetro. Dividamos en 2 cadenas:
cad1 = cadenas_par[1:5000, ]
cad2 = cadenas_par[5001:10000, ]

# Convertimos a objetos mcmc
mc1 = mcmc(cad1)
mc2 = mcmc(cad2)

# Creamos la lista de cadenas
mcmc_chains = mcmc.list(mc1, mc2)

gelman.diag(mcmc_chains)
gelman.plot(mcmc_chains)

heidel.diag(mc1)
heidel.diag(mc2)

raftery.diag(mcmc_chains)

# ==========================================
# 3. Construcción del Modelo BSTS
# ==========================================

# Serie de tiempo de la demanda
y_ts <- ts(day$cnt, frequency = 7)  # frecuencia semanal

# Matriz de regresores (la misma que usamos en la regresión bayesiana)
X_bsts <- day[, c("registered", "temp", "atemp", "hum", "windspeed", "workingday")]

# ==========================================
# Definir los componentes del modelo BSTS
# ==========================================

ss <- list()

# 1. Componente de tendencia local ("local level + slope")
ss <- AddLocalLinearTrend(ss, y_ts)

# 2. Componente de estacionalidad semanal (7 días)
ss <- AddSeasonal(ss, y_ts, nseasons = 7)

# ==========================================
# Priors automáticos del paquete BSTS
# (incluyen varianzas de tendencia y estacionalidad)
# ==========================================

modelo_bsts <- bsts(
  formula = cnt ~ registered + temp + atemp + hum + windspeed + workingday,
  state.specification = ss,
  data = day,
  niter = 5000,
  ping = 0,
  seed = 123
)

# Resumen del modelo
summary(modelo_bsts)

# ==========================================
# 7. Modelo BSTS para la demanda diaria
# ==========================================


# ================================
# 1. Serie de tiempo
# ================================
y_ts <- ts(day$cnt, frequency = 7)

# ================================
# 2. Matriz de predictores
# ================================
X_bsts <- day[, c("registered", "temp", "atemp", "hum", "windspeed", "workingday")]

# ================================
# 3. Especificación del estado
# ================================
ss <- list()

# 3.1 Tendencia local (Local Linear Trend)
ss <- AddLocalLinearTrend(ss, y_ts)

# 3.2 Estacionalidad semanal (7 días)
ss <- AddSeasonal(ss, y_ts, nseasons = 7)

# ================================
# 4. Modelo BSTS completo
# ================================
modelo_bsts <- bsts(
  formula = cnt ~ registered + temp + atemp + hum + windspeed + workingday,
  state.specification = ss,
  data = day,
  niter = 5000,
  seed = 123,
  ping = 0
)

# ================================
# 5. Resultados
# ================================
summary(modelo_bsts)
plot(modelo_bsts)
plot(modelo_bsts, "components")
plot(modelo_bsts, "coefficients")
plot(modelo_bsts, "residuals")
plot(modelo_bsts, "predictors")



# ================================
# 4. Pronóstico a 30 días
# ================================

# ------------------------------
# 1. Serie de tiempo
# ------------------------------
y_ts <- ts(day$cnt, frequency = 7)

# ------------------------------
# 2. Componentes estructurales
# ------------------------------
ss <- list()
ss <- AddLocalLinearTrend(ss, y_ts)
ss <- AddSeasonal(ss, y_ts, nseasons = 7)

# ------------------------------
# 3. Ajustar modelo SIN regresores
# ------------------------------
modelo_bsts_estructural <- bsts(
  y_ts ~ 1,
  state.specification = ss,
  niter = 5000,
  ping = 0,
  seed = 123
)

# ------------------------------
# 4. Pronóstico a 30 días
# ------------------------------
# Horizonte de pronóstico
h <- 30  

# Crear matriz de predictores para los 30 días futuros
newdata_bsts <- day %>% 
  tail(1) %>%                      # tomar el último día real
  slice(rep(1:n(), each = h)) %>%  # replicarlo 30 veces
  select(registered, temp, atemp, hum, windspeed, workingday)

# Pronóstico con BSTS
pred <- predict(
  modelo_bsts_estructural,
  horizon = h,
  newdata = newdata_bsts,
  burn = 1000
)

# ------------------------------
# 5. Data frame de pronóstico
# ------------------------------
forecast_df <- data.frame(
  day = max(day$fecha) + 1:h,
  mean = pred$mean,
  low  = pred$interval[, 1],
  high = pred$interval[, 2]
)

# ------------------------------
# 6. Gráfico corregido
# ------------------------------
library(ggplot2)

ggplot(forecast_df, aes(x = day, y = mean)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "skyblue2", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  labs(
    title = "Pronóstico BSTS (30 días) – Modelo Estructural",
    subtitle = "Tendencia + Estacionalidad semanal",
    x = "Fecha",
    y = "Demanda estimada (cnt)"
  ) +
  theme_minimal()

# ============================
# Gráfico de descomposición del BSTS
# ============================

plot(modelo_bsts, "components")


# ============================
# MODELO 1: Estacional
# ============================

# Modelo 1: solo tendencia + estacionalidad
y_ts <- ts(day$cnt, frequency = 7)
ss1 <- AddLocalLinearTrend(list(), y_ts)
ss1 <- AddSeasonal(ss1, y_ts, nseasons = 7)

modelo_estacional <- bsts(
  cnt ~ 1,
  state.specification = ss1,
  data = day,
  niter = 3000,
  ping = 0,
  seed = 123
)

# Modelo 2: completo con regresores
ss2 <- AddLocalLinearTrend(list(), y_ts)
ss2 <- AddSeasonal(ss2, y_ts, nseasons = 7)

modelo_completo <- bsts(
  cnt ~ registered + temp + atemp + hum + windspeed + workingday,
  state.specification = ss2,
  data = day,
  niter = 3000,
  ping = 0,
  seed = 123
)
comparacion <- CompareBstsModels(list(
  Estacional   = modelo_estacional,
  Completo     = modelo_completo
))
comparacion
loglik_est  <- logLik(modelo_estacional)
loglik_comp <- logLik(modelo_completo)
library(bsts)

comparacion <- CompareBstsModels(
  list(
    Estacional = modelo_estacional,
    Completo = modelo_completo
  )
)

print(comparacion)
plot(comparacion)


