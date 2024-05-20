#install.packages("devtools")
#devtools::install_github("DougLuke/UserNetR")
library(UserNetR)
library(network)
library(igraph)
library(intergraph)
library(RColorBrewer)
data("lhds")
class(lhds)
lhds<-asIgraph(lhds)
summary(lhds)

# Cambiar los valores NA en "NA" para hacer la binarizacion y que no se presenten problemas
V(lhds)$hivscreen[is.na(V(lhds)$hivscreen)] <- "NA"
V(lhds)$nutrition[is.na(V(lhds)$nutrition)] <- "NA"

table(V(lhds)$hivscreen)
table(V(lhds)$nutrition)


hivscreen.numeric<-ifelse(V(lhds)$hivscreen == "Y", 1, 0)
nutrition.numeric<-ifelse(V(lhds)$nutrition == "Y", 1, 0)

V(lhds)$hivscreen.numeric<-hivscreen.numeric
V(lhds)$nutrition.numeric<-nutrition.numeric

table(V(lhds.gc)$nutrition.numeric)

V(lhds)[hivscreen.numeric == 1]$color <- adjustcolor("pink",0.5)
V(lhds)[hivscreen.numeric == 0]$color <- adjustcolor("lightblue",0.5)



 


# Asegurarse de que el atributo hivscreen es de tipo carácter

# Verificar los valores únicos de hivscreen después del cambio
print(unique(V(lhds)$hivscreen))

# Calcular la asortatividad con respecto al atributo hivscreen
# Convertir los atributos a factores para asegurar que la función los maneje correctamente
hivscreen_factor <- as.factor(V(lhds)$hivscreen)

# Calcular la asortatividad
assortativity_nominal(lhds, as.numeric(as.factor(V(lhds)$hivscreen)), directed = FALSE)

clu <- components(lhds)
# grafo de la componente gigante 
lhds.gc <- induced_subgraph(lhds, 
                              clu$membership==which.max(clu$csize))
set.seed(42)
l<-layout_as_tree(lhds.gc)
plot(lhds.gc,vertex.size=10, vertex.label=NA,layout=l)


# cálculo promedio del vecino más cercano
nn.ave <- sapply(V(lhds.gc),
                 function(x) mean(V(lhds.gc)[.nei(x)]$hivscreen.numeric))

nn.pred <- as.numeric(nn.ave > 0.5)

# cálculo de la tasa de error

print(mean(as.numeric(nn.pred != V(lhds.gc)$hivscreen.numeric)))

library(ROCR)


plot(NA, NA, xlim = c(0,1), ylim = c(0,1), xlab = "Tasa Falsos Positivos", ylab = "Tasa verdaderos positivos")
pred <- ROCR::prediction(predictions = nn.ave, labels = V(lhds.gc)$hivscreen.numeric)
perf <- ROCR::performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
# ROC
lines(x = perf@x.values[[1]], y = perf@y.values[[1]], type = "l", col = 4)
# AUC
perf <- ROCR::performance(prediction.obj = pred, measure = "auc")
auc <- perf@y.values[[1]]
text(x = 0.7, y = 0.05, labels = paste0("AUC = ", round(mean(auc), 4)), cex = 1.5)
abline(a = 0, b = 1, col = "gray", lwd = 2)
###################################################################################################
summary(lhds.gc)
table(V(lhds.gc)$hivscreen)
suppressMessages(suppressWarnings(library(ngspatial)))

# se toma el atributos nodal que se desea predecir
X <- V(lhds.gc)$hivscreen.numeric
# matriz de adjacencia
A <- as_adjacency_matrix(lhds.gc, sparse=FALSE)

# se toman caracteristica de la red asociadas con cada nodo
nutrition<- V(lhds.gc)$nutrition.numeric # vectores binarios
# se asigna la correspondiente forma a un modelo mrf que tiene atributos nodales
formula <- X ~ nutrition
formula1 <- X~1

length(X)
length(nutrition)



# modelo autologistico sin atributos nodales
m1.mrf <- autologistic(formula1, A=A,
                       control=list(confint="none"))




# estimaciones 
m1.mrf$coefficients



# vector booleano que marca verdadero si se supera el umbral
mrf1.pred <- as.numeric((m1.mrf$fitted.values > 0.5))



# cálculo de la tasa de error
mean(as.numeric(mrf1.pred != V(lhds.gc)$hivscreen.numeric))







# modelo autologistico con atributos nodales
m2.mrf <- autologistic(X~nutrition, A=A,
                       control=list(confint="none"))





# estimaciones 
m2.mrf$coefficients



# vector booleano que marca verdadero si se supera el umbral
# cálculo de la tasa de error
mrf.pred2 <- as.numeric((m2.mrf$fitted.values > 0.5))
mean(as.numeric(mrf.pred2 != V(lhds.gc)$hivscreen.numeric))




### 2.2 Bondad Y ajuste



set.seed(42)  
# numero de intentos
ntrials <- 100
a1.mrf <- numeric(ntrials)
a2.mrf <- numeric(ntrials)
# intercepto
Z1 <- rep(1,length(X))
# matriz que tiene intercepto y atributos nodales
Z2 <- cbind(Z1, nutrition)
# simulacion
for(i in 1:ntrials){
  # simular red unicamente con intercepto
  X1.mrf <- rautologistic(as.matrix(Z1), A=A,
                          theta=m1.mrf$coefficients)
  # simular red que tiene intercepto y atributos nodales
  X2.mrf<- rautologistic(as.matrix(Z2), A=A,
                         theta=m2.mrf$coefficients)
  #asortatividad 
  a1.mrf[i] <- assortativity(lhds.gc, X1.mrf+1,
                             directed=FALSE)
  a2.mrf[i] <- assortativity(lhds.gc, X2.mrf+1,
                             directed=FALSE)
}





# asortatividad de la red original
assortativity(lhds.gc, X+1, directed=FALSE)



hist(x = a1.mrf, freq = FALSE, col = "gray95", border = "gray95", 
     xlab = 'Asortatividad', ylab = "Frecuencia", 
     main = paste('ppp =', mean(a1.mrf >= assortativity(lhds.gc, X + 1, directed = FALSE))))

abline(v = mean( assortativity(lhds.gc, X + 1, directed = FALSE)), col = 2, lty = 3)




summary(a1.mrf)


hist(x = a2.mrf, freq = FALSE, col = "gray95", border = "gray95", 
     xlab = 'Asortatividad', ylab = "Frecuencia", 
     main = paste('ppp =', mean(a2.mrf >= assortativity(lhds.gc, X + 1, directed = FALSE))))

abline(v = mean( assortativity(lhds.gc, X + 1, directed = FALSE)), col = 2, lty = 3)


plot(NA, NA, xlim = c(0,1), ylim = c(0,1), xlab = "Tasa Falsos Positivos", ylab = "Tasa verdaderos positivos")
pred <- ROCR::prediction(predictions = nn.ave, labels = V(lhds.gc)$hivscreen.numeric)
perf <- ROCR::performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
# ROC
lines(x = perf@x.values[[1]], y = perf@y.values[[1]], type = "l", col = 4)
# AUC
perf <- ROCR::performance(prediction.obj = pred, measure = "auc")
auc <- perf@y.values[[1]]
text(x = 0.7, y = 0.05, labels = paste0("AUC = ", round(mean(auc), 4)), cex = 1.5)
abline(a = 0, b = 1, col = "gray", lwd = 2)

plot(NA, NA, xlim = c(0,1), ylim = c(0,1), xlab = "Tasa Falsos Positivos", ylab = "Tasa verdaderos positivos")
pred <- ROCR::prediction(predictions = m2.mrf$fitted.values, labels = V(lhds.gc)$hivscreen.numeric)
perf <- ROCR::performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
# ROC
lines(x = perf@x.values[[1]], y = perf@y.values[[1]], type = "l", col = 4)
# AUC
perf <- ROCR::performance(prediction.obj = pred, measure = "auc")
auc <- perf@y.values[[1]]
text(x = 0.7, y = 0.05, labels = paste0("AUC = ", round(mean(auc), 4)), cex = 1.5)
abline(a = 0, b = 1, col = "gray", lwd = 2)


plot(NA, NA, xlim = c(0,1), ylim = c(0,1), xlab = "Tasa Falsos Positivos", ylab = "Tasa verdaderos positivos")
pred <- ROCR::prediction(predictions =m1.mrf$fitted.values, labels = V(lhds.gc)$hivscreen.numeric)
perf <- ROCR::performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
# ROC
lines(x = perf@x.values[[1]], y = perf@y.values[[1]], type = "l", col = 4)
# AUC
perf <- ROCR::performance(prediction.obj = pred, measure = "auc")
auc <- perf@y.values[[1]]
text(x = 0.7, y = 0.05, labels = paste0("AUC = ", round(mean(auc), 4)), cex = 1.5)
abline(a = 0, b = 1, col = "gray", lwd = 2)

################################################################################################

suppressMessages(suppressWarnings(library(kernlab)))


# matriz laplaciana
L <- as.matrix(laplacian_matrix(lhds.gc))
# descomposicion propia de la matriz laplaciana
e.L <- eigen(L)
# numero de vertices en la componente gigante
nv <- vcount(lhds.gc)
# como es conectada unicamente el ultimo valor es cero
e.vals <- e.L$values[1:(nv-1)]
# se hace la inversa generalizada 
f.e.vals <- c((e.vals)^(-1), 0)
# ultimo valor propio
e.vals[nv]



# se costruye el kernel
K1.tmp <- e.L$vectors %*% diag(f.e.vals) %*%t(e.L$vectors)
K1 <- as.kernelMatrix(K1.tmp)



# se genera otro kernel con los motifs
K.nutrition <- V(lhds.gc)$nutrition.numeric %*% t(V(lhds.gc)$nutrition.numeric)


# se suman ambos kernel para crear uno solo que tenga en cuenta tanto los motif como la matriz L
K2.tmp <- 0.5 * K1.tmp + 0.5 * K.nutrition
K2 <- as.kernelMatrix(K2.tmp)




# indica el modelo, se utiliza C-svc porque es clasificacion y el tipo de minimizacion que se hace coincida con 
# la propuesta 
summary(lhds.gc)
m1.svm <- ksvm(K1, X, type="C-svc")
# valores ajustados
m1.svm.fitted <- fitted(m1.svm)

# tasa de error 
mean(as.numeric(m1.svm.fitted != V(lhds.gc)$hivscreen.numeric))




# tasa de error 
m2.svm <- ksvm(K2, X, type="C-svc")



m2.svm.fitted <- fitted(m2.svm)
mean(as.numeric(m2.svm.fitted != V(lhds.gc)$hivscreen.numeric))


plot(NA, NA, xlim = c(0,1), ylim = c(0,1), xlab = "Tasa Falsos Positivos", ylab = "Tasa verdaderos positivos")
pred <- ROCR::prediction(predictions =m1.svm.fitted , labels = V(lhds.gc)$hivscreen.numeric)
perf <- ROCR::performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
# ROC
lines(x = perf@x.values[[1]], y = perf@y.values[[1]], type = "l", col = 4)
# AUC
perf <- ROCR::performance(prediction.obj = pred, measure = "auc")
auc <- perf@y.values[[1]]
text(x = 0.7, y = 0.05, labels = paste0("AUC = ", round(mean(auc), 4)), cex = 1.5)
abline(a = 0, b = 1, col = "gray", lwd = 2)


plot(NA, NA, xlim = c(0,1), ylim = c(0,1), xlab = "Tasa Falsos Positivos", ylab = "Tasa verdaderos positivos")
pred <- ROCR::prediction(predictions =m2.svm.fitted , labels = V(lhds.gc)$hivscreen.numeric)
perf <- ROCR::performance(prediction.obj = pred, measure = "tpr", x.measure = "fpr")
# ROC
lines(x = perf@x.values[[1]], y = perf@y.values[[1]], type = "l", col = 4)
# AUC
perf <- ROCR::performance(prediction.obj = pred, measure = "auc")
auc <- perf@y.values[[1]]
text(x = 0.7, y = 0.05, labels = paste0("AUC = ", round(mean(auc), 4)), cex = 1.5)
abline(a = 0, b = 1, col = "gray", lwd = 2)

