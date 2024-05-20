
library(igraph)
library(sand)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(update = FALSE)
BiocManager::install(c("GOstats","GO.db"), update = FALSE)
library(GOstats)
library(GO.db)
BiocManager::install("org.Sc.sgd.db", update = FALSE)
library(org.Sc.sgd.db)

library(ngspatial)
library(kernlab)

##################################################################################

set.seed(42)
data(ppi.CC)

summary(ppi.CC)

V(ppi.CC)$ICSC[1:10]

V(ppi.CC)[ICSC == 1]$color <- "yellow"
V(ppi.CC)[ICSC == 0]$color <- "blue"
plot(ppi.CC, vertex.size=5, vertex.label=NA)

clu <- components(ppi.CC)
ppi.CC.gc <- induced_subgraph(ppi.CC, clu$membership==which.max(clu$csize))
nn.ave <- sapply(V(ppi.CC.gc), function(x) mean(V(ppi.CC.gc)[nei(x)]$ICSC))

par(mfrow=c(2,1))
hist(nn.ave[V(ppi.CC.gc)$ICSC == 1], col="yellow",ylim=c(0, 30),
     xlab="Proportion Neighbors w/ ICSC",main="Egos w/ ICSC")
hist(nn.ave[V(ppi.CC.gc)$ICSC == 0], col="blue",ylim=c(0, 30),
     xlab="Proportion Neighbors w/ ICSC",main="Egos w/out ICSC")
par(mfrow=c(1,1))

nn.pred <- as.numeric(nn.ave > 0.5)
mean(as.numeric(nn.pred != V(ppi.CC.gc)$ICSC))

##################################################################################
x <- as.list(org.Sc.sgdGO2ALLORFS)

current.icst <- x[names(x) == "GO:0035556"]
ev.code <- names(current.icst[[1]])
icst.ida <- current.icst[[1]][ev.code == "IDA"]


orig.icsc <- V(ppi.CC.gc)[ICSC == 1]$name

candidates <- intersect(icst.ida, V(ppi.CC.gc)$name)

new.icsc <- setdiff(candidates, orig.icsc)

new.icsc

nn.ave[V(ppi.CC.gc)$name %in% new.icsc]

#################################################################################

X <- V(ppi.CC.gc)$ICSC
A <- as_adjacency_matrix(ppi.CC.gc, sparse=FALSE)

formula1 <- X~1

gene.motifs <- cbind(V(ppi.CC.gc)$IPR000198, V(ppi.CC.gc)$IPR000403, V(ppi.CC.gc)$IPR001806,
                    V(ppi.CC.gc)$IPR001849, V(ppi.CC.gc)$IPR002041, V(ppi.CC.gc)$IPR003527)
formula2 <- X ~ gene.motifs

m1.mrf <- autologistic(formula1, A=A, control=list(confint="none"))

m1.mrf$coefficients

mrf1.pred <- as.numeric((m1.mrf$fitted.values > 0.5))

mean(as.numeric(mrf1.pred != V(ppi.CC.gc)$ICSC))

m1.mrf$fitted.values[V(ppi.CC.gc)$name %in% new.icsc]


m2.mrf <- autologistic(formula2, A=A, control=list(confint="none"))
m2.mrf$coefficients

mrf.pred2 <- as.numeric((m2.mrf$fitted.values > 0.5))
mean(as.numeric(mrf.pred2 != V(ppi.CC.gc)$ICSC))

m2.mrf$fitted.values[V(ppi.CC.gc)$name %in% new.icsc]

################################################################################

set.seed(42)
ntrials <- 100

a1.mrf <- numeric(ntrials)
a2.mrf <- numeric(ntrials)

Z1 <- rep(1,length(X))
Z2 <- cbind(Z1, gene.motifs)

for(i in 1:ntrials){
  X1.mrf <- rautologistic(as.matrix(Z1), A=A,theta=m1.mrf$coefficients)
  X2.mrf<- rautologistic(as.matrix(Z2), A=A,theta=m2.mrf$coefficients)
  a1.mrf[i] <- assortativity(ppi.CC.gc, X1.mrf+1,directed=FALSE)
  a2.mrf[i] <- assortativity(ppi.CC.gc, X2.mrf+1,directed=FALSE)
}

assortativity(ppi.CC.gc, X+1, directed=FALSE)

summary(a1.mrf)
summary(a2.mrf)

#############################################################################

par(mfrow=c(1,1))
L <- as.matrix(laplacian_matrix(ppi.CC.gc))
e.L <- eigen(L)
nv <- vcount(ppi.CC.gc)
e.vals <- e.L$values[1:(nv-1)]
f.e.vals <- c((e.vals)^(-1), 0)
plot(f.e.vals, col="magenta", xlim=c(1, nv),
     xlab=c("Index i"), ylab=expression(f(gamma[i])))

e.vec <- e.L$vectors[, (nv-1)]
v.colors <- character(nv)
v.colors[e.vec >= 0] <- "red"
v.colors[e.vec < 0] <- "blue"
v.size <- 15 * sqrt(abs(e.vec))
l <- layout_with_fr(ppi.CC.gc)
plot(ppi.CC.gc, layout=l, vertex.color=v.colors,
     vertex.size=v.size, vertex.label=NA)

################################################################################

K1.tmp <- e.L$vectors %*% 
  diag(f.e.vals) %*%
  t(e.L$vectors)
K1 <- as.kernelMatrix(K1.tmp)

K.motifs <- gene.motifs %*%
  t(gene.motifs)

K2.tmp <- 0.5 * K1.tmp + 0.5 * K.motifs
K2 <- as.kernelMatrix(K2.tmp)

################################################################################

m1.svm <- ksvm(K1, X, type="C-svc")
m1.svm.fitted <- fitted(m1.svm)
mean(as.numeric(m1.svm.fitted != V(ppi.CC.gc)$ICSC))
m1.svm.fitted[V(ppi.CC.gc)$name %in% new.icsc]

m2.svm <- ksvm(K2, X, type="C-svc")
m2.svm.fitted <- fitted(m2.svm)
mean(as.numeric(m2.svm.fitted != V(ppi.CC.gc)$ICSC))
m2.svm.fitted[V(ppi.CC.gc)$name %in% new.icsc]

################################################################################

set.seed(42)

gl <- list()
gl$ba <- sample_pa(250, m=5, directed=FALSE)
gl$er <- sample_gnm(250, 1250)
gl$ws <- sample_smallworld(1, 250, 5, 0.01)

beta <- 0.5
gamma <- 1

ntrials <- 100

sim <- lapply(gl, sir, beta=beta, gamma=gamma,no.sim=ntrials)
 
plot(sim$er)
plot(sim$ba, color="palegoldenrod",
     median_color="gold", quantile_color="gold")
plot(sim$ws, color="pink", median_color="red",
     quantile_color="red")

x.max <- max(sapply(sapply(sim, time_bins), max))
y.max <- 1.05 * max(sapply(sapply(sim, function(x)
  median(x)[["NI"]]), max, na.rm=TRUE))

plot(time_bins(sim$er), median(sim$er)[["NI"]],
     type="l", lwd=2, col="blue", xlim=c(0, x.max),
     ylim=c(0, y.max), xlab="Time",
     ylab=expression(N[I](t)))
lines(time_bins(sim$ba), median(sim$ba)[["NI"]],
      lwd=2, col="gold")
lines(time_bins(sim$ws), median(sim$ws)[["NI"]],
      lwd=2, col="red")
legend("topright", c("ER", "BA", "WS"),
       col=c("blue", "gold", "red"), lty=1)


