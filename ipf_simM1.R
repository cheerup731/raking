
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)


#library(mipfp);library(sampling);library(INLA);library(splines);library(data.table) 

set.seed(108)

N=5000

z1 = sample(c(1:4), size = N, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.4))
z2 = sample(c(1:3), size = N, replace = TRUE, prob = c(0.3, 0.3, 0.4))
x1 = rnorm(N, mean = 5)
x2 = rnorm(N, mean = 5)
p = 1/(1+exp(0.1*z1 + 0.1*z2 - 0.1*x1 - 0.1*x2))
y = rbinom(N,1,p)

pop = data.frame(y = as.factor(y),z1 = as.factor(z1), z2 = as.factor(z2), x1, x2)
  
N.P=c(500, 1000, 1500, 2000); n.P = 500; true.p = table(pop$y)[2]/N
true.table = table(pop[, c(1:3)])/N # target table (true values)

M = 500

scen.fun <- function(i){
  
  n.NP = N.P[i]
  error.a = sq.error.a = error.y = sq.error.y = matrix(nrow = M, ncol = 8)
  
for (m in 1:M) {

# srs probability sample 
indx = sample.int(N, size = n.P, replace = F)
S = rep(F,N)
S[indx] = T
sample.P = pop[S, ]

# nonprobability sample 
f = function(theta) sum(exp(theta+pop$x1+pop$x2) / (1 + exp(theta+pop$x1+pop$x2))) - n.NP
theta = uniroot(f, c(-100, 100))$root
includeNP = exp(theta+pop$x1+pop$x2) / (1 + exp(theta+pop$x1+pop$x2))
Sstar = as.logical(UPrandomsystematic(includeNP))
sample.NP = pop[Sstar,]
d = N/n.P

# B set
B = pop[Sstar+S == 1, ]
B$Z = Sstar[Sstar+S == 1]

# estimate O
glmO = glm(Z ~ x1+x2, data = B, family = "binomial")
O = exp(predict(glmO, newdata = sample.NP))
w = 1+(d-1)/O


#### design based ####

# unweighted seed 
seed.table = table(sample.NP[, c(1:3)]) # checking 0
tgt.list.dims.23 = list(c(2, 3)) # 1=y, 2=z1, 3=z2
tgt.data.23 = list(table(pop$z1, pop$z2)) # known table
ipf.v1 = Ipfp(seed = seed.table, target.list = tgt.list.dims.23, target.data = tgt.data.23)

# weighted seed 
sample.NP = cbind(sample.NP, w)
xx = aggregate(sample.NP$w, list(sample.NP$y, sample.NP$z1, sample.NP$z2), sum)
colnames(xx) = c( "y", "z1", 'z2', 'w')
seed.table.w = merge(as.data.frame(true.table), xx, by = c('z2','z1','y'), all = TRUE )
seed.table.w[is.na(seed.table.w$w), 5] <- 0
seed.table.w = array(seed.table.w[,5], c(2, 4, 3))
ipf.w1 = Ipfp(seed = seed.table.w, target.list = tgt.list.dims.23, target.data = tgt.data.23)

wy = aggregate(sample.NP$w, list(sample.NP$y), sum)$x
tgt.data.w = list(wy, table(pop$z1, pop$z2)) 
tgt.list.dims = list(1, c(2, 3))
ipf.w2 = Ipfp(seed = seed.table.w, target.list = tgt.list.dims, target.data = tgt.data.w)

#unweighted with known y
ipf.v2 = Ipfp(seed = seed.table, target.list = tgt.list.dims, target.data = tgt.data.w)



#### model based ####

# post stratified weights and merge it into number count
sample.NP = sample.NP[, -c(4,5)]
AA = as.data.frame(table(sample.NP))
BB = as.data.frame(table(sample.NP[,c(2,3,4)]))
BB = BB[BB$Freq != 0,]
ABBA = merge(AA, BB, by = c("z1", "z2", "w"))
ABBA = ABBA[ABBA$y == 1, -4]
colnames(ABBA) = c("z1", 'z2', 'w', 'y', 'n')
ABBA$w = as.numeric(ABBA$w)

setDT(ABBA)
CC = ABBA[ ,list(n.z = sum(n), w.z = sum(w)), by= c('z1','z2')]
CC$regular = CC$n.z/CC$w.z
CC = subset(CC, select = c('z1','z2', 'regular'))
ABBA = merge(ABBA, CC, by= c('z1','z2')) 
ABBA$w.reg = ABBA$w*ABBA$regular
ABBA$idnum = 1:nrow(ABBA)

intKnots <- quantile(unique(ABBA$w.reg), seq(0,1,length=(19))[-c(1,19)])
intKnots = as.vector(intKnots)
B.basis = bs(ABBA$w.reg, knots=intKnots, degree=2, Boundary.knots=c(min(ABBA$w.reg),max(ABBA$w.reg)), intercept=TRUE)
P = diff(diag(20), diff = 2)
K = t(P)%*% P
eigOmega <- eigen(K)
indsZ <- 1:18
UZ <- eigOmega$vectors[,indsZ]
LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))     
ZSpline <- B.basis%*%LZ


formula = (y ~ z1 + z2 + w.reg + f(idnum, model="z", Z=ZSpline, initial=5, param=c(1,0.01)))

mod.pred.spline <- inla(formula, family = "binomial", data = ABBA, Ntrials=n, 
                        control.predictor=list(compute=TRUE,link=1),
                        control.compute=list(config = TRUE), num.threads=1, safe = FALSE)                       

ABBA$pred.spline = mod.pred.spline$summary.fitted.values[,"0.5quant"]

# population size predictions
ABBA$x = 1
formula = n ~ -1 +
  f(z1, model="iid", hyper = list(prec = list(initial = log(0.0001),fixed = TRUE))) +
  f(z2, model="iid", hyper = list(prec = list(initial = log(0.0001),fixed = TRUE))) +
  f(idnum, x, model="iid", constr = FALSE, hyper = list(prec = list(initial = log(0.0001),fixed = TRUE)))

mod.N = inla(formula, family = "poisson", data = ABBA, 
             control.predictor=list(compute=TRUE,link=1), control.compute = list(config=TRUE))

ran.age = mod.N$summary.random$z1[,c(1,2)]
colnames(ran.age) = c('z1', 'n.z1')
abba = merge(ABBA, ran.age, by = 'z1')

ran.z2 = mod.N$summary.random$z2[,c(1,2)]
colnames(ran.z2) = c('z2', 'n.z2')
abba = merge(abba, ran.z2, by = 'z2')

abba$est.n = exp(abba$n.z1+abba$n.z2)*abba$w

setDT(abba)
est.sum = abba[ ,list(sum=sum(est.n)), by= c('z1','z2')]
known.N = as.data.frame(tgt.data.23)
colnames(known.N) = c('z1','z2', 'N.z')
est.sum = merge(est.sum, known.N, by = c('z1','z2'))

abba = merge(abba, est.sum, by = c('z1','z2'))
abba$N.hat = abba$N.z * abba$est.n / abba$sum

abba$est = (abba$y + (abba$N.hat-abba$n)*abba$pred.spline)/N
est.inla = abba[ ,list(est =sum(est)), by = c('z1','z2')]
est.inla = merge(est.inla, known.N, by = c('z1','z2'))
est.inla = rbind(est.inla, est.inla)
est.inla$y = as.factor(c(rep(1, 12),rep(0, 12)))
est.inla$est[c(13:24)] = est.inla$N.z[c(1:12)]/N - est.inla$est[c(1:12)]
est.inla = est.inla[,-4]
est.inla = merge(est.inla, true.table, by = c('z1','z2', 'y') )

###

formula = (y ~ z1 + z2)

mod.pred.spline <- inla(formula, family = "binomial", data = ABBA, Ntrials=n, 
                        control.predictor=list(compute=TRUE,link=1),
                        control.compute=list(config = TRUE), num.threads=1, safe = FALSE)                       

abba$pred.spline2 = mod.pred.spline$summary.fitted.values[,"0.5quant"]
est.inla2 = abba[ ,list(pred.spline2 = mean(pred.spline2),y = sum(y), n = sum(n)), by = c('z1','z2')]
est.inla2 = merge(est.inla2, known.N, by = c('z1','z2'))
est.inla2$est = (est.inla2$y + (est.inla2$N.z-est.inla2$n)*est.inla2$pred.spline2)/N
est.inla2 = rbind(est.inla2, est.inla2)
est.inla2$y = as.factor(c(rep(1, 12),rep(0, 12)))
est.inla2$est[c(13:24)] = est.inla2$N.z[c(1:12)]/N - est.inla2$est[c(1:12)]
est.inla2 = merge(est.inla2, true.table, by = c('z1','z2', 'y') )


#### performance ####

error.a[m,] = c(sum((seed.table/nrow(sample.NP) - true.table)^2/true.table),
                sum((seed.table.w/sum(w) - true.table)^2/true.table),
                sum((ipf.v1$p.hat - true.table)^2/true.table),
                sum((ipf.v2$p.hat - true.table)^2/true.table),
                sum((ipf.w1$p.hat - true.table)^2/true.table),
                sum((ipf.w2$p.hat - true.table)^2/true.table),
                sum((est.inla$est - est.inla$N)^2/est.inla$N),
                sum((est.inla2$est - est.inla2$N)^2/est.inla$N))

        
sq.error.a[m,] = c(sum((seed.table/nrow(sample.NP) - true.table)^2),
                   sum((seed.table.w/sum(w) - true.table)^2),
                   sum((ipf.v1$p.hat - true.table)^2),
                   sum((ipf.v2$p.hat - true.table)^2),
                   sum((ipf.w1$p.hat - true.table)^2),
                   sum((ipf.w2$p.hat - true.table)^2),
                   sum((est.inla$est - est.inla$N)^2),
                   sum((est.inla2$est - est.inla2$N)^2))


error.y[m,] = c(sum(seed.table[2, ,])/nrow(sample.NP),
                wy[2]/sum(w),
                sum(ipf.v1$p.hat[2, ,]),  
                sum(ipf.v2$p.hat[2, ,]),  
                sum(ipf.w1$p.hat[2, ,]), 
                sum(ipf.w2$p.hat[2, ,]),
                sum(est.inla$est[est.inla$y == 1]),
                sum(est.inla2$est[est.inla2$y == 1]))  

sq.error.y[m,] = c((sum(seed.table[2, ,])/nrow(sample.NP) - true.p)^2,
                   (wy[2]/sum(w) - true.p)^2,
                   (sum(ipf.v1$p.hat[2, ,]) - true.p)^2,  
                   (sum(ipf.v2$p.hat[2, ,]) - true.p)^2,  
                   (sum(ipf.w1$p.hat[2, ,]) - true.p)^2, 
                   (sum(ipf.w2$p.hat[2, ,]) - true.p)^2,
                   (sum(est.inla$est[est.inla$y == 1])- true.p)^2,
                   (sum(est.inla2$est[est.inla2$y == 1])- true.p)^2)  

}

  return(list(nNP = n.NP, 
              RB.a = colMeans(error.a), 
              MSE.a = colMeans(sq.error.a),
              RB.y  = colMeans(error.y) - true.p, 
              MSE.y = colMeans(sq.error.y)))
  

}

library(doRNG); library(doParallel)

registerDoParallel(cores = min(c(length(N.P), (detectCores()-1))))
registerDoRNG(13)

result <- foreach(i = c(1:length(N.P)),
                  .packages = c('sampling', 'mipfp', 'INLA', 'splines', 'data.table'),
                  .combine = rbind) %dopar% scen.fun(i)

save.image("ipf_simM1.RData")

table = matrix(data = unlist(result[,1]), nrow = length(N.P))
table = cbind(table, matrix(data = unlist(result[,2])*1e3, nrow = length(N.P), byrow = T))
table = cbind(table, matrix(data = unlist(result[,3])*1e3, nrow = length(N.P), byrow = T))
table = cbind(table, matrix(data = unlist(result[,4])*1e3, nrow = length(N.P), byrow = T))
table = cbind(table, matrix(data = unlist(result[,5])*1e3, nrow = length(N.P), byrow = T))


library("xtable")

Latex = xtable(table[, c(1:9)], digits = c(0,0, c(rep(2,8))), caption = 'RB.A')
names(Latex) <-  c('NP','naive', 'weighted', 'ipf.uw', 'ipf.uw+wy',
                   'ipf.w', 'ipf.w+wy','Bays')
Latex
(Latex = xtable(table[, c(1,10:17)], digits = c(0,0, c(rep(2,8))), caption = 'MSE.A'))
(Latex = xtable(table[, c(1,18:25)], digits = c(0,0, c(rep(2,8))), caption = 'RB.Y'))
(Latex = xtable(table[, c(1,26:33)], digits = c(0,0, c(rep(2,8))), caption = 'MSE.Y'))

#t.table = round(true.table*100,2)
#t.table = as.data.frame(t.table)
#(Latex = xtable(cbind(t.table[t.table$y == 1,],t.table[t.table$y == 0,])))
