library(rugarch)
library(ccgarch)
library(rmgarch)
library(gogarch)
library(MTS)

#**************************************************************#
######################## PENDENCIAS ############################
#**************************************************************#

#### Conferir processo simulado no ccgarch - ou o processo simulado nao eh o especificado ou o ajuste do ccgarch eh ruim mesmo.
#### Rodar o dccfit do pacote ccgarch para ajustar o ccgarch simulado e comparar os resultados
### Simular e ajustar dcc-garch usando o rmgarch para comparar com os resultados do ccgarch

### Replicar este processo para o dccgarch
### SIMULAR PROCESSO COM CARGAS FATORIAIS!!!!!

#**************************************************************#
########################## M-GARCH #############################
#**************************************************************#
### EXEMPLOS DE CCC-GARCH E DCC-GARCH USANDO O PACOTE 'ccgarch'
### CCC-GARCH
### DCCt-GARCH (Tse e Tsui)
### DCCe-GARCH (Engle)
### CP-GARCH
### F-GARCH
### (G)O-GARCH
### (G)ICA-GARCH

h11.spec<-garchSpec(model = list(omega=0.1,alpha=0.05,beta=0.7),cond.dist = "norm")
h22.spec<-garchSpec(model = list(omega=0.2,alpha=0.5,beta=0.5),cond.dist = "norm")
h11<-garchSim(spec = h11.spec, n = 1000)
h22<-garchSim(spec = h22.spec, n = 1000)
R<-matrix(c(1,0.5,-0.5,1),ncol=2,byrow=T)

#**************************************************************#
######################### 1-CCC-GARCH ##########################
#**************************************************************#

########## 1.1-EXEMPLOS USANDO O PACOTE 'ccgarch' ##############

# *eccc.sim*
# Args:
# 'nobs'  a number of observations to be simulated (T)
# 'a'     a vector of constants in the GARCH equation (N×1)
# 'A'     an ARCH parameter matrix in the GARCH equation. A can be a diagonal matrix for the
#         original CCC-GARCH model or a full matrix for the extended model (N×N)
# 'B'     a GARCH parameter matrix in the GARCH equation. B can be a diagonal matrix for the
#         original CCC-GARCH model  or a full matrix for the extended model (N×N)
# 'R'     a constant conditional correlation matrix (N×N)
# 'd.f'   the degrees of freedom parameter for the t-distribution
# 'model' "diagonal" for the diagonal model and "extended" for the extended model

c.a1=c(0.05,0.02)
c.A1=matrix(c(0.1,0,0,0.05),ncol=2)
c.B1=matrix(c(0.6,0,0,0.75),ncol=2)
c.R1=matrix(c(1,0.3,0.3,1),ncol=2)
c.H1<-eccc.sim(nobs=1000, c.a1, c.A1, c.B1, c.R1, d.f=5, model="diagonal")

#'h'      a matrix of the simulated conditional variances (T × N )
#'eps'    a matrix of the simulated time series with (E)CCC-GARCH process (T × N )

par(mfrow=c(1,2))
plot.ts(H1$eps)
plot.ts(H1$h)

# > a1
# [1] 0.05 0.02
# > A1
# [,1] [,2]
# [1,]  0.5  0.0
# [2,]  0.0  0.7
# > B1
# [,1] [,2]
# [1,]  0.2  0.0
# [2,]  0.0  0.1
# > R1
# [,1] [,2]
# [1,]  1.0  0.8
# [2,]  0.8  1.0

# > fit$para.mat
# $a
# [1] 0 0
# 
# $A
# [,1]   [,2]
# [1,] 0.95362 0.0000
# [2,] 0.00000 1.2606
# 
# $B
# [,1]    [,2]
# [1,] 0.54404 0.00000
# [2,] 0.00000 0.44626
# 
# $R
# [,1]    [,2]
# [1,] 1.00000 0.80935
# [2,] 0.80935 1.00000


w0=c(1,1)
A0=diag(2)
B0=diag(2)
R0=diag(2)
fit1<-eccc.estimation(a=w0,A=A0,B=B0,R=R0,model="diagonal",dvar=H1$eps,method="CG" )
fit2<-dcc.estimation(inia=w0,iniA=A0,iniB=B0,ini.dcc=c(0.1,0.1),model="diagonal",dvar=H1$eps,method="CG" ) ### Note que para termos um ccc-garch teriramos que ter ini.dcc=c(0,0), o que não é possivel

names(fit1)
# "out"       "h"         "std.resid" "opt"       "para.mat" 

fit1$opt
# NOTE: $convergence == 0 ## indica que o processo iterativo de estimação convergiu 

fit1$out
fit1$para.mat

#**************************************************************#
######################### 2-DCC-GARCH ##########################
#**************************************************************#

########### 2.1-SIMULACAO USANDO O PACOTE ccgarch ##############

## Usando o pacote ccgarch tambem é possivel simular processos dcc garch. abaixo simularemos um processo dcc-garch com o ccgarch
#'dcc.para' a vector of the DCC parameters (2 × 1) -- corresponds to (theta1,theta2) in 12.33
# 'd.f'     the degrees of freedom parameter for the t-distribution
# 'model'   "diagonal" for the diagonal model and "extended" for the extended model

a <- c(0.003, 0.005, 0.001)
A <- diag(c(0.2,0.3,0.15))
B <- diag(c(0.4, 0.1, 0.2))
R <- matrix(c(1.0, -0.4, 0.3, -0.4, 1.0, -0.12, 0.3, -0.12, 1.0),3,3)
theta <- c(0.05,0.9)
dcc.data.t <- dcc.sim(1000, a, A, B, R, theta, d.f=3, model="diagonal")
plot.ts(dcc.data.t$eps)
plot.ts(dcc.data.t$h)

#Results:
# 'z'      a matrix of random draws from N (0, I). (TxN)
# 'std.z'  a matrix of the standardised residuals. std.z t ∼ N (0, R t ) where R t is the DCC
#          matrix at t. If d.f is set to a finite positive real number, z t ∼ t d.f (0, R t ) (TxN)
# 'dcc'    a matrix of the simulated dynamic conditional correlations (TxN)
# 'h'      a matrix of the simulated conditional variances (TxN)
# 'eps'    a matrix of the simulated time series with (E)CCC-GARCH process (TxN)

########## 2.2-ESTIMACAO USANDO O PACOTE 'ccgarch' ###############

# dcc.sim Arguments:

# 'nobs'     a number of observations to be simulated (T )
# 'a'        a vector of constants in the vector GARCH equation (N × 1)
# 'A'        an ARCH parameter matrix in the vector GARCH equation (N × N )
# 'B'        a GARCH parameter matrix in the vector GARCH equation (N × N )
# 'R'        an unconditional correlation matrix (N × N )
# 'dcc.para' a vector of the DCC parameters (2 × 1)
# 'd.f'      the degrees of freedom parameter for the t-distribution
# 'model'    a character string describing the model. 
#            "diagonal" for the diagonal model  and "extended" for the extended 
#            (full ARCH and GARCH parameter matrices) model

d.a1=c(0.05,0.02)
d.A1=matrix(c(0.1,0,0,0.05),ncol=2)
d.B1=matrix(c(0.6,0,0,0.75),ncol=2)
d.R1=matrix(c(1,0.3,0.3,1),ncol=2)
d.alpha1 = 0.1
d.beta1  = 0.2
d.H1<-dcc.sim(nobs=1000, d.a1, d.A1, d.B1, d.R1, dcc.para=c(d.alpha1,d.beta1), d.f=5, model="diagonal")


# 'dcc'      a matrix of the simulated dynamic conditional correlations (T × N 2 )
# 'h'        a matrix of the simulated conditional variances (T × N )
# 'eps'      a matrix of the simulated time series with DCC-GARCH process (T × N )

plot.ts(d.H1$eps, main = "Processos simulados")
plot.ts(d.H1$h, main="Volatilidade observada nos processos simulados")
plot.ts(d.H1$dcc, main="Dt")

d.w0=c(0.5,0.5)
d.A0=diag(2)
d.B0=diag(2)
d.alpha0 = 0.5
d.beta0  = 0.5
d.fit1<-dcc.estimation(inia=d.w0,iniA=d.A0,iniB=d.B0,ini.dcc=c(0.2,0.2),model="diagonal",dvar=d.H1$eps)

############ 2.3-ESTIMACAO USANDO O PACOTE 'rmgarch' ###############

spec1 = ugarchspec(distribution = "norm")
mspec = multispec(c(spec1,spec1))
fitspec=dccspec(mspec,VAR=TRUE,lag=1,dccOrder=c(1,1),model="DCC",distribution="mvnorm")

fit2=dccfit(fitspec,H2$eps)

############ 2.4-ESTIMACAO USANDO O PACOTE 'MTS' ###############


mtspre.fit3=dccPre(d.H1$eps, include.mean = T, p = 0, cond.dist = "std")
d.fit3=dccFit(mtspre.fit3$sresi, type = "Engle", theta = c(d.alpha0, d.beta0),
              ub = c(0.92, 0.92), lb = c( 1e-04, 1e-04), cond.dist = "std", df = 5)

names(d.fit3$estimates)<-c("theta1","theta2","df.t")
mtspre.fit3$est
d.fit3$estimates

plot.ts(d.fit3$rho.t,cex.lab=0.5,cex.axis=0.8, main = "Dt")














































##################### GO-GARCH ######################
### EXEMPLOS USANDO O PACOTE 
gogarchsim(fit, n.sim = 1, n.start = 0, m.sim = 1,
           startMethod = c("unconditional", "sample"), prereturns = NA, preresiduals = NA,
           presigma = NA, mexsimdata = NULL, rseed = NULL, cluster = NULL, ...)
# Arguments
# fit A GO-GARCH fit object of class goGARCHfit.
# n.sim The simulation horizon.
# n.start The burn-in sample.
# m.sim The number of simulations.
# startMethod Starting values for the simulation. Valid methods are “unconditional” for the
# expected values given the density, and “sample” for the ending values of the
# actual data from the fit object.
# prereturns Allows the starting return data to be provided by the user.
# preresiduals Allows the starting factor residuals to be provided by the user.
# presigma Allows the starting conditional factor sigma to be provided by the user.
# mexsimdata A list of matrices with the simulated lagged external variables (if any). The list
# should be of size m.sim and the matrices each have n.sim + n.start rows.
# rseed Optional seeding value(s) for the random number generator.
# cluster A cluster object created by calling makeCluster from the parallel package. If
# it is not NULL, then this will be used for parallel estimation (remember to stop
#                                                                 the cluster on completion).