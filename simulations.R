library(rugarch)
library(ccgarch)

library(rmgarch)
library(gogarch)

################################################################
########################## M-GARCH #############################
################################################################
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

##################### CCC-GARCH ######################
### EXEMPLOS USANDO O PACOTE 'ccgarch'

# 'nobs'  a number of observations to be simulated (T)
# 'a'     a vector of constants in the GARCH equation (N×1)
# 'A'     an ARCH parameter matrix in the GARCH equation. A can be a diagonal matrix for the original CCC-GARCH model 
#         or a full matrix for the extended model (N×N)
# 'B'     a GARCH parameter matrix in the GARCH equation. B can be a diagonal matrix for the original CCC-GARCH model 
#         or a full matrix for the extended model (N×N)
# 'R'     a constant conditional correlation matrix (N×N)
# 'd.f'   the degrees of freedom parameter for the t-distribution
# 'model' "diagonal" for the diagonal model and "extended" for the extended model

a1=c(0.05,0.02)
A1=matrix(c(0.5,0,0,0.7),ncol=2)
B1=matrix(c(0.2,0,0,0.1),ncol=2)
R1=matrix(c(1,0.8,0.8,1),ncol=2)
H1<-eccc.sim(nobs=1000, a=a1, A=A1, B=B1, R=R1,  model="diagonal")

#'h'      a matrix of the simulated conditional variances (T × N )
#'eps'    a matrix of the simulated time series with (E)CCC-GARCH process (T × N )

par(mfrow=c(1,2))
plot.ts(H$eps)
plot.ts(H$h)

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
fit<-eccc.estimation(a=w0,A=A0,B=B0,R=R0,model="diagonal",dvar=H1$eps)
#fit$out
fit$para.mat

DCCtest(H1$eps)

H2<-dccsim(fitORspec, n.sim = 1000, n.start = 0, m.sim = 2,
           startMethod = c("unconditional"),Qbar = NULL, Nbar = NULL,
           rseed = NULL, mexsimdata = NULL, vexsimdata = NULL, cluster = NULL,
           VAR.fit = NULL, prerealized = NULL, ...)

#### Conferir processo simulado no ccgarch - ou o processo simulado nao eh o especificado ou o ajuste do ccgarch eh ruim mesmo.
### Simular e ajustar usando o rmgarch para comparar com os resultados do ccgarch

### REplicar este processo para o dccgarch
### SIMULAR PROCESSO COM CARGAS FATORIAIS!!!!!
















##################### DCC-GARCH #################
### EXEMPLOS USANDO O PACOTE 'ccgarch' 
#Args:
# 'nobs'    a number of observations to be simulated (T)
# 'a'       a vector of constants in the GARCH equation (N×1)
# 'A'       an ARCH parameter matrix in the GARCH equation. A can be a diagonal matrix for the original CCC-GARCH model 
#           or a full matrix for the extended model (N×N)
# 'B'       a GARCH parameter matrix in the GARCH equation. B can be a diagonal matrix for the original CCC-GARCH model 
#           or a full matrix for the extended model (N×N)
# 'R'       a constant conditional correlation matrix (N×N)
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


















































##################### GO-GARCH ######################
### EXEMPLOS USANDO O PACOTE 
gogarchsim(fit, n.sim = 1, n.start = 0, m.sim = 1,
           startMethod = c("unconditional", "sample"), prereturns = NA, preresiduals = NA,
           presigma = NA, mexsimdata = NULL, rseed = NULL, cluster = NULL, ...)
Arguments
fit A GO-GARCH fit object of class goGARCHfit.
n.sim The simulation horizon.
n.start The burn-in sample.
m.sim The number of simulations.
startMethod Starting values for the simulation. Valid methods are “unconditional” for the
expected values given the density, and “sample” for the ending values of the
actual data from the fit object.
prereturns Allows the starting return data to be provided by the user.
preresiduals Allows the starting factor residuals to be provided by the user.
presigma Allows the starting conditional factor sigma to be provided by the user.
mexsimdata A list of matrices with the simulated lagged external variables (if any). The list
should be of size m.sim and the matrices each have n.sim + n.start rows.
rseed Optional seeding value(s) for the random number generator.
cluster A cluster object created by calling makeCluster from the parallel package. If
it is not NULL, then this will be used for parallel estimation (remember to stop
                                                                the cluster on completion).