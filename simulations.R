library(rugarch)
library(ccgarch)
library(rmgarch)
library(gogarch)
library(MTS)
library(mgarch)
save.image("~/Documents/sheetR/MGARCH-Examples/Simulations.RData")
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



#**************************************************************#
######################### 3-BEKK-GARCH ##########################
#**************************************************************#

bekk.spec = garchSpec(model = list(alpha = 0.2, beta = 0.7))
bekk.A0=matrix(c(0.1,0,0,-0.1,0.2,0,0.3,-0.1,0.4),ncol=3,byrow=T)
bekk.A1=matrix(c(0.7,0,0.1,0,0.2,-0.1,0.2,-0.1,0.4),ncol=3,byrow=T)
bekk.B1=matrix(c(0.2,0,-0.5,0,0,0.2,0,-0.1,0.1),ncol=3,byrow=T)
bekk.S0=diag(3)
bekk.epsilon1=garchSim(bekk.spec, n = 1000)
bekk.epsilon2=garchSim(bekk.spec, n = 1000)
bekk.epsilon3=garchSim(bekk.spec, n = 1000)
bekk.epsilon<-cbind(bekk.epsilon1,bekk.epsilon2,bekk.epsilon3)

bekk.S<-array(0,dim=c(1000,3,3))
bekk.S[1,,]=diag(3)

for(t in 2:1000){  
  bekk.S[t,,]=bekk.A0%*%t(bekk.A0) + bekk.A1 %*%  t(bekk.epsilon[t,]) %*% bekk.epsilon[t,] %*% t(bekk.A1) + bekk.B1%*%bekk.S[t-1,,]
}
MTSplot(bekk.S)

###### 3.1-Simulacao utilizando o pacote mgarch ##### 
bekk.sim <- mvBEKK.sim(series.count = 3, T = 2500)
names(bekk.sim)
# > names(bekk.sim)
# [1] "length"            "series.count"      "order"            
# [4] "params"            "true.params"       "eigenvalues"      
# [7] "uncond.cov.matrix" "white.noise"       "eps"              
# [10] "cor"               "sd" 
bekk.sim$length
bekk.sim$series.count
bekk.sim$order
bekk.sim$params
bekk.sim$true.params
bekk.sim$eigenvalues
bekk.sim$uncond.cov.matrix
matrix.bekk.sim<-cbind(bekk.sim$eps[[1]],bekk.sim$eps[[2]],bekk.sim$eps[[3]])
plot.ts(matrix.bekk.sim)

##### 3.1.a-Estimacao do BEKK via mgarch #####
fit.bekk.mgarch<-mvBEKK.est(matrix.bekk.sim)
##### 3.2-Estimacao do BEKK via MTS #####
system.time(fit.bekk.mts<-BEKK11(matrix.bekk.sim))
names(fit.bekk.mts)
# [1] "estimates"  "HessianMtx" "Sigma.t"

fit.bekk.mts$estimates

# Coefficient(s):
#   Estimate  Std. Error   t value   Pr(>|t|)    
# mu1  -0.04128586  0.03738236  -1.10442 0.26941048    
# mu2   0.03661509  0.04147022   0.88292 0.37727688    
# mu3   0.05237963  0.06317117   0.82917 0.40700829    
# A011  1.08667629          NA        NA         NA    
# A021  0.17363956  0.04793755   3.62220 0.00029210 ***
# A031  0.39239361  0.06920195   5.67027 1.4257e-08 ***
# A022  1.21375768          NA        NA         NA    
# A032 -0.03629572  0.07315632  -0.49614 0.61979618    
# A033  1.84875451          NA        NA         NA    
# A11   0.10000000  0.02602777   3.84205 0.00012201 ***
# A21   0.02000000  0.02938531   0.68061 0.49611697    
# A31   0.02000000  0.04494247   0.44501 0.65631007    
# A12   0.02000000  0.02258799   0.88543 0.37592668    
# A22   0.10000000  0.02645073   3.78061 0.00015644 ***
# A32   0.02000000  0.03941846   0.50738 0.61189063    
# A13   0.02000000  0.01464603   1.36556 0.17207776    
# A23   0.02000000  0.01723125   1.16068 0.24577126    
# A33   0.10000000  0.02909439   3.43709 0.00058800 ***
# B11   0.80000000  0.00212753 376.02264 < 2.22e-16 ***
# B21   0.02000000  0.01159774   1.72447 0.08462227 .  
# B31   0.02000000  0.01685708   1.18645 0.23544657    
# B12   0.02000000  0.00921145   2.17121 0.02991531 *  
# B22   0.80000000  0.00175136 456.78913 < 2.22e-16 ***
# B32   0.02000000  0.01700638   1.17603 0.23958312    
# B13   0.02000000  0.00587621   3.40356 0.00066515 ***
# B23   0.02000000  0.00699098   2.86083 0.00422532 ** 
# B33   0.80000000  0.00285770 279.94580 < 2.22e-16 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# user   system  elapsed 
# 3394.252   37.831 3443.343 