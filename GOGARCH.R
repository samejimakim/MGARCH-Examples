library(fGarch)
library(gogarch)
library(rmgarch)

######### SIMULACAO DO GO-GARCH MANUALMENTE ###############

#M=matrix(c(-0.009,0.06,-0.03,0.03,0.03,-0.01,0.04,-0.001,-0.05),ncol=3,byrow=T)

M=matrix(c(-0.4,0.1,-0.3,0.3,0.3,-0.1,0.4,-0.1,-0.5),ncol=3,byrow=T)

# > M
#        [,1]   [,2]  [,3]
# [1,] -0.009  0.060 -0.03
# [2,]  0.030  0.030 -0.01
# [3,]  0.040 -0.001 -0.05

########### MATRIZ DE FATORES ORTOGONAIS: ############### 
##   TRES COMPONENTES ORTOGONAIS AO LONGO DO TEMPO     ##
##   Pacote fGarch para simular os fatores             ##

# GARCH(1,1) - use default omega and specify alpha/beta

gog.spec = garchSpec(model = list(alpha = 0.2, beta = 0.7))
gog.b1<-garchSim(gog.spec, n = 1000)
gog.b2<-garchSim(gog.spec, n = 1000)
gog.b3<-garchSim(gog.spec, n = 1000)

plot(gog.b1)
hist(gog.b1,breaks=30)
qqnorm(gog.b1,pch=20,col=8)
qqline(gog.b1)

gog.bt<-cbind(gog.b1,gog.b2,gog.b3)
gog.rt<-t(M%*%t(bt))

#### Estimacao a partir do pacote GOGARCH ####

system.time(fit1_gogarch<-gogarch(gog.rt*10^3,~garch(1,1),estby="ml"))

# Slots:
#   opt: Object of class "list": List returned by nlminb.
#   estby: Object of class "character": Estimation method.
#   models: Object of class "list": List of univariate GARCH model fits.
#   garchf: Object of class "formula": Garch formula used for uncorrelated component GARCH models.
#   name: Object of class "character": The name of the original data object.
#   X: Object of class "matrix": The data matrix.
#   V: Object of class "matrix": Covariance matrix of X.
#   Z: Object of class "matrix": Transformation matrix.
#   H: Object of class "list": List of conditional variance/covariance matrices.

#   U: Object of class "matrix": Orthogonal matrix.
#   Y: Object of class "matrix": Extra  cted component matrix.

#   P: Object of class "matrix": Left singular values of Var/Cov matrix of X.
#   Dsqr: Object of class "matrix": Square roots of eigenvalues on diagonal, else zero.

temp1<-fit1_gogarch@Z%*%t(fit1_gogarch@Y)
fit1_gogarch@Z%*%t(fit1_gogarch@Z)

t(temp1[,1:10])
at[1:10,]
plot(t(temp1[,1:10]), at[1:10,])


fit1_gogarch@P%*%fit1_gogarch@Dsqr%*%fit1_gogarch@U
fit1_gogarch@Z
fit1_gogarch@U

## APENAS O METODO DE MAXIMA VEROSSIMILHANCA ESTA IMPLEMENTADO ## 

fit1_gogarch
fit1_gogarch@Z # Decomoposicao da matriz V = Z %*% t(Z)
fit1_gogarch@V # Matriz de covariancias de X: cov(X)=V
fit1_gogarch@U
fit1_gogarch@P
fit1_gogarch@Dsqr
fit1_gogarch@Z%*%fit1_gogarch@U


#### Estimacao a partir do pacote RMGARCH ####

gogarchspec2<-gogarchspec(mean.model=list(model="constant"),distribution.model="mvnorm")
gogarchfit(gogarchspec2,gog.rt)
gogarchfit2<-gogarchfit(gogarchspec2,gog.rt*10^3)
coef(gogarchfit2)

nisurface(gogarchfit2,pair=c(1,2))
nisurface(gogarchfit2,pair=c(1,3))
nisurface(gogarchfit2,pair=c(2,3))
rcov(gogarchfit2)
