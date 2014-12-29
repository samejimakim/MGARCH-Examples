library(rmgarch)
data(dji30retw)
library(parallel)
cl = makePSOCKcluster(10)
X = as.matrix(dji30retw)

ugarchspec(mean.model = list(armaOrder = c(1,0)))


uspec.n = multispec(replicate(2, ugarchspec(mean.model = list(armaOrder = c(0,0)))))
# DCC (MVN)
spec.dccn = dccspec(uspec.n, dccOrder = c(1, 1), distribution = 'mvnorm')
# DCC (MVT) with QML first stage
spec.dcct = dccspec(uspec.n, dccOrder = c(1, 1), distribution = 'mvt')
# GO-GARCH (MVN)
spec.ggn = gogarchspec(mean.model = list(model = 'AR', lag = 1), distribution.model = 'mvnorm', ica = 'fastica')
# GO-GARCH (maNIG)
spec.ggg = gogarchspec(mean.model = list(model = 'AR', lag = 1), distribution.model = 'manig', ica = 'fastica')

fit.1 = dccfit(spec.dccn, data = H1$eps, solver = 'solnp', fit.control = list(eval.se = FALSE))


fit.2 = dccfit(spec.dcct, data = X, solver = 'solnp', cluster = cl, fit.control = list(eval.se = FALSE))



fit.3 = gogarchfit(spec.ggn, data = X, solver = 'hybrid', cluster = cl, gfun = 'tanh', maxiter1 = 40000, epsilon = 1e-08, rseed = 100)
fit.4 = gogarchfit(spec.ggg, data = X, solver = 'hybrid', cluster = cl, gfun = 'tanh', maxiter1 = 40000, epsilon = 1e-08, rseed = 100)
