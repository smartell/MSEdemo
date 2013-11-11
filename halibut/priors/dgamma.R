require(MCMCpack) #for inverse gamma distribution
require(pscl)
x <- c(0.04^2  , 0.15^2) # settings for observation error prior
x <- c(0.10^2  , 0.60^2) # settings for processs error prior
y <- c(0.10  , 0.90)

theta <- c(alpha=3.78, beta=0.0102)

fn <- function(theta)
{
	ss = sum( (y - pgamma(1/rev(x),theta[1],theta[2]) )^2 )
}

fit <- optim(theta,fn)

xx  <- seq(0.0,(1.2*max(x))^0.5,length=100)
yy  <- dinvgamma(xx^2,fit$par[1],fit$par[2])
plot(xx^2,yy,type="l")

