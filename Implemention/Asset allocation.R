rm(list = ls())
options(stringsAsFactors=FALSE)
Sys.setenv (TZ="Asia/Taipei")

library(CVXR)
## loading R packages and environment (not all packages are used in the script) ####
## loading R packages and environment (not all packages are used in the script) ####
package.list <- c("xts", "quantmod", "lubridate", "glmnet", "PerformanceAnalytics", "dplyr", "tidyr", "lazyeval",
                  "quadprog", "Rglpk")


if (sum(!is.element(package.list, installed.packages()[,1]))>0) {
        install.packages(package.list[!is.element(package.list, installed.packages()[,1])], dependencies=TRUE)
}

sapply(package.list, require, character.only=TRUE)
##
## IMPLEMENTATION ####
##

##
## CCQP: solves the following non-negative QP problem:
##      min { 1/2 * t(w) %*% (t(X) %*% X + lambda*diag(N)) %*% w - t(w) %*% t(X) %*% y == 1/2*L2(y - X %*% w) + lambda * 1/2*L2(w) }, s.t. w > 0, sum(w)=1
##      where w a [Nx1] vector, R is a [TxN] matrix and y is a [Tx1] vector
##  - constraints is a vector with values {0,1,-1} indicating coefficients that are unconstrained, or constrained to be positive or negative, resp.
##
CCQP <- function (
        y,
        X,
        lambda = 0,
        constraints = NULL
)
{
        #browser()##
        N <- ncol(X)
        T <- nrow(X)
        ones.N <- matrix(1, N)
        zeros.N <- matrix(0, N)
        
        ## setting up the QP
        Dmat <- t(X) %*% X + lambda*diag(N)
        dvec <- t(X) %*% y
        
        ## defining constraints
        if (is.null(constraints))
                constraints <- rep(0, N)
        
        Amat <- diag(constraints)
        Amat <- cbind(matrix (1, nrow=N), Amat)
        
        bvec <- c(1, zeros.N)
        meq <- 1
        
        ## solving the problem
        qp <- solve.QP (Dmat, dvec, Amat, bvec, meq)
        w <- matrix (qp$solution, dimnames=list(colnames(X), NULL))
        
        list(coeff=w)
}





##
## Constrained Maximum Sharpe Ratio Portfolio: QP formulation delivering the portfolio maximizing the SR over a history of sample returns
##
CMSRP = function (
        returns,
        alpha = 0 ## alpha is a covariance regularization parameter
)
{  
        num.assets <- ncol(returns)
        num.samples <- nrow (returns)
        ones.T <- as.matrix (rep (1, num.samples))
        ones.M <- as.matrix (rep (1, num.assets))
        
        mu <- matrix (colMeans (returns))
        Sigma <- cov(returns) + alpha*diag(ncol(returns))
        
        Dmat <- 2*Sigma
        dvec <- matrix (0, nrow=num.assets)
        
        Amat <- cbind (mu, ones.M, diag(num.assets))
        bvec <- matrix (c(1,0,rep(0,num.assets)))
        meq <- 1
        
        ## solving the problem
        qp <- solve.QP (Dmat, dvec, Amat, bvec, meq)
        
        t <- as.numeric(t(ones.M) %*% matrix (qp$solution))
        w <- matrix (qp$solution / t)
        
        list (w=w, qp=qp)
}



##
## CVaROP: Conditional VaR Optimal Portfolio
##
CVaROP = function (
        returns,
        beta = 0.95,
        rmin = 0,
        wmin = 0,
        wmax = 1
)
{
        num.assets <- ncol(returns)
        num.samples <- nrow (returns) ## number of periods/samples/scenarios
        mu <- colMeans (returns)
        
        # creat objective vector, constraint matrix, constraint rhs
        Amat <- rbind (cbind(rbind(1, mu), matrix(data=0, nrow=2, ncol=num.samples+1)),
                       cbind(coredata(returns), diag(num.samples), 1))
        
        obj <- c(rep(0, num.assets), rep(-1/((1-beta)*num.samples), num.samples), -1)
        bvec <- c(1, rmin, rep(0, num.samples))
        
        # direction vector
        dir.vec <- c("==", ">=", rep(">=", num.samples))
        
        # bounds on weights
        bounds <- list (lower=list(ind=1:num.assets, val=rep(wmin, num.assets)),
                        upper=list(ind=1:num.assets, val=rep(wmax, num.assets)))
        
        res <- Rglpk_solve_LP (obj=obj,
                               mat=Amat,
                               dir=dir.vec,
                               rhs=bvec,
                               types=rep("C",length(obj)),
                               max=TRUE,
                               bounds=bounds)
        
        ## NB: res$solution = [w1, w2, ..., wN, d1, d2, ..., dS, RVaR], for N = #assets and S = #scenarios
        w <- matrix(as.numeric(res$solution[1:num.assets])) ## only extract the optimal portfolio
        
        list (w=w, status=list (solver.status=res$status))
}




##
## PerformanceSummary: returns some relevant performance and risk metrics for a series of returns
##
PerformanceSummary <- function (
        rets,
        periods = 252
)
{
        ER <- function (x, periods = periods) { mean(x, na.rm=TRUE)*periods }
        SD <- function (x, periods = periods) { sd(x, na.rm=TRUE)*sqrt(periods) }
        SR <- function (x, periods = periods) { ER(x, periods) / SD(x, periods) }
        
        perf <- rbind ("Mean Ret." = apply (rets, 2, ER, periods=periods),
                       "Std." = apply (rets, 2, SD, periods=periods),
                       "Sharpe" = apply (rets, 2, SR, periods=periods),
                       "MaxDD" = PerformanceAnalytics::maxDrawdown(rets),
                       "ES" = PerformanceAnalytics::ES(rets, p=0.95, method="historical"),
                       "Calmar" = as.vector(PerformanceAnalytics::CalmarRatio(rets)),
                       "Skewness" = PerformanceAnalytics::skewness(rets),
                       "Kurtosis" = PerformanceAnalytics::kurtosis(rets))
        
        round(perf, digits=4)
}




Data.df <- read.csv("E:/Users/lenovo/GBI1data.csv")
Data <- xts(Data.df[,-1], as.Date(Data.df$date)) 
prices <- Data 
## computing returns
rets <- exp(diff(log(prices))) - 1
rets <- rets[-1,]
rets[is.na(rets)] <- 0
Rets <- rets

######################################################################
# Case 1

assetc <- cbind(Rets$IEF,Rets$TLT,Rets$IVV,Rets$GLD ,Rets$EMB)
PerformanceSummary(assetc, periods=52)
cormat<-signif(cor(assetc),2)
cormat



#######################################################################3
# Case 2

library(CVXR)

# create function for MVP
MVP <- function(mu, Sigma, lmd = 0.5) {
        w <- Variable(nrow(Sigma))
        prob <- Problem(Maximize(t(mu) %*% w - lmd*quad_form(w, Sigma)),
                        constraints = list(w >= 0, sum(w) == 1))
        result <- solve(prob)
        w <- as.vector(result$getValue(w))
        names(w) <- colnames(Sigma)
        return(w)
}

# this function can now be used as
# w_MVP <- MVP(mu, Sigma, lmd = 2)
mu.MVP <- colMeans(Rets)
sigma.MVP <- cov(Rets)
w_MVP2 <- MVP(mu.MVP, sigma.MVP, lmd = 2)
w_MVP0.5 <- MVP(mu.MVP, sigma.MVP, lmd = 0.5)
w_MVP4 <- MVP(mu.MVP, sigma.MVP, lmd = 4)
w_MVP8 <- MVP(mu.MVP, sigma.MVP, lmd = 8)
w_MVP16 <- MVP(mu.MVP, sigma.MVP, lmd = 16)
w_MVP100 <- MVP(mu.MVP, sigma.MVP, lmd = 100)

rw_MVP2 <- xts(Rets %*% w_MVP2, index(Rets))
rw_MVP0.5 <- xts(Rets %*% w_MVP0.5, index(Rets))
rw_MVP4 <- xts(Rets %*% w_MVP4, index(Rets))
rw_MVP8 <- xts(Rets %*% w_MVP8, index(Rets))
rw_MVP16 <- xts(Rets %*% w_MVP16, index(Rets))
rw_MVP100 <- xts(Rets %*% w_MVP100, index(Rets))



PerformanceSummary(cbind(rw_MVP0.5,rw_MVP2,rw_MVP4,rw_MVP8,rw_MVP16,rw_MVP100 ), periods=52) ## 52 periods if weekly returns
round(cbind(w_MVP0.5,w_MVP2,w_MVP4,w_MVP8,w_MVP16,w_MVP100 ), digits=2)


#Mean-downside risk portfolio
portfolioDR <- function(X, lmd = 0.5, alpha = 2) {
        T <- nrow(X)
        N <- ncol(X)
        X <- as.matrix(X)
        mu <- colMeans(X)
        w <- Variable(N)
        
        prob <- Problem(Maximize(t(w) %*% mu - (lmd/T) * sum(pos(t(mu) %*% w - X %*% w))^alpha),
                        constraints = list(w >= 0, sum(w) == 1))
        result <- solve(prob)
        return(as.vector(result$getValue(w)))
}


w_DR_alpha3 <- portfolioDR(Rets, alpha = 3)

###############################################
# Risk parity vanilla convex case
# initial point
N<- ncol(Rets)
Sigma<- cov(Rets)
x0 <- rep(1/N, N)

# function definition
fn_convex <- function(x, Sigma) {
        N <- nrow(Sigma)
        return(0.5 * t(x) %*% Sigma %*% x - (1/N)*sum(log(x)))
}

# optimize with general-purpose solver
result <- optim(par = x0, fn = fn_convex, Sigma = Sigma, method = "BFGS")
x_convex <- result$par
w_RPP_convex <- x_convex/sum(x_convex)
rw_RPP <- xts(Rets %*% w_RPP_convex, index(Rets))
rw_DR_alpha3  <- xts(Rets %*% w_DR_alpha3 , index(Rets))
###########################################################
# SP CVar



## (CVaR)
wcvarp <- CVaROP (Rets, beta=0.95, rmin=0, wmin=0, wmax=1)$w
dimnames(wcvarp) <- list(colnames(Rets), "CVaRP")
rcvarp <- xts(Rets %*% wcvarp, index(Rets))


## (MSR)
wmsrp <- CMSRP(Rets)$w
dimnames(wmsrp) <- list(colnames(Rets), "MSRP")
rmsrp <- xts(Rets %*% wmsrp, index(Rets))


## showing the CVaR and MSR optimal portfolios
round(cbind(wcvarp, wmsrp,w_DR_alpha3,w_RPP_convex ), digits=2)

## getting some performance and risk metrics for each optimal portfolio
PerformanceSummary(cbind(rmsrp, rcvarp,rw_DR_alpha3,rw_RPP), periods=52) ## 52 periods if weekly returns

#############################################3
#Inverse volatility portfolio (IVP) 
IVP <- function(Sigma) {
        sigma <- sqrt(diag(Sigma))
        w <- 1/sigma
        w <- w/sum(w)
        return(w)
}

# this function can now be used as
Sigma<- cov(Rets)
w_IVP <- IVP(Sigma)
rw_IVP <- xts(Rets %*% w_IVP, index(Rets))
PerformanceSummary(cbind(rmsrp, rcvarp,rw_DR_alpha3,rw_RPP,rw_IVP ), periods=52) ## 52 periods if weekly returns
round(cbind(wcvarp, wmsrp,w_DR_alpha3,w_RPP_convex,w_IVP ), digits=2)
#########################################################
# create Most diversified portfolio (MDP) 
MSRP <- function(mu, Sigma) {
        w_ <- Variable(nrow(Sigma))
        prob <- Problem(Minimize(quad_form(w_, Sigma)),
                        constraints = list(w_ >= 0, t(mu) %*% w_ == 1))
        result <- solve(prob)
        w <- as.vector(result$getValue(w_)/sum(result$getValue(w_)))
        names(w) <- colnames(Sigma)
        return(w)
}
mu <- colMeans(Rets)
# this function can now be used as
w_MSRP <- MSRP(mu, Sigma)
w_MDP <- MSRP(mu = sqrt(diag(Sigma)), Sigma)
rw_MDP <- xts(Rets %*% w_MDP, index(Rets))
PerformanceSummary(cbind(rmsrp, rcvarp,rw_DR_alpha3,rw_RPP,rw_IVP,rw_MDP  ), periods=52) ## 52 periods if weekly returns
round(cbind(wcvarp, wmsrp,w_DR_alpha3,w_RPP_convex,w_IVP,w_MDP), digits=2)
############################################

# create function for maximum decorrelation portfolio (MDCP)  based on GMVP()

# create function for GMVP
GMVP <- function(Sigma) {
        w <- Variable(nrow(Sigma))
        prob <- Problem(Minimize(quad_form(w, Sigma)), 
                        constraints = list(w >= 0, sum(w) == 1))
        result <- solve(prob)
        w <- as.vector(result$getValue(w))
        names(w) <- colnames(Sigma)
        return(w)
}
MDCP <- function(Sigma) {
        C <- diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
        colnames(C) <- colnames(Sigma)
        return(GMVP(Sigma = C))
}
# this function can now be used as
w_MDCP <- MDCP(Sigma)
rw_MDCP <- xts(Rets %*% w_MDCP, index(Rets))
PerformanceSummary(cbind(rmsrp, rcvarp,rw_DR_alpha3,rw_RPP,rw_IVP,rw_MDP, rw_MDCP ), periods=52) ## 52 periods if weekly returns
round(cbind(wcvarp, wmsrp,w_DR_alpha3,w_RPP_convex,w_IVP,w_MDP,w_MDCP), digits=2)