##
## HOMEWORK 03



## loading R packages and environment (not all packages are used in the script) ####
package.list <- c("xts", "quantmod", "lubridate", "glmnet", "PerformanceAnalytics", "dplyr", "tidyr", "lazyeval", "tseries")


if (sum(!is.element(package.list, installed.packages()[,1]))>0) {
        install.packages(package.list[!is.element(package.list, installed.packages()[,1])], dependencies=TRUE)
}

sapply(package.list, require, character.only=TRUE)



rm(list = ls())
options(stringsAsFactors=FALSE)
Sys.setenv (TZ="Asia/Taipei")


##
GetMarketData <- function (src = "yahoo", instruments, Tstart, Tend)
{
        market.data <- new.env()
        
        symbols <- getSymbols (instruments, src=src, env=market.data, auto.assign=TRUE,
                               from=as.Date(Tstart), to=as.Date(Tend))
        
        pairs.list <- mget (ls(market.data), market.data)
        output <- do.call ("cbind", pairs.list)
        if (src != "yahoo")
                output <- output[sprintf("%s::%s", Tstart, Tend)]
        
        output
}




#¡¡Credit (short, intermediate, long, broad)
instruments = c("IVV","VEA","IEF","TLT","TIP","EMB","HYG","GLD","RWX","SHY")
Data2 <- GetMarketData(instruments=instruments, Tstart="2010-01-01", Tend="2019-12-01") ## dataset from Yahoo Finance!
prices0 <- Data2
prices <- Data2[,grep("Adjusted", colnames(prices0))]
colnames(prices) <- gsub(".Adjusted", "", colnames(prices))
Data2 <- prices[,instruments]


Data2.week <- Data2[endpoints(Data2, "weeks"),]

data.all <- cbind(Data2.week)

names(data.all)
## dumping data
write.csv(data.frame(date=as.character(index(data.all )), coredata(data.all )), "GBI1data.csv", row.names=FALSE)
