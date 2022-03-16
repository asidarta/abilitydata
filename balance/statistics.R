

rm(list=ls()) 

require(ez)
require(Rmisc)

# First load the data into a DATA FRAME
mdir <- "C:/Users/ananda.sidarta/Documents/MATLAB/balance/"
mydata <- read.csv(paste(mdir,"balancestatic2.csv",sep=""),header=T,sep=",")

# Convert multiple columns to be factor!!
mydata[,1:4] <- lapply(mydata[,1:4], factor)

ezANOVA( data=mydata
         , dv=area
         , wid=subj
         , between=.(age)
         , within=.(side,task))

summarySE(mydata, "area", c("task"))
summarySE(mydata, "mFreqx", "task")
summarySE(mydata, "p50_Y", "task")
