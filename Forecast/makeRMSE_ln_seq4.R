#=======================================================================
# Loading Packages
#=======================================================================
library(tidyverse)
library(rjags)
library(coda)
library(bayesplot)
library(MCMCvis)
library(runjags, warn.conflicts = F)
library(DescTools)
library(magrittr, warn.conflicts = F)
library(grid)
library(gridExtra, warn.conflicts = F)
library(corrplot)
library(matrixStats)

library(compiler)
compiler::loadcmp("predict_cases_fn4.cmpR")

work_data <- pull_cases('iEN', 'icases', extended=TRUE, fname='canotify_augmented_2022-1-11.tsv') # changes the variable her
work_data[,"cases"] <- log(work_data[,"cases"])

fore.at = c(150, 151)
fore.at = c(150:360)
fore.at = c(150:(360+35))
fore.at = c(150:(360+37))
seq(along=fore.at)

#load("fore_en.RData")
#fore <- fore[,-c(1,9)]
load("SavedResults/fore.noen_ln_10x4.RData")
dim(fore)
fore_noen <- fore[seq(along=fore.at),]
load("SavedResults/fore.en4_ln_10x4.RData")
dim(fore)
fore_en   <- fore[seq(along=fore.at),]
head(fore_en)

cor(work_data$cases[fore.at-7+1], fore_en[,1+1])
cor(work_data$cases[fore.at-7+3], fore_en[,1+3])
cor(work_data$cases[fore.at-7+7], fore_en[,1+7])

a <- cbind(work_data$cases[fore.at], fore_en[seq(along=fore.at),1+1])
a

a <- cbind(work_data$cases[fore.at-7], fore_en[,1+(1:7)])
cor(a, use="complete")

a <- cbind(work_data$cases[fore.at-7], fore_noen[,1+(1:7)])
cor(a, use="complete")



pdf("makeRMSE_ln4_seq.pdf")
plot(  x=fore_en[,1]-7+1, y=work_data$cases[fore.at-7+1],pch=16, xlab="day", ylab="log-cases",
       main="Recorded cases verses EN 1-day forecast")
points(x=fore_en[,1]-7+1, y=fore_en[,1+1],pch=16,col=3)
plot(  x=fore_en[,1]-7+3, y=work_data$cases[fore.at-7+3],pch=16, xlab="day", ylab="log-cases",
       main="Recorded cases verses EN 3-day forecast")
points(x=fore_en[,1]-7+3, y=fore_en[,3+1],pch=16,col=3)
plot(  x=fore_en[,1]-7+7, y=work_data$cases[fore.at-7+7],pch=16, xlab="day", ylab="log-cases",
       main="Recorded cases verses EN 7-day forecast")
points(x=fore_en[,1]-7+7, y=fore_en[,7+1],pch=16,col=3)

plot(  x=fore_noen[,1]-7+1, y=work_data$cases[fore.at-7+1],pch=16, xlab="day", ylab="log-cases",
       main="Recorded cases verses no EN 1-day forecast")
points(x=fore_noen[,1]-7+1, y=fore_noen[,1+1],pch=16,col=3)
plot(  x=fore_noen[,1]-7+3, y=work_data$cases[fore.at-7+3],pch=16, xlab="day", ylab="log-cases",
       main="Recorded cases verses no EN 3-day forecast")
points(x=fore_noen[,1]-7+3, y=fore_noen[,3+1],pch=16,col=3)
plot(  x=fore_noen[,1]-7+7, y=work_data$cases[fore.at-7+7],pch=16, xlab="day", ylab="log-cases",
       main="Recorded cases verses no EN 7-day forecast")
points(x=fore_noen[,1]-7+7, y=fore_noen[,7+1],pch=16,col=3)

mad <- NULL
for(ahead in 1:7){
 a <- c(mean(abs(work_data$cases[fore.at-7+ahead]-fore_en[,1+ahead])),
        mean(abs(work_data$cases[fore.at-7+ahead]-fore_noen[,1+ahead])) )
 print(c(ahead, format(a,3), format(a[1]/a[2],3)))
 mad <- rbind(mad, c(ahead, format(a,digits=3), format(100*(1-a[1]/a[2]),digits=0)))
}
colnames(mad) <- c("Forecast days ahead","MAD EN", "MAD baseline", "Percent Reduction in MAD by EN")

mse <- NULL
for(ahead in 1:7){
 a <- c(mean(((work_data$cases[fore.at-7+ahead]-fore_en[,1+ahead])^2)[201:nrow(fore_en)]),
        mean(((work_data$cases[fore.at-7+ahead]-fore_noen[,1+ahead])^2)[201:nrow(fore_en)]) )
 print(c(ahead, a, a[1]/a[2]))
 mse <- rbind(mse, c(ahead, format(a,digits=3), format(100*(1-a[1]/a[2]),digits=0)))
}
colnames(mse) <- c("Forecast days ahead","MSE EN", "MSE baseline", "Percent Reduction in MSE by EN")
mse

work_data$cases[fore.at-7+ahead][119]
fore_en[,1+ahead][119]
fore_noen[,1+ahead][119]

mse <- NULL
for(ahead in 1:7){
 a <- c((mean((work_data$cases[fore.at-7+ahead]-fore_en[,1+ahead])^2)),
        (mean((work_data$cases[fore.at-7+ahead]-fore_noen[,1+ahead])^2)) )
 print(c(ahead, a, a[1]/a[2]))
 mse <- rbind(mse, c(ahead, format(a,digits=3), format(100*(1-a[1]/a[2]),digits=0)))
}
colnames(mse) <- c("Forecast days ahead","MSE EN", "MSE baseline", "Percent Reduction in MSE by EN")

# Export to pdf
pdf(sprintf("compare_ln4_seq.pdf"))
grid.arrange(tableGrob(mad), tableGrob(mse), top=textGrob("Mean Squared Error of Forecasting with and without EN"))

q()

df_rmse_noen <- df_rmse
cbind(df_rmse_en, df_rmse_noen, df_rmse_en[,2]/df_rmse_noen[,2])
plot(df_rmse_en[,2],type="l",log="y")
lines(df_rmse_noen[,2],col=3)
plot(df_rmse_en[,2]/df_rmse_noen[,2],type="l",log="y")
abline(h=1, lty=2)
mean(df_rmse_en[,2] < df_rmse_noen[,2])
