install.packages('mlbench')
library(mlbench)
library(corrplot)
data(PimaIndiansDiabetes)
plot1 <- corrplot(cor(PimaIndiansDiabetes[,!(names(PimaIndiansDiabetes) %in% c("diabetes"))]), 
                  method="square",
                  order="hclust", tl.cex=0.7, cl.cex=0.5, tl.col="black", addrect=2)

?dist()
