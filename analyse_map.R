require(MASS)

setwd("/Users/MAEL/Documents/M2_BI/Genomique/omiques_floobits/achaz_bridier/")

data = read.table("results_file2.txt", header = FALSE)

length(which(data$V2 == -1))
length(which(data$V2 != -1))

length(which(data$V2 != -1)) / (length(which(data$V2 == -1))+length(which(data$V2 != -1)))


truehist(data$V2[which(data$V2 != -1)], nbins = 1000)
