val <- read.delim(file = '/Users/amoral70/Dropbox (ASU)/Output_cluster_IAV/testvalues/fullmuts.txt', header = FALSE,sep = "\t")
class(val)
nrow(val)
names(val) <-c("base")
head(val)
uniq <- unique(val)
nrow(uniq)

n_occur <- data.frame(table(val$base))
n_occur[n_occur$Freq > 1,]

library(plyr)
z <- count(val, vars=c("base"))
library(dplyr)
ft <- val %>% group_by(base) %>% summarise(freq=n())
ft

library(data.table)
setDT(val)[,freq := .N, by = c("base")]
val[order(freq, decreasing = T),]
setDT(val)[, freq := .N, by = .(base)][order(-freq)]