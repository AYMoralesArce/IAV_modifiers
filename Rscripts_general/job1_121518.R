#rm(list = ls())
library(PopGenome)
slim_path <- "/usr/local/bin/slim"
script_path <- "/Users/anamoralesarae/test.slim"

num=3

mu <- runif(num, min = 1e-9, max = 1e-6)
mudf <- data.frame(mu)

p <- runif(num, min = 0, max = 0.2)
pdf <- data.frame(p)

b <- floor(runif(num, min = 1, max = 100))
bdf <- data.frame(b)

master <- cbind(mudf, pdf, bdf)
row.master <- dim(master)[1]
tajd <- NULL
#seed <- set.seed(seed = NULL)

for (i in (1:row.master)){
 
  #print(i)
  #seed <- sample(1:1000000, 1)
  for (j in 1:1000){
  #system2(slim_path, c("-s", seed, "-d", paste0("mu=", master$mu[i]), "-d", paste0("p=",master$p[i]), "-d", paste0("b=", master$b[i]),script_path),stdout=T)

  system2(slim_path, c("-d", paste0("mu=", master$mu[i]), "-d", paste0("p=",master$p[i]), "-d", paste0("b=", master$b[i]),script_path),stdout=T)
    
  header <- "ms 10 1"
  onefile <- "/Users/anamoralesarae/Desktop/slim_outputs/slim_output/ms.txt"
  tmp <- c(header, readLines(onefile))
  write(tmp, "tmp.msout")
  ms <- readMS("tmp.msout")
  n <- neutrality.stats(ms)
  gn <- get.neutrality(n)
  res <- c(master$mu[i], master$p[i], master$b[i], gn[[1]])
  tajd <- rbind(tajd, res)
  }
  
  if ((i %% 3) == 0) {
    saveRDS(tajd, file = "tajs_out.rds")
  }
}

readRDS("tajs_out.rds")
