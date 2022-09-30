####Generate data frame

mu <- runif(10000, min = 0.000000001, max = 0.000001 )
mudf <- data.frame(mu)

p <- runif(10000, min = 0, max = 0.2)
pdf <- data.frame(p)

b <- floor(runif(10000, min = 1, max = 100))
bdf <- data.frame(b)

master <- cbind(mudf, pdf, bdf)

# In case I want to see how mu and others have distributed in the space
plot(master$mu)
plot(pdf$p)
plot(mudf$mu, pdf$p)

View(pdf)
View(mudf)
View(master)


doOneSLiMRun <- function(x)
{
  seed <- x[1]
  mu <- x[2]
  p <- x[3]
  b <- x[4]
  
  system2(slim_path, c("-s", seed, "-d", paste0("mu=", mu), "-d", paste0("p=",p), "-d", paste0("b=", b),script_directory),stdout=T)

  header <- "ms 10 1"
  onefile <- "/Users/anamoralesarae/Desktop/slim_outputs/slim_output/ms.txt"
  tmp <- c(header, readLines(onefile))
  write(tmp, "tmp.msout")
  ms <- readMS("tmp.msout")
  n <- neutrality.stats(ms)
  gn <- get.neutrality(n)
  return(gn[[1]][c(1,2)])
  save(gn, file = "gn_out_temp.Rdata", append= TRUE)
  
}

#to test if runs
doOneSLiMRun(c(100030, 1e-7 , 0.1, 1))


doOneSLiMRun2 <- function(x)
{

  mu <- master[,1]
  p <- master[,2]
  b <- master[,3]
  
  system2(slim_path, c("-s", seed, "-d", paste0("mu=", mu), "-d", paste0("p=",p), "-d", paste0("b=", b),script_directory),stdout=T)
  
  header <- "ms 10 1"
  onefile <- "/Users/anamoralesarae/Desktop/slim_outputs/slim_output/ms.txt"
  tmp <- c(header, readLines(onefile))
  write(tmp, "tmp.msout")
  ms <- readMS("tmp.msout")
  n <- neutrality.stats(ms)
  gn <- get.neutrality(n)
  return(gn[[1]][c(1,2)])
  save(gn, file = "gn_out_temp.Rdata", append= TRUE)
  
}
#lapply?

lapply(master [1,], doOneSLiMRun2)

#?
for (row in 1:master)
{ 
  x= script_directory(master$mu, master$p, master$b)
  doOneSLiMRun(x)
} 

reps <- 1000

  list.append(mean(na.omit(gn2[,1])), sd(na.omit(gn2[,1])), mean(na.omit(gn2[,2])), sd(na.omit(gn2[,2]))) 
  


