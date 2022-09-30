# Arguments
#args <- commandArgs(trailingOnly=TRUE)
#arg1 <- args[1]
#arg2 <- args[2]
#
# For trouble shooting
arg1 <- "slim_output"
arg2 <- "slim_script.txt"

# slim location
slim_path <- "/usr/local/bin/slim"

# set working directories
slim_directory <- "/Users/anamoralesarae/Dropbox\ (ASU)/slim_outputs"
script_directory <- "/Users/anamoralesarae/Dropbox\ (ASU)/slim_scripts"

run_directory <- paste(slim_directory, arg1, sep = "/")
if (!dir.exists(run_directory)) {
  dir.create(run_directory)
}

# specify the path for input and output
script_path <- paste(script_directory, arg2, sep = "/")

script_text <- readLines(script_path)
code <- data.frame(gsub('filePath="/Users/anamoralesarae/Dropbox (ASU)/slim_outputs/slim_output/ms.txt"', paste0('filePath="', run_directory, '/ms.txt"'), script_text))
names(code) <- "//"

# run slim script function
doOneSLiMRun <- function(seed, script)
{
  system2(slim_path, args = c("-s", seed, shQuote(script)), stdout=T, stderr=T)
}

# run one simulation function
doOneSimulation <- function(seed, script, nsim, nrep)
{
  # collect output from running the script once
  output <- doOneSLiMRun(seed, script)
  
  # Rename every output
  outfile <- paste0(run_directory, "/ms.txt")
  new <- paste0(run_directory, "/ms_simulation_", nsim, "_replicate_", nrep, ".txt")
  file.rename(outfile, new)
  
  header <- "ms 10 1"
  onefile <- new
  tmp <- c(header, readLines(onefile))
  write(tmp, "tmp.msout")
  ms <- readMS("tmp.msout")
  n <- neutrality.stats(ms)
  gn <- get.neutrality(n)
  tajd <- rbind(tajd, gn[[1]])
  colnames(tajd) <- colnames(gn[[1]])
  
  # Sometimes I get some simulations with zero segregating sites so I have to add "na.omit" for R to  omit them from the average
  mean(na.omit(tajd[,1]))
  sd(na.omit(tajd[,1]))
  
  # mean for segregating sites
  mean(na.omit(tajd[,2]))
  sd(na.omit(tajd[,2]))
  
  # calculate summary statistics, this is done by PopGenome package
  ## where is sumstats declared?
  write(sumstats, "masterfile", append = TRUE)
  
  # print message
  print(paste0('simulation #', nsim, ' replicate #', nrep, ' done!'))
}

#should I add here the loop?? so it runs this 1 million of times?

# number of simulations
nsims <- 1000000

seed <- 123
mu <- runif(n = nsims, min = 0.000000001, max = 0.000001)
p <- runif(n = nsims, min = 0, max = 0.2)
bn <- runif(n = nsims, min = 1, max = 100)

#number of replicates per simulation
nreps <- 1000

for (i in c(1:nsims))
{
  codeSim <- mgsub(c("initializeMutationRate(mu)", 'defineConstant("psi",p)', "p1.setSubpopulationSize(bn);}"), c(paste0("initializeMutationRate(", mu[i], ")"), paste0('defineConstant("psi",', p[i], ')'), paste0("p1.setSubpopulationSize(", bn[i], ");}"), code)
  
  script_path_sim <- gsub(".txt", paste0("_sim_", i, ".txt"), script_path)
  
  script_file <- write.csv2(codeSim, file = script_path_sim, quote = FALSE, row.names = FALSE)
  
  for (j in c(1:nreps))
  {
    doOneSimulation(seed, script_path_sim, i, j)
  }
}