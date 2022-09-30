slim_path <- "/usr/local/bin/slim"

#aqui le digo a R a cual directorio enviar los resultados, que es donde estoy trabajando

slim_directory <- "/Users/ana/Dropbox (ASU)/slim_outputs"
script_directory <- "/Users/ana/Dropbox (ASU)/slim_scripts"

run_directory <- paste(slim_directory, arg1, sep = "/")
if (!dir.exists(run_directory))
  dir.create(run_directory)

#aqui le digo a R cual script usar, y le doy la direccion de donde encontrarlo

script_path <- paste(script_directory, arg2, sep = "/")

script_text <- readLines(script_path)
code <- data.frame(gsub('filePath="ms.txt"', paste0('filePath="', run_directory, '/ms.txt"'), script_text))
names(code) <- "//"

script_file <- write.csv2(code, file = script_path, quote = FALSE, row.names = FALSE)
#values
mu <- seq(.000000001,.000001, by=.0000000000001)
p <- seq(0,0.2, by = 0.0000000001)
bn <- runif(max = 100, min = 1, n = 1)

data.frame(mu)



#aqui le digo a R cual funcion buscar

doOneSLiMRun <- function(seed, script)
{
  system2(slim_path, args = c("-s", seed, "-d mu=", mu, "-d psi", p, "-d bn", bn, shQuote(script)), stdout=T, stderr=T)
}


#aqui creo el objeto reps que significa replications y en este caso son 1000
reps <- 1000

#aqui hago un for loop para que haga de 1 a 1000 replications
iter <- 1
for (iter in c(1:reps))
{
  # collect output from running the script once
  output <- doOneSLiMRun(iter, script_path)
  
  #Rename every output
  outfile <- paste0(run_directory, "/ms.txt")
  file.rename(outfile, paste0(run_directory, "/ms_replicate_", iter, ".txt"))
  new <- paste0(run_directory, "/ms_replicate_", iter, ".txt")
  
  header <- "ms 10 1"
  onefile <- new
  tmp <- c(header, readLines(onefile))
  write(tmp, "tmp.msout")
  ms <- readMS("tmp.msout")
  n <- neutrality.stats(ms)
  gn <- get.neutrality(n)
  tajd <- rbind(tajd, gn[[1]])
  
  
  # calculate summary statistics
  write(sumstats, "masterfile", append = TRUE)
  
  # print message
  print(paste0('replicate #', iter, ' done!'))
}