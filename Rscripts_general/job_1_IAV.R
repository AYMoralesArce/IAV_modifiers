slim_path <- "/usr/local/bin/slim"
script_path <- "/Users/ana/Dropbox\\ \\(ASU\\)/Output_JensenLab_computer1/virus_sim1/job1_iav.slim"

f0 <- 20
f0df <- data.frame(f0)

f1 <- 40
f1df <- data.frame(f1)

f2 <- 20
f2df <- data.frame(f2)

f3 <- 20
f3df <- data.frame(f3)

master <- cbind(f0df, f1df, f2df, f3df)
row.master <- dim(master)[1]
write.csv(master, "./row_parameters.csv")

for (i in (1:row.master)){
 
for (j in 1:100){
 system2(slim_path, c("-d", paste0("f0=", master$f0[i]), "-d", paste0("f1=",master$f1[i]), "-d", paste0("f2=", master$f2[i]), "-d", paste0("f3=", master$f3[i]),script_path),stdout=T)

 onefile <- "./ms.ms"
 newfilename <- paste0("./job1_row", i, "_", "rep", j, ".ms");
 file.copy(onefile, newfilename)

 secondfile <- "./muts.txt"
 newfilename2 <- paste0("./job1_row", i, "_", "rep", j, ".txt");
 file.copy(secondfile, newfilename2)
 
zipName <- paste0("row", i, ".zip")
msfiles <- dir(path = ".", pattern = ".ms$")
zip(zipName, msfiles)
file.remove(msfiles)

zipName2 <- paste0("row", i, ".zip")
mutfiles2 <- dir(path = ".", pattern = ".txt")
zip(zipName2, mutfiles2)
file.remove(mutfiles2)

print(paste0('replicate #',j,'done!'))

}}