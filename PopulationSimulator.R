{
initial_allele_freq <- 0.5

Ne_marine <- 50000
Ne_Gosling <- 10000
Ne_Roberts <- 10000
Ne_vector <- c(Ne_marine, Ne_Gosling, Ne_Roberts)

T_Robcolonization <- 2000
T_Goscolonization <- 4000


N_generations <- 12000
N_simulations <- 50
record <- matrix(nrow = N_generations, ncol = 3)
allele_fre_result <- matrix(nrow = N_simulations, ncol = 4)

colnames(allele_fre_result) <- c("Initial","Marine", "Gosling", "Roberts")

initial_allele_freq <- runif(N_simulations, 0,1)
system.time(
for (sim in 1:N_simulations){
  record[1,1] <- initial_allele_freq[sim]
  for(gen in 2:N_generations){
	  # get allele freq next generation in Marine
  	if(gen == T_Robcolonization){record[gen-1,3] <- record[gen-1,1]}
	  if(gen == T_Goscolonization){record[gen-1,2] <- record[gen-1,1]}

  	record[gen,] <- rbinom(3, size = Ne_vector, record[gen-1,])/Ne_vector
  }
  allele_fre_result[sim,] <- append(record[1,1],record[N_generations,])

})
}
time <-1: N_generations
plot(record[,1] ~ time, type = "l", ylim = c(0,1))
lines(x = time, y = record[,2], col = "green")
lines(x = time, y = record[,3], col = "red")


# Calculate Fst