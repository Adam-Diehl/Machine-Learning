Elitism = numeric(1000)
for(i in 1:5000) {
  print(i)
  x = SSAGA_LS_EL(fitness_function, 1000, 10, 500, 10)
  for(j in 1:length(x$FitnessGrowth)) {
    Elitism[j] = Elitism[j] + (1/i)*x$FitnessGrowth[j]
  }
}
RandomImmigrants = numeric(1000)
for(i in 1:5000) {
  print(i)
  x = SSAGA_LS_RI(fitness_function, 1000, 10, 500, 10)
  for(j in 1:length(x$FitnessGrowth)) {
    RandomImmigrants[j] = RandomImmigrants[j] + (1/i)*x$FitnessGrowth[j]
  }
}
CrossOver = numeric(1000)
for(i in 1:5000) {
  print(i)
  x = SSAGA_LS_CO(fitness_function, 1000, 10, 500, 10)
  for(j in 1:length(x$FitnessGrowth)) {
    CrossOver[j] = CrossOver[j] + (1/i)*x$FitnessGrowth[j]
  }
}
Mutation = numeric(1000)
for(i in 1:5000) {
  print(i)
  x = SSAGA_LS_MU(fitness_function, 1000, 10, 500, 10)
  for(j in 1:length(x$FitnessGrowth)) {
    Mutation[j] = Mutation[j] + (1/i)*x$FitnessGrowth[j]
  }
}

plot(CrossOver, type = "l", col = "red", ylab = "Fitness", xlab = "Generation", main = "Fitness Gain as a Function of State")
lines(Elitism, col = "black")
lines(Mutation, col = "blue")
lines(RandomImmigrants, col = "green")