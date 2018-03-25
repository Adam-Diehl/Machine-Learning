NoSwitchingFitness = numeric(1000)
NoSwitchingIterations = numeric(1000)
NoSwitchingRunTime = numeric(1000)
SwitchingFitness = numeric(1000)
SwitchingIterations = numeric(1000)
SwitchingRunTime = numeric(1000)

for(i in 1:length(NoSwitchingFitness)) {
  print(i)
  x = SSAGA(fitness_function_SPHCOS, 1000, 10, 100, 10)
  NoSwitchingFitness[i] = x$FinalFitness
  NoSwitchingIterations[i] = x$Iterations
  NoSwitchingRunTime[i] = x$RunTime
  y = SSAGA11(fitness_function_SPHCOS, 1000, 10, 100, 10)
  SwitchingFitness[i] = y$FinalFitness
  SwitchingIterations[i] = y$Iterations
  SwitchingRunTime[i] = y$RunTime
}

summary(as.data.frame(NoSwitchingRunTime) - as.data.frame(SwitchingRunTime))
summary((as.data.frame(SwitchingFitness) - as.data.frame(NoSwitchingFitness))/mean(SwitchingFitness))
summary(as.data.frame(NoSwitchingIterations) - as.data.frame(SwitchingIterations))
