SSAGA = function(fitness_func, max_iterations, max_runtime, max_stagnation, input_width) {
  
  ################ TOP LEVEL DOCUMENTATION ################
  #Input documentation
  #Fitness function: objective function (READ IN THE INPUTS)
  #Max iterations: maximum number of generations to allow
  #Max runtime: number of seconds the algorithm can run (runtime will exceed this by the amount it takes to finish the current loop)
  #Max stagnation: number of sequential generations without an improvement in fitness
  #Sample input: SSAGA(fitness_func, 1000, 30, 500, 10)  
  #input_width: how many inputs/dims the input takes 
  
  #Genetic operators
  #Elistism (E) - preserve best solutions across generations
  #Cross-Over/Cross-Breeding (CO) - mix the best solutions
  #Mutation (M) - randomly change the best solutions
  #Random Immigrants (RI) - generate totally new solutions
  
  ################ CODE SETUP ################
  
  #Capture the entry time
  CodeEnterTime = Sys.time()
  
  #Set algorithm parameters
  PopulationSize = 64
  GeneticOperator_Primary = (1/2)
  GeneticOperator_Secondary = (1/4)
  GeneticOperator_Tertiary = (3/16)
  GeneticOperator_Quaternary = (1/16)
  
  #Set states as vector - position 1 = E, position 2 = CO, position 3 = M, position 4 = RI
  StateBias_Elitism = c(PopulationSize*GeneticOperator_Primary, PopulationSize*GeneticOperator_Tertiary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Quaternary)
  StateBias_CrossOver = c(PopulationSize*GeneticOperator_Quaternary, PopulationSize*GeneticOperator_Primary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Tertiary)
  StateBias_Mutation = c(PopulationSize*GeneticOperator_Quaternary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Primary, PopulationSize*GeneticOperator_Tertiary)
  StateBias_RandomImmigrants = c(PopulationSize*GeneticOperator_Quaternary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Tertiary, PopulationSize*GeneticOperator_Primary)
  
  ################ ALGORITHM INTERFACE ################
  
  #Define parameters
  BurnInPeriod = 100 #Number of iterations to allow before checking stagnation
  
  MutationRate = 1 #Number of genes to mutate per chromosome
  
  #Initialize the algorithm
  #Create initial population
  OldPopulation = apply(matrix(0, input_width, PopulationSize), MARGIN = 1, runif)
  
  #Normalize the population
  OldPopulation = sweep(OldPopulation, 1, rowSums(OldPopulation), "/")
  
  #Instantiate a blank matrix to hold the new solutions
  NewPopulation = matrix(0, PopulationSize, input_width)
  
  #Create a blank list to store the best solutions
  SolutionFitness = c(0)
  
  #Record the time entering the optimization procedure
  CurrentTime = Sys.time()
  
  #Initialize iterations count
  iterations = 0
  
  #Initialize population stagnation index
  populationStagnationIndex = 0
  
  #Set initial state
  CurrentState = StateBias_CrossOver
  
  #Algorithm iteration
  #STEPS CORRESPOND TO THE ORDER THE NEW MATRIX IS FILLED
  #Step 1: Elitism
  #Step 2: Cross-Over
  #Step 3: Mutation
  #Step 4: Random Immigrants (totally new solutions)
  
  #Alert user to entrance
  print("Beginning optimization...")
  
  while(iterations < max_iterations && (CurrentTime - CodeEnterTime) < max_runtime && populationStagnationIndex < max_stagnation) {
    
    #Create a score vector
    scoreVect = matrix(0, PopulationSize, 1)
    
    #Score the entries
    for(i in 1:PopulationSize) {
      scoreVect[i] = fitness_func(OldPopulation[i,])
    }
    
    #Bind the score column to the solutions matrix
    OldPopulation = cbind(OldPopulation, scoreVect)
    
    #Sort by score
    OldPopulation = OldPopulation[order(OldPopulation[,input_width+1]),]
    
    #Record the best score
    SolutionFitness = append(SolutionFitness, log(1 + 1/OldPopulation[1, input_width+1]))
    
    #Delete the scores
    OldPopulation = OldPopulation[,-(input_width+1)]
    
    #Carry over the top individuals to the next population
    for(i in 1:CurrentState[1]) {
      NewPopulation[i,] = OldPopulation[i,]
    }
    
    #Cross-Breed the top solutions using two-point cross-over
    for(i in seq(1, CurrentState[2], 2)) {
      points = sample(1:input_width,2)
      #Set new genes for Chromosome i
      NewPopulation[(i+CurrentState[1]),] = OldPopulation[i,]
      NewPopulation[(i+CurrentState[1]), points[1]:points[2]] = OldPopulation[i+1,points[1]:points[2]]
      #Set new genes for Chromosome i+1
      NewPopulation[(i+CurrentState[1]+1),] = OldPopulation[i+1,]
      NewPopulation[(i+CurrentState[1]+1), points[1]:points[2]] = OldPopulation[i,points[1]:points[2]]
    }
    
    #Mutate the top half of solutions
    for(i in 1:CurrentState[3]) {
      mutationIndex = sample(1:input_width, MutationRate) #Randomly select the gene(s) which will be mutated
      OldPopulation[i,mutationIndex] = runif(1,-1,1)
      NewPopulation[(i + CurrentState[1] + CurrentState[2]),] = OldPopulation[i,]
    }
    
    #Receive random immigrants
    for(i in 1:CurrentState[4]) {
      NewPopulation[(i + CurrentState[1] + CurrentState[2] + CurrentState[3]),] = runif(input_width, -1, 1)
    }
    
    #Assign new matrix to the old matrix
    OldPopulation = NewPopulation
    
    #Update exit conditions
    iterations = iterations + 1
    CurrentTime = Sys.time()
    if(iterations > BurnInPeriod && SolutionFitness[iterations] == SolutionFitness[iterations-1]) {
      populationStagnationIndex = populationStagnationIndex + 1
    } else {
      populationStagnationIndex = 0
    }
  }
  
  #Delete the dummy entry from the solution matrix
  SolutionFitness = SolutionFitness[-1]
  
  #Algorithm wrap up
  CodeExitTime = Sys.time()
  print("Solution convergence.")
  #plot(SolutionFitness, main = "Solution Fitness", xlab = "Generation", ylab = "Relative Fitness", type = "l", ylim = c(0,25), col = "black")
  #plot(diff(SolutionFitness), main = "Solution Growth", xlab = "Generation", ylab = "Growth in Relative Fitness", type = "l", col = "black")
  print(paste("Fitness:", fitness_func(OldPopulation[1,])))
  print("Best solution:")
  print(OldPopulation[1,])
  print(paste("Iterations: ", iterations))
  print(paste("Time elapsed: ", CodeExitTime - CodeEnterTime))
  
  #Algorithm output
  ReturnObject = c(0)
  ReturnObject$Solution = OldPopulation[1,]
  ReturnObject$FinalFitness = fitness_func(OldPopulation[1,])
  ReturnObject$Iterations = iterations
  ReturnObject$RunTime = CodeExitTime - CodeEnterTime
  ReturnObject$FitnessHistory = SolutionFitness
  ReturnObject$FitnessGrowth = diff(SolutionFitness)
  ReturnObject = ReturnObject[-1]
  return(ReturnObject)
  
}