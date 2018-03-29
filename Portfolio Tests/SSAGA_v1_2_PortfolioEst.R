SSAGA12_PE = function(fitness_func, max_iterations, max_runtime, max_stagnation, input_width, returns_mtx, risk_mtx, lambda) {
  
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
  
  #Capture the entry time
  CodeEnterTime = Sys.time()
  
  ################ STATE INSTATIATION ################
  
  #Set algorithm parameters
  PopulationSize = 100
  GeneticOperator_Primary = (1/2)
  GeneticOperator_Secondary = (1/4)
  GeneticOperator_Tertiary = (1/5)
  GeneticOperator_Quaternary = (1/20)
  
  #Set states as vector - position 1 = E, position 2 = CO, position 3 = M, position 4 = RI
  StateBias_Elitism = c(PopulationSize*GeneticOperator_Primary, PopulationSize*GeneticOperator_Tertiary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Quaternary)
  StateBias_CrossOver = c(PopulationSize*GeneticOperator_Quaternary, PopulationSize*GeneticOperator_Primary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Tertiary)
  StateBias_Mutation = c(PopulationSize*GeneticOperator_Quaternary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Primary, PopulationSize*GeneticOperator_Tertiary)
  StateBias_RandomImmigrants = c(PopulationSize*GeneticOperator_Quaternary, PopulationSize*GeneticOperator_Secondary, PopulationSize*GeneticOperator_Tertiary, PopulationSize*GeneticOperator_Primary)
  
  #Generate array of states for easy indexing
  StatesList = cbind(StateBias_Elitism, StateBias_CrossOver, StateBias_Mutation, StateBias_RandomImmigrants)
  
  ################ ALGORITHM PARAMETERS ################
  
  #Define parameters
  BurnInPeriod = 5 #Number of iterations to allow before checking stagnation
  
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
  
  #Set initial state & record the index of that state in the State List
  CurrentState = StateBias_RandomImmigrants
  CurrentStateBit = 4
  
  #Number of generations to wait before switching state
  StateSwitchDelay = 5
  
  #Fitness target 
  FitnessTarget = 10
  CurrentFitness = 0
  
  #Array to record average 
  SumDistance = numeric(max_iterations)
  
  ################ ALGORITHM LOOP ################
  
  #Algorithm iteration
  #STEPS CORRESPOND TO THE ORDER THE NEW MATRIX IS FILLED
  #Step 1: Elitism
  #Step 2: Cross-Over
  #Step 3: Mutation
  #Step 4: Random Immigrants (totally new solutions)
  
  #Alert user to entrance
  print("Beginning optimization...")
  
  while(iterations < max_iterations && (CurrentTime - CodeEnterTime) < max_runtime && populationStagnationIndex < max_stagnation && CurrentFitness < FitnessTarget) {
    
    #Create a score vector
    scoreVect = matrix(0, PopulationSize, 1)
    
    #Score the entries
    for(i in 1:PopulationSize) {
      scoreVect[i] = fitness_func(OldPopulation[i,], alpha, omega, lambda)
    }
    
    #Bind the score column to the solutions matrix
    OldPopulation = cbind(OldPopulation, scoreVect)
    
    #Sort by score
    OldPopulation = OldPopulation[order(-OldPopulation[,input_width+1]),]

    #Record the best score (relative)
    SolutionFitness = append(SolutionFitness, log(1/(1+OldPopulation[1, input_width+1])))
    
    #Delete the scores
    OldPopulation = OldPopulation[,-(input_width+1)]
    
    #Calculate the current fitness of the best score (objective)
    CurrentFitness = fitness_func(OldPopulation[1,], alpha, omega, lambda)
    
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
      OldPopulation[i,mutationIndex] = runif(1,0,1)
      NewPopulation[(i + CurrentState[1] + CurrentState[2]),] = OldPopulation[i,]
    }
    
    #Receive random immigrants
    for(i in 1:CurrentState[4]) {
      NewPopulation[(i + CurrentState[1] + CurrentState[2] + CurrentState[3]),] = runif(input_width, 0, 1)
    }
    
    #Assign new matrix to the old matrix and constrain portfolio weights
    OldPopulation = NewPopulation
    OldPopulation = sweep(OldPopulation, 1, rowSums(OldPopulation), "/")
    
    #Update exit conditions
    iterations = iterations + 1
    CurrentTime = Sys.time()
    if(iterations > BurnInPeriod && SolutionFitness[iterations] == SolutionFitness[iterations-1]) {
      populationStagnationIndex = populationStagnationIndex + 1
    } else {
      populationStagnationIndex = 0
    }
    
    #Determine genetic distance 
    Distance = numeric(PopulationSize)
    j = 1
    for(i in 1:(PopulationSize-1)) {
      Distance[j] = sum(abs(OldPopulation[i,] - OldPopulation[i+1,]))
    }
    SumDistance[iterations] = sum(Distance)
    
    #Switch states
    if((populationStagnationIndex %% StateSwitchDelay == 0) && (populationStagnationIndex != 0)) {
      
      #Select new state
      NewStateBit = sample(1:4, 1)
      #Set state based on genetic diversity -- if low genetic diversity, then 50/50 chance state = RI, 25% M, 25% CO
      if (SumDistance < mean(SumDistance) && NewStateBit == 1) {
        CurrentState = StatesList[,4] #Random Immigrants
      } else {
        CurrentState = StatesList[,CurrentStateBit]
      }
      
    }
    
  }
  
  #Delete the dummy entry from the solution matrix
  SolutionFitness = SolutionFitness[-1]
  
  #Algorithm wrap up
  CodeExitTime = Sys.time()
  print("Solution convergence.")
  print(paste("Annual returns: ", (1 + as.numeric((t(OldPopulation[1,]) %*% alpha)))^252 - 1))
  print(paste("Annual std.dev.: ", as.numeric(sqrt((t(OldPopulation[1,]) %*% omega %*% OldPopulation[1,])*sqrt(252)))))
  #plot(SolutionFitness, main = "Solution Fitness", xlab = "Generation", ylab = "Relative Fitness", type = "l", ylim = c(0,25), col = "black")
  #plot(diff(SolutionFitness), main = "Solution Growth", xlab = "Generation", ylab = "Growth in Relative Fitness", type = "l", col = "black")
  print(paste("Fitness:", fitness_func(OldPopulation[1,], alpha, omega, lambda)))
  #print("Best solution:")
  #print(OldPopulation[1,])
  print(paste("Iterations: ", iterations))
  print(paste("Time elapsed: ", CodeExitTime - CodeEnterTime))
  
  #Algorithm output
  ReturnObject = c(0)
  ReturnObject$Solution = OldPopulation[1,]
  ReturnObject$FinalFitness = as.numeric(fitness_func(OldPopulation[1,], alpha, omega, lambda))
  ReturnObject$Iterations = iterations
  ReturnObject$RunTime = CodeExitTime - CodeEnterTime
  ReturnObject$FitnessHistory = SolutionFitness
  ReturnObject$FitnessGrowth = diff(SolutionFitness)
  ReturnObject$GeneticDiversity = SumDistance[1:iterations]
  ReturnObject$Returns = as.numeric(1 + (t(OldPopulation[1,]) %*% alpha)^252 - 1)
  ReturnObject$Risk = as.numeric((t(OldPopulation[1,]) %*% omega %*% OldPopulation[1,])*sqrt(252))
  ReturnObject = ReturnObject[-1]
  return(ReturnObject)
  
}