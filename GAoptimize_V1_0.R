GAoptimize = function(risk_matrix, return_vector) {
  #DEVMODE - time the script
  codeEnterTime = Sys.time()
  
  #Define the fitness function:
  #d = portfolio weights, risk_matrix = omega, return_vector = alpha
  fitnessFunction = function(d, risk_matrix, return_vector){
    return(t(d) %*% risk_matrix %*% d - t(d) %*% return_vector)
  }

  #Set constraints
  BuyInThreshold = 0.01
    
  #Define parameters
  #Runtime controls - 3 different exit conditions for the optimization routine
  maxIterations = 1500
  maxTime = 5 #30 seconds
  populationStagnationCap = 1000 #Number of iterations without a new best solution
  #stagnationReductionFactor = 0.75 #Factor to reduce the stagnation cap when triggered
  burnInPeriod = 100 #Number of iterations to allow before checking stagnation
  
  #Portfolio parameters - a gene corresponds to a weight for a single stock
  numberOfGenes = as.numeric(dim(omega)[2])
  
  #Genetic algorithm parameters
  numberOfIndividuals = 64
  
  numberCrossOvers = (3/16) * numberOfIndividuals
  numberElitists = (1/16) * numberOfIndividuals
  numberMutatedChromosomes = (1/2) * numberOfIndividuals
  numberRandomImmigrants = (1/4) * numberOfIndividuals
  
  mutationRate = 1 #Number of genes to mutate per chromosome
  
  #Initialize the algorithm
  #Create initial weight matrix
  populationOld = apply(matrix(0, numberOfGenes, numberOfIndividuals), MARGIN = 1, runif)
  
  #Normalize the Weights Matrix
  populationOld = sweep(populationOld, 1, rowSums(populationOld), "/")
  
  #Instantiate a blank matrix to hold the new solutions
  populationNew = matrix(0, numberOfIndividuals, numberOfGenes)
  
  #Create a blank list to store the best solutions
  solutionFitness = c(0)
  
  #Record the time entering the optimization procedure
  startTime = Sys.time()
  currentTime = Sys.time()
  
  #Initialize iterations count
  iterations = 0
  
  #Initialize population stagnation index
  populationStagnationIndex = 0
  
  #Algorithm iteration
  #STEPS CORRESPOND TO THE ORDER THE NEW MATRIX IS FILLED
  #Step 1: Elitism
  #Step 2: Cross-Over
  #Step 3: Mutation
  #Step 4: Random Immigrants (totally new solutions)
  
  while(iterations < maxIterations && (currentTime - startTime) < maxTime && populationStagnationIndex < populationStagnationCap) {
    
    #Create a score vector
    scoreVect = matrix(0, numberOfIndividuals, 1)
    
    #Score the entries
    for(i in 1:numberOfIndividuals) {
      scoreVect[i] = fitnessFunction(populationOld[i,], omega, alpha)
    }
    
    #Bind the score column to the solutions matrix
    populationOld = cbind(populationOld, scoreVect)
    
    #Sort by score
    populationOld = populationOld[order(populationOld[,numberOfGenes+1]),]
    
    #Record the best score
    solutionFitness = append(solutionFitness, populationOld[1, numberOfGenes+1])
    
    #Delete the scores
    populationOld = populationOld[,-(numberOfGenes+1)]
    
    #Carry over the top individuals to the next population
    for(i in 1:numberElitists) {
      populationNew[i,] = populationOld[i,]
    }
    
    #Cross-Breed the top solutions using two-point cross-over
    for(i in seq(1, numberCrossOvers, 2)) {
      points = sample(1:numberOfGenes,2)
      #Set new genes for Chromosome i
      populationNew[(i+numberElitists),] = populationOld[i,]
      populationNew[(i+numberElitists), points[1]:points[2]] = populationOld[i+1,points[1]:points[2]]
      #Set new genes for Chromosome i+1
      populationNew[(i+numberElitists+1),] = populationOld[i+1,]
      populationNew[(i+numberElitists+1), points[1]:points[2]] = populationOld[i,points[1]:points[2]]
    }
    
    #Mutate the top half of solutions
    for(i in 1:numberMutatedChromosomes) {
      mutationIndex = sample(1:numberOfGenes, mutationRate) #Randomly select the gene(s) which will be mutated
      populationOld[i,mutationIndex] = runif(1,0,1)
      populationNew[(i + numberElitists + numberCrossOvers),] = populationOld[i,]
    }
    
    #Receive random immigrants
    for(i in 1:numberRandomImmigrants) {
      populationNew[(i + numberElitists + numberCrossOvers + numberMutatedChromosomes),] = runif(numberOfGenes, 0, 1)
    }
    
    #Assign new matrix to the old matrix
    populationOld = populationNew
    
    #Normalize the Weights Matrix
    populationOld = sweep(populationOld, 1, rowSums(populationOld), "/")
    
    #Update exit conditions
    iterations = iterations + 1
    currentTime = Sys.time()
    if(iterations > burnInPeriod && solutionFitness[iterations] == solutionFitness[iterations-1]) {
      populationStagnationIndex = populationStagnationIndex + 1
    } else {
      populationStagnationIndex = 0
    }
  }
  
  #Delete the dummy entry from the solution matrix
  solutionFitness = solutionFitness[-1]
  
  #DEVTOOLS - code exit time
  codeExitTime = Sys.time()
  
  #Notify exit status
  print("Solution convergence")
  
  #Record col names
  colnames(populationOld) = colnames(NASDAQ_No_Date)
  
  #Set buy-in threshold
  for(i in 1:dim(populationOld)[2]) {
    if(populationOld[1,i] < BuyInThreshold) {
      populationOld[1,i] = 0
    }
  }
  
  #Print Non-Zero stock weights
  print("Stock Weights (1% Buy-In Threshold)")
  for(i in 1:dim(populationOld)[2]) {
    if(populationOld[1,i] > 0.01) {
      print(populationOld[1,i])
    }
  }
  
  #Algorithm output
  #print("Score:")
  #print(fitnessFunction(populationOld[1,], omega, alpha))
  #print("Best solution:")
  #print(populationOld[1,])
  print(paste("Iterations:",iterations))
  print(codeExitTime - codeEnterTime)
  return(populationOld[1,])
  #par(mfrow=c(1,2),oma=c(0,0,2,0))
  #plot(populationOld[1,]*100, main = "Stock Weights", xlab = "Stock", ylab = "Weight (Percent)")
  #plot((1/(solutionFitness))-min(1/solutionFitness), type = "l", main = "Solution Convergence", xlab = "Generation", ylab = "Fitness")
  #title("Optimal Portfolio for the NASDAQ 100", outer=TRUE)
}