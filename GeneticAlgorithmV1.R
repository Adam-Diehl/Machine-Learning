#Terminology
  #The fitness function is what we're trying to optimize
  #A particular solution is an individual
  #The data of a solution is the Chromosome
  #The Chromosome is made up of Genes
  #Cross-Over is new individuals created by 'cross breeding' top solutions
  #Elitism is keeping the top X solutions from the previous iteration
  #Mutation is changing one(+) gene(s) in a Chromosome
  #Random immigrants are totally new solutions to the population

#DEVMODE - time the script
codeEnterTime = Sys.time()

#Define the fitness function - this is like an objective function but can include penalty functions and other useful characteristics
  #Our algorithm finds minimums, so we don't need to modify our fitness function
  #This function is the sphere function in 5 dimensions (https://en.wikipedia.org/wiki/Test_functions_for_optimization)
fitnessFunction = function(x, y, z, a, b) {
  #return(abs(x*cos(2*x) + y*cos(2*y) + z*cos(2*z) + a*cos(2*a) + b*cos(2*b)))
  return (x^2 + y^2 + z^2 + a^2 + b^2)
}

#Define parameters
  #Runtime controls - 3 different exit conditions for the optimization routine
  maxIterations = 5000
  maxTime = 170 #170 seconds, 10 seconds short of 3 minutes just for safety's sake
  populationStagnationCap = 1000 #Number of iterations without a new best solution
  burnInPeriod = 100 #Number of iterations to allow before checking stagnation

  #Portfolio parameters - a gene corresponds to a weight for a single stock
  numberOfGenes = 5
  
  #Genetic algorithm parameters
  numberOfIndividuals = 64
  
  numberCrossOvers = (3/16) * numberOfIndividuals
  numberElitists = (1/16) * numberOfIndividuals
  numberMutatedChromosomes = (1/2) * numberOfIndividuals
  numberRandomImmigrants = (1/4) * numberOfIndividuals
  
  mutationRate = 1 #Number of genes to mutate per chromosome
  
#Initialize the algorithm
  #Create initial weight matrix
  weightsMatrixOld = apply(matrix(0, numberOfGenes, numberOfIndividuals), MARGIN = 1, runif)
  
  #Normalize the Weights Matrix
  weightsMatrixOld = sweep(weightsMatrixOld, 1, rowSums(weightsMatrixOld), "/")

  #Instantiate a blank matrix to hold the new solutions
  weightsMatrixNew = matrix(0, numberOfIndividuals, numberOfGenes)
  
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
    scoreVect[i] = fitnessFunction(weightsMatrixOld[i,1], weightsMatrixOld[i,2], weightsMatrixOld[i,3], weightsMatrixOld[i,4], weightsMatrixOld[i,5])
  }
  
  #Bind the score column to the solutions matrix
  weightsMatrixOld = cbind(weightsMatrixOld, scoreVect)
  
  #Sort by score
  weightsMatrixOld = weightsMatrixOld[order(weightsMatrixOld[,numberOfGenes+1]),]
  
  #Record the best score
  solutionFitness = append(solutionFitness, log(1/weightsMatrixOld[1, numberOfGenes+1]))
  
  #Delete the scores
  weightsMatrixOld = weightsMatrixOld[,-(numberOfGenes+1)]
  
  #Carry over the top individuals to the next population
  for(i in 1:numberElitists) {
    weightsMatrixNew[i,] = weightsMatrixOld[i,]
  }
  
  #Cross-Breed the top solutions using two-point cross-over
  for(i in seq(1, numberCrossOvers, 2)) {
    points = sample(1:numberOfGenes,2)
    #Set new genes for Chromosome i
    weightsMatrixNew[(i+numberElitists),] = weightsMatrixOld[i,]
    weightsMatrixNew[(i+numberElitists), points[1]:points[2]] = weightsMatrixOld[i+1,points[1]:points[2]]
    #Set new genes for Chromosome i+1
    weightsMatrixNew[(i+numberElitists+1),] = weightsMatrixOld[i+1,]
    weightsMatrixNew[(i+numberElitists+1), points[1]:points[2]] = weightsMatrixOld[i,points[1]:points[2]]
  }
  
  #Mutate the top half of solutions
  for(i in 1:numberMutatedChromosomes) {
    mutationIndex = sample(1:numberOfGenes, mutationRate) #Randomly select the gene(s) which will be mutated
    weightsMatrixOld[i,mutationIndex] = runif(1,-1,1)
    weightsMatrixNew[(i + numberElitists + numberCrossOvers),] = weightsMatrixOld[i,]
  }

  #Receive random immigrants
  for(i in 1:numberRandomImmigrants) {
    weightsMatrixNew[(i + numberElitists + numberCrossOvers + numberMutatedChromosomes),] = runif(numberOfGenes, -1, 1)
  }
  
  #Assign new matrix to the old matrix
  weightsMatrixOld = weightsMatrixNew
  
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
  
#Algorithm output
  print("Score:")
  print(fitnessFunction(weightsMatrixOld[i,1], weightsMatrixOld[i,2], weightsMatrixOld[i,3], weightsMatrixOld[i,4], weightsMatrixOld[i,5]))
  print("Best solution:")
  print(weightsMatrixOld[1,])
  print("Iterations:")
  print(iterations)
  print("Time elapsed:")
  print(codeExitTime - codeEnterTime)
  print("Solution convergence")
  plot(solutionFitness, type = "l", main = "Convergence of Optimal Solutions", xlab = "Iteration", ylab = "Log Inverse of Fitness")
  