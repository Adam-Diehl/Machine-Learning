#Function to be optimized
fitness_function = function(coordinates) {
  return (abs(sum(coordinates*cos(2*coordinates))))
}