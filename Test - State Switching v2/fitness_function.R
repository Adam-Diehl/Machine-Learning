#Functions to be optimized
fitness_function_SPH = function(coordinates) {
  return (abs(sum(coordinates^2)))
}

fitness_function_SPHCOS = function(coordinates) {
  return (abs(coordinates*cos(2*coordinates)))
}
