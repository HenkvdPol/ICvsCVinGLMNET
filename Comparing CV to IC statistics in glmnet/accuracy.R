accuracy <- function(actual, prediction){
  # Variable selection by looking at the Beta accuracy (number of TP) method:
  # %TP = (correctly predicted number of non-zeros) / (total actual number of non-zeros).
  
  actual_number_of_non_zero <- length(which(actual != 0))
  actual_number_of_zero <- length(which(actual == 0))
  
  # Get correct predictions.
  correctly_predicted_non_zero <- length(intersect(which(actual != 0), 
                                                   which(prediction != 0)))
  correctly_predicted_zero <- length(intersect(which(actual == 0), 
                                               which(prediction == 0)))
  
  # TP is correctly predicted non zero elements. So that TN is correctly predicted zero elements.
  TP <- correctly_predicted_non_zero / actual_number_of_non_zero
  TN <- correctly_predicted_zero / actual_number_of_zero
  
  # We also want to see their Total sum.
  T_sum <- (correctly_predicted_non_zero + correctly_predicted_zero) / length(actual)
  
  # Now also want to know the percentage of the prediction values have
  # been shrunk to zero.
  shrunk <- length(which(prediction == 0)) / length(prediction)
  
  # Return result as a list.
  accuracy <- list(TP, TN, T_sum, shrunk)
  names(accuracy) = c('TP', 'TN', 'T_sum', 'shrunk')
  return(accuracy)
}
