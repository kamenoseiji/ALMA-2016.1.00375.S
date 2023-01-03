gaussConv <- function(vec_x, y_vec, vec_sample, width){
    # IN: vec_x : input vector of x axis
    # IN: vec_y : input vector of y value
    # IN: vec_sample : re-sampling points
    # IN: width : Gaussian convolution width
    outVec <- vec_sample
    for(index in 1:length(vec_sample)){
        x <- vec_sample[index]
        weight <- dnorm(vec_x - x, 0, width)
        weight <- weight / sum(weight)
        outVec[index] <- y_vec %*% weight
    }
    return( outVec )
}
