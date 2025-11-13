#' Regression of data with Equal Uncertanties
#' 
#' This code implements a function to perform linear regression on x and y data which as equal measurement 
#' uncertainties. The function depend on 'ggplot2' for plotting.
#' The function RegEqualSigma() require the following parameters.
#' 
#' @param Xi Numeric vector containing the x variable measured values
#' @param Yi Numeric vector containing the y variable measured values 
#' @param sXi Numeric vector containing the uncertainties of the x measured values
#' @param sYi Numeric vector containing the uncertainties of the y measured values
#' @param plot Default False. Will return plots if True of: i) iterativly computed b values against iteration number,
#'                                              ii) Y vs X plot with the OLS regression in blue and York linear regression in red,
#'                                              iii) Chi-Squared density function for the corresponding DF, with the S statistic
#'                                                   as a vertical red line.
#'
#' The functions returns:
#' @return York slope and York intercept of York regression,Xi-xi,Yi-yi,sqrt(Sb2),sqrt(Sa2),S,Q,length(Xi)-2, bVect
#' @return Expected xi and Expected yi - Expected x and x values, 
#' @return ResXi,ResYi - Residuals of all Xi and Yi measurement respect to the corresponding expected xi and yi,
#' @return Rsquared - Coefficient of determinaton on the computed linear model,
#' @return MAE_X,MAE_y,MAEdiag - Mean Absolute Error on X and Y directions and Diagonal Mean Absolute Error,
#' @return sqrt(Sb2),sqrt(Sa2) - Uncertainty on slope and intercept,
#' @return S - the minimized sum of squared weighted residuals,
#' @return Q - probability of obtaining the computed S value or higher for the relative Chi-Squared density function for the corresponding DF,
#' @return length(Xi)-2 - the DF degrees of freedom,
#' 
#' @export
RegEqualSigma <- function(Xi, Yi, sXi, sYi, plot = F){
  
  sigmaEqualityCheck <- sum((sXi[1] != sXi[2:length(sXi)]), (sYi[1] != sYi[2:length(sYi)]))
  warnsigmaEqualityCheck <- character()
  if(sigmaEqualityCheck == 0){
    
    warnsigmaEqualityCheck <- " ARE EQUAL. Equal sigma regression is appropriate."
    
    sXi <- sXi[1]
    sYi <- sYi[1]
    
    # Computation of the terms of the second degree equation b^2+B*b+C=0
    B <- ((-N*sum(Yi^2)+(sum(Yi))^2)*(sXi)^2-(-N*sum(Xi^2)+(sum(Xi))^2)*(sYi)^2)/((N*sum(Xi*Yi)-sum(Xi)*sum(Yi))*(sXi)^2)
    C <- -(sYi/sXi)^2
    
    # Compute Delta of the second degree equation and the slope of the linear regression (b) the positive solution
    delta=B^2-4*C;
    b <- (-B+sqrt(delta))/2 
    
    # Compute the intercept (a), the wheight (W), the minimized sum of squared weighted regression (S), the Q probability 
    a <- (sum(Yi)-b*sum(Xi))/N
    W <- 1/(sYi^2+(b)^2*sXi^2)
    S <- W*sum((Yi-(a+b*Xi))^2)
    Q <- 1-pchisq(S,N-2)
    
    # Use the observed points (Xi ,Yi) and Wi to calculate  ̄X  and  ̄Y , from which Ui and Vi , and hence bi can be evaluated for each point
    X <- (sum(Wi*Xi))/(sum(Wi))
    Y <- (sum(Wi*Yi))/(sum(Wi))
    Ui <- Xi - X
    Vi <- Yi - Y
    BETAi <- Wi*((Ui/wYi)+((bVect[i]*Vi)/wXi)-(bVect[i]*Ui+Vi)*(ri/alphai))
    xi <- X + BETAi # expectation for the Xi values 
    yi <- Y + b*BETAi # expectation for the Yi values 
    x <- (sum(Wi*xi))/(sum(Wi))
    y <- (sum(Wi*yi))/(sum(Wi))
    ui <- xi - x
    vi <- yi - y
    Sb2 <- (1/sum(Wi*ui^2)) # error on estimated slope 
    Sa2 <- (1/sum(Wi))+(x^2)*Sb2 # error on estimated intercept

    Sb <- sqrt(Sb2)
    Sa <- sqrt(Sa2)
    
    ResXi <- Xi-xi # residuals
    ResYi <- Yi-yi
    
    MAE_X <- mean(abs(Xi-xi)) # mean absolute error on X
    MAE_Y <- mean(abs(Yi-yi)) # mean absolute error on Y
    MAEdiag <- sqrt(MAE_X^2 + MAE_Y^2) # diagonal mean absolute error
    
    SSres <- sum((ResYi)^2) # The sum of squares of residuals
    SStot <- sum((Yi-mean(Yi))^2) # The total sum of squares
    
    Rsquared <- 1 - SSres/SStot # Coefficient of determination
    
    res <- list(warnsigmaEqualityCheck,
                b,a,xi,yi,ResXi,ResYi,
                Rsquared,MAE_X,MAE_y,MAEdiag,
                Sb,Sa,S,Q,length(Xi)-2) # results to return
    
  }else{
    warnsigmaEqualityCheck <- " DIFFERS. You MUST use YorkRegression()"
    b <- NULL
    a <- NULL
    xi <- NULL
    yi <- NULL
    ResXi <- NULL
    ResYi  <- NULL
    Rsquared <- NULL
    MAE_X <- NULL
    MAE_y <- NULL
    MAEavg <- NULL
    Sb <- NULL
    Sa <- NULL
    S <- NULL
    Q <- NULL
    res <- list(warnsigmaEqualityCheck,
                b,a,xi,yi,ResXi,ResYi,
                Rsquared,MAE_X,MAE_y,MAEdiag,
                Sb,Sa,S,Q,length(Xi)-2) # results to return
  }
  
  names(res) <- c("Xi and Yi provided uncertanties:", 
                  "Slope", 
                  "Intercept", 
                  "Expected xi", 
                  "Expected xi", 
                  "X-residuals", 
                  "Y-residuals", 
                  "Rsquared",
                  "Mean Absolute Error on X",
                  "Mean Absolute Error on Y",
                  "Average Mean Absolute Error",
                  "Slope uncertainty", 
                  "Intercept uncertainty", 
                  "S", 
                  "Q", 
                  "df")
  return(res)
  
}