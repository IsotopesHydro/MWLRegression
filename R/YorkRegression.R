#' York linear regression function
#'
#'
#' This code implements a function to perform York linear regression which accounts for
#' uncertainties in both y and y variable, with uncertainties also potentially correlated.
#' The function depend on 'ggplot2' for plotting.
#' The function YorkRegression() require the following parameters.
#' 
#' @param Xi Numeric vector containing the x variable measured values
#' @param Yi Numeric vector containing the y variable measured values 
#' @param sXi Numeric vector containing the uncertainties of the x measured values
#' @param sYi Numeric vector containing the uncertainties of the y measured values
#' @param iter Number of iterations to be performed, default is 100
#' @param plot Default False. Will return plots if True of: i) iterativly computed b values against iteration number,
#'                                              ii) Y vs X plot with the OLS regression in blue and York linear regression in red,
#'                                              iii) Chi-Squared density function for the corresponding DF, with the S statistic
#'                                                   as a vertical red line.
#'
#' The functions returns:
#' @return York slope and York intercept of York regression,Xi-xi,Yi-yi,sqrt(Sb2),sqrt(Sa2),S,Q,length(Xi)-2, bVect
#' @return Expected xi and Expected yi - Expected x and x values, 
#' @return Xi-xi,Yi-yi - Residuals of all Xi and Yi measurement respect to the corresponding expected xi and yi,
#' @return sqrt(Sb2),sqrt(Sa2) - Uncertainty on slope and intercept,
#' @return Rsquared - Coefficient of determinaton on the computed linear model,
#' @return MAE_X,MAE_y,MAEavg - Mean Absolute Error on X and Y directions and Average Mean Absolute Error,
#' @return S - the minimized sum of squared weighted residuals,
#' @return Q - probability of obtaining the computed S value or higher for the relative Chi-Squared density function for the corresponding DF,
#' @return length(Xi)-2 - the DF degrees of freedom,
#' @return bVect - Iterations of b value computation.
#' 
#' @export
YorkRegression <- function(Xi, Yi, sXi, sYi, iter = 100, plot = F){ ##### Initialize the function to perform LSE York regression
  
  sigmaEqualityCheck <- sum((sXi[1] != sXi[2:length(sXi)]), (sYi[1] != sYi[2:length(sYi)]))
  warnsigmaEqualityCheck <- character()
  if(sigmaEqualityCheck == 0){
    warnsigmaEqualityCheck <- " ARE EQUAL. Consider the use of RegEqualSigma()"
  }else{
    warnsigmaEqualityCheck <- " DIFFERS. York regression is appropriate"
  }
  
  
  bVect <-c(NA) # Vector to store all estimated slope values
  OLS <- summary(lm(Yi~Xi)) # OLS to obtain the starting slope value
  bOLS <- as.numeric(OLS$coefficients[2,1])  #starting OLS slope value stored for computation and subsequent plotting
  aOLS <- as.numeric(OLS$coefficients[1,1])  #OLS intercept stored for subsequent plotting
  bVect[1] <- bOLS # insert the starting OLS slope value in the first position of the Vector to store all estimated slope values
  
  wXi <- 1/(sXi^2) # Xi weights equal to 1 / Xi error 
  wYi <- 1/(sYi^2) # Yi weights equal to 1 / Yi error 
  
  options(warn=-1)
  ifelse(is.na(cor(sXi, sYi)), yes = ri <- 10^-20, no = ri <- cor(sXi, sYi)) # estimate correlation between errors in X and Y. If the errors are equal for all Xi and Yi respectively,
  options(warn=0)                                                                    # then function corr() gives an 'NA' result, in this case the 'ri' is set = 10^-20
  
  alphai <- sqrt(wXi*wYi) # compute parameter alpha
  
  
  for(i in 1:iter){ # iteration for computing the slope value
    
    Wi <- (wXi*wYi)/(wXi+(bVect[i]^2)*wYi-2*bVect[i]*ri*alphai) # Use the weights, with the value of last computed slope and the correlations ri, and alpha values, to compute Wi for each point
    
    # Use the observed points (Xi ,Yi) and Wi to calculate  ̄X  and  ̄Y , from which Ui and Vi , and hence bi can be evaluated for each point
    X <- (sum(Wi*Xi))/(sum(Wi))
    Y <- (sum(Wi*Yi))/(sum(Wi))
    Ui <- Xi - X
    Vi <- Yi - Y
    BETAi <- Wi*((Ui/wYi)+((bVect[i]*Vi)/wXi)-(bVect[i]*Ui+Vi)*(ri/alphai))
    
    # Use Wi , Ui , Vi , and bi in the expression for b in Eq. ~13b! to calculate an improved estimate of b
    bVect[i+1] <- sum(Wi*BETAi*Vi)/sum(Wi*BETAi*Ui)
  }
  
  b <- bVect[length(bVect)] # final slope value
  a <- Y - b*X # intercept given the final slope value
  xi <- X + BETAi # expectation for the Xi values 
  yi <- Y + b*BETAi # expectation for the Yi values 
  x <- (sum(Wi*xi))/(sum(Wi))
  y <- (sum(Wi*yi))/(sum(Wi))
  ui <- xi - x
  vi <- yi - y
  Sb2 <- (1/sum(Wi*ui^2)) # error on estimated slope 
  Sa2 <- (1/sum(Wi))+(x^2)*Sb2 # error on estimated intercept
  
  S <- sum(Wi*(Yi-b*Xi-a)^2) # sum of weighted squared residuals
  Q <- 1-pchisq(S,length(Xi)-2) # probability 
   
  ResXi <- Xi-xi # residuals
  ResYi <- Yi-yi
  
  MAE_X <- mean(abs(ResXi)) # mean absolute error on X
  MAE_Y <- mean(abs(ResYi)) # mean absolute error on Y
  MAEdiag <- sqrt(MAE_X^2 + MAE_Y^2) # diagonal mean absolute error
  
  SSres <- sum((Yi-yi)^2) # The sum of squares of residuals
  SStot <- sum((Yi-mean(Yi))^2) # The total sum of squares
  Rsquared <- 1 - SSres/SStot # Coefficient of determination
  
  if(plot){ # plotting the Y vs X biplot with OLS and York LSE regression lines 
    
    df <- data.frame(Xi, Yi, sXi, sYi)
    
    reg_lines <- data.frame(
      intercept = c(aOLS, a),
      slope = c(bOLS, b),
      line = c("OLS", "York LSE")
    )
    
    p1 <- ggplot(df, aes(x = Xi, y = Yi)) +
      geom_point() +
      geom_errorbar(aes(ymin = Yi - sYi, ymax = Yi + sYi), width = 0.2, color = "gray", linewidth = 0.8) +
      geom_errorbarh(aes(xmin = Xi - sXi, xmax = Xi + sXi), width = 0.2, color = "gray", linewidth = 0.8) +
      geom_abline(data = reg_lines, aes(intercept = intercept, slope = slope, color = line), linewidth = 0.5) +
      scale_color_manual(values = c("red", "blue")) +
      labs(color = "Regression Lines", x = "Xi", y = "Yi")
    
    df3=data.frame(X=rchisq(1:100000, df=length(Xi)-2))
    p2 <- ggplot(df3,aes(x=X,y=after_stat(density))) + geom_density(fill='blue') +
      geom_vline(aes(xintercept= S), col = "red")
    
    dfbreps <- data.frame(Iterations = c(1:length(bVect)), bValues = bVect)
    p3 <- ggplot(dfbreps, aes(x = Iterations, y = bValues)) + 
      geom_point() + 
      geom_line() +
      labs(
        x = "Iterations",
        y = "Iterated b values")
    
    print(p3)
    print(p1)
    print(p2)
  }
  
  res <- list(warnsigmaEqualityCheck,
              b,a,xi,yi,ResXi,ResYi,
              Rsquared,MAE_X,MAE_Y,MAEavg,
              sqrt(Sb2),sqrt(Sa2),
              S,Q,length(Xi)-2,
              bVect) # results to return
  names(res) <- c("Xi and Yi provided uncertanties:", 
                  "York Slope", "York Intercept", "Expected xi", "Expected xi", "X-residuals", "Y-residuals", 
                  "Rsquared", "Mean Absolute Error on X", "Mean Absolute Error on Y", "Diagonal Mean Absolute Error",
                  "Slope uncertainty", "Intercept uncertainty", "S", "Q", "df",
                  "Iterated slope values")
  return(res)
}



