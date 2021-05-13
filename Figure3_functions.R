###############################################################################
######## AUXILIARY FUNCTIONS FOR COMPARISON OF OMICS DATA TYPE PROPERTIES #####
########     Tarazona et al. 2021 Nature Computational Sciences           #####
###############################################################################



#### FUNCTION 1: getPower() ####

getPower = function (parameters, power = NULL, n = NULL, fdr = 0.05, alpha = 0.05) {
  
  if (is.null(fdr)) fdr = alpha*(parameters$m - parameters$m1)/(parameters$m1 + alpha*(parameters$m- parameters$m1))
  
  if (is.null(power)) { # Compute power for given n
    if (is.null(n)) stop("Please, indicate a value for either power or n arguments. \n")
    
    if (parameters$type == 1) { # COUNT DATA
      potencia = est_power(n = n, w = parameters$w, rho = parameters$minFC, lambda0 = parameters$meanCounts, 
                           phi0 = parameters$dispersion, f = fdr, m = parameters$m, m1 = parameters$m1)
      
    } 
    
    if (parameters$type == 2) {  # NORMAL DATA
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      potencia = power.t.test(n = n, delta = parameters$delta, sd = parameters$dispersion, 
                              sig.level = myAlpha, type = "two.sample", alternative = "two.sided")$power
      # print(potencia); print(parameters$dispersion); print(n); print(parameters$delta); print(myAlpha)
    }
    
    
    if (parameters$type == 3) {  # BINARY DATA
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      potencia = power.prop.test(n = n, p1 = parameters$delta[1], p2 = parameters$delta[2], 
                                 sig.level = myAlpha, alternative = "two.sided")$power
    }
    
    
    return(potencia)
    
  } else {  # Compute n for given power
    
    if (parameters$type == 1) { # COUNT DATA
      tamany = sample_size(power = power, m = parameters$m, m1 = parameters$m1, f = fdr,
                           k = 1, w = parameters$w, rho = parameters$minFC, lambda0 = parameters$meanCounts, 
                           phi0 = parameters$dispersion)
      
    } 
    
    if (parameters$type == 2) {  # NORMAL DATA
      # Alpha estimation
      # r1 = power * parameters$m1
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      
      tamany = try(power.t.test(power = power, delta = parameters$delta, sd = parameters$dispersion, 
                                sig.level = myAlpha, type = "two.sample", alternative = "two.sided")$n, silent = TRUE)
      if (inherits(tamany, "try-error")) tamany = 2
    }
    
    
    if (parameters$type == 3) {  # BINARY DATA
      
      r1 = parameters$m1
      myAlpha = min(r1 * fdr / ((parameters$m - parameters$m1) * (1 - fdr)), 0.05)  ## NEW
      
      tamany = try(power.prop.test(power = power, p1 = parameters$delta[1], p2 = parameters$delta[2],
                                   sig.level = myAlpha, alternative = "two.sided")$n, silent = TRUE)
      if (inherits(tamany, "try-error")) tamany = 2
    }
    
    return(max(ceiling(tamany), 2))
    
  }
  
}


#### FUNCTION 2: powerPlot() ####

powerPlot = function(parameters, optimalSampleSize, omicCol = NULL) {
  
  if (is.null(omicCol)) {
    if (length(parameters) > 12) {
      stop("Too many omics to be plotted. Please, select a lower number of omics to plot. \n")
    }
    omicCol = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]
    omicCol = omicCol[1:length(parameters)]
  } 
  
  omicShape = 1:length(parameters)
  names(omicCol) = names(omicShape) = names(parameters)
  
  
  ## Power versus Sample Size
  
  # Sample Sizes
  nmax = max(optimalSampleSize$n)
  xMin = 2
  xMax = round(max(nmax+20, (3*nmax - xMin)/2),0)
  xValues = sort(unique(c(round(seq(xMin, xMax, (xMax - xMin)/10)), optimalSampleSize$n)))
  
  # Powers
  yValues = matrix(NA, ncol = length(parameters), nrow = length(xValues))
  rownames(yValues) = xValues
  colnames(yValues) = names(parameters)
  
  for (i in 1:nrow(yValues)) {
    for (j in 1:ncol(yValues)) {
      yValues[i,j] = getPower(parameters[[j]], power = NULL, n = xValues[i], fdr = optimalSampleSize$fdr, alpha = optimalSampleSize$alpha)
    }
  }
  
  matplot(xValues, yValues, type = "l", lwd = 2, xlab = "Sample size", 
          ylab = "Statistical power",
          main = "Power vs Sample Size", col = omicCol, lty = omicShape)
  
  legend("bottomright", names(parameters), lwd = 2, col = omicCol, lty = omicShape, bty = "n")
  
  return(list("PowerVsSsampleSize" = yValues))
}