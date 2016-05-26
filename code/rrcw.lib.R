require("MASS")
require("glmnet")

predictAllPhenotypes <- function(x, y, xnew, method = "ridge"){
  
  if (method == "ridge"){
    Alpha = 0
  }else{
    Alpha = 1
  }
  
  Yhat <- matrix(0, nrow = dim(x)[1], ncol = dim(y)[2])
  Yhatnew <- matrix(0, nrow = dim(xnew)[1], ncol = dim(y)[2])
  
  for (i in 1:dim(y)[2]){
    glm.fit <- glmnet(x = as.matrix(x), y = y[,i], alpha = Alpha, family = "gaussian", lambda = .1)
    Yhat[,i] <- predict(object = glm.fit, as.matrix(x), type = "link")
    Yhatnew[,i] <- predict(object = glm.fit, as.matrix(xnew), type = "link")
  }
  return(list(Yhattrn=Yhat, Yhattst=Yhatnew))
}

rrCWGCV <- function(y, yhat, r){
  
  y <- as.matrix(y)
  yhat <- as.matrix(yhat)
  
  Q <- t(y) %*% yhat %*% solve( t(y) %*% y ) %*% t(yhat) %*% y %*% solve( t(yhat) %*% yhat)
  TC2T <- eigen(Q)
  C2 <- TC2T$values
  ThatInv <- TC2T$vectors
  D <- ((1-r)*(C2-r))/((1-r)^2*C2 + r^2*(1 - C2)) 
  
  D[D<0] <- 0
  D[D>1] <- 1
  
  TDT <- ThatInv %*% diag(x = D, nrow = length(D)) %*% solve(ThatInv)
  
  return(TDT)
  }


rrCWFCV <- function(x, y, yhat, nfold = 5){
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  yhat <- as.matrix(yhat)
  
  #Q <- t(y) %*% yhat %*% solve( t(y) %*% y ) %*% t(yhat) %*% y %*% solve( t(yhat) %*% yhat)
  #TC2T <- eigen(Q)
  #C2 <- TC2T$values
  #ThatInv <- TC2T$vectors
  
  k <- as.integer(dim(y)[1]/nfold)
  foldLabels <- sample(rep(1:nfold, k)[1:dim(y)[1]])
  
  d <- rep(0,dim(y)[2] )
  for(i in 1:nfold){
    
    IXtstfold <- which(foldLabels ==i)
    IXtrnfold <- setdiff(1:dim(y)[1], IXtstfold)
    
    xtrnfd <- x[IXtrnfold, ]
    ytrnfd <- y[IXtrnfold, ]
    xtstfd <- x[IXtstfold, ]
    ytstfd <- y[IXtstfold, ]
    
    Yfd <- predictAllPhenotypes(x = xtrnfd, y = ytrnfd, xnew = xtstfd)
    
    Qfd <- t(ytstfd) %*% Yfd$Yhattst %*% solve( t(ytstfd) %*% ytstfd ) %*% t(Yfd$Yhattst) %*% ytstfd %*% solve( t(Yfd$Yhattst) %*% Yfd$Yhattst)
    TC2T <- eigen(Qfd)
    C2 <- TC2T$values
    That <- solve(TC2T$vectors)
      
    tHatyTilde <- ytstfd %*% t(That)
    tHatyHat <- Yfd$Yhattst %*% t(That)
    
    for (j in 1:dim(y)[2]){
      d[j] <- d[j] + solve(t(tHatyHat[,j]) %*% tHatyHat[,j]) %*% t(tHatyHat[,j]) %*% tHatyTilde[,j]
    }
      
  }
  
  d <- d/nfold
  
  d[d<0] <- 0
  d[d>1] <- 1
  
  TDT <- solve(That) %*% diag(x = d, nrow = length(d)) %*% That
  
  return(TDT)
}





