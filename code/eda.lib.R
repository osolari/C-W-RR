panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "blue", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


phenoTypeTest <- function(data){
  
  tTest.df <- data.frame()
  for (i in 1:(dim(data)[2] - 1)){
    
    x <- data[data$fish == "fish", i] 
    y <- data[data$fish == "NoFish", i] 
    
    tTest <- t.test(x = x, y = y)
    
    tTest.df <- rbind(tTest.df, c(tTest$statistic,tTest$p.value) )
  }
  colnames(tTest.df) <- c("statistic", "pVal")
  rownames(tTest.df) <- colnames(data)[1:dim(data)[2]- 1]
  
  
  return(tTest.df)
}