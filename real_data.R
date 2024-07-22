################################################################################
#               Code to generate real and benchmark data sets                  #
################################################################################
rm(list=ls())

# Coffee data set
{
  require(pgmm)
  data(coffee)
  data <- coffee
  data1 <- data[,-c(1, 2)]
  # calculate labels assignment according to bean variety
  dat <- t(data)
  part1 <- dat[1,]
  part1=as.numeric(part1)
  
  T= dim(data1)[1]
  p= dim(data1)[2]
  Y=t(as.matrix(data1))
  # de-meaning and standardising the variables
  Ym=apply(Y,1,mean)
  Ys=sqrt(apply(Y,1,var))
  Y1=Y
  for (i in 1:p) {
    Y[i,]=(Y[i,]-Ym[i])/Ys[i]
  }
}

# Italian wines data set
{
  require(pgmm)
  data(wine)
  data <- wine
  data1 <- data[,-1]
  # calculate labels assignment according to wine sort
  dat <- t(data)
  part1 <- dat[1,]
  
  T= dim(data1)[1]
  p= dim(data1)[2]
  Y=t(as.matrix(data1))
  # de-meaning and standardising the variables
  Ym=apply(Y,1,mean)
  Ys=sqrt(apply(Y,1,var))
  Y1=Y
  for (i in 1:p) {
    Y[i,]=(Y[i,]-Ym[i])/Ys[i]
  }
  # just de-meaning
  #Ym=apply(Y,1,mean)
  #for (i in 1:p) {
  #Y[i,]=(Y[i,]-Ym[i])
  #}
}

# Olive oils data set
{
  require(FlexDir)
  data <- oliveoil
  data1 <- data[,-c(1, 2)]
  # calculate labels assignment according to areas
  dat <- t(data)
  part1 <- dat[1,]
  
  # calculate labels assignment according to regions
  dat <- t(data)
  part2 <- dat[2,]
  
  T= dim(data1)[1]
  p= dim(data1)[2]
  Y=t(as.matrix(data1))
  # visualise plots
  par(mfrow=c(1,1))
  plot(Y[1,], Y[2,])
  # de-meaning and standardising the variables
  Ym=apply(Y,1,mean)
  Ys=sqrt(apply(Y,1,var))
  Y1=Y
  for (i in 1:p) {
    Y[i,]=(Y[i,]-Ym[i])/Ys[i]
  }
}