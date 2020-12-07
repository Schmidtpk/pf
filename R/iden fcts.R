
############################ Define the spline identification function

#' @export
spline.median <- function(x,y,stateVariable,theta,model)
{

  e <- y-x

  model_variable<-y
  parameter <- theta

  node <-  c(-2,0,2)

  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }

  a <- cbind( gamm(1,parameter,model_variable),
              gamm(2,parameter,model_variable),
              gamm(3,parameter,model_variable),
              gamm(4,parameter,model_variable))

  b1 <- 2* a[,2] * node[1]-2*a[,1]*node[1]
  c1 <-  a[,2]*node[1]^2-a[,1]*node[1]^2-b1 * node[1]

  b4 <- 2* a[,3] * node[3]-2*a[,4]*node[3]
  c4 <-  a[,3]*node[3]^2-a[,4]*node[3]^2-b4 * node[3]

  ff <- (
    (2* a[,1] * e + b1)       *(e <= node[1])
    + 2* a[,2] * e       *(e > node[1])*(e <= node[2])
    + 2* a[,3] * e       *(e > node[2])*(e <= node[3])
    + (2* a[,4] * e + b4)       *(e > node[3])
  )


  return(ff)
}

#' @export
fixed <- function (x, y)
{
  (y <= x) - stats::pnorm(-0.25)
}

#' @export
fixed_exp <- function (x, y)
{
  ((y <= x) - 0.3508) * abs(x-y)
}

#' @export
fixed_mean <- function (x, y)
{
  ((y <= x) - .5)*abs(x-y)
}

#' @export
fixed_median <- function (x, y)
{
  (y <= x) - .5
}

#' @export
MZ <- function(x,y)
  NULL

#' @export
MZ_PointFore <- function(x, y, stateVariable, theta, model, ...)
{
  (y <= (theta[1]+x)) - stats::pnorm(-0.25)
}

#' @export
probit_constant <- function (stateVariable, theta, ...)
{
  return(stats::pnorm(theta))
}



#' @export
probit_constant2 <- function (stateVariable, theta, ...)
{
  return(theta)
}


#' @export
spline.flex.median <- function(x,y,stateVariable,theta,model)
{

  e <- y-x

  model_variable<-y
  parameter <- theta

  node <-  c(median(e[which(e<0)]),0,median(e[which(e>0)]))

  if(sum(is.na(node))!=0)
  {
    warning("One node is NA. Subistitue -1,0,1.")
    node <- c(-1,0,1)
  }

  a <- cbind( gamm(1,parameter,model_variable),
              gamm(2,parameter,model_variable),
              gamm(3,parameter,model_variable),
              gamm(4,parameter,model_variable))

  b1 <- 2* a[,2] * node[1]-2*a[,1]*node[1]
  c1 <-  a[,2]*node[1]^2-a[,1]*node[1]^2-b1 * node[1]

  b4 <- 2* a[,3] * node[3]-2*a[,4]*node[3]
  c4 <-  a[,3]*node[3]^2-a[,4]*node[3]^2-b4 * node[3]

  ff <- (
    (2* a[,1] * e + b1)       *(e <= node[1])
    + 2* a[,2] * e       *(e > node[1])*(e <= node[2])
    + 2* a[,3] * e       *(e > node[2])*(e <= node[3])
    + (2* a[,4] * e + b4)       *(e > node[3])
  )


  return(ff)
}

Tau <- function(y)
{
  return((1+exp(-y))^(-1))
}


gamm <- function(i,teta,y)
{
  ifelse(i==1 , f <- Tau(teta[1] + teta[2]*y - log (3)),
         ifelse(i==2 , f <- (1- gamm(1,teta,y)) * Tau(teta[3]+teta[4]*y-log(2)) ,
                ifelse(i==3 , f <- (1- gamm(1,teta,y)-gamm(2,teta,y)) * Tau(teta[5]+teta[6]*y-log(1)) ,
                       f <- 1-gamm(1,teta,y)-gamm(2,teta,y)-gamm(3,teta,y) )))
  return(f)
}


#' @export
linearprobit <- function (stateVariable, theta)
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for linear probit model")
  }
  return(pnorm(stateVariable*theta[2] + theta[1]))
}

#' @export
breakprobit <- function (stateVariable, theta)
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for break model")
  }
  return(pnorm((stateVariable>0)*theta[2] + theta[1]))
}

#' @export
periodprobit <- function (stateVariable, theta)
{
  if (length(theta) != 2) {
    stop("Wrong dimension of parameter theta for periodic model")
  }
  return(pnorm((sin(stateVariable*2*pi/4))*theta[2] + theta[1]))
}
