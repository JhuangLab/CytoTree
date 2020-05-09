#################################################################
# This script contains a set of mathematical computing equations
#################################################################

# Mean of positive values
#
mean.pos <- function(x) {
  y <- which(x > 0)
  mean(x[y])
}

# Arithmetic mean of log-transformed values
#
mean.of.logs <- function(x, base=2) {
  log(mean((base ^ x) - 1) + 1, base = base)
}

# Sum of log-transformed values
#
sum.of.logs <- function(x, base=2) {
  log(sum((base ^ x) - 1) + 1, base = base)
}

# Arithmetic mean of positive log-transformed values
#
# Un-logs, takes the arithmetic mean of those values greater than 0, and re-logs
mean.of.logs.pos <- function(x, base=2) {
  y <- x[x>0]
  if (length(y) > 0) {
    return(log(mean((base ^ y) - 1) + 1, base = base))
  } else {
    return(0)
  }
}

# Return the signed value with greater magnitude
#
# Determines which value has greater magnitude (i.e. by absolute value) and returns the original signed value.
# Compares vectors x and y pairwise.
pmax.abs <- function(x, y) {
  z <- y
  x.bigger <- (abs(x) > abs(y))
  z[x.bigger] <- x [x.bigger]
  return(z)
}

# Returns proportion of expressing cells
#
prop.exp <- function(x) {
  return(length(which(x>0)) / length (x))
}

# Returns proportion of non-expressing cells
#
# Determines proportion of input values <= 0.
#
prop.nonexp <- function(x) {
  return(length(which(x<=0)) / length (x))
}

# Logistic function
#
logistic <- function(x, x0, k, c=1) {
  c / (1 + exp(-1*k*(x-x0)))
}

# Inverse logistic function
#
inv.logistic <- function(x, x0, k, c=1) {
  log(c/x-1)/(-1*k)+x0
}

# Determining preference between a pair of values
#
preference <- function(x, y, signed = FALSE) {
  z <- as.data.frame(cbind(x, y))
  if (signed) {
    z$p <- apply(z, 1, function(q) (q[1]-q[2])/(q[1]+q[2]))
  } else {
    z$p <- apply(z, 1, function(q) abs(q[1]-q[2])/(q[1]+q[2]))
  }
  # If both x and y are 0, preference will be NA, but want to return 0.
  z[is.na(z$p),"p"] <- 0
  return(z$p)
}

# Arithmetic Mean of As Numeric
#
num.mean <- function(x) mean(as.numeric(x))

sigmoid <- function(x, a, b, c, d) {
  c / (1 + exp(-a*(x-b))) + d
}

d.dx.sigmoid <- function(x, a, b, c, d=NULL) {
  a*c*exp(a*(x+b)) / ((exp(a*x) + exp(a*b))^2)
}

inv.sigmoid <- function(y, a, b, c, d) {
  -log(c/(y-d) - 1)/a + b
}

# Combine probabilities
#
# Calculates cumulative probability
#
combine.probs <- function(x) {
  1-prod(1-x)
}
