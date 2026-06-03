# read in a bunch of Lx's, Fx's and Kx's from a few countries
# They're 5-year age groups so let's say our projection step is 5 years
# I'm assuming the radix is 1. If not, I have to adjust everything by the radix. 

dat =read.csv("projection-rates.csv")
head(dat)

# Define some functions to use 
TFR <- function(f, n=5) return(n * sum(f))   # default n = 5-year rates
GRR <- function(f, n=5, ffab=0.4886) return(n * sum(f) * ffab)
NRR <- function(l,f,ffab=0.4886) return( ffab*sum(l*f))
mu <- function(l,f, n=5) {
   x=seq_along(l)*n - n/2
   return( sum(x*l*f)/sum(l*f)) }

leslie <- function(l,f, ffab=0.4886) {
  n = length(l)-1
  A = matrix(0, n, n)
  A[1,] = ffab * l[1]*(f[1:n]+f[2:(n+1)]*l[2:(n+1)]/l[1:n])/2  # top row
  diag(A[2:n, 1:(n-1)]) = l[2:n] / l[1:(n-1)]  # subdiagonal
return(A)
}

attach(dat)

# Examples of functions   =========
TFR(Brazil.Fx) 
TFR(Niger.Fx)
TFR(USA1934)
plot(x,Niger.Fx,type="l",ylab="5Fx")
lines(x,USA1934.Fx,col="red")
lines(x,Brazil.Fx,col="blue")

NRR(USA1934.Lx,USA1934.Fx, 0.4877) # in this case, we know the Ffab = 0.4877
mu(USA1934.Lx,USA1934.Fx)

NRR(Canada.Lx,Brazil.Fx)
mu(Kenya.Lx,Niger.Fx)

# Let's create some Leslie matrices
CB = leslie(Canada.Lx,Brazil.Fx)  # We can create Leslie matrices using a variety of rates
CB
UU= leslie(USA1934.Lx,USA1934.Fx, 0.4877) #These are the US rates
UU

# ================ Projection code here, using UU and initial Russian Kx
A = UU                       # For this projection, use the US 1934 rates
K = matrix(0,nr=10,nc=51)    # Create a big matrix to hold our projection results
K[1:10,1] = Russia.Kx[1:10]  # Here's where we initialize the population in year 0.

for(i in 1:50) K[,i+1] = A %*% K[,i]  # Here's where we do the entire projection
K                            # That was quick

totalpop <- colSums(K)
year=seq(0,250,by=5)
plot(year,totalpop)            # plot the total pop
plot(year, totalpop, log="y")  # plot the total on log scale

# Looks pretty linear on log scale, i.e., exponential growth. Let's calculate 
# growth rate for the last five years of the projection period

log(totalpop[51]/totalpop[50])/5 


# Hmmm. We've seen that number before

# We do the projection by age group, so let's look at the age distribution
# We could either look at the raw totals, or the percentage age distribution
# K already tells us the raw totals in each 5-year age group. To get the 
# age distribution in percentage terms, we'll divide each column by the totalpop.

K.distrib <- sweep(K, 2, totalpop, "/")

K.distrib

# Woah. It looks like the proportions in each age group converge. 
# Let's check that with a plot

plot(year, K.distrib[1,], type="l", col=1)
lines(year, K.distrib[2,], col=2)
lines(year, K.distrib[3,], col=3)

lines(year, K.distrib[10,], col=1)

# Hmmm. The proportions in each age group have stabilized and are in equilibrium

# We could repeat this with a different initial population. (Do this for diff K(0) )



# ================ Ignore the part below this for now
eigen.UU = eigen(UU)
r.UU = log(Re(eigen.UU$value[1]))/5
# Hmmm. Where have we seen that before?

v = Re(eigen.UU$vectors[,1])
k.UU = v/sum(v)
# Hmmm^2. Where have we seen that before?