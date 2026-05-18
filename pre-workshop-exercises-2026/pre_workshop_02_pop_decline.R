## ============================================================================
##  Population Decline in a Brazil-like Country
## ----------------------------------------------------------------------------
##  Berkeley Formal Demography Workshop, 2026
##  Pre-workshop exercise
##  Author: Joshua R. Goldstein
## ============================================================================
##
##  Motivation
##  ----------
##  We simulate population decline in a Brazil-like country.  Mortality is
##  held fixed at 1970 levels; fertility declines continuously, reaching
##  replacement (NRR = 1) around the year 2000.
##
##  Two things to watch:
##    (1) The birth stream, or B(t)    -- when do births peak?
##    (2) Total population, or K(t)    -- when does the population peak?
##
##  What you'll learn
##  -----------------
##    - How to implement a time-varying projection by updating the
##      Leslie matrix at each step.
##    - Some stylized facts about populations on the way to decline.
##    - Motivation for learning how to do useful mathematical analysis of
##      peak pop
##
##  Reference Material:
##    "The Formal Demography of Peak Population" by Goldstein and Cassidy (2024).
##    - Google Drive link to reading: http://bit.ly/4uQnTfe
##
##  Two surprises to look out for
##  -----------------------------
##    - Births start to decline *before* fertility reaches replacement.
##    - Population starts to decline *after* fertility reaches replacement.
## ============================================================================

library(data.table)

fraction_female <- 0.4886    # Wachter standard fraction female at birth 
n    <- 5         # width of age interval (years)


## ============================================================================
##  1. Read 1970 demographic rates
## ============================================================================

brazil.dt <- fread("resources/brazil_like_country_data.csv")
x   <- brazil.dt$x      # age at start of interval
nLx <- brazil.dt$nLx    # life-table person-years (gets us portion surviving to each age interval)
nFx <- brazil.dt$nFx    # age-specific fertility rates

m <- length(x)          # number of age groups

par(mfrow = c(1, 2))
plot(x, nLx, type = 'o', main = "nLx, Brazil-like 1970",
     xlab = "Age", ylab = "nLx")
plot(x, nFx, type = 'o', main = "nFx, Brazil-like 1970",
     xlab = "Age", ylab = "nFx")

## ----------------------------------------------------------------------------
##  Q1.  What is the upper age interval?
##  Q2.  How wide are the age intervals?
##  Q3.  How many age groups are there?
## ----------------------------------------------------------------------------


## ============================================================================
##  2. Summary measures of the 1970 schedules
## ============================================================================

##  Q4.  What is e0 in 1970?
e0 <- sum(nLx)

##  Q5.  What is TFR in 1970?
tfr <- n * sum(nFx)

##  Q6.  What is NRR in 1970?  (And what does the 0.4886 stand for?)
nrr_1970 <- sum(nFx * nLx * fraction_female)


## ============================================================================
##  3. A fertility trajectory: exponential decline to replacement
## ============================================================================
##
##  We anchor the trajectory at two points,
##      NRR(1950) = 2     # rough pre-transition level
##      NRR(2000) = 1     # replacement reached around 2000
##  and let NRR decline exponentially in between (and continue past 2000).

nrr_1950 <- 2
nrr_2000 <- 1
k        <- log(nrr_2000 / nrr_1950) / 50    # exponential rate of decline

## Rescale 1970 fertility so that, at the rescaled levels, NRR = 1.
nFx_2000 <- nFx / nrr_1970

## sanity check: should be 1 if we recalculate the NRR at t=2000 manually
sum(nFx_2000 * nLx * fraction_female)

## NRR trajectory on a 5-year grid
t     <- seq(1950, 2100, by = n)
nrr_t <- nrr_2000 * exp(k * (t - 2000))

par(mfrow = c(1, 1))
plot(t, nrr_t, type = 'o',
     main = "Simulated NRR(t), reaching replacement in 2000",
     xlab = "Year", ylab = "NRR")
abline(v = 2000, lty = 2)
abline(h = 1.0,  lty = 2)


## ============================================================================
##  4. Replacement-level Leslie matrix A0
## ============================================================================
##
##  Standard discrete Leslie construction:
##    - Sub-diagonal:  nLx(x+n) / nLx(x)   (survival from one age interval to the next)
##    - First row:     period-averaged birth rate, weighted by female share

nFx   <- nFx_2000                      # use the rescaled (NRR=1) schedule
nFxpn <- c(nFx[-1], 0)
nLxpn <- c(nLx[-1], 0)

firstrow <- (nLx[1] / 2) * (nFx + nFxpn * nLxpn / nLx) * fraction_female
subdi    <- (nLxpn / nLx)[-m]

A0 <- matrix(0, m, m)
A0[1, ] <- firstrow
for (i in 1:(m - 1)) { 
    A0[i + 1, i] <- subdi[i] ## loop thru subdiagonal
}

## sanity check: leading eigenvalue should be 1
Re(eigen(A0)$values[1])


## ============================================================================
##  5. Starting population: stable at 1950 fertility (and 1970 life table)
## ============================================================================

brazil_pop_1950 <- 50    # million

A_1950      <- A0
A_1950[1, ] <- A0[1, ] * nrr_1950

K_1950_stable <- Re(eigen(A_1950)$vectors[, 1])                          # unnormalized stable pop
K_1950_stable_proportions <- (K_1950_stable / sum(K_1950_stable))        # normalize
K_1950_stable_counts <- K_1950_stable_proportions * brazil_pop_1950      # scale to correct size

## sanity check: should be 50 (million)
sum(K_1950_stable_counts)


## ============================================================================
##  6. Time-varying projection
## ============================================================================
##
##  At each step we rebuild the Leslie matrix by scaling the first row of
##  A0 by NRR(t) -- equivalent to scaling all age-specific fertility rates
##  uniformly.

n_steps    <- length(t) - 1
K.mat      <- matrix(NA, m, n_steps + 1)
K.mat[, 1] <- K_1950_stable_counts               # same initial population
lambda.vec <- numeric(n_steps)

for (i in 1:n_steps) {
    At             <- A0
    At[1, ]        <- A0[1, ] * nrr_t[i]         # update first row of leslie matrix with exponential decline from earlier
    lambda.vec[i]  <- Re(eigen(At)$values[1])
    K.mat[, i + 1] <- At %*% K.mat[, i]
}


## ============================================================================
##  7. Visualize: when do births peak?  When does population peak?
## ============================================================================

Kt <- colSums(K.mat)     # total pop
Bt <- K.mat[1, ]         # youngest age group, proxy for births

tB <- t[which.max(Bt)]   # year of peak births
tN <- t[which.max(Kt)]   # year of peak population

par(mfrow = c(3, 1))

plot(t, nrr_t, type = 'o',
     main = "NRR(t)",
     xlab = "Year", ylab = "NRR")
abline(v = 2000, lty = 2)
abline(h = 1, lty = 2)

plot(t, Bt, type = 'o', col = 'red',
     main = "Youngest age group (0-4), millions",
     xlab = "Year", ylab = "Bt")
abline(v = 2000, lty = 2)
abline(v = tB, col = 'red')
axis(1, at = tB, col = 'red')

plot(t, Kt, type = 'o', col = 'blue',
     main = "Total population K(t), millions",
     xlab = "Year", ylab = "Kt")
abline(v = 2000, lty = 2)
abline(v = tN, col = 'blue')
axis(1, at = tN, col = 'blue')


## ============================================================================
##  Discussion questions
## ============================================================================

## ----------------------------------------------------------------------------
##  Q7.  (Review) When does NRR = 1?
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
##  Q8.  When does population peak?  Before or after NRR = 1?
##       Sketch your intuition in a sentence or two.
##       (Wrong answers fine :)
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
##  Q9.  When do births peak?  Again, sketch your intuition.
## ----------------------------------------------------------------------------

## ----------------------------------------------------------------------------
##  Q10. What is the mean age of the stationary population, and what
##       might it have to do with the timing of the birth peak and pop peak?
## ----------------------------------------------------------------------------

## mean age of stationary population
sum((x + n / 2) * nLx) / sum(nLx)
## [1] 36.46233

## ----------------------------------------------------------------------------
##  Q11. What if we had a country with faster or slower rate of NRR decline? Would
##  births still peak about the same time before t0 and pop peak about the same time
##  after? Or ...?
## ----------------------------------------------------------------------------


## ----------------------------------------------------------------------------
##  Q12. What might be unrealistic about this simulation?
## ----------------------------------------------------------------------------

## Hint: Compare your results with Figure 2 of Goldstein and Cassidy paper. 
