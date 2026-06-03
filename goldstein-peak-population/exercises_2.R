## Berkeley Formal Demography Workshop, 2026
## Peak Population, Exercise 2
## Timing of peaks using time-varying Leslie matrices
## declining fertility

## features: 
## -- Leslie matrices, rather than renewal equation
## -- Realistic demography from Brazil 2000
## -- Get to explore different speeds of fertility decline k
## -- Get to play with shiny app
## -- Get to how "hump" of births moves thru age-structure to give
##    peak pop.

## We PROJECT a female population forward with a time-varying Leslie
## matrix (Brazil 2000 schedules, Coale decline NRR(t) = exp(k t),
## k<0) and read births and population size straight off the state
## vector.
##
## We compare the numerical peak times to the analytic predictions:
##   births peak at      t_B = -mu / 2
##   population peaks at  t_N = t_B + A_0
##
## STARTING STATE (assumption): we begin the projection at t = -150 in
## the STATIONARY age structure (K0 proportional to L_x).  At t = -150
## fertility is very high, so this start is deliberately
## "wrong", but the 150-year burn-in lets the transient decay.

## ---- schedules: Brazil 2000 (UN WPP2024), single year of age ----------
csv <- "brazil_2000_1_year_rates.csv" ## pre-extracted for you
if (!file.exists(csv)) csv <- file.path("..", csv)   # script lives in slides/
br  <- read.csv(csv)

x  <- br$x                             # ages 0..100
m  <- length(x)                        # number of age classes
n  <- 1                                # age interval
Lx <- br$Lx                            # Person-year survivorship (l0 = 1)
Fx_2000 <- br$Fx                       # ASFR
e0 <- sum(Lx)                          # life expectancy at birth
A0 <- sum((x + n/2) * Lx) / sum(Lx) # mean age of the stationary population
mu <- sum((x + n/2) * Fx_2000) / sum(Fx_2000)  # mean age of childbearing

## ---- build the Leslie matrix (NRR normalized to 1) --------------------
nrr_2000 <- sum(Lx * Fx_2000 * 0.4886)  ## .4886 is fraction female births
Fx <- Fx_2000 / nrr_2000                ## normalized so nrr = 1 at time 0
print(sum(Fx * Lx * 0.4886))            ## check: should be 1

## now get L(x+n) and F(x+n) for building the Leslie matrix
Lxpn <- c(Lx[-1], 0)
Fxpn <- c(Fx[-1], 0)

## subdiagonal survival: s_a = L_{a+1} / L_a
subdi <- (Lxpn / Lx)[-m]
## first-row fertility (birth-flow): avg of this & next age's net fertility
firstrow <- Lx[1] / 2 * (Fx + Fxpn * Lxpn / Lx) * 0.4886

## build the stationary Leslie matrix
L0 <- matrix(0, m, m)                   # the t = 0 Leslie matrix
L0[1, ] <- firstrow
L0[cbind(2:m, 1:(m - 1))] <- subdi

## check the leading eigenvalue is ~1 (since NRR was normalized to 1)
cat(sprintf("dominant eigenvalue of L0 = %s\n", eigen(L0)$values[1]))

## ---- project with declining fertility ---------------------------------
## K(t+1) = L_t %*% K(t);
## L_t scales the fertility row by NRR(t) = exp(k t), k<0.
k     <- -0.01    # k < 0. >>>>> CHANGE "k" HERE <<<<<
t.vec <- -150:100
K.mat <- matrix(NA, m, length(t.vec))
K.mat[, 1] <- Lx                        # start in stationary age structure
for (i in 1:(length(t.vec) - 1)) {
    t  <- t.vec[i]
    Lt <- L0
    Lt[1, ] <- L0[1, ] * exp(k * t)     # this year's fertility row
    K.mat[, i + 1] <- Lt %*% K.mat[, i]
}
Bt <- K.mat[1, ]                        # age-0 entries = births
Kt <- colSums(K.mat)                    # population size

## get peak births and peak pop (denoting "N" for "K")
t_B <- t.vec[which.max(Bt)]
t_N <- t.vec[which.max(Kt)]

## ---- compare numerical peaks to analytic predictions ------------------
cat(sprintf("mu = %.1f,  A_0 = %.1f,  e_0 = %.1f\n", mu, A0, e0))
cat(sprintf("peak births: Leslie t_B = %d  analytic t_B = -mu/2 = %.1f\n",
            t_B, -mu/2))
cat(sprintf("peak pop: Leslie t_N = %d  analytic t_N = t_B + A_0 = %.1f\n",
            t_N, -mu/2 + A0))


## ---- plot -------------------------------------------------------------
phi <- Lx * Fx * 0.4886                 # net maternity (sums to NRR = 1)
par(mfrow = c(3, 2))

plot(x, Lx, type = 'l', lwd = 2, main = "survivorship l(x): Brazil 2000",
     xlab = "age", ylab = "l(x)", ylim = c(0, 1))
abline(v = A0, col = "grey60", lty = 2)
legend("topleft", bty = "n", sprintf("e0 = %.0f, A0 = %.0f", e0, A0))

plot(x, phi, type = 'l', lwd = 2, main = "net maternity: Brazil 2000",
     xlab = "age", ylab = "l(x) f(x)")

plot(t.vec, exp(k * t.vec), type = 'l', lwd = 2, main = "NRR(t) = exp(k t), k<0",
     xlab = "t", ylab = "NRR"); abline(h = 1, v = 0, col = "grey60")

## Births
plot(t.vec, Bt, type = 'l', lwd = 2, col = "red",
     main = "births B(t)  [Leslie]", xlab = "t", ylab = "B(t)")
abline(v = t_B,   col = "red", lty = 2)
## abline(v = -mu/2, col = "black",   lty = 3)
grid()
legend("topleft", bty = "n", sprintf("obs. t_B = %.0f", t_B))
       
## blank plot to align B(t) and K(t)
plot(0,0, type = "n", axes = FALSE, ylab = "", xlab = "") 
title("Blank space, ...")

## Pop size
plot(t.vec, Kt, type = 'l', lwd = 2, col = "blue",
     main = "population K(t)  [Leslie]", xlab = "t", ylab = "K(t)")
abline(v = t_N,        col = "blue", lty = 2)
## abline(v = -mu/2 + A0, col = "black",   lty = 3)
grid()
legend("topright", bty = "n", sprintf("obs. t_N = %.0f", t_N))

## Questions 1

## 1. Do things look right? How did you check?

## 2. Fill in the following table

## k         observed t_B       observed t_N
## -.005
## -.01
## -.02
## -.05

## Do the peaks vary with k, or are they pretty much invariant?

## 3. Play with the shiny app.

## 4. With k = -.01, visualize peak births moving thru age structure
## by uncommenting  code below

## dimnames(K.mat) <- list(x, t.vec)
## my_t.vec = seq(-20, 50, 10)
## par(mfcol = c(4,2))
## for (i in 1:length(my_t.vec))
## {
##      t = my_t.vec[i]
##      plot(x, K.mat[,paste(t)], main = t, ylim = c(0, 100),
##           xlab = "K(x)",
##           cex = .5)
##      lines(x, K.mat[, "25"], col = 'red') ## K(a, t = 25) -- about peak pop
##      legend("bottomleft", c(paste0("Pop size = ",
##                                  round(sum(K.mat[,paste(t)]))),
##                           paste0("Peak size = ",
##                                  round(max(Kt)))),
##             text.col = c("black", "red"),
##             bty = 'n')
##      grid()
## }

## Can you see that before peak, more young don't make up for fewer old; and after peak, vice-versa?

## Congratulations! You are done with Exercises 2!




