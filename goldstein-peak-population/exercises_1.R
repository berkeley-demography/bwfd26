## Berkeley Formal Demography Workshop, 2026
## Peak Population, Exercise 1
## Verifying analytical solution for renewal equation with
## declining fertility


## Full (time-varying) renewal:
##     B(t) = sum_x B(t-x) * phi(x,0) * exp(k*t)        [k < 0]
##
## Concentrated-age ("difference equation") solution from the slides:
##     log B(t) = log B(0) + (k/2) t + k/(2 mu) t^2     [k < 0]
##
## We build a simple net-maternity schedule phi, run the full renewal
## numerically, and compare it to the log-quadratic analytical solution.

## note: net maternity phi(x) = f(x) * l(x)

## ---- net maternity schedule -------------------------------------------
## A smooth bell centered at age 30; normalize so NRR = sum(phi) = 1.
x       <- 0:50 ## age
phi.x   <- dnorm(x, mean = 30, sd = 4) ## <-- ADJUST VARIANCE HERE
phi.0   <- phi.x / sum(phi.x)        # NRR(0) = 1
names(phi.0) <- x                    # name by age: phi.0["30"] is age 30
mu      <- round(sum(phi.0 * x))     # mean age of childbearing = 30

## ---- time grid and decline rate ---------------------------------------
k     <- -0.01    # k < 0: rate of exponential fertility change 
t.vec <- -200:200 ## time, relative to t0, when NRR = 1
B0    <- 1 ## arbitrary
log.B0 = log(B0)

## ages used in the recurrence: exclude x = 0 so B(t) never references
## itself on the right-hand side when we calculate (t-x). (A hack!)
x.rec <- 1:50


## ---- analytical (log-quadratic) solution ------------------------------
Bt_analytic <- exp(log.B0 + (k/2) * t.vec + k / (2 * mu) * t.vec^2)

## ---- numerical renewal ------------------------------------------------
## seed the first years of history (t = -200 .. -151) with the
## analytical values, then march recurrence forward.
Bt_simu        <- rep(NA, length(t.vec))
names(Bt_simu) <- t.vec

t.pre        <- -(200:151)                       # 50 years of history
Bt_simu[paste(t.pre)] <- Bt_analytic[t.vec %in% t.pre]
## alternative: Bt_simu[paste(t.pre)] <- rep(1, length(t.pre)
Bt_simu[paste(t.pre)] <- rep(.1, length(t.pre))

start.i <- which(t.vec == -150)                  # first t we project for
for (i in start.i:length(t.vec)) {
    t      <- t.vec[i] ## this value of t

    ## time-varying fertility  phi(x,t)
    phi.xt <- phi.0[paste(x.rec)] * exp(k * t)   

    ## RENEWAL EQUATION
    Bt_simu[paste(t)] <- sum(Bt_simu[paste(t - x.rec)] * phi.xt) 
}

## ---- plot the results
par(mfrow = c(2,2))
plot(x, phi.0, type = 'l', main = "Net maternity phi(x)")
plot(t.vec, Bt_simu, type = 'l', lty = 2, lwd = 3,
     xlab = "Time t", ylab = "Births B(t)",
     main = "Projected vs. analytic births") ## renewal-based projection
grid()
lines(t.pre, Bt_simu[paste(t.pre)], lwd = 5, col = 'blue') ## init conditions
lines(t.vec, Bt_analytic, col = 'red', lwd = 2) ## analytic


legend("topright",
       c("initial conditions",
         "projection",
         "analytic"),
       col = c("blue", "black", "red"),
       lwd = c(5, 3, 2),
       lty = c(1, 2, 1))

## Qs

## 1. Does your projection come close to the analytic solution?
## What are differences/similarities?

## 2. See what happens if you make the phi have very little variance,
## e.g. 0.4 instead of 4. And what happens if you have phi have big
## variance: e.g., 14 instead of 4.

## 3. Add some code to if timing of peak births is ok, even when SD is
## high hint: try putting "abline(v = t.vec[which.max(B)], lty = 2)"
## after your plot to add a line at maximum of projected births.

## 4. Plot log(births). Does it look quadratic?

## 5. If you have time, try different initial conditions, e.g.,
## a constant set of births to start:
##         Bt_simu[paste(t.pre)] <- rep(1, length(t.pre))
## And if you want try varying SD. Describe your results. Do you still
## get a bell-shaped curve? Does it still have roughly the same peak?
## What if anything changes?

## Congratulations! You are done with Exercises 1!
