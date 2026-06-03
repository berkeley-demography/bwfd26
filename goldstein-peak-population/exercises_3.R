## Berkeley Formal Demography Workshop, 2026
## Peak Population, Exercise 3
## Analysis of the UN projection for a country (default Brazil)

## features:
## -- we'll analyze a real UN projection
## -- we'll fit curvature in log fertility decline
## -- we'll estimate the longevity effect
## -- we'll omit migration
## -- we'll see how accurate Coale + Goldstein/Cassidy is
## -- a baseline for doing your project on other countries

## ---- choose a country (change this for the group project) -------------
loc <- "Brazil"
## loc <- "China"
## loc <- "Germany"
## loc <- "Italy"
## loc <- "Iran (Islamic Republic of)"
## loc <- "United States of America"
## for others, try table(dt_full$Location) after you've loaded data

## read indicators and restrict to this country + selected variables
## (optional: look at names(dt_full) to list all variables)
library(data.table)
dt_full <- fread("~/Data/wpp/WPP2024_Demographic_Indicators_Medium.csv")
dt <- dt_full[Location == loc & !is.na(NRR),
              .(Location, Time, NRR, LEx, NetMigrations, Births, TPopulation1July)]
setorder(dt, Time)

## ---- 1. Choose t0 BY EYE, then fit log NRR centered at t0 -------------
## Automatic crossing-finders are unreliable when NRR is wiggly or hovers
## near 1 (Germany, Italy, US...), so we read t0 off the plot by hand.
## First, look at the NRR series:
par(mfrow = c(1, 1))
dt[, plot(Time, NRR, type = 'l',
          main = paste(loc, "NRR -- pick t0 where it crosses 1"))]
abline(h = 1, lty = 2); grid()

## eyeballed t0 (year NRR crosses replacement) for the sample countries --
## LOOK at the plot and adjust; for a new country, just set t0 by hand.
## t0_guess <- c("Brazil" = 2002, "China" = 1991, "Germany" = 1970,
##               "Italy" = 1977, "Iran (Islamic Republic of)" = 2010,
##               "United States of America" = 1972)
t0 <- 2000 ## >>>>>>> ADJUST <<<<<<
abline(v = t0, lty = 2)

## Fit log NRR as a quadratic in time CENTERED at t0, with NO intercept:
## t0 is (by definition) where NRR = 1, so log NRR = 0.
##  k1, k2 are then the linear and quadratic coefficients.
fit_dt <- dt[(Time - t0) %in% -20:30] ## ADJUSTABLE
fit_dt[, t := Time - t0]
m.nrr  <- lm(log(NRR) ~ -1 + t + I(t^2), data = fit_dt)
k1 <- coef(m.nrr)["t"]
k2 <- coef(m.nrr)["I(t^2)"]

## overlay the fitted log-quadratic curve
lines(fit_dt$Time, exp(fitted(m.nrr)), col = "red", lwd = 2)

## ---- 2. Longevity: exponential rate of e0 improvement near t0 ---------
## rho = slope of log(e0); fit over t0 .. t0+30
m.e0 <- lm(log(LEx) ~ I(Time - t0), data = dt[Time >= t0 & Time <= t0 + 30])
rho  <- coef(m.e0)[2]

## ---- 3. Standard stationary-population constants (as in the paper) ----
## Paper uses Sweden's 2015 stationary population: A0 = 42, sigma^2_l = 630,
## and mean age of childbearing mu0 = 30. (Country-specific values would come
## from that country's own life table + fertility schedule.)
mu0    <- 30
A0     <- 42
sig2_l <- 630

## ---- 4. Decompose the peak-population lag t_N - t0  (paper Eq. 13) ----
coale_lag  <- -mu0 / 2 + A0           # Coale baseline    (-mu/2 + A0)
curv_lag   <- (k2 / -k1) * sig2_l     # slowing fertility decline (k2 term)
longev_lag <- (rho * mu0) / -k1       # rising longevity  (rho term)
mig_lag    <- 0                       # migration omitted here
tN_hat_lag <- coale_lag + curv_lag + longev_lag + mig_lag

## ---- 5. Observed peak population (UN medium projection) ---------------
tN_obs     <- dt$Time[which.max(dt$TPopulation1July)]
tN_obs_lag <- tN_obs - t0

## ---- summary row ------------------------------------------------------
cat(sprintf("\n%s   (k1 = %.4f, k2 = %.5f, rho = %.4f, t0 = %.1f)\n",
            loc, k1, k2, rho, t0))
cat(sprintf("   -mu/2 + A0  (Coale) : %5.1f\n", coale_lag))
cat(sprintf("   + curvature (k2)    : %5.1f\n", curv_lag))
cat(sprintf("   + longevity (rho)   : %5.1f\n", longev_lag))
cat(sprintf("   + migration         : %5.1f   (omitted)\n", mig_lag))
cat(sprintf("   = t_N_hat lag       : %5.1f   (year %.0f)\n", tN_hat_lag, t0 + tN_hat_lag))
cat(sprintf("   t_N_obs lag (UN)    : %5.1f   (year %d)\n", tN_obs_lag, tN_obs))


## ---- 6. Four-panel summary (paper Fig. 2 idea, one country) -----------
## See how well the model does: NRR (+ fit), e0 (+ fit), births and pop with
## PREDICTED peaks (dashed orange) vs OBSERVED UN peaks (solid, panel color).
tB_hat <- t0 - mu0 / 2                        # predicted peak births
tN_hat <- t0 + tN_hat_lag                     # predicted peak population
tB_obs <- dt$Time[which.max(dt$Births)]       # observed peak births
## tN_obs (observed peak population) computed in step 5

nrrgreen <- "#1B9E77"; birthred <- "#C0392B"; popblue <- "#1F6FB2"; pred <- "orange"
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))

## NRR with log-quadratic fit
dt[, plot(Time, NRR, type = 'l', col = nrrgreen, lwd = 2,
          main = paste(loc, ": NRR"), xlab = "year", ylab = "NRR")]
lines(fit_dt$Time, exp(fitted(m.nrr)), col = "red", lwd = 2, lty = 2)
abline(h = 1, col = "grey60"); abline(v = t0, col = "grey60", lty = 3); grid()
abline(v = t0, lty = 2)
legend("topright", bty = "n", cex = 0.8, c("UN", "log-quad fit", "eye'd t0"),
       col = c(nrrgreen, "red", "black"), lwd = c(2,2,1), lty = c(1, 2,2))

## life expectancy e0 with fitted exponential (longevity, rate rho)
dt[, plot(Time, LEx, type = 'l', lwd = 2,
          main = paste(loc, ": life expectancy"), xlab = "year", ylab = "e0")]
lines(dt[Time >= t0 & Time <= t0 + 30]$Time, exp(fitted(m.e0)),
      col = "red", lwd = 2, lty = 2)
abline(v = t0, lty = 2); grid()

## births: observed peak vs predicted (t0 - mu/2)
dt[, plot(Time, Births / 1e3, type = 'o', col = birthred, lwd = 2,
          cex = .5,
          main = paste(loc, ": births"), xlab = "year", ylab = "births (millions)")]
abline(v = tB_obs, col = birthred, lwd = 1)            # observed peak
abline(v = tB_hat, col = pred, lwd = 2, lty = 2)       # predicted
abline(v = t0, col = "grey60", lty = 3); grid()
abline(v = t0, lty = 2)

## population: observed peak vs predicted (t0 + lag)
dt[, plot(Time, TPopulation1July / 1e3, type = 'l', col = popblue, lwd = 2,
          main = paste(loc, ": population"), xlab = "year", ylab = "population (millions)")]
abline(v = tN_obs, col = popblue, lwd = 2)             # observed peak
abline(v = tN_hat, col = pred, lwd = 2, lty = 2)       # predicted
abline(v = t0, lty = 2); grid()
legend("topright", bty = "n", cex = 0.8, c("observed peak", "predicted peak"),
       col = c(popblue, pred), lwd = 2, lty = c(1, 2))
