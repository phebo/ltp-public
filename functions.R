# Copyright (C) 2024  Phebo D. Wibbens
#   
#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

ddln <- function(x, mu1, mu2, sig1, sig2, p) {
  # Density for the double-log-normal distribution
  (1-p)*dlnorm(-x, mu1, sig1) + p*dlnorm(x, mu2, sig2)
}

dst <- function(x, m, s, df) {
  # Density for the scaled-t distribution
  dt((x-m)/s, df)/s
}

daep <- function(x, m, al, bl, ar, br, log = FALSE) {
  A0 <- function(x) x^(1/x-1) * gamma(1/x)
  C <- al * A0(bl) + ar * A0(br)
  out <- -log(C) - ifelse(x < m, 1/bl * abs((x-m)/al)^bl, 1/br * abs((x-m)/ar)^br)
  if(log) return(out)
  else return(exp(out))
}

dstable2 <- function(x, alpha, beta, gamma, delta, log = FALSE) {
  # Density of stable distribution that gives back NA instead of an error
  tryCatch(dstable(x, alpha, beta, gamma, delta, log = log), error = function(e) NA)
}

dgnorm2 <- function(x, mu, lalpha, beta, log = FALSE) {
  # Use logarithm of scale parameter alpha for numerical stability in mle calculation
  dgnorm(x, mu, exp(lalpha), beta, log = log)
}

mle.nln <- function(u, n = NA, ...) {
  # Calculate NLN parameters based on vector of observations u
  # When specifying n, do MLE regression on that number of bins to speed up, otherwise ordinary MLE
  if(is.na(n)) {
    nll <- function(mux=0, muy=0, sig=1, rho=0) -sum(log(dnln(u, mux, muy, sig, rho)))
  } else {
    df <- tibble(lu = log(abs(u)), s = ifelse(u > 0, 1, -1)) %>%
      group_by(cut(lu, n), s) %>% summarize(
        k = n(), 
        lu = mean(lu)) %>% ungroup() %>%
      mutate(
        u = s*exp(lu))
    nll <- function(mux=0, muy=0, sig=1, rho=0) -sum(df$k * log(dnln(df$u, mux, muy, sig, rho)))
  }
  fit <- mle(nll, ...)
  muy = coef(fit)["muy"]
  sig = coef(fit)["sig"]
  rho = coef(fit)["rho"]
  out <- tibble(
    par = c(names(coef(fit)), "r", "rho.sig", "p"),
    est = c(coef(fit), muy+sig*rho, rho*sig, pnorm(muy)),
    se = c(sqrt(diag(vcov(fit))),
           sqrt(t(c(0,1,rho,sig)) %*% vcov(fit) %*% c(0,1,rho,sig)),
           sqrt(t(c(rho,sig)) %*% vcov(fit)[c("sig","rho"), c("sig","rho")] %*% c(rho,sig)),
           dnorm(muy) * sqrt(vcov(fit)["muy","muy"])))
  names(out$est) <- names(out$se) <- out$par
  out
}

calc.hist <- function(u, n=sqrt(length(u)/20), psLog = NA, lim = c(min(u), max(u))) {
  # Make a histogram with equal bins on pseudoLog / asinh or linear (if psLog = NA) scale
  nobs <- length(u)
  u <- u[between(u,lim[1],lim[2])]
  frac <- length(u) / nobs # Fraction of included observations
  if(is.na(psLog)) {
    hst <- hist(u, plot=F, breaks=seq(lim[1], lim[2], length.out=n))
    out <- tibble(u = hst$mids)
  } else {
    u2 <- asinh(u / 2 /psLog)
    lim2 <- asinh(lim / 2 /psLog)
    hst <- hist(u2, plot=F, breaks=seq(lim2[1], lim2[2], length.out=n))
    out <- tibble(u2 = hst$mids,
                  u = 2 * psLog * sinh(u2))
  }
  out %>%
    mutate(
      k = hst$counts,
      d = hst$dens * frac,
      dmin = qgamma(0.025, k) * d / k,
      dmax = qgamma(0.975, k) * d / k) %>%
    filter(k != 0)
}

calc.dens <- function(df, lim, n, psLog) {
  # Caculate the density function on the pseudoLog / asinh scale, based on parameter estiamtes df
  lim2 <- asinh(lim / 2/ psLog)
  tibble(u2 = seq(lim2[1], lim2[2], length.out=n),
         u = 2 * psLog * sinh(u2),
         d = do.call(dnln, c(list(x=u), as.list(df$est[1:4]))),
         d2 = d * sqrt(u^2 + (2*psLog)^2))
}

