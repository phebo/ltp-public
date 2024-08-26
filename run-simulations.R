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

library(Rcpp)
library(mvtnorm)
library(gnorm)
library(stats4)
library(tidyverse)
sourceCpp("dnln.cpp")
source("functions.R")

writeOutput <- T

#### Simulation check of NLN density and MLE ####

psLog <- 0.02 
mux <- -1; muy <- -0.3; sig <- 2; rho <- 0.5
n <- 1e4
set.seed(497031) # Fixed seed for code replicability
xy <- rmvnorm(n, mean = c(mux, muy),
              sig = rbind(c(sig^2, rho*sig), c(rho*sig, 1)))
z <- xy[,2] * exp(xy[,1])
pars <- mle.nln(z)
pars %>% mutate(lower = est-1.96*se, upper = est+1.96*se) #95% intervals

calc.hist(z, n=50, psLog=psLog) %>%
  ggplot() + 
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), size=0.02, linewidth = 0.3) +
  geom_line(data = calc.dens(pars, c(-100,2000), 500, psLog),
            mapping = aes(x=u, y=d2), linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-4,1e-2,1e-0),
                labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  xlab("Simulated LTP (pseudo-log)") + ylab("Density (log)") + theme_bw() +
  guides(linetype = guide_legend(byrow = TRUE))
if(writeOutput) ggsave("fig-sim-check.pdf", width = 3, height = 2.5)

mean(z)
(muy + rho*sig) * exp(mux + sig^2/2)

#### Dynamic process ####

psLog <- 1

sim.process <- function(n=10, Tt=20, alpha = 0.1, beta = 1, rho = 0, X0 = 1) {
  aPl = sqrt(.5 + .5*sqrt(1-rho^2))
  aMin = sqrt(.5 - .5*sqrt(1-rho^2))
  expand_grid(i=1:n, t=0:Tt) %>%
    mutate(
      z1 = rgnorm(n * (Tt+1), alpha=alpha, beta=beta),
      z2 = rgnorm(n * (Tt+1), alpha=alpha, beta=beta),
      dx = aPl * z1 + aMin * z2,
      dy = aMin * z1 + aPl * z2) %>%
    group_by(i) %>% mutate(X = X0 * exp(cumsum(c(0,dx[-length(dx)])))) %>% ungroup() %>%
    mutate(dY = X * dy)
}

set.seed(4938749) # Set random seed to allow for exact replication
dfSim <- expand_grid(beta = c(0.7, 0.8, 0.9), rho = c(0,.4,.8)) %>%
  rowwise() %>% mutate(
    panel = list(sim.process(n=1e4, Tt=40, beta=beta, rho=rho)),
    Y = list(panel %>% group_by(i) %>% summarize(Y=sum(dY)) %>% pull(Y)),
    #par = list(mle.nln(Y, n=100)),
    par = list(mle.nln(Y)),
    hist = list(calc.hist(Y, n=40, psLog = psLog)),
    fit = list(calc.dens(par, c(min(Y),max(Y)), 500, psLog))) %>%
  ungroup()
dfSimHist <- dfSim %>% select(beta, rho, hist) %>% unnest(hist) %>%
  mutate(beta = paste("beta ==", beta), rho = paste("rho[0] ==", rho))
dfSimFit <- dfSim %>% select(beta, rho, fit) %>% unnest(fit) %>%
  mutate(beta = paste("beta ==", beta), rho = paste("rho[0] ==", rho))

ggplot() + facet_grid(rho ~ beta, scales="free_x", space="free_x", labeller = label_parsed) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), dfSimHist, size=0.02, linewidth = 0.3) +
  geom_line(aes(x=u, y=d2), dfSimFit, linewidth = 0.3) +
  scale_y_log10(limits=c(5e-6,2), breaks = c(1e-4,1e-2,1), labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(3,1)),0,10^c(1,3)),
                     minor_breaks = c(-100,-1,1,100)) +
  xlab("Simulated LTP (pseudo-log)") + ylab("Density (log)") + theme_bw() +
  guides(linetype = guide_legend(byrow = TRUE))
if(writeOutput) ggsave("fig-sim-dynamic.pdf", width = 6, height = 5)

