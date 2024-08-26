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

library(gnorm)
library(stabledist)
library(stats4)
library(MASS, include.only = "fitdistr")
library(xtable)
library(Rcpp)
library(tidyverse)
library(mvtnorm)
sourceCpp("dnln.cpp")
source("functions.R")
options("scipen"=100)

dataDir <- "."
writeOutput <- T 
calcFit <- F 
yrs <- seq(2002, 2022, 5) # Analysis of subperiods
psLog <- 0.02 # Parameter used for pseudo-log transformation x -> asinh(x / 2 / psLog)

#### Read data ####

dfGics2 <- read_csv("db-gics2.csv")
dfP <- read_csv(file.path(dataDir, "ltep-panel-02-21.csv")) %>%
  mutate(giccd2 = giccd %/% 100) %>% left_join(dfGics2) %>% mutate(ind2 = paste(giccd2, ind2)) %>%
  mutate(geo = case_when(
    curcd == "USD" ~ "USA",
    curcd == "EUR" ~ "Euro",
    curcd == "GBP" ~ "UK"))
dfLt <- dfP %>% group_by(gvkey) %>%
  summarize(ltep = sum(epd), giccd2 = last(giccd) %/% 100, conm=last(conm), geo=last(geo),
            nyr = n(), start=min(fyear), end=max(fyear)) %>% ungroup() %>%
  mutate(s = factor(ltep > 0, levels=c(T, F), labels=c("Positive", "Negative")))
dfWealth <- read_csv(file.path(dataDir, "wealthcreation16_0.csv"))
  # Downloaded from https://wpcarey.asu.edu/department-finance/faculty-research/do-stocks-outperform-treasury-bills in November 2023
dfLiva <- read_csv(file.path(dataDir, "db-liva.csv")) 
  # Downloaded from https://github.com/phebo/liva-public in November 2023

ltepUS <- dfLt %>% filter(geo == "USA") %>% pull(ltep)
altp <- dfLt %>% filter(geo == "USA") %>% mutate(altp = ltep / nyr) %>% pull(altp) # Annual LTP
liva <- dfLiva %>% filter(between(year, 2002, 2021)) %>% group_by(gvkey) %>%
  summarize(liva = sum(liva)) %>% arrange(-liva) %>% pull(liva) 
wealth <- dfWealth  %>% pull(wealth)

dfLt %>% group_by(geo) %>% summarize(n=n())

# VISA example caluclation missing data
dfP %>% filter(gvkey == 179534) %>% select(fyear, epn, epd) %>%
  xtable(digits = c(0,0,1,1)) %>% print(include.rownames = F)
dfP %>% filter(gvkey == 179534) %>% 
  group_by(period = cut(fyear, breaks = yrs, labels = paste(yrs[-5], yrs[-5]+4, sep = "-"), right = F)) %>%
  summarize(ltep = sum(epd))


#### US plots ####

# Cumulative chart a la Bessembinder (2018, Fig 2)
dfLt %>% filter(geo == "USA") %>% arrange(-ltep) %>%
  mutate(n=row_number(), cum = cumsum(ltep)/1e3) %>%
  bind_rows(tibble(n=0, cum=0, s="Positive"), .) %>%
  ggplot(aes(x=n, y=cum)) + geom_line() + geom_vline(xintercept = 4398, linetype=2) +
  xlab("Number of firms") + ylab("Cumulative LTP ($T)") +
  scale_color_brewer(palette = "Set1", direction = -1) + theme_bw() +
  scale_x_continuous(labels = function(x) scales::comma(x)) + ylim(c(0,27.5)) +
  annotate("text", x=0, y=27, label = expression("Positive LTP " %<-% " "), hjust = "left", vjust = "top") +
  annotate("text", x=4300, y=27, label = expression(" " %->% " Negative LTP"), hjust = "left", vjust = "top") +
  annotate("label", x = 150, y = 3.3, label = "Top 10: $3.3T", hjust = "left") +
  annotate("point", x = 10, y = 3.3) +
  annotate("label", x = 250, y = 13.3, label = "Top 100: $13.3T", hjust = "left") +
  annotate("point", x = 100, y = 13.3) +
  annotate("label", x = 12000, y = 17.8, label = "All 11,656: $18.2T", hjust = "right", vjust = "top") +
  annotate("point", x = 11656, y = 18.2) +
  theme(plot.margin = unit(c(5,15,5,5), "pt"))
if(writeOutput) ggsave("fig-us-cum.pdf", width = 4, height = 3)
(dfLt %>% filter(geo == "USA") %>% arrange(-ltep) %>%
    mutate(cum = cumsum(ltep)/1e3))[c(1,10,100,4398,4399,length(ltepUS)),]
length(ltepUS)
1-4398/11656
13.3/18.2

dfLt %>% filter(geo == "USA") %>%
  ggplot(aes(x=ltep)) + geom_histogram(bins=50) +
  geom_vline(xintercept = 0, linetype=2) +
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(3:-1),0,10^(-1:3))) +
  xlab("LTP ($B, pseudo-log)") + ylab("#Firms") + theme_bw() +
  theme(plot.title.position = "plot")
if(writeOutput) ggsave("fig-us-hist.pdf", width = 6, height = 3)

#### US fits ####

# Calculate 5 fits for 4 data sets
dfDataUS <-with(dfP %>% filter(geo == "USA", fyear==2021),
                tibble(name=c("roa", "eroc", "epn", "ltp"),
                       u=list(roa, roc-cc, epn, ltepUS)))
dfUS <- expand_grid(
  dfDataUS,
  tibble(dist = c("st", "dln", "nln", "aep", "stable"),
         expr = c(expression(fitdistr(u, "t")),
                  expression(fitdistr(u[u!=0], ddln, start = list(mu1=0, mu2=0, sig1=1, sig2=1, p=0.5))),
                  expression(fitdistr(u, dnln, start = list(mux=0, muy=0, sig=1, rho=0))),
                  expression(fitdistr(u, daep, list(bl=1,br=1,al=1,ar=1,m=-0.1))),
                  expression(fitdistr(u, dstable2, list(alpha=0.5, beta=-0.1, gamma=0.1, delta=-0.02))))))
dfUS[dfUS$dist == "aep" & dfUS$name == "roa", "expr"][[1]][[1]] <- 
  expression(fitdistr(u, daep, list(bl=1,br=1,al=1,ar=1,m=0.2)))
dfUS[dfUS$dist == "aep" & dfUS$name == "eroc", "expr"][[1]][[1]] <-
  expression(fitdistr(u, daep, list(bl=1,br=1,al=1,ar=1,m=0.3)))

if(calcFit) {
  dfUS <- dfUS %>% rowwise() %>% mutate(fit = list(eval(expr))) %>% ungroup()
  dfParUS <- dfUS %>% rowwise() %>% mutate(out = list({
    vars <- names(coef(fit))
    bind_rows(
      tibble(par = vars, type = "est", value = coef(fit)),
      tibble(par = vars, type = "se", value = sqrt(diag(vcov(fit)))),
      tibble(par = c("n", "loglik", "AIC"), type = "stat",
             value = c(fit$n, logLik(fit), AIC(fit))))})) %>%
    ungroup() %>% select(name, dist, out) %>% unnest(out)
  if(writeOutput) write_csv(dfParUS, "fits-us.csv")
} else {
  dfParUS <- read_csv("fits-us.csv")
}
dfDist <- tibble(dist = c("st", "stable", "aep", "dln", "nln"),
                 name = c("Scaled t", "Stable", "Asymmetric ex-\nponential power",
                          "Double\nlog-normal", "Normal-\nlog-normal"),
                 name2 = c("Scaled t", "Stable", "Asymmetric exponential power",
                          "Double log-normal", "Normal-log-normal"))

# LTP & economic profit charts
psLog <- 0.02 
lim <- c(-1e3, 1e3)
dfFit <- dfParUS %>% filter(name %in% c("ltp", "epn"), type=="est") %>%
  select(-type) %>% group_by(name, dist) %>% nest() %>%
  mutate(out = list({
    pars <- as.list(data[[1]]$value)
    names(pars) <- data[[1]]$par
    lim2 <- asinh(lim / 2/ psLog)
    tibble(u2 = seq(lim2[1], lim2[2], length.out=501),
           u = 2 * psLog * sinh(u2),
           d = do.call(paste0("d", dist), c(list(x=u), pars)),
           d2 = d * sqrt((u)^2 + (2*psLog)^2))})) %>%
  ungroup() %>% select(-data) %>% unnest(out) %>%
  filter(d2 > 1e-5) %>% mutate(
    dist = factor(dist, levels = dfDist$dist, labels = dfDist$name),
    name = factor(name, levels = c("ltp", "epn"), labels = c("20-year LTP", "1-year EP")))
dfHist <- dfDataUS %>% filter(name %in% c("ltp", "epn")) %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=50, psLog=psLog))) %>% ungroup() %>% select(-u) %>% unnest(hist) %>%
  mutate(name = factor(name, levels = c("ltp", "epn"), labels = c("20-year LTP", "1-year EP")))
ggplot() + facet_wrap(~name, ncol = 1) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), dfHist, size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2, color = dist, linetype = dist), dfFit, linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-4,1e-2,1e0), labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-4,1)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("Profit ($B, pseudo-log)") + ylab("Density (log)") + theme_bw() +
  guides(linetype = guide_legend(byrow = TRUE)) +
  theme(legend.title = element_blank(), legend.spacing.y = unit(10,"pt"))
if(writeOutput) ggsave("fig-fits-us.pdf", width = 6.5, height = 7)

# LTP only
ggplot() + 
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), filter(dfHist, name == "20-year LTP"), size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2, color = dist, linetype = dist), filter(dfFit, name == "20-year LTP"), linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-4,1e-2,1e0), labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-4,1)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("LTP ($B, pseudo-log)") + ylab("Density (log)") + theme_bw() +
  guides(linetype = guide_legend(byrow = TRUE)) +
  theme(legend.title = element_blank(), legend.spacing.y = unit(10,"pt"))

# LTP & EP charts: zoom around zero
lim <- c(-0.01,0.01)
dfFit <- dfParUS %>% filter(name %in% c("ltp", "epn"), type=="est") %>% select(-type) %>%
  group_by(name, dist) %>% nest() %>%
  mutate(out = list({
    pars <- as.list(data[[1]]$value)
    names(pars) <- data[[1]]$par
    lim2 <- asinh(lim / 2/ psLog)
    tibble(x = seq(lim[1], lim[2], length.out=2e3),
           dens = do.call(paste0("d", dist), c(list(x=x), pars)))})) %>%
  ungroup() %>% select(-data) %>% unnest(out) %>% mutate(
    dist = factor(dist, levels = dfDist$dist, labels = dfDist$name),
    name = factor(name, levels = c("ltp", "epn"), labels = c("20-year LTP", "1-year EP")))
dfHist <- dfDataUS %>% filter(name %in% c("ltp", "epn")) %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=30, lim=lim))) %>% ungroup() %>% select(-u) %>% unnest(hist) %>%
  mutate(name = factor(name, levels = c("ltp", "epn"), labels = c("20-year LTP", "1-year EP")))
ggplot() + facet_wrap(~ name, ncol = 1) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), dfHist, size=0.1, linewidth = 0.3) +
  geom_line(data = dfFit, mapping = aes(x=x, y=dens, linetype=dist, color=dist), linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0,25)) +
  xlab("Profit ($B)") + ylab("Density") + theme_bw() +
  guides(linetype = guide_legend(byrow = TRUE)) +
  theme(legend.title = element_blank(), legend.spacing.y = unit(10,"pt"))
if(writeOutput) ggsave("fig-fits-us-zoom.pdf", width = 6.5, height = 5)

# ROA & EROC charts
psLog <- 0.5
lim <- c(-1e2, 200)
dfFit <- dfParUS %>% filter(name %in% c("roa", "eroc"), type=="est") %>%
  select(-type) %>% group_by(name, dist) %>% nest() %>%
  mutate(out = list({
    pars <- as.list(data[[1]]$value)
    names(pars) <- data[[1]]$par
    lim2 <- asinh(lim / 2/ psLog)
    tibble(u2 = seq(lim2[1], lim2[2], length.out=501),
           u = 2 * psLog * sinh(u2),
           d = do.call(paste0("d", dist), c(list(x=u), pars)),
           d2 = d * sqrt((u)^2 + (2*psLog)^2))})) %>%
  ungroup() %>% select(-data) %>% unnest(out) %>%
  filter(d2 > 1e-5) %>% mutate(
    dist = factor(dist, levels = dfDist$dist, labels = dfDist$name),
    name = factor(name, levels = c("roa", "eroc"), labels = c("Return on assets", "Excess return on capital")))
dfHist <- dfDataUS %>% filter(name %in% c("roa", "eroc")) %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=70, psLog=psLog))) %>% ungroup() %>% select(-u) %>% unnest(hist) %>%
  mutate(name = factor(name, levels = c("roa", "eroc"), labels = c("Return on assets", "Excess return on capital")))
ggplot() + facet_wrap(~name, ncol = 1) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), dfHist, size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2, color = dist, linetype = dist), dfFit, linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-3,1e-1,1e1), labels=c(expression(10^-3), expression(10^-1), expression(10^1))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,1,0)),0,10^c(0,1,2)),
                     minor_breaks = NULL) +
  coord_cartesian(ylim = c(5e-4,10), xlim = c(-100,100)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("Profitability (pseudo-log)") + ylab("Density (log)") + theme_bw() +
  guides(linetype = guide_legend(byrow = TRUE)) +
  theme(legend.title = element_blank(), legend.spacing.y = unit(10,"pt"))
if(writeOutput) ggsave("fig-fits-roa.pdf", width = 6.5, height = 7)

# ROA & EROC charts: zoom around zero
lim <- c(-0.1,0.1)
dfFit <- dfParUS %>% filter(name %in% c("roa", "eroc"), type=="est") %>% select(-type) %>%
  group_by(name, dist) %>% nest() %>%
  mutate(out = list({
    pars <- as.list(data[[1]]$value)
    names(pars) <- data[[1]]$par
    lim2 <- asinh(lim / 2/ psLog)
    tibble(x = seq(lim[1], lim[2], length.out=2e3),
           dens = do.call(paste0("d", dist), c(list(x=x), pars)))})) %>%
  ungroup() %>% select(-data) %>% unnest(out) %>% mutate(
    dist = factor(dist, levels = dfDist$dist, labels = dfDist$name),
    name = factor(name, levels = c("roa", "eroc"), labels = c("Return on assets", "Excess return on capital")))
dfHist <- dfDataUS %>% filter(name %in% c("roa", "eroc")) %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=50, lim=lim))) %>% ungroup() %>% select(-u) %>% unnest(hist) %>%
  mutate(name = factor(name, levels = c("roa", "eroc"), labels = c("Return on assets", "Excess return on capital")))
ggplot() + facet_wrap(~name, ncol = 1) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), dfHist, size=0.1, linewidth = 0.3) +
  geom_line(data = dfFit, mapping = aes(x=x, y=dens, linetype=dist, color=dist), linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("Profitability") + ylab("Density") + theme_bw() +
  guides(linetype = guide_legend(byrow = TRUE)) +
  theme(legend.title = element_blank(), legend.spacing.y = unit(10,"pt"))
if(writeOutput) ggsave("fig-fits-roa-zoom.pdf", width = 6.5, height = 5)

# Tables
dfParUS %>% filter(par == "AIC") %>% select(-c(par, type)) %>% mutate(
  dist = factor(dist, levels = dfDist$dist, labels = dfDist$name2),
  name = factor(name, levels = c("ltp", "epn", "roa", "eroc"), labels = c("LTP", "EP", "ROA", "EROC"))) %>%
  arrange(dist, name) %>% pivot_wider() %>%
  xtable(digits = 0) %>% print(format.args = list(big.mark = ","), include.rownames = F)

dfParUS %>% filter(type != "stat") %>% mutate(
  dist = factor(dist, levels = dfDist$dist, labels = dfDist$name2),
  name = factor(name, levels = c("ltp", "epn", "roa", "eroc"), labels = c("LTP", "EP", "ROA", "EROC"))) %>%
  arrange(dist, name) %>%
  group_by(dist, name, par) %>% summarize(value = sprintf("%.3f (%.3f)", value[1], value[2])) %>%
  ungroup() %>% pivot_wider() %>% select(-dist) %>% xtable() %>% print(include.rownames = F)
dfParUS %>% filter(par == "n")
dfParUS %>% filter(name == "ltp", dist == "nln", type != "stat") %>% pivot_wider(names_from = type)

# Contour plot of NLN fit
pars <- dfParUS %>% filter(name == "ltp", dist == "nln", type == "est") %>% pull(var=value, name=par)
dfLab <- expand_grid(expx=10^(-2:2),y=c(-1,1)) %>% add_row(expx=0.003,y=0) %>% mutate(ltp=y*expx)
expand_grid(x=seq(log(1e-3),log(200), length.out=101), y=seq(-3,2.5, length.out=101)) %>%
  mutate(expx = exp(x), ltp=y*expx,
         p = with(as.list(pars), dmvnorm(cbind(x,y), mean=c(mux,muy),
                                         sigma= rbind(c(sig^2, rho*sig), c(rho*sig, 1))))) %>%
  ggplot(aes(x=expx, y=y)) +
  geom_contour_filled(aes(z=p), bins=20) +
  geom_contour(aes(z=ltp), breaks=c(-10^(2:-2),0,10^(-2:2)), color="red") +
  geom_label(data = dfLab, aes(label=ltp)) +
  scale_x_log10(name = "exp(X)") + ylab("Y") +
  scale_fill_grey(guide="none", start=1, end=0) +
  theme_bw()
if(writeOutput) ggsave("fig-nln-contour.pdf", width = 4, height = 3)

#### NLN fit of various data sets ####

psLog <- 0.02
dfData <- bind_rows(
  dfLt %>% group_by(geo) %>% summarize(u = list(ltep)) %>%
    mutate(type = "geo") %>% rename(name = geo),
  dfP %>%
    mutate(period = cut(fyear, breaks = yrs, labels = paste(yrs[-5], yrs[-5]+4, sep = "-"), right = F)) %>%
    group_by(geo, period, gvkey) %>% summarize(ltep = sum(epd)) %>% summarize(u = list(ltep)) %>%
    ungroup() %>% mutate(name = paste(geo, period), type = paste0("geo5",geo)) %>% 
    select(-c(geo, period)),
  dfLt %>% filter(geo == "USA") %>% group_by(giccd2) %>% summarize(u = list(ltep)) %>%
    left_join(dfGics2) %>% mutate(name = ind2, type = "ind") %>% select(-c(giccd2, ind2)),
  tibble(type = "measure", name = "LIVA", u = list(liva)),
  tibble(type = "measure", name = "Wealth", u = list(wealth/1e3)),
  tibble(type = "sup", name = "altp", u = list(altp)))
if(calcFit) {
  dfData <- dfData %>% rowwise() %>% mutate(par = list(mle.nln(u))) %>% ungroup() %>%
      mutate(name = fct_inorder(name))
  dfPar <- dfData %>% select(-u) %>% unnest(par)
  if(writeOutput) write_csv(dfPar, "fits-nln.csv")
} else {
  dfPar <- read_csv("fits-nln.csv")
  dfData <- dfPar %>% group_by(name, type) %>% nest(par = par:se) %>%
    ungroup() %>% inner_join(dfData)
}
names <- dfPar %>% filter(par == "rho.sig") %>% arrange(est) %>% pull(name)
dfHist <- dfData %>% rowwise() %>% mutate(hist = list(calc.hist(u, psLog=psLog))) %>%
  ungroup() %>% select(-c(u,par)) %>% unnest(hist)
dfFit <- dfData %>% rowwise() %>% mutate(fit = list(calc.dens(par, c(-500,1000), 500, psLog))) %>%
  ungroup() %>% select(-c(u,par)) %>% unnest(fit)

ggplot() + facet_wrap(~name) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), dfHist, size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2), dfFit, linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-5,1e-3,1e-1), labels=c(expression(10^-5), expression(10^-3), expression(10^-1))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-5,2)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("LTP ($B, pseudo-log)") + ylab("Density (log)") + theme_bw()
if(writeOutput) ggsave("fig-fits-all.pdf", width = 16, height = 8)

# Charts by geography

dfData %>% filter(type == "geo") %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=40, psLog=psLog))) %>%
  ungroup() %>% select(-c(u,par)) %>% unnest(hist) %>%
  ggplot() + facet_wrap(~name, ncol = 1) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2), filter(dfFit, type == "geo"), linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-4,1e-2,1e-0),
                labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-4,1), xlim = c(-200,400)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("LTP (€B/£B/$B, pseudo-log)") + ylab("Density (log)") + theme_bw()
if(writeOutput) ggsave("fig-fits-geo.pdf", width = 6, height = 8)

dfData %>% filter(type == "geo", name == "USA") %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=40, psLog=psLog))) %>%
  ungroup() %>% select(-c(u,par)) %>% unnest(hist) %>%
  ggplot() +
  geom_pointrange(aes(x = u, y = d / sqrt(u^2 + (2*psLog)^2), ymin = dmin/ sqrt(u^2 + (2*psLog)^2), ymax = dmax/ sqrt(u^2 + (2*psLog)^2)), size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d), filter(dfFit, type == "geo", name == "USA"), linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-6, 1e-4,1e-2,1e-0),
                labels=c(expression(10^-6), expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-6,1e1), xlim = c(-200,400)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("LTP ($B, pseudo-log)") +
  ylab(expression(paste("Untransformed density (",`$B`^-1,", log)"))) +
  theme_bw()
if(writeOutput) ggsave("fig-fit-us-untransformed.pdf", width = 6, height = 3)

dfData %>% filter(substr(type, 1, 4) == "geo5") %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=20, psLog=psLog))) %>%
  ungroup() %>% select(-c(u,par)) %>% unnest(hist) %>%
  ggplot() + facet_wrap(~name, ncol = 3, dir = "v") +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2), filter(dfFit, substr(type, 1, 4) == "geo5"), linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-4,1e-2,1e-0),
                labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(xlim = c(-200,500), ylim = c(1e-4,2)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("LTP (B, pseudo-log)") + ylab("Density (log)") + theme_bw()
if(writeOutput) ggsave("fig-fits-geo5yr.pdf", width = 6, height = 8)

# Charts by industry

# Cumulative all industries
dfLt %>% filter(geo == "USA") %>% arrange(giccd2, -ltep) %>%
  group_by(giccd2) %>% 
  mutate(n=row_number(), cum = cumsum(ltep)/1e3) %>% ungroup() %>%
  bind_rows(tibble(n=0, cum=0, giccd2=unique(dfLt$giccd2)), .) %>%
  left_join(dfGics2) %>%
  mutate(ind2 = factor(ind2, levels=names),
         firstNeg = ifelse((ltep < 0) & (lag(ltep) > 0), n, NA)) %>%
  ggplot(aes(x=n, y=cum)) + facet_wrap(~ind2, ncol = 3) +
  geom_line() + geom_vline(aes(xintercept = firstNeg), linetype=2) +
  xlab("Number of firms") + ylab("Cumulative LTP ($T)") + theme_bw() +
  scale_x_continuous(labels = function(x) scales::comma(x), breaks = c(0,1000,2000))

# Cumulative selected industries
dfLt %>% filter(geo == "USA") %>% arrange(giccd2, -ltep) %>%
  group_by(giccd2) %>% 
  mutate(n=row_number(), cum = cumsum(ltep)/1e3,
         nrel = n/last(n), cumrel = cum/last(cum)) %>% ungroup() %>%
  bind_rows(tibble(n=0, cum=0, nrel=0, cumrel=0, giccd2=unique(dfLt$giccd2)), .) %>%
  filter(giccd2 %in% c(45,30,25)) %>% left_join(dfGics2) %>%
  mutate(ind2 = factor(ind2, levels=names)) %>%
  ggplot(aes(x=nrel, y=cumrel, color=ind2, linetype=ind2)) +
  geom_line() + 
  xlab("Number of firms (% of total)") + ylab("Cumulative LTP (% of total)") +
  scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) +
  scale_color_brewer(palette = "Set1", direction = -1) + theme_bw() +
  theme(legend.title = element_blank())
if(writeOutput) ggsave("fig-us-cum-ind.pdf", width = 5, height = 2.5)

dfData %>% filter(type == "ind") %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=20, psLog=psLog))) %>%
  ungroup() %>% select(-c(u,par)) %>% unnest(hist) %>%
  mutate(name = factor(name, levels=names)) %>%
  ggplot() + facet_wrap(~name, ncol = 3) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2), linewidth = 0.3, 
            data = dfFit %>% filter(type == "ind") %>% mutate(name = factor(name, levels=names))) +
  scale_y_log10(breaks = c(1e-4,1e-2,1e-0),
                labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-4,1), xlim = c(-200,500)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("LTP ($B, pseudo-log)") + ylab("Density (log)") + theme_bw()
if(writeOutput) ggsave("fig-fits-ind.pdf", width = 6, height = 5)

# LIVA and Wealth measure charts

dfData %>% filter(type == "measure") %>% rowwise() %>%
  mutate(hist = list(calc.hist(u, n=50, psLog=psLog))) %>%
  ungroup() %>% select(-c(u,par)) %>% unnest(hist) %>%
  ggplot() + facet_wrap(~name, ncol = 1) +
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2), filter(dfFit, type == "measure"), linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-5,1e-3,1e-1), labels=c(expression(10^-5), expression(10^-3), expression(10^-1))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-5,0.3), xlim = c(-1e3,1e3)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("Shareholder value ($B, pseudo-log)") + ylab("Density (log)") + theme_bw()
if(writeOutput) ggsave("fig-fits-measures.pdf", width = 6, height = 6)

# Overall parameter chart

dfVars <- tibble(levels = c("mux", "sig", "muy", "p", "rho", "r", "rho.sig"),
                 labels = c(
                   expression(paste("Mean log size ", (mu[italic(X)]))),
                   expression(paste("St dev log size ", (sigma))),
                   expression(paste("Conditional ", (mu[italic(Y)]))),
                   expression(paste("Probability positive profit ", (italic(p)))),
                   expression(paste("Correlation size-return "(rho))),
                   expression(paste("Unconditional ", (mu[italic(Y)]+rho*sigma))),
                   expression(paste("Asymmetry ", (rho*sigma)))))

dfPar %>% 
  filter(par %in% c("muy", "r", "rho.sig"), name != "altp") %>%
  mutate(par = factor(par, levels = dfVars$levels, labels = dfVars$labels),
         name = fct_rev(factor(name, levels = names))) %>%
  ggplot(aes(x=name, y=est, ymin = est-1.96*se, ymax = est+1.96*se)) +
  facet_grid(type~par, scales = "free", space = "free", labeller = label_parsed) +
  geom_hline(yintercept = 0, color = "grey50") + geom_pointrange(size=0.1, linewidth = 0.3) + 
  ylab("Mean returns (estimate and 95% interval)") + xlab(element_blank()) +
  coord_flip() + theme_bw() +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank())
if(writeOutput) ggsave("fig-est.pdf", width = 6.5, height = 6)

dfPar %>% 
  filter(par %in% c("mux", "sig", "rho"), name != "altp") %>%
  mutate(par = factor(par, levels = dfVars$levels, labels = dfVars$labels),
         name = fct_rev(factor(name, levels = names))) %>%
  ggplot(aes(x=name, y=est, ymin = est-1.96*se, ymax = est+1.96*se)) +
  facet_grid(type~par, scales = "free", space = "free_y", labeller = label_parsed) +
  geom_pointrange(size=0.1, linewidth = 0.3) + 
  ylab("Mean returns (estimate and 95% interval)") + xlab(element_blank()) +
  coord_flip() + theme_bw() +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank())
if(writeOutput) ggsave("fig-est2.pdf", width = 6.5, height = 6)

dfPar %>% filter(type == "geo", par %in% c("muy", "r", "rho.sig")) %>% mutate(pm = se*1.96)
dfPar %>% filter(substr(type, 1, 4) == "geo5", par == "rho.sig")
dfPar %>% filter(substr(type, 1, 4) == "geo5", par == "rho.sig") %>%
  group_by(type) %>% summarize(est = mean(est), se = mean(se))

# Industry table
dfInd <- dfLt %>% filter(geo == "USA") %>% left_join(dfGics2) %>%
  mutate(ind2 = factor(ind2, levels=names)) %>% select(ind2, conm, ltep) %>%
  group_by(ind2)
dfIndTable <- bind_cols(
    dfInd %>% slice_max(ltep, n=5) %>% arrange(ind2, -ltep),
    dfInd %>% slice_min(ltep, n=5) %>% arrange(ind2, ltep) %>%
      ungroup() %>% select(-ind2)) %>% ungroup()
dfIndTable %>% select(-ind2) %>% xtable(digits = 1) %>%
  print(include.rownames = F)
dfInd %>% arrange(ind2, -ltep) %>% group_by(ind2) %>% 
  summarize(ltep.tot = sum(ltep), ltep5 = sum(ltep[1:5]))

# Supplementary analysis annual LTP

dfLt %>% filter(geo == "USA") %>%
  ggplot(aes(x = nyr)) + geom_histogram(bins = 20) +
  scale_x_continuous(name = "Number of years in sample", breaks = c(0,10,20)) +
  scale_y_continuous(name = "Number of firms", minor_breaks = NULL) +
  theme_bw()

altp <- dfLt %>% filter(geo == "USA") %>% mutate(altp = ltep / nyr) %>% pull(altp)
altpFit <- dfPar %>% filter(name == "altp") %>% calc.dens(c(-100,100), 500, psLog)
altpHist <- calc.hist(altp, n=50, psLog=psLog)

ggplot() + 
  geom_pointrange(aes(x = u, y = d, ymin = dmin, ymax = dmax), altpHist, size=0.1, linewidth = 0.3) +
  geom_line(aes(x = u, y = d2), altpFit, linewidth = 0.3) +
  scale_y_log10(breaks = c(1e-4,1e-2,1), labels=c(expression(10^-4), expression(10^-2), expression(10^0))) + 
  scale_x_continuous(trans=scales::pseudo_log_trans(sigma=psLog), breaks=c(-10^(c(2,0)),0,10^c(0,2)),
                     minor_breaks = c(-1000,-10,-0.1,0.1,10,1000)) +
  coord_cartesian(ylim = c(1e-4,1), xlim = c(-100,100)) +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "dashed", "solid")) +
  scale_color_brewer(palette = "Set1") +
  xlab("LTP per year ($B, pseudo-log)") + ylab("Density (log)") + theme_bw()
if(writeOutput) ggsave("fig-altp.pdf", width = 3, height = 2)
