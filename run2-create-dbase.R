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

library(tidyverse)
options("scipen"=100)

dataDir <- "."
years <- c(2002,2021)
excludeGvkeys <- c("015520","011644","034290") # Companies with data integrity issues
excludeGics <- c("40", "60")  # Exclude financials and real estate from analysis
minN <- 10 # minimum number of companies in an industry for beta calculation
mrp <- 0.05 # Market risk premium; see Koller et al. Ch 15, p. 313
betaDebt <- 0.15 # Beta of debt; see Koller et al. Ch 15, p. 320
writeOutput <- T

dfFedRaw <- read_csv(file.path(dataDir, "rf-fed.csv"))
dfBundRaw <- read_csv(file.path(dataDir, "rf-bund.csv"))
dfGiltRaw <- read_csv(file.path(dataDir, "rf-gilt.csv"))
dfFundRaw <- read_csv(file.path(dataDir, "compu-fund-90-22.csv"))
dfEurRaw <- read_csv(file.path(dataDir, "compu-eur-00-22.csv"))
dfSecRaw <- read_csv(file.path(dataDir, "compu-secm-00-22.csv"))
dfGics <- read_csv(file.path(dataDir, "compu-gics.csv"))
dfGicsG <- read_csv(file.path(dataDir, "compu-g-gics.csv"))
dfGnames <- read_csv("db-gics.csv")

# Basic cleaning of databases:
dfGicsC <- bind_rows(dfGics, dfGicsG) %>%
  arrange(gvkey,indfrom) %>% group_by(gvkey) %>% summarize(giccd = last(ggroup)) %>%
  mutate(giccd = ifelse(giccd == 2540, 5020, giccd))
# GICS 5020 is a reclassification from 2540, see https://www.msci.com/documents/1296102/1339060/GICS+Structure+Revisions+in+2018.pdf/8c297da5-8160-4c2b-98e1-36998f6d2117

dfSec <- dfSecRaw %>%
  filter(
    !is.na(prccm + ajexm + trfm),
    pmin(prccm, ajexm, trfm) >= 0.01, # To prevent rounding issues
    !gvkey %in% excludeGvkeys) %>%  
  mutate(
    mcend = prccm * cshom / 1e9,
    month = as.numeric(format(datadate,"%Y%m")),
    year = as.numeric(format(datadate,"%Y"))) %>%
  inner_join(dfGicsC) 

# Firms in both NA and global database, to be excluded from EUR/GBP
dups <- intersect(unique(dfFundRaw$gvkey), unique(dfEurRaw$gvkey)) 

dfFund <- dfFundRaw %>%
  mutate(curcd = "USD") %>%
  bind_rows(dfEurRaw) %>%
  filter(icapt > 0, !is.na(icapt + ebit + txt + dltt + dlc + at),
         !(curcd != "USD" & gvkey %in% dups)) %>%
  arrange(gvkey, fyear) %>%
  mutate(month = as.numeric(format(datadate,"%Y%m"))) %>%
  left_join(dfSec %>% select(gvkey, month, mcend)) %>%
  inner_join(dfGicsC) %>%
  mutate(
    nopat = (ebit - txt) / 1e3,
    capend = icapt / 1e3,
    atend = at / 1e3,
    cons = c(F, ((gvkey == lag(gvkey)) & (fyear == lag(fyear) + 1)) [-1]),
    capbeg = ifelse(cons, lag(capend), NA),
    atbeg = ifelse(cons, lag(atend), NA),
    roc = nopat / capbeg,
    roa = nopat / atbeg) %>%
  filter(cons)
stopifnot(!anyDuplicated(dfFund %>% select(gvkey, fyear)))
xtabs( ~ fyear + factor(curcd), dfFund)

# Calculate TSRs
dfSec <- dfSec %>%
  inner_join(dfFund %>% select(gvkey, year = fyear)) %>%
  arrange(gvkey, month) %>%
  inner_join(dfSec %>% group_by(month) %>% summarize() %>% arrange(month) %>% mutate(monthno = row_number() - 1)) %>%
  group_by(gvkey) %>% mutate(
    cons = (monthno == lag(monthno) + 1),
    mcbeg = lag(mcend),
    tsr = (prccm * trfm / ajexm) / lag(prccm * trfm / ajexm) - 1) %>%
  ungroup() %>%
  filter(
    cons,
    !is.na(tsr + mcbeg),
    year >= years[1],
    year <= years[2]) %>%
  mutate(usdtsr = tsr * mcbeg)
stopifnot(!anyDuplicated(dfSec %>% select(gvkey, month)))

# Calculate market index
dfRf <-
  bind_rows(
    dfFedRaw %>% mutate(curcd = "USD"),
    dfBundRaw %>% mutate(curcd = "EUR"),
    dfGiltRaw %>% mutate(curcd = "GBP")) %>%
  filter(!is.na(yield)) %>%
  mutate(
    rf = yield / 100,
    month = as.numeric(format(date,"%Y%m"))) %>%
  group_by(curcd, month) %>% summarize(rfa = mean(rf)) %>% ungroup() %>%
  mutate(lrf = log1p(rfa) / 12, rf = expm1(lrf))

dfMkt <- dfSec %>% group_by(month) %>%
  summarize(
    mcbeg = sum(mcbeg),
    usdtsr = sum(usdtsr)) %>%
  mutate(
    mkttsr = usdtsr / mcbeg,
    lmkttsr = log1p(mkttsr),
    lmktex = cumsum(lmkttsr),
    index = exp(lmktex),
    mcindex = lead(mcbeg)/first(mcbeg),
    delta = index[month == years[2] * 100 + 12] / index) %>%
  left_join(dfRf %>% filter(curcd == "USD") %>% select(-curcd)) %>%
  mutate(lrmrf = lmkttsr - lrf)
dfMkt %>% summarize(across(c(lmkttsr, lrf, lrmrf), function(x) expm1(mean(x)*12)))

# Calculate LIVA (beta = 1)
dfSec <- dfSec %>%
  inner_join(dfMkt %>% select(month, mkttsr, delta)) %>%
  mutate(
    ler = log1p(tsr) - log1p(mkttsr),
    liva = (tsr - mkttsr) * mcbeg * delta,
    cf = (1 + tsr) * mcbeg - mcend,
    livac = cf * delta,  # Cash contribution to LIVA
    mcadj = mcbeg * delta)
stopifnot(abs(sum(dfSec$liva)) < 1e-10)

# LIVAs for checking; note this differs from global database, b/c calculating discount rate only using US companies
dfCoYr <- dfSec %>% group_by(gvkey, year) %>% summarize(mcend = last(mcend),mcbeg=first(mcbeg), liva = sum(liva)) %>% 
  left_join(dfFund %>% group_by(gvkey) %>% summarize(conm = last(conm)))
dfCo <- dfSec %>% group_by(gvkey) %>% summarize(mcend = last(mcend), liva = sum(liva)) %>%
  left_join(dfFund %>% group_by(gvkey) %>% summarize(conm = last(conm)))
dfInd <- dfSec %>% group_by(giccd) %>% summarize(liva = sum(liva)) %>% left_join(dfGnames)

# Calculate industry betas
dfIndP <- dfSec %>%
  group_by(giccd, month) %>% summarize(
    mcbeg = sum(mcbeg),
    usdtsr = sum(usdtsr),
    n = n()) %>%
  ungroup() %>%
  mutate(
    indtsr = usdtsr / mcbeg,
    lindtsr = log1p(indtsr),
    date = as.Date(paste0(month,"01"), format = "%Y%m%d")) %>%
  left_join(dfMkt %>% select(month, lrf, lrmrf)) %>%
  mutate(lrindrf = lindtsr - lrf) %>%
  filter(n >= minN)
dfInd <- dfIndP %>%
  group_by(giccd) %>%
  summarize(
    data = list({
      fit <- lm(lrindrf ~ lrmrf)
      tibble(alpha = coef(fit)[1], alpha.se = sqrt(vcov(fit)[1,1]), betaL = coef(fit)[2],  betaL.se = sqrt(vcov(fit)[2,2]))}),
    nCoMo = sum(n), nMo = n()) %>%
  unnest(data) %>%
  left_join(dfFund %>% filter(!is.na(mcend)) %>% group_by(giccd) %>%
              summarize(debt = sum((dltt+dlc)/1000), equity = sum(mcend))) %>%
  mutate(beta = (equity*betaL + debt*betaDebt) / (equity + debt)) %>% # Rearrangement of Koller et al, Ch15, p. 320
  left_join(dfGnames) 

dfRfYr <- dfRf %>% mutate(year = month %/% 100) %>%
  group_by(curcd, year) %>% summarize(rf = mean(rfa)) %>% ungroup() %>%
  filter(between(year, years[1], years[2]))
# Calculate cost of capital and discount factor by year
dfCc <- expand_grid(
  dfInd %>% select(giccd, beta) %>% filter(!substr(giccd,1,2) %in% excludeGics),
  dfRfYr
) %>%
  mutate(cc = rf + beta * mrp) %>%
  arrange(giccd, curcd, year) %>%
  group_by(giccd, curcd) %>% mutate(delta = prod(1+cc) / cumprod(1+cc)) %>% ungroup()

# Calculate long-term economic profit
dfP <- dfFund %>%
  inner_join(dfCc %>% select(giccd, curcd, fyear = year, cc, delta)) %>% # Removes excluded industries
  left_join(dfGnames) %>%
  mutate(
    epn = nopat - cc * capbeg, # Nominal
    epd = epn * delta,         # Discounted (w. cost of capital)
    capadj = capbeg * delta) %>%
  select(gvkey, curcd, conm, datadate, fyear, month, giccd, nopat, atbeg, atend, roa, capbeg, capend, capadj, roc, cc, delta, epn, epd)
if(writeOutput) write_csv(dfP, file.path(dataDir, "ltep-panel-02-21.csv"))

