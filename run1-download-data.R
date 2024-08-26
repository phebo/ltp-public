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
library(RPostgres)
# Setting up Postgres: https://wrds-www.wharton.upenn.edu/pages/support/programming-wrds/programming-r/r-from-your-computer/
# Also note how to work with 2FA and Postgres: https://wrds-www.wharton.upenn.edu/pages/about/log-in-to-wrds-using-two-factor-authentication/multi-factor-authentication-for-sasconnect-and-postgresql/

username <- "YOUR_WRDS_USERNAME" # Need to update this with your WRDS details
password <- "YOUR_WRDS_PASSWORD"
dataDir <- "."
writeOutput <- T

# Risk-free rates
dfFedRaw <- read_csv("https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=DGS10&scale=left&cosd=2000-01-01&coed=2022-01-01&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Daily&fam=avg&fgst=lin&fgsnd=2022-01-01&line_index=1&transformation=lin&vintage_date=2022-01-01&revision_date=2022-01-01&nd=2000-01-01",
                     skip = 1, col_names = c("date", "yield"), col_types = "Dd")
if(writeOutput) write_csv(dfFedRaw, file.path(dataDir, "rf-fed.csv"))

dfBundRaw <- read_csv("https://api.statistiken.bundesbank.de/rest/download/BBK01/WT1010?format=csv&lang=en",
                      skip = 8, col_names = c("date", "yield"), col_types = "Dd-")
if(writeOutput) write_csv(dfBundRaw, file.path(dataDir, "rf-bund.csv"))

# Downloaded from https://www.bankofengland.co.uk/boeapps/database/fromshowcolumns.asp?Travel=NIxIRxSUx&FromSeries=1&ToSeries=50&DAT=RNG&FD=1&FM=Jan&FY=2000&TD=14&TM=Oct&TY=2022&FNY=&CSVF=TT&html.x=139&html.y=36&C=C6S&Filter=N#
dfGiltRaw <- read_csv(file.path(dataDir, "BoE-Database_export.csv"),
                      skip = 1, col_names = c("date", "yield")) %>%
  mutate(date = as.Date(date, format = "%d-%m-%Y"))
if(writeOutput) write_csv(dfGiltRaw, file.path(dataDir, "rf-gilt.csv"))

# WRDS session
wrds <- dbConnect(Postgres(), host='wrds-pgdata.wharton.upenn.edu', port=9737, dbname='wrds',
                  sslmode='require', user=username, password=password)

# Compustat annual fundamental data
res <- dbSendQuery(wrds, "select gvkey,conm,datadate,fyear,icapt,at,dltt,dlc,ebit,pi,xint,txt from comp.funda where
                     fyear between '1990' and '2022' and
                     indfmt = 'INDL' and
                     consol = 'C' and
                     popsrc = 'D' and
                     datafmt = 'STD' and
                     curcd = 'USD'")
dfFundRaw <- dbFetch(res, n=-1)
dbClearResult(res)
if(writeOutput) write_csv(dfFundRaw, file.path(dataDir, "compu-fund-90-22.csv"))

# Ibid Euro and GBP
res <- dbSendQuery(wrds, "select gvkey,conm,datadate,fyear,icapt,at,dltt,dlc,ebit,pi,xint,txt,curcd from comp.g_funda where
                     fyear between '2000' and '2022' and
                     indfmt = 'INDL' and
                     consol = 'C' and
                     popsrc = 'I' and
                     datafmt = 'HIST_STD' and
                     (curcd = 'EUR' or curcd = 'GBP')")
dfEurRaw <- dbFetch(res, n=-1)
dbClearResult(res)
if(writeOutput) write_csv(dfEurRaw, file.path(dataDir, "compu-eur-00-22.csv"))

# Compustat monthly security data
res <- dbSendQuery(wrds, "select gvkey,iid,datadate,ajexm,curcdm,prccm,trfm,cshom,exchg,fic,tpci
                          from comp.secm where
                          cyear between '2000' and '2022' and
                          primiss = 'P' and tpci = '0' and cshom > 0 and curcdm = 'USD'")
dfSecRaw <- dbFetch(res, n=-1)
dbClearResult(res)
if(writeOutput) write_csv(dfSecRaw, file.path(dataDir, "compu-secm-00-22.csv"))

# GICS codes by gvkey
res <- dbSendQuery(wrds, "select * from comp.co_hgic")
dfGics <- dbFetch(res, n=-1)
dbClearResult(res)
if(writeOutput) write_csv(dfGics, file.path(dataDir, "compu-gics.csv"))

res <- dbSendQuery(wrds, "select * from comp.g_co_hgic")
dfGicsG <- dbFetch(res, n=-1)
dbClearResult(res)
if(writeOutput) write_csv(dfGicsG, file.path(dataDir, "compu-g-gics.csv"))

# GICS names
res <- dbSendQuery(wrds, "select * from comp.r_giccd")
dfGnames <- dbFetch(res, n=-1)
dbClearResult(res)
if(writeOutput) write_csv(dfGnames, file.path(dataDir, "gics-names.csv"))
