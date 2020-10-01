## ----setup, include = FALSE, cache=FALSE--------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.asp = 0.618
)

## ----packages, message=FALSE, warning=TRUE------------------------------------
library(ggplot2)
library(dplyr)
library(sdmTMB)

## ----glimpse-pcod-------------------------------------------------------------
glimpse(pcod)

## ----spde, fig.asp=0.8--------------------------------------------------------
pcod_spde <- make_spde(pcod$X, pcod$Y, n_knots = 100)

## ----model--------------------------------------------------------------------
m <- sdmTMB(
  data = pcod, 
  formula = density ~ 0 + as.factor(year),
  time = "year", spde = pcod_spde, 
  family = tweedie(link = "log"),
  threshold_parameter = "depth_scaled",
  threshold_function = "logistic")

## ----output-------------------------------------------------------------------
df = data.frame("par"=names(m$sd_report$value), "value"=as.numeric(m$sd_report$value),
  "sd"=as.numeric(m$sd_report$sd))
df

