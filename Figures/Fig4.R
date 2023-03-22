# Fig 4
rm(list = ls())
library('pacman')
pacman::p_load(survival,ranger,ggplot2,dplyr,ggfortify)

tab <- as.tibble(read.csv(here("Figures","MI.csv")))

tab <- tab %>% mutate(percent = if_else(suma < quantile(suma, 0.25), 1, 0) )

survobj <- Surv(time=tab$age, event=tab$state, type="right")
cox_fit <- survfit(survobj ~ as.factor(tab$pr) )
cox_fit


autoplot(cox_fit)
