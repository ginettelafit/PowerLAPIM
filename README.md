# PowerLAPIM
A Shiny application and R package to conduct power analysis for Longitudinal Actor-Partner Interdependence Models (L-APIMs) that include quadratic effects.

The repository contains functions used in the following [paper](https://psyarxiv.com/mnce4/):

Lafit, Ginette, Sels, Laura, Adolf, Janne K., Loeys, Tom, and Eva Ceulemans. (2021). “PowerLAPIM: An Application to Conduct Power Analysis for Longitudinal Actor-Partner Interdependence Models that Include Quadratic Effects” 

## Shiny app and R package to perform power analysis to select the number of dyads for the L-APIM in intensive longitudinal dyadic studies

Users can download the app and run locally on their computer by executing the following commands in R or Rstudio. 

```
# Check if R packages are installed

list.of.packages = c("nlme","MASS","tidyverse","future.apply","gridExtra","formattable","htmltools",
"shiny","DT","ggplot2","gridExtra","data.table","plyr","dplyr","tidyr","shinyjs")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(nlme)
library(MASS)
library(tidyverse)
library(future.apply)
library(gridExtra)
library(formattable)
library(htmltools)
library(shiny)
library(DT)
library(ggplot2)
library(gridExtra)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(shinyjs)

library(devtools)
devtools::install_github("ginettelafit/PowerLAPIM", force = T)

library(PowerLAPIM)

# Using Gist: users can launch this app with:
shiny::runGist('1d186b6d9bc76f5e41871ce40e5cee47')

```

## How the app works in a nutshell

The shiny app focuses on a set of research questions regarding the L-APIM. The application allows conducting power analysis to select the number of dyads for 32 L-APIM variants. We distinguish five different model categories: 

* L-APIMs with linear effects only
* L-APIMs with additional quadratic effects
* L-APIMs with group differences in the linear and quadratic effects
* L-APIMs including a continuous or dichotomous time-varying moderator
* L-APIMs including autoregressive effects

These categories are not mutually exclusive in that the models included in the application sometimes combine multiple category features (e.g., models with linear, quadratic, and autoregressive effects).  All models are implemented in a multilevel regression framework, which means that they distinguish differences between dyads from differences within dyads and partners. We restrict the application to models with random intercepts, capturing stable inter-dyadic differences, and fixed actor, partner, and moderation effects. A more detailed description of the models can be found in the supplementary material available on [OSF](https://osf.io/vtb9e/).

A hands-on tutorial for conducting simulation-based power analyses for the L-APIM can be fuond [here](https://psyarxiv.com/mnce4/).



