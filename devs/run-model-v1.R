# This is a development script for a delay-difference model in rtmb
# The first draft attempts to reproduce the iscam delay difference model from
#   the 2018 Pacific Cod stock assessment for Haida Gwaii/Queen Charlotte
#   Sound stock in Pacific Canada (Area 5ABCD) [published in 2020]

# **For this first version, this is a ONE group, ONE area, ONE sex model
# Therefore, we can strip out a lot of the indexing in the iscam inputs**

# Package name: DDRTMB

# Authors: Robyn Forrest and Catarina Wor (Pacific Biological Station)

# Date created: May 15, 2024
# Last Modified: May 15, 2024

# Notes:
# - The iscam input files are already loaded into the package:
#   1. pcod2020dat (the input data)
#   2. pcod2020ctl (controls for running the model,
#      including prior descriptions and misc settings)
#   3. pcod2020pfc (controls for projections and decision tables)
#
# - The input files include some inputs for the iscam age structured model
#       that are not used in the delay-difference implementation. The input
#       files will eventually be customised
# - Iscam uses an errors-in-variables approach to partition observation error,
#       with multiplicative weighting of each index observation using annual CVs
# - Want to change this to additive weightings (as per SS3), then also explore
#       state space options. But first try to reproduce the iscam results!

# TODO:
# CHECK how to set dat$alloc if more than one commercial gear
# Can probably relax requirement of setting nmeanwt to 1 when no mean weight data
# Put the dat, ctl and pfc stripping below into functions
# Change parameter values to log
# TIDY UP the three recruitment parameters - currently set to all be the same as per Paul Starr's request
# Probably don't need all the counters

# Document and build package (these are also buttons in Rstudio)
#    this will incorporate new functions saved in the R and data folders
devtools::document()
devtools::load_all()

library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up data and parameter controls. Better to make this a function, or even just
#  make it part of the package.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Parameters and controls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


