########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Module for Data Input of the Profound Naloxone distribution model: 
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: Dec 30, 2020
# Last update: Jan 11, 2020
#
#############################################################################
###################        DATA INPUT function      #########################
#############################################################################

## Define input excel file ##
WB        <- loadWorkbook(paste0("Inputs/MasterTable.xlsx"))

## Parameters for initial cohort
prev.oud       <- read.xlsx(WB, sheet="PrevOUD")$pe               #prevalence of OUD (ar at risk for OUD)
rgn            <- colnames(read.xlsx(WB, sheet="InitialPop"))[-1] #region names (city/town)
ini.pop.p      <- read.xlsx(WB, sheet="InitialPop")               #demographics for initial population
prop.preb      <- read.xlsx(WB, sheet="DrugPattern")   
prop.illicit.m <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "illicit" & sex == "m")$pe
prop.illicit.f <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "illicit" & sex == "f")$pe
prop.il.hr.m   <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "il.hr" & sex == "m")$pe
prop.il.hr.f   <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "il.hr" & sex == "f")$pe
prop.inact     <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "inact")$pe
prop.fx.op     <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "fx")$pe
prop.ever.od.preb   <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "ever.od")$pe[1]
prop.ever.od.il.lr  <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "ever.od")$pe[1]
prop.ever.od.il.hr  <- subset(read.xlsx(WB, sheet="OpioidPattern"), Group == "ever.od")$pe[1]


prev.NODU.m    <- subset(read.xlsx(WB, sheet="StimulantPattern"), Group == "NODU" & sex == "m")$pe
prev.NODU.f    <- subset(read.xlsx(WB, sheet="StimulantPattern"), Group == "NODU" & sex == "f")$pe
prop.fx.op     <- subset(read.xlsx(WB, sheet="StimulantPattern"), Group == "fx")$pe

