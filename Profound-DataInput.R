########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Module for Data Input of the Profound Naloxone distribution model: 
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: Dec 30, 2020
# Last update: Jan 31, 2021
#
#############################################################################
###################        DATA INPUT function      #########################
#############################################################################

## Define input excel file ##
WB        <- loadWorkbook(paste0("Inputs/MasterTable.xlsx"))

## Parameters for initial cohort ##
InitialPop    <- read.xlsx(WB, sheet="InitialPop")
pop.size      <- round(with(InitialPop, pe[par == "pop.size"]) * with(InitialPop, pe[par == "prop.12older"]), 0)   #size of initial population 12 and older
prev.oud      <- with(InitialPop, pe[par == "prev.oud"])                 #prevalence of OUD (or at risk for OUD)
prev.NODU.m   <- with(InitialPop, pe[par == "prev.NODU" & sex == "m"])   #prevalence of non-opioid drug use among males, in addition to OUD
prev.NODU.f   <- with(InitialPop, pe[par == "prev.NODU" & sex == "f"])   #prevalence of non-opioid drug use among females, in addition to OUD

Demographic   <- read.xlsx(WB, sheet="Demographic")
demo.mx       <- data.matrix(Demographic[ , 4:ncol(Demographic)])
v.rgn         <- colnames(Demographic)[-c(1:3)]                          #region names (city/town)
v.demo.sex    <- Demographic$sex
v.demo.race   <- Demographic$race
v.demo.age    <- Demographic$age

OUDDemo       <- read.xlsx(WB, sheet="OUDPrevNSDUH")$pe
StimDemo      <- read.xlsx(WB, sheet="StimPrevNSDUH")$pe

OpioidPattern <- read.xlsx(WB, sheet="OpioidPattern")
ini.il.m      <- with(OpioidPattern, pe[par == "ini.il" & sex == "m"])   # % of illicite opioid use among OUD pop
ini.il.f      <- with(OpioidPattern, pe[par == "ini.il" & sex == "f"])   # % of illicite opioid use among OUD pop
ini.il.hr.m   <- with(OpioidPattern, pe[par == "ini.il.hr" & sex == "m"])# % of high-risk among illicit opioid pop
ini.il.hr.f   <- with(OpioidPattern, pe[par == "ini.il.hr" & sex == "f"])# % of high-risk among illicit opioid pop
ini.inact     <- with(OpioidPattern, pe[par == "ini.inact"])
ini.OUD.fx    <- with(OpioidPattern, pe[par == "ini.OUD.fx"])
ini.everod.preb  <- with(OpioidPattern, pe[par == "ini.everod" & group == "preb"])
ini.everod.il.lr <- with(OpioidPattern, pe[par == "ini.everod" & group == "il.lr"])
ini.everod.il.hr <- with(OpioidPattern, pe[par == "ini.everod" & group == "il.hr"])

ini.NOUD.fx      <- read.xlsx(WB, sheet="StimulantPattern")$pe

initials      <- list(pop.size = pop.size, prev.oud = prev.oud, prev.NODU.m = prev.NODU.m, prev.NODU.f = prev.NODU.f, demo.mx = demo.mx, v.rgn = v.rgn, OUDDemo = OUDDemo, StimDemo = StimDemo, 
                      ini.il.m = ini.il.m, ini.il.f = ini.il.f, ini.il.hr.m = ini.il.hr.m, ini.il.hr.f = ini.il.hr.f, ini.inact = ini.inact, ini.OUD.fx = ini.OUD.fx,
                      ini.everod.preb = ini.everod.preb, ini.everod.il.lr = ini.everod.il.lr, ini.everod.il.hr = ini.everod.il.hr, ini.NOUD.fx = ini.NOUD.fx)


## Parameters for microsimulation ##
#life table: for mortality
mor.bg.y      <- read.xlsx(WB, sheet="LifeTable")$pe
mor.bg        <- 1-(1-mor.bg.y/1000000)^(1/12)
mor.drug.y    <- read.xlsx(WB, sheet="LifeTable")$drug
mor.drug      <- 1-(1-mor.drug.y/1000000)^(1/12)
mor.gp        <- read.xlsx(WB, sheet="LifeTable")$age
#risk of overdose
OverdoseRisk  <- read.xlsx(WB, sheet="OverdoseRisk")
od.preb.sub   <- with(OverdoseRisk, pe[par == "od.preb.sub"])
od.il.lr.sub  <- with(OverdoseRisk, pe[par == "od.il.lr.sub"])
od.NODU.sub   <- with(OverdoseRisk, pe[par == "od.NODU.sub"])
multi.hr      <- with(OverdoseRisk, pe[par == "multi.hr"])
multi.fx      <- with(OverdoseRisk, pe[par == "multi.fx"])
multi.relap   <- with(OverdoseRisk, pe[par == "multi.relap"])
multi.sub     <- with(OverdoseRisk, pe[par == "multi.sub"])
multi.NODU.fx <- with(OverdoseRisk, pe[par == "multi.NODU.fx"])
#transition probability
TransProb     <- read.xlsx(WB, sheet="TransProb")
p.preb2il.lr  <- with(TransProb, pe[par == "p.preb2il.lr"])
p.preb2inact  <- with(TransProb, pe[par == "p.preb2inact"])
p.il.lr2il.hr <- with(TransProb, pe[par == "p.il.lr2il.hr"])
p.il.lr2inact <- with(TransProb, pe[par == "p.il.lr2inact"])
p.il.hr2il.lr <- with(TransProb, pe[par == "p.il.hr2il.lr"])
p.il.hr2inact <- with(TransProb, pe[par == "p.il.hr2inact"])
p.inact2relap <- with(TransProb, pe[par == "p.inact2relap"])


## Parameters for decision tree ##
OD_loc_priv   <- read.xlsx(WB, sheet="ODSettingEMS")$private
OD_loc_pub    <- read.xlsx(WB, sheet="ODSettingEMS")$public
OD_loc        <- rbind(OD_loc_priv, OD_loc_pub)
rownames(OD_loc) <- c("priv", "pub")
colnames(OD_loc) <- v.rgn
DecisionTree  <- read.xlsx(WB, sheet="DecisionTree")
OD_wit_priv   <- with(DecisionTree, pe[par == "OD_wit" & group == "priv"])
OD_wit_pub    <- with(DecisionTree, pe[par == "OD_wit" & group == "pub"])
OD_911_priv   <- with(DecisionTree, pe[par == "OD_911" & group == "priv"])
OD_911_pub    <- with(DecisionTree, pe[par == "OD_911" & group == "pub"])
OD_hosp       <- with(DecisionTree, pe[par == "OD_hosp"])
OD_cess       <- with(DecisionTree, pe[par == "OD_cess"])

Survival      <- read.xlsx(WB, sheet="Survival")
sur_bl        <- with(Survival, pe[par == "sur_bl"])
rr_sur_Nx     <- with(Survival, pe[par == "rr_sur_Nx"])
rr_sur_EMS    <- with(Survival, pe[par == "rr_sur_EMS"])

## Parameters for naloxone kits ##
r.LossExp     <- 1/with(read.xlsx(WB, sheet="NxKit"), pe[par == "LossExp"])
Low2Priv      <- with(read.xlsx(WB, sheet="NxKit"), pe[par == "Low2Priv"])
nlx.adj       <- with(read.xlsx(WB, sheet="NxKit"), pe[par == "nlx.adj"])
NxData        <- read.xlsx(WB, sheet="NxData")
array.Nx.full <- array(0, dim = c(length(unique(NxData$year)), length(unique(NxData$risk)), length(v.rgn)))
dimnames(array.Nx.full)[[1]] <- unique(NxData$year)
dimnames(array.Nx.full)[[2]] <- unique(NxData$risk)
dimnames(array.Nx.full)[[3]] <- v.rgn
for (i in 1:length(unique(NxData$year))){
  array.Nx.full[i, , ] <- data.matrix(subset(NxData, year == unique(NxData$year)[i])[-c(1,2)])
}
NxMvt         <- data.matrix(read.xlsx(WB, sheet="NxMvt")[,-1])
row.names(NxMvt) <- v.rgn

## Parameters for cost ##
Cost          <- read.xlsx(WB, sheet="Cost")
c.preb        <- with(Cost, pe[par == "c.preb"])
c.il.lr       <- with(Cost, pe[par == "c.il.lr"])
c.il.hr       <- with(Cost, pe[par == "c.il.hr"])
c.inact       <- with(Cost, pe[par == "c.inact"])
c.NODU        <- with(Cost, pe[par == "c.NODU"])
c.nlx.kit     <- with(Cost, pe[par == "c.nlx.kit"])
c.nlx.dtb     <- with(Cost, pe[par == "c.nlx.dtb"])
c.relap.v          <- numeric(0)                      # cost of remaining one cycle: relapsed, as the average of inactive and prior state
c.relap.v["preb"]  <- (c.preb + c.inact)/2
c.relap.v["il.lr"] <- (c.il.lr + c.inact)/2
c.relap.v["il.hr"] <- (c.il.hr + c.inact)/2

c.EMS         <- with(Cost, pe[par == "c.EMS"])
c.hospcare    <- with(Cost, pe[par == "c.hospcare"])