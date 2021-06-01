## Model parameter updates for calibration process ##
if(file.exists(paste0("Inputs/Calib_par_table.rds"))){
  calib.par  <- readRDS(paste0("Inputs/Calib_par_table.rds"))
} else if (!file.exists(paste0("Inputs/Calib_par_table.rds"))){
  library(FME)
  CalibPar <- read.xlsx(WB, sheet="CalibPar")
  parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  row.names(parRange) <- CalibPar$par
  set.seed(4122021)
  calib.par <- Latinhyper(parRange, 10000)
  calib.par <- data.frame(calib.par)
  saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))
}

nm.calp <- names(calib.par)

calib.rs.table <- matrix(0, nrow = nrow(calib.par), ncol = 15+length(nm.calp))
colnames(calib.rs.table) <- c("index", "seed", nm.calp, 
                              "od.death16", "od.death17", "od.death18", "od.death19",
                              "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                              "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19", "gof")
calib.rs.table[ , "index"] <- c(1:nrow(calib.par))
calib.rs.table[ , 3:(2+length(nm.calp))] <- as.matrix(calib.par)
for (ss in 1:nrow(calib.par)){
  print(paste0("Parameter set: ", ss))
  calib.seed = seed + ss
  calib.rs.table[ss, "seed"] = calib.seed
  for (pp in 1:length(nm.calp)){
    assign(nm.calp[pp], calib.par[ss,pp])
  }
  ## Fentanyl use status for initial population determined externally (allow to vary) ##
  n.opioid <- sum(init.pop$curr.state != "NODU")
  n.noud   <- sum(init.pop$curr.state == "NODU")
  # determine fentanyl use among initial population who use opioids 
  set.seed(calib.seed)
  fx         <- sample(0:1, size = n.opioid, prob = c(1-ini.OUD.fx, ini.OUD.fx), replace = T)
  init.pop$fx[init.pop$curr.state != "NODU"] <- fx
  # determine fentanyl use among initial population who use stimulants (non-opioid)
  set.seed(calib.seed*2)
  fx         <- sample(0:1, size = n.noud, prob = c(1-ini.NOUD.fx, ini.NOUD.fx), replace = T)
  init.pop$fx[init.pop$curr.state == "NODU"] <- fx
  
  # Overdose probability matrix (per month)
  od.matrix             <- matrix(0, nrow = 4, ncol = 2)
  rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
  colnames(od.matrix)   <- c("first", "subs")
  od.matrix["preb", "subs"]   <- od.preb.sub
  od.matrix["il.lr", "subs"]  <- od.il.lr.sub
  od.matrix["il.hr", "subs"]  <- od.il.lr.sub * multi.hr
  od.matrix["NODU", "subs"]   <- od.NODU.sub
  od.matrix[ , "first"]       <- od.matrix[ , "subs"] / multi.sub
  
  sim_sq    <- MicroSim(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = calib.seed)        # run for no treatment
  
  calib.rs.table[ss, c("od.death16", "od.death17", "od.death18", "od.death19")] <- colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5))[-5]
  calib.rs.table[ss, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")] <- (colSums(matrix(sim_sq$m.oddeath.fx, nrow=12, ncol=5)) / colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5)))[-5]
  calib.rs.table[ss, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")] <- colSums(matrix(sim_sq$m.EDvisits, nrow=12, ncol=5))[-5]
}

Target   <- read.xlsx(WB, sheet="Target")
tar.data <- Target$pe

for (ss in 1:nrow(calib.rs.table)){
  model.pred <- calib.rs.table[ss, c("od.death16", "od.death17", "od.death18", "od.death19",
                                     "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                                     "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")]
  gof <- 0
  for (j in 1:length(tar.data)){
    if (j %in% c(5:8)){
      gof <- gof + (abs(model.pred[j]*100 - tar.data[j]*100)/(tar.data[j]*100))/length(tar.data)
    } else {
      gof <- gof + (abs(model.pred[j] - tar.data[j])/tar.data[j])/length(tar.data)
    }
  }
  calib.rs.table[ss, "gof"] <- gof
}

sorted.mx <- calib.rs.table[order(calib.rs.table[ , "gof"], decreasing = F), ]

library(xlsx)
write.xlsx(sorted.mx, file="Calibration_preliminary.xlsx", 
           col.names = T, row.names = F)

calib.sp <- 20
calib.result.mx <- sorted.mx[1:calib.sp, ]

par(mfrow=c(1,3))

md.oddeath <- apply(calib.result.mx[ , c("od.death16", "od.death17", "od.death18", "od.death19")], 2, median)
ymax       <- max(calib.result.mx[ , c("od.death16", "od.death17", "od.death18", "od.death19")])
plot(x = 2016:2019,  tar.data[1:4], col = "black" , pch = 18, xlab="Year", ylab="Number of opioid overdose deaths", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', ylim = c(0, 500), frame.plot = FALSE)
# title("A", adj = 0, line =-0.5)
for (i in 1:calib.sp){
  lines(x = 2016:2019, y = calib.result.mx[i, c("od.death16", "od.death17", "od.death18", "od.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd=2)
}
lines(x = 2016:2019, y = md.oddeath, col = "red", lwd=3)
points(x = 2016:2019,  tar.data[1:4], col = "black" , pch = 16, cex=1.2, cex.axis=0.95)
axis(1, at = 2016:2019, pos=0, lwd.ticks=0, cex.axis=1.2)
# axis(side=1, pos=0, lwd.ticks=0)
abline(h=0)
# legend(x="bottomright",
#        col=c("dodgerblue","firebrick2", "lightgrey"), 
#        lwd=c(1.2, 1, 0.5), lty = c(1, NA, NA), pch=c(16, 16,16), pt.cex = c(1,1.1,1), cex = 0.8, bty = "n")


md.fxoddeath <- apply(calib.result.mx[ ,c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], 2, median)
ymax       <- max(calib.result.mx[ , c("fx.death16", "fx.death17", "fx.death18", "fx.death19")])
plot(x = 2016:2019,  tar.data[5:8], col = "black" , pch = 18, xlab="Year", ylab="Proportion of OOD deaths with fentanyl present", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', yaxt='n', ylim = c(0, 1), frame.plot = FALSE)
# title("A", adj = 0, line =-0.5)
for (i in 1:calib.sp){
  lines(x = 2016:2019, y = calib.result.mx[i, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd=2)
}
lines(x = 2016:2019, y = md.fxoddeath, col = "red", lwd=3)
points(x = 2016:2019,  tar.data[5:8], col = "black" , pch = 16, cex=1.2, cex.axis=0.95)
axis(1, at = 2016:2019, pos=0, lwd.ticks=0, cex.axis=1.2)
axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c("0", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.2)
# axis(side=1, pos=0, lwd.ticks=0)
abline(h=0)


md.edvisits <- apply(calib.result.mx[ ,c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], 2, median)
ymax       <- max(calib.result.mx[ , c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")])
plot(x = 2016:2019,  tar.data[9:12], col = "black" , pch = 18, xlab="Year", ylab="Number of ED visists for opioid overdose", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', ylim = c(0, 2500), frame.plot = FALSE)
# title("A", adj = 0, line =-0.5)
for (i in 1:calib.sp){
  lines(x = 2016:2019, y = calib.result.mx[i, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], col = adjustcolor("indianred1", alpha = 0.2), lwd=2)
}
lines(x = 2016:2019, y = md.edvisits, col = "red", lwd=3)
points(x = 2016:2019,  tar.data[9:12], col = "black" , pch = 16, cex=1.2, cex.axis=1.2)
axis(1, at = 2016:2019, pos=0, lwd.ticks=0, cex.axis=1.2)
# axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c("0", "20%", "40%", "60%", "80%", "100%"))
# axis(side=1, pos=0, lwd.ticks=0)
abline(h=0)

par(mfrow=c(1,1))
legend("top", 
       legend = c("Target", "Model: median", "Model: simulation"), 
       col = c("black",  "red", adjustcolor("indianred1", alpha = 0.2)), 
       pch = c(16, NA, NA), 
       lty = c(NA, 1, 1), 
       lwd=3,
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = T)


library(ggplot2)
calib.post  <- calib.result.mx[ , 3:(2+length(nm.calp))]
par  <- rep(colnames(calib.post), 2)
case <- c(rep("prior", length(nm.calp)), rep("posterior", length(nm.calp)))
pe   <- rep(0, length(nm.calp)*2)
lend <- rep(0, length(nm.calp)*2)
uend <- rep(0, length(nm.calp)*2)
ggplot.data <- data.frame(par = par, case = case, pe = pe, lend = lend, uend = uend)
CalibPar <- read.xlsx(file = "Inputs/MasterTable.xlsx", sheetIndex = 20)
ggplot.data$pe[ggplot.data$case == "prior"]   <- CalibPar$pe
ggplot.data$lend[ggplot.data$case == "prior"] <- CalibPar$pe - CalibPar$lower
ggplot.data$uend[ggplot.data$case == "prior"] <- CalibPar$upper - CalibPar$pe

ggplot.data$pe[ggplot.data$case == "posterior"]   <- apply(calib.post, 2, median)
ggplot.data$lend[ggplot.data$case == "posterior"] <- apply(calib.post, 2, median) - apply(calib.post, 2, quantile, probs = 0.025)
ggplot.data$uend[ggplot.data$case == "posterior"] <- apply(calib.post, 2, quantile, probs = 0.975) - apply(calib.post, 2, median)

ggplot.data$case <- factor(ggplot.data$case, levels = c("prior", "posterior"))

p<- ggplot(ggplot.data, aes(x=case, y=pe, color=case)) + 
  geom_point()+
  geom_errorbar(aes(ymin=pe-lend, ymax=pe+uend), width=.2, position=position_dodge(0.05)) +
  facet_wrap(~ par, scales = "free_y", ncol = 5) +
  labs(y="Value", x="") + theme_bw()
