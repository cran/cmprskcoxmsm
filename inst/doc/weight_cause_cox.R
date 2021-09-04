### R code from vignette source 'weight_cause_cox.rnw'

###################################################
### code chunk number 1: weight_cause_cox.rnw:74-75
###################################################
options(width=80)


###################################################
### code chunk number 2: weight_cause_cox.rnw:134-136
###################################################
library(cmprskcoxmsm)
load("follic.RData")


###################################################
### code chunk number 3: weight_cause_cox.rnw:141-157
###################################################
## Change the treatment name
follic$treatment <- ifelse(follic$ch=="Y","Combination treatment","Radiation alone")

## Distribution of the treatment
table(follic$treatment)

## make the stage as the character variable:
follic$clinstg <- ifelse(follic$clinstg==1,"Stage I","Stage II")

## Generate the weight:
OUT1 <- doPS(dat = follic,
             Trt = "treatment",
             Trt.name = "Combination treatment",
             VARS. = c("age","hgb","clinstg"))

follic1 <- OUT1[["Data"]]


###################################################
### code chunk number 4: weight_cause_cox.rnw:162-163 (eval = FALSE)
###################################################
## plot(OUT1)


###################################################
### code chunk number 5: weight_cause_cox.rnw:181-185
###################################################
follic1$status.1 <- NA
follic1$status.1[which(follic$status==0)] <- "No response"
follic1$status.1[which(follic$status==1)] <- "Disease relapse"
follic1$status.1[which(follic$status==2)] <- "Death"


###################################################
### code chunk number 6: weight_cause_cox.rnw:188-196
###################################################
tab1 <- weight_cause_cox(follic1,
                         time = "time",
                         time2 = NULL,
                         Event.var = "status.1",
                         Event = "Disease relapse",
                         weight.type = "Stabilized",
                         ties = NULL)
tab1


###################################################
### code chunk number 7: weight_cause_cox.rnw:205-217 (eval = FALSE)
###################################################
## cif.dr <- cif_est(follic1,
##                   time = "time",
##                   time2 = NULL,
##                   Event.var = "status.1",
##                   Events = c("Disease relapse","Death"),
##                   cif.event = "Disease relapse",
##                   weight.type = "Stabilized",
##                   ties = NULL,
##                   risktab = TRUE,
##                   risk.time = 10)
## cif_dr <- cif.dr$cif_data
## risk_dr10 <- cif.dr$risk_tab


###################################################
### code chunk number 8: weight_cause_cox.rnw:220-225
###################################################
cif_dr <- read.csv("cif_dr.csv")[,-1]
risk_dr10 <- read.csv("risk_dr10.csv")[,-1]
colnames(risk_dr10) <- c("Risk Difference (95\\% CI)",
                         "Risk Ratio (95\\% CI)")
rownames(risk_dr10) <- c("time: 10")


###################################################
### code chunk number 9: weight_cause_cox.rnw:229-232
###################################################
plot_est_cif(cif.dat = cif_dr,
             color = c("#1c9099","#756bb1"),
             ci.cif = TRUE)


###################################################
### code chunk number 10: weight_cause_cox.rnw:236-237
###################################################
risk_dr10


###################################################
### code chunk number 11: weight_cause_cox.rnw:241-253 (eval = FALSE)
###################################################
## cif.death <- cif_est(follic1,
##                   time = "time",
##                   time2 = NULL,
##                   Event.var = "status.1",
##                   Events = c("Disease relapse","Death"),
##                   cif.event = "Death",
##                   weight.type = "Stabilized",
##                   ties = NULL,
##                   risktab = TRUE,
##                   risk.time = 10)
## cif_death <- cif.death$cif_data
## risk_death10 <- cif.death$risk_tab


###################################################
### code chunk number 12: weight_cause_cox.rnw:256-261
###################################################
cif_death <- read.csv("cif_death.csv")[,-1]
risk_death10<- read.csv("risk_death10.csv")[,-1]
colnames(risk_death10) <- c("Risk Difference (95\\% CI)",
                         "Risk Ratio (95\\% CI)")
rownames(risk_death10) <- c("Time: 10")


###################################################
### code chunk number 13: weight_cause_cox.rnw:265-268
###################################################
plot_est_cif(cif.dat = cif_death,
             color = c("#2c7fb8","#f03b20"),
             ci.cif = TRUE)


###################################################
### code chunk number 14: weight_cause_cox.rnw:273-274
###################################################
risk_death10


