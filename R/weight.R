##################################################################
######################### Generate PS  ###########################
##################################################################

doPS <- function(data,
                 Trt,
                 Trt.name,
                 VARS.) {

  ## exclude all the missing value for the covariates
  tDat <- data[complete.cases(data[, c(Trt, VARS.)]), ]

  ## Make all the categorical covariates as factor
  for (i in 1:length(VARS.)) {
    if (paste0(class(tDat[, VARS.[i]]), collapse = '') == 'character') {
      tDat[, VARS.[i]] <- factor(tDat[, VARS.[i]])
      }
  }

  ### Make the group as numeric variable
  if(paste0(class(tDat[, Trt]), collapse = '') == 'factor'){
    tDat[, Trt] <- as.character(tDat[, Trt])
    tDat[, Trt] <- factor(tDat[, Trt],
                             levels = c(unique(tDat[, Trt])[which(unique(tDat[, Trt])!=Trt.name)],
                                        Trt.name))
  }else{
    tDat[, Trt] <- factor(tDat[, Trt],
                            levels = c(unique(tDat[, Trt])[which(unique(tDat[, Trt])!=Trt.name)],
                                       Trt.name))
  }

  GROUP <- paste(Trt,".1",sep = "")
  tDat[, GROUP] <- as.numeric(tDat[, Trt]) - 1
  tDat[, "Trt"] <- tDat[, Trt]

  ### Generate the propensity score
  ps1 <- ps(as.formula(paste(GROUP, ' ~', paste0(VARS., collapse = '+'))),
            data = tDat,
            n.trees = 30000, interaction.depth = 2,
            shrinkage = 0.01, perm.test.iters = 0,
            stop.method = c('es.mean'),
            estimand = 'ATE', verbose = FALSE)

  ps.1 <- as.numeric(unlist(ps1$ps))

  ## ATE
  tDat$ps_ate <- ps.1
  tDat$ipw_ate_unstab <- ( tDat[,GROUP] / ps.1 ) + ( (1 - tDat[,GROUP])/(1 - ps.1) )
  tDat$ipw_ate_stab <-   mean(tDat[, GROUP]) * tDat[, GROUP] / ps.1 + (1 - mean(tDat[, GROUP])) * (1 - tDat[, GROUP])/(1 - ps.1)

  res <- list(tDat,ps1)
  names(res) <- c("Data","PS")
  class(res) <- "PS"
  return(res)
}

##################################################################
#########################  Plot  #################################
##################################################################

#### Generate the SMD plot for the covariates in the PS model and the histogram of PS
plot.PS <- function(x,...){

  #### Histogram of the propensity score
  c1 <- "#2C7FB85A"
  c2 <- "#F03B205A"
  ax <- pretty(0:1,40)
  a <- x$Data
  t <- unique(a$Trt)
  hgA <- hist(a$ps_ate[which(a$Trt==t[1])], breaks = ax, plot = FALSE)
  hgB <- hist(a$ps_ate[which(a$Trt==t[2])], breaks = ax, plot = FALSE)
  par(xpd=TRUE)
  plot(hgA, col = c1,xlim = c(0,1),xlab = "Propensity Score",ylab = "count",
       cex.axis=2,cex.lab=1.5,xaxt="n",main = "Histogram")
  axis(1, at=seq(0,1,0.1), labels=seq(0,1,0.1),cex.axis=2)
  plot(hgB, col = c2,add = TRUE)
  legend("topright",
         col = c("#2c7fb8","#f03b20"),pch = 15,pt.cex = 3,
         legend = t,bty = 'n',text.width = 0.5,cex=0.8)

  #### SMD plot for the covariates in the PS model
  PS <- x$PS
  btps1 <- bal.table(PS)
  old_pars <- par(mar = c(4, 12, 1,1 ) + 0.1)
  on.exit(par(old_pars), add = TRUE)
  x1 <- btps1[[1]]$std.eff.sz
  tOrd <- order(x1)
  x2 <- btps1[[2]]$std.eff.sz
  tNames <- rownames(btps1[[1]])[tOrd]
  par(xpd=FALSE)
  plot(sort(x1), 1:length(x1), pch = 20, col = 'firebrick2', xlab = 'Standardized Mean Difference', ylab = '',
       yaxt = 'n',cex.axis=2,cex.lab=2)
  rect(-0.1, -100, 0.1, 100,
       col = 'gray93', border = NA)
  abline(v = 0, lty = 'dashed', col = 'gray')
  points(sort(x1), 1:length(x1), pch = 20, col = 'firebrick2',cex=2)
  points(x2[tOrd], 1:length(x2), pch = 18, col = 'darkslategray',cex=2)
  axis(side = 2, at = 1:length(x2), labels = tNames, las = 1, cex.axis = 2)
  legend('bottomright',
         legend = c('Weighted', 'Unweighted'),
         fill = c('darkslategray', 'firebrick3'),
         cex = 0.7, bty = 'n')
  box(lwd = 1.5)
  invisible()
}
