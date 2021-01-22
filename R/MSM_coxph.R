###################################################################
##### Weighted Cause-Specific Cox Proportional Hazard model #######
###################################################################
weight_cause_cox <- function(data,
                             time,
                             time2 = NULL,
                             Event.var,
                             Event,
                             weight.type,
                             ties = NULL){

  ## treating both censoring and competing risk as censoring
  Event. <- ifelse(data[,Event.var]==Event,1,0)

  if(is.null(time2)){
    s1 <- Surv(time = data[,time], event = Event.)
  }else{
    s1 <- Surv(time = data[,time], time2 = data[,time2], event = Event.)
  }

  ## Weighted Cox ph:
  if(weight.type=="Unstabilized"){

    if(is.null(ties)){
      fit1 <- coxph(s1 ~ data[,"Trt"],
                    weights = data[,"ipw_ate_unstab"])
    }else{
      fit1 <- coxph(s1 ~ data[,"Trt"],
                    weights = data[,"ipw_ate_unstab"],ties = ties)
    }

  }else if(weight.type=="Stabilized"){

    if(is.null(ties)){
      fit1 <- coxph(s1 ~ data[,"Trt"],
                    weights = data[,"ipw_ate_stab"])
    }else{
      fit1 <- coxph(s1 ~ data[,"Trt"],
                    weights = data[,"ipw_ate_stab"],ties = ties)
    }

  }else{
    stop("Weights type misspecified")
  }


  sand_var <- as.numeric(sandwich(fit1))

  out <- summary(fit1)$coefficients[, c(1, 4, 5, 6), drop = FALSE]
  out[,2] <- sqrt(sand_var)
  out[,3] <- out[,1]/out[,2]
  out[,4] <- 2 * pnorm(abs(out[,3]),lower.tail=FALSE)
  CIs <- summary(fit1)$conf.int[, c(1, 3:4), drop = FALSE]
  CIs[, 1] <- paste('$',round(CIs[,1],3),'$',sep = "")
  CIs[, 2] <- paste('($', round(exp(as.numeric(out[,1] - qnorm(0.975) * out[,2])),3), '$, $',
                    round(exp(as.numeric(out[,1] + qnorm(0.975) * out[,2])),3),
                    '$)', sep = '')

  out <- cbind(out, CIs[, 1:2, drop = FALSE])
  out[, 1:3] <- paste('$', round(as.numeric(out[, 1:3]),3), '$', sep = '')
  out[, 4] <- pvalFormat(out[, 4])
  colnames(out) <- c('Estimate', 'Robust SE', 'z-value',
                     'p-value', 'Hazard Ratio', '95\\% CI')

  rownames(out) <- levels(data[,"Trt"])[2]
  return(out)

}


