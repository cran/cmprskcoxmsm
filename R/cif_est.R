###### Estimate the CIF

cif_est <- function(data,
                    time,
                    time2 = NULL,
                    Event.var,
                    Events,
                    cif.event,
                    weight.type,
                    ties = NULL,
                    risktab = TRUE,
                    risk.time = NULL){

  ### Create the event indicator for all the events:
  event.int <- matrix(NA,nrow = nrow(data),ncol = length(Events))
  for (i in 1:ncol(event.int)) {
    event.int[,i] <- ifelse(data[,Event.var]==Events[i],1,0)
  }

  ### fit the regression model for all event:
  cox.res <- list()
  for (i in 1:ncol(event.int)) {

    ## Generate the survival outcome for each event
    if(is.null(time2)){
      s1 <- Surv(time = data[,time], event = event.int[,i])
    }else{
      s1 <- Surv(time = data[,time], time2 = data[,time2], event = event.int[,i])
    }

    ## Weighted Cox ph for each event:
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

    }

    cox.res[[i]] <- fit1

   }

  ### Estimated the basline cumulative hazard of each event at each time
  ## \Lambda_j(t;a) = \Lambda_0j(t)*exp(\beta * a)
  bashaz.1 <- list()
  for (i in 1:ncol(event.int)) {
      coef.1 <- cox.res[[i]]$coefficients
      a <- basehaz(cox.res[[i]],centered = FALSE)
      a <- a[,c("time","hazard")]
      a$hazard.A <- a$hazard
      a$hazard.B <- a$hazard*exp(coef.1)
      bashaz.1[[i]] <- a
  }

  ### Estimated Survival function:
  ## S(t;a) = exp^(-(\Lambda_1(t;a) + \Lambda_1(t;a)))
  S.all <- matrix(NA,nrow = nrow(bashaz.1[[1]]),ncol = 3)
  S.all <- as.data.frame(S.all)
  colnames(S.all) <- c("time","S.all.A","S.all.B")
  S.all$time <- bashaz.1[[1]]$time
  S.all$S.all.A <- exp( - (Reduce("+",bashaz.1)$hazard.A) )
  S.all$S.all.B <- exp( - (Reduce("+",bashaz.1)$hazard.B) )

  ### Calculate the CIF: F_j(t;a) = sum( S(t;a)*\delta(Lambda_j(t;a)) )

  cif.1 <- matrix(NA,nrow = nrow(S.all),ncol = 3)
  cif.1 <- as.data.frame(cif.1)
  colnames(cif.1) <- c("time","cif.A","cif.B")
  cif.1$time <- S.all$time
  bashaz.2 <- bashaz.1[[which(Events==cif.event)]]
  for (i in 1:nrow(cif.1)) {
    if(i==1){
      cif.1$cif.A[i] <- S.all$S.all.A[i] * bashaz.2$hazard.A[i]
      cif.1$cif.B[i] <- S.all$S.all.B[i] * bashaz.2$hazard.B[i]
    }else{
      temp.A <- temp.B <- NULL
      for (j in 2:i) {
        temp.A<- c(temp.A,
                   S.all$S.all.A[j] * (bashaz.2$hazard.A[j] - bashaz.2$hazard.A[j-1]))
        temp.B<- c(temp.B,
                   S.all$S.all.B[j] * (bashaz.2$hazard.B[j] - bashaz.2$hazard.B[j-1]))
      }
      cif.1$cif.A[i] <- cif.1$cif.A[1] + sum(temp.A)
      cif.1$cif.B[i] <- cif.1$cif.B[1] + sum(temp.B)
    }
  }

  ### Get the confidence interval of CIF:
  ## using bootstrap
    ## bootstrap data
    boot.A <- matrix(NA,nrow = nrow(cif.1),ncol = 200)
    boot.B <- matrix(NA,nrow = nrow(cif.1),ncol = 200)
    rdiff.1<- matrix(NA,nrow = nrow(cif.1),ncol = 200)
    rratio.1 <- matrix(NA,nrow = nrow(cif.1),ncol = 200)

    ## do the bootstrap:
    for (k in 1:200) {
      id <- sample(1:nrow(data),nrow(data),replace = TRUE)
      data.boot <- data[id,]
      ### fit the regression model for all event:
      cox.res.boot <- list()
      for (i in 1:ncol(event.int)) {

        ## Generate the survival outcome for each event
        if(is.null(time2)){
          s1 <- Surv(time = data.boot[,time], event = event.int[,i][id])
        }else{
          s1 <- Surv(time = data.boot[,time], time2 = data.boot[,time2], event = event.int[,i][id])
        }

        ## Weighted Cox ph for each event:
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

        }

        cox.res.boot[[i]] <- fit1

      }


      ### Estimated the basline cumulative hazard of each event at each time
      ## \Lambda_j(t;a) = \Lambda_0j(t)*exp(\beta * a)
      bashaz.1.boot <- list()
      for (i in 1:ncol(event.int)) {
        coef.1 <- cox.res.boot[[i]]$coefficients
        a <- basehaz(cox.res.boot[[i]],centered = FALSE)
        a <- a[,c("time","hazard")]
        a$hazard.A <- a$hazard
        a$hazard.B <- a$hazard*exp(coef.1)
        bashaz.1.boot[[i]] <- a
      }

      ### Estimated Survival function:
      ## S(t;a) = exp^(-(\Lambda_1(t;a) + \Lambda_1(t;a)))
      S.all.boot <- matrix(NA,nrow = nrow(bashaz.1.boot[[1]]),ncol = 3)
      S.all.boot <- as.data.frame(S.all.boot)
      colnames(S.all.boot) <- c("time","S.all.boot.A","S.all.boot.B")
      S.all.boot$time <- bashaz.1.boot[[1]]$time
      S.all.boot$S.all.boot.A <- exp( - (Reduce("+",bashaz.1.boot)$hazard.A) )
      S.all.boot$S.all.boot.B <- exp( - (Reduce("+",bashaz.1.boot)$hazard.B) )

      ### Calculate the CIF: F_j(t;a) = sum( S(t;a)*\delta(Lambda_j(t;a)) )

      cif.1.boot <- matrix(NA,nrow = nrow(S.all.boot),ncol = 3)
      cif.1.boot <- as.data.frame(cif.1.boot)
      colnames(cif.1.boot) <- c("time","cif.A","cif.B")
      cif.1.boot$time <- S.all.boot$time
      bashaz.2 <- bashaz.1.boot[[which(Events==cif.event)]]
      for (i in 1:nrow(cif.1.boot)) {
        if(i==1){
          cif.1.boot$cif.A[i] <- S.all.boot$S.all.boot.A[i] * bashaz.2$hazard.A[i]
          cif.1.boot$cif.B[i] <- S.all.boot$S.all.boot.B[i] * bashaz.2$hazard.B[i]
        }else{
          temp.A <- temp.B <- temp.C <- temp.D <- NULL
          for (j in 2:i) {
            temp.A<- c(temp.A,
                       S.all.boot$S.all.boot.A[j] * (bashaz.2$hazard.A[j] - bashaz.2$hazard.A[j-1]))
            temp.B<- c(temp.B,
                       S.all.boot$S.all.boot.B[j] * (bashaz.2$hazard.B[j] - bashaz.2$hazard.B[j-1]))
          }
          cif.1.boot$cif.A[i] <- cif.1.boot$cif.A[1] + sum(temp.A)
          cif.1.boot$cif.B[i] <- cif.1.boot$cif.B[1] + sum(temp.B)
        }
      }

      ## use step function to estimate the cif of the event time
      a1 <- stepfun(cif.1.boot$time,c(0,cif.1.boot$cif.A))
      a2 <- stepfun(cif.1.boot$time,c(0,cif.1.boot$cif.B))

      boot.A[,k] <- a1(cif.1$time)
      boot.B[,k] <- a2(cif.1$time)

      rdiff.1[,k] <- a2(cif.1$time) - a1(cif.1$time)
      rratio.1[,k] <- a2(cif.1$time) / a1(cif.1$time)
    }

    ci.cif.A <- apply(boot.A,1,function(x){quantile(x,probs = c(0.025,0.975))})
    ci.cif.B <- apply(boot.B,1,function(x){quantile(x,probs = c(0.025,0.975))})

    cif.data <- cbind(cif.1$time,cif.1$cif.A,t(ci.cif.A),cif.1$cif.B,t(ci.cif.B))
    cif.data <- as.data.frame(cif.data)
    colnames(cif.data) <- c("time","cif.control","cif.control.lower","cif.control.upper",
                           "cif.exposure","cif.exposure.lower","cif.exposure.upper")

    if(risktab){
      tab1  <- matrix(NA,nrow = nrow(cif.1),ncol = 6)
      tab1  <- as.data.frame(tab1)
      tab1$V1 <- cif.1$cif.B - cif.1$cif.A
      tab1$V2 <- apply(rdiff.1,1,function(x){quantile(x,probs = 0.025,na.rm = TRUE)})
      tab1$V3 <- apply(rdiff.1,1,function(x){quantile(x,probs = 0.975,na.rm = TRUE)})
      tab1$V4 <- cif.1$cif.B / cif.1$cif.A
      tab1$V5 <- apply(rratio.1,1,function(x){quantile(x,probs = 0.025,na.rm = TRUE)})
      tab1$V6 <- apply(rratio.1,1,function(x){quantile(x,probs = 0.975,na.rm = TRUE)})

      tab1.fun <- list()
      for (j in 1:ncol(tab1)) {
        if(is.na(tab1[1,j])){
          tab1.fun[[j]] <- stepfun(cif.1$time,c(NA,tab1[,j]))
        }else{
          tab1.fun[[j]] <- stepfun(cif.1$time,c(0,tab1[,j]))
        }
      }

      risk.est <- NULL
      for (j in 1:ncol(tab1)) {
        risk.est[j] <- round(tab1.fun[[j]](risk.time),3)
      }

      table.1 <- matrix(NA,nrow = 1,ncol = 2)

      table.1[1,1] <- paste("$",risk.est[1],"$ ($",risk.est[2],",",risk.est[3],"$)",sep = "")
      table.1[1,2] <- paste("$",risk.est[4],"$ ($",risk.est[5],",",risk.est[6],"$)",sep = "")

      colnames(table.1) <- c("Risk Difference (95\\% CI)","Risk Ratio (95\\% CI)")
      rownames(table.1) <- paste("time: ",risk.time,sep = "")

      res <- list(cif.data,table.1)
      names(res) <- c("cif_data","risk_tab")
    }else{
      res <- list(cif.data)
      names(res) <- c("cif_data")
    }
    return(res)
}


plot_est_cif <- function(cif.data,
                         color = color,
                         ci.cif = FALSE){
  if(ci.cif){
    ggplot(cif.data, aes_(x = ~time, y = ~cif.control)) +
    geom_line(aes(colour=color[1]),size=1.4)+
    geom_line(aes_(x=~time, y=~cif.exposure,colour=color[2]),size=1.4)+
    ylim(c(0,1))+
    geom_ribbon(aes_(ymin=~cif.control.lower, ymax=~cif.control.upper), fill=color[1], alpha=.3, show.legend = F)+
    geom_ribbon(aes_(ymin=~cif.exposure.lower, ymax=~cif.exposure.upper), fill=color[2], alpha=.3, show.legend = F)+
    xlab("Time") +
    ylab("Cumulative incidence function")+
    scale_color_manual(name = "Group", values = color, labels = c("Control","Exposure"))+
    scale_linetype_manual(name = "Group", values = color, labels = c("Control","Exposure"))+
    theme(
      plot.title =   element_text(size=20, face="bold",hjust = 0.5),
      axis.title.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.x = element_text(face="bold",size=18),
      axis.text.y = element_text(face="bold",size=18),
      legend.title = element_text(size=18),
      legend.text = element_text(size=18),
      legend.position= c(0.85,0.85))
  }else{
    ggplot(cif.data, aes_(x = ~time, y = ~cif.control)) +
      geom_line(aes(colour=color[1]),size=1.4)+
      geom_line(aes_(x=~time, y=~cif.exposure,colour=color[2]),size=1.4)+
      ylim(c(0,1))+
      xlab("Time") +
      ylab("Cumulative incidence function")+
      scale_color_manual(name = "Group", values = color, labels = c("Control","Exposure"))+
      scale_linetype_manual(name = "Group", values = color, labels = c("Control","Exposure"))+
      theme(
        plot.title =   element_text(size=28, face="bold",hjust = 0.5),
        axis.title.x = element_text(size=28, face="bold"),
        axis.title.y = element_text(size=28, face="bold"),
        axis.text.x = element_text(face="bold",size=22),
        axis.text.y = element_text(face="bold",size=22),
        legend.title = element_text(size=22),
        legend.text = element_text(size=22),
        legend.position= c(0.8,0.8))
   }
}
