# R codes for epidemiology dynamics for one stage model for LOW & HIGH Infection and Recovery rates: Saikat Batabyal #

library(deSolve)
library(rootSolve)
library(bayestestR)

################################################################################
## ODE FRAMEWORK 
################################################################################

ode.solve=function(fit, fix, modelname, init.vector, timepoints, timename){ 
  allpars <- c(fit, fix)
  init.conds <- c()
  for(i in 1:length(init.vector)){
    init.conds <- c(init.conds, allpars[init.vector[[i]]])
  }
  soln=as.data.frame(lsoda(y=init.conds, times=timepoints, func=modelname, parms=allpars, rtol = 1e-04, atol = 1e-04))
  if (is.null(timename)) timename="years" 
  names(soln)[1]=timename
  return(soln)
}

###############################################################################
### ODE Model System
###############################################################################

mtb_epdm <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    {
      dS.dt = birth - mu*S - beta*S + gamma*L + delta_L*L + delta_L*I #Susceptible
      dI.dt = eta*beta*S + beta*(1-epsilon)*L - theta*I - beta*I - lambda*I - mu*I - delta_L*I #Progressive_Stage_1
      dL.dt = (1-eta)*beta*S + lambda*I - beta*(1-epsilon)*L - gamma*L - mu*L - delta_L*L #Non-progressive_Stage_1
      dA.dt = theta*I + beta*I - mu*A - delta_R*A #Final Progressive Stage
    }
    
    return(list(c(dS.dt, dI.dt, dL.dt, dA.dt)))})}


################################################################################
### EPDM dynamics for one stage model for LOW Infection and Recovery rates
################################################################################

################################################################################
### Part 1 (Before treatment starts)
### Parameters set of steady state for different epsilons
### In the file name S6 suggests SET no. 6, just a random pick
################################################################################

pars_S6_ss_eps_1 <- list(birth = 0.013,
                         mu = 0.013,
                         beta = 0.01,
                         gamma = 0.05,
                         eta = 0.05,
                         theta = 2,
                         epsilon = 1,
                         lambda = 1,
                         delta_R = 0.3,
                         delta_L = 0
)


# Repeat the same process for epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


################################################################################
### Steady state solutions
################################################################################

### epsilon=1 ###

init.guess_S6_ss_eps_1<-c(S = 4.1, I = 0.01, L = 0.01, A = 0.05)

fit_S6_ss_value_eps_1<- as.data.frame(steady(y = init.guess_S6_ss_eps_1, time = 0,
                                             func = mtb_epdm, parms = pars_S6_ss_eps_1, pos = TRUE))

summary(fit_S6_ss_value_eps_1)

fit_S6_ss_value_eps_1

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


###########################################################################################
### Sum over Progressive, Non-progressive and Active steady-state (Before treatment starts)
###########################################################################################

### epsilon=1 ###

sumP_S6_ss_eps_1 <- fit_S6_ss_value_eps_1[2,]
sumNP_S6_ss_eps_1 <- fit_S6_ss_value_eps_1[3,]
Active_S6_ss_eps_1 <- fit_S6_ss_value_eps_1[4,]

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


################################################################################
### Fixed Values (Before Treatment Starts)
################################################################################

init.vector <- c("S", "I", "L", "A")

fix_S6_ss_eps_1 <- c(S = fit_S6_ss_value_eps_1[1, ], 
                     I = fit_S6_ss_value_eps_1[2, ], 
                     L = fit_S6_ss_value_eps_1[3, ], 
                     A = fit_S6_ss_value_eps_1[4, ])

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


########################################################################################
### Fitted Parameters set of ODE system for different epsilons (Before Treatment Starts)
########################################################################################

fit_S6_ss_eps_1 <- c(birth = 0.013,
                     mu = 0.013,
                     beta = 0.01,
                     gamma = 0.05,
                     eta = 0.05,
                     theta = 2,
                     epsilon = 1,
                     lambda = 1,
                     delta_R = 0.3,
                     delta_L = 0)


### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


################################################################################
### Time-period (Before Treatment Starts)
################################################################################

mytimes_S6_ss_eps_1<- seq(0,100,0.01) #data taken from years 1-10

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### Solutions of ODE system (Before treatment starts)
################################################################################

sol_S6_ss_eps_1<-ode.solve(fit_S6_ss_eps_1,fix_S6_ss_eps_1,mtb_epdm,init.vector,mytimes_S6_ss_eps_1,"years")

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### Part 2 (During treatment)
### Parameters set for different epsilons
################################################################################

pars_S6_mstarts_eps_1 <- list(birth = 0.013,
                              mu = 0.013,
                              beta = 0.01,
                              gamma = 0.05,
                              eta = 0.05,
                              theta = 2,
                              epsilon = 1,
                              lambda = 1,
                              delta_R = 0.3,
                              delta_L = 1.5
)

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


################################################################################
### Steady state solutions (During treatment)
################################################################################

### epsilon=1 ##################################################################

init.guess_S6_mstarts_eps_1<-c(S = fit_S6_ss_value_eps_1[1, ], 
                               I = fit_S6_ss_value_eps_1[2, ], 
                               L = fit_S6_ss_value_eps_1[3, ], 
                               A = fit_S6_ss_value_eps_1[4, ])


fit_S6_ss_mstarts_eps_1<- as.data.frame(steady(y = init.guess_S6_mstarts_eps_1, time = 10,
                                               func = mtb_epdm, parms = pars_S6_mstarts_eps_1, pos = TRUE))

summary(fit_S6_ss_mstarts_eps_1)


### Dynamic run to steady-state (During Treatment Effect) ###################


fit_S6_mstarts_value_eps_1 <- as.data.frame(steady(y = c(S = fit_S6_ss_value_eps_1[1, ], 
                                                         I = fit_S6_ss_value_eps_1[2, ], 
                                                         L = fit_S6_ss_value_eps_1[3, ], 
                                                         A = fit_S6_ss_value_eps_1[4, ]),
                                                   time = c(0,10),
                                                   func = mtb_epdm, parms = pars_S6_mstarts_eps_1, method = "runsteady"))

summary(fit_S6_mstarts_value_eps_1)
fit_S6_mstarts_value_eps_1


### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


###########################################################################################
### Sum over Progressive, Non-progressive and Active steady-state (During Treatment Effect)
###########################################################################################

### epsilon=1 ### 

sumP_S6_mstarts_eps_1 <- fit_S6_mstarts_value_eps_1[2,]
sumNP_S6_mstarts_eps_1 <- fit_S6_mstarts_value_eps_1[3,]
Active_S6_mstarts_eps_1 <- fit_S6_mstarts_value_eps_1[4,]

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### Fixed Values (During treatment)
################################################################################

init.vector <- c("S", "I", "L", "A")

fix_S6_mstarts_eps_1 <- c(S = fit_S6_ss_value_eps_1[1, ], 
                          I = fit_S6_ss_value_eps_1[2, ], 
                          L = fit_S6_ss_value_eps_1[3, ], 
                          A = fit_S6_ss_value_eps_1[4, ])


### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

#################################################################################
### Fitted Parameters set of ODE system for different epsilons (During treatment)
#################################################################################

fit_S6_mstarts_eps_1 <- c(birth = 0.013,
                          mu = 0.013,
                          beta = 0.01,
                          gamma = 0.05,
                          eta = 0.05,
                          theta = 2,
                          epsilon = 1,
                          lambda = 1,
                          delta_R = 0.3,
                          delta_L = 1.5
)

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### Time-period (During treatment)
################################################################################

mytimes_S6_mstarts_eps_1<- seq(0,10,0.01) #data taken from years 1-10

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### SOLUTION (During treatment)
################################################################################

sol_S6_mstarts_eps_1<-ode.solve(fit_S6_mstarts_eps_1,fix_S6_mstarts_eps_1,mtb_epdm,init.vector,mytimes_S6_mstarts_eps_1,"years")

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### Part 3
### Fixed Values of ODE system for different epsilons (After stopping treatment)
################################################################################

init.vector <- c("S", "I", "L", "A") 

fix_S6_mstops_eps_1 <- c(S = fit_S6_mstarts_value_eps_1[1,], 
                         I = fit_S6_mstarts_value_eps_1[2,], 
                         L = fit_S6_mstarts_value_eps_1[3,], 
                         A = fit_S6_mstarts_value_eps_1[4,])


### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

#########################################################################################
### Fitted Parameters set of ODE system for different epsilons (After stopping treatment)
#########################################################################################

fit_S6_mstops_eps_1 <- c(birth = 0.013,
                         mu = 0.013,
                         beta = 0.01,
                         gamma = 0.05,
                         eta = 0.05,
                         theta = 2,
                         epsilon = 1,
                         lambda = 1,
                         delta_R = 0.3,
                         delta_L = 0
)

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###


################################################################################
### Time-period (After stopping treatment)
################################################################################

mytimes_S6_mstops_eps_1 <- seq(10,100,0.01) #data taken from years 1-10

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### SOLUTION (After stopping treatment)
################################################################################

sol_S6_mstops_eps_1 <-ode.solve(fit_S6_mstops_eps_1,fix_S6_mstops_eps_1,mtb_epdm,init.vector,mytimes_S6_mstops_eps_1,"years")

### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

#################################################################################
### Normalised Active Cases & CALCULATE area OVER (AOC) and UNDER (AUC) the curve
#################################################################################

################################################################################
### Epsilon = 1 ################################################################
################################################################################

norm_active_S6_ss_eps_1 <- (sol_S6_ss_eps_1[,5])/(Active_S6_ss_eps_1)
norm_active_S6_mstarts_eps_1 <- (sol_S6_mstarts_eps_1[,5])/(Active_S6_ss_eps_1)
norm_active_S6_mstops_eps_1 <- (sol_S6_mstops_eps_1[,5])/(Active_S6_ss_eps_1)


### Calculate AOC & AUC ########################################################

auc_starts_S6_eps_1 <- area_under_curve(mytimes_S6_mstarts_eps_1,norm_active_S6_mstarts_eps_1, method = "spline")

auc_starts_S6_eps_1

final_auc_1_S6_eps_1 <- 10 - auc_starts_S6_eps_1

final_auc_1_S6_eps_1

auc_mstops_bind_S6_eps_1 <-cbind.data.frame(norm_active_S6_mstops_eps_1,mytimes_S6_mstops_eps_1)
attach(auc_mstops_bind_S6_eps_1)

### PART 1 | delta_L=3, f=0.5

condition_1_S6_eps_1 <- norm_active_S6_mstops_eps_1 < 1 
auc_mstops_dataset_S6_condition_1_eps_1 <- auc_mstops_bind_S6_eps_1[1:(min(which(condition_1_S6_eps_1 == FALSE))-1),]

auc_stops_1_S6_eps_1 <- area_under_curve(auc_mstops_dataset_S6_condition_1_eps_1$mytimes_S6_mstops_eps_1,
                                         auc_mstops_dataset_S6_condition_1_eps_1$norm_active_S6_mstops_eps_1, method = "spline")
auc_stops_1_S6_eps_1

final_auc_2_S6_eps_1 <- (max(auc_mstops_dataset_S6_condition_1_eps_1$mytimes_S6_mstops_eps_1)-10) - auc_stops_1_S6_eps_1
final_auc_2_S6_eps_1

### PART 2 | delta_L=3, f=0.5

condition_2_S6_eps_1 <- norm_active_S6_mstops_eps_1 > 1 
auc_mstops_dataset_S6_condition_2_eps_1 <- auc_mstops_bind_S6_eps_1[(min(which(condition_2_S6_eps_1 == TRUE))):nrow(auc_mstops_bind_S6_eps_1),]


auc_stops_2_S6_eps_1 <- area_under_curve(auc_mstops_dataset_S6_condition_2_eps_1$mytimes_S6_mstops_eps_1,
                                         auc_mstops_dataset_S6_condition_2_eps_1$norm_active_S6_mstops_eps_1, method = "spline")
auc_stops_2_S6_eps_1
final_auc_3_S6_eps_1 <- auc_stops_2_S6_eps_1 - (100-min(auc_mstops_dataset_S6_condition_2_eps_1$mytimes_S6_mstops_eps_1))
final_auc_3_S6_eps_1

### Calculate Total Area

auc_decrease_S6_eps_1 <- final_auc_1_S6_eps_1 + final_auc_2_S6_eps_1
auc_increase_S6_eps_1 <- final_auc_3_S6_eps_1

total_auc_S6_eps_1 <- auc_decrease_S6_eps_1 - auc_increase_S6_eps_1 
total_auc_S6_eps_1


### Repeat again the same with epsilon=0.99/0.98/0.97/0.96/0.95/0.94/0.9/0.85/0.75/0.5/0.25/0 ###

################################################################################
### Percentage of Incidence Rate
################################################################################

### Epsilon=1 ##################################################################

incidence_rate_ss_set2_eps_1 <- data.frame(Time = sol_S6_ss_eps_1$years,
                                           S = sol_S6_ss_eps_1$S, 
                                           NP = sol_S6_ss_eps_1$L,
                                           P = sol_S6_ss_eps_1$I,
                                           A = sol_S6_ss_eps_1$A,
                                           percentage = (sol_S6_ss_eps_1$A/(sol_S6_ss_eps_1$S
                                                                            +sol_S6_ss_eps_1$I
                                                                            +sol_S6_ss_eps_1$L
                                                                            +sol_S6_ss_eps_1$A))*100)


min(incidence_rate_ss_set2_eps_1$percentage)


incidence_rate_mstarts_set2_eps_1 <- data.frame(Time = sol_S6_mstarts_eps_1$years,
                                                S = sol_S6_mstarts_eps_1$S, 
                                                NP = sol_S6_mstarts_eps_1$L,
                                                P = sol_S6_mstarts_eps_1$I,
                                                A = sol_S6_mstarts_eps_1$A,
                                                percentage = (sol_S6_mstarts_eps_1$A/(sol_S6_mstarts_eps_1$S
                                                                                      +sol_S6_mstarts_eps_1$I
                                                                                      +sol_S6_mstarts_eps_1$L
                                                                                      +sol_S6_mstarts_eps_1$A))*100)

min(incidence_rate_mstarts_set2_eps_1$percentage)


incidence_rate_mstops_set2_eps_1 <- data.frame(Time = sol_S6_mstops_eps_1$years,
                                               S = sol_S6_mstops_eps_1$S, 
                                               NP = sol_S6_mstops_eps_1$L,
                                               P = sol_S6_mstops_eps_1$I,
                                               A = sol_S6_mstops_eps_1$A,
                                               percentage = (sol_S6_mstops_eps_1$A/(sol_S6_mstops_eps_1$S
                                                                                    +sol_S6_mstops_eps_1$I
                                                                                    +sol_S6_mstops_eps_1$L
                                                                                    +sol_S6_mstops_eps_1$A))*100)

max(incidence_rate_mstops_set2_eps_1$percentage)

incidence_rate_set2_eps_1 <- data.frame(SS = min(incidence_rate_ss_set2_eps_1$percentage),
                                        Tstarts = min(incidence_rate_mstarts_set2_eps_1$percentage), 
                                        Tstops = max(incidence_rate_mstops_set2_eps_1$percentage))


### Repeat again the same for epsilon = 0 ###

################################################################################
### FIGURE for Low Transmission Rate (epsilon=1)
### ### Repeat again the same with epsilon = 0
################################################################################

{
par(mfrow = c(1,2), mar = c(3,4.1,1.8,1), mgp=c(2.0,0.7,0), 
      xaxs = "i", yaxs = "i", cex.lab = 1.2, cex.axis = 1.2)
  
################################################################################
### FIGURE 3B | Treatment Starts | Epsilon = 1
################################################################################
  
  plot(mytimes_S6_mstarts_eps_1,log10(sol_S6_mstarts_eps_1[,5]), 
       xlab='', ylab='',
       ylim = c(-5,2), xlim = c(-10,100),
       pch=1, tck = 0.02,
       type= 'l', col = "black", cex = 0.9, lwd = 3, lty = 1,
       yaxt="n", xaxt="n", bty = "l")
  
  mtext("Years after treating Mtb infection", side=1, line=1.9, cex=1.2)
  mtext("Relative no. of individual", side=2, line=2.8, cex=1.2)
  
  lines(mytimes_S6_mstarts_eps_1,log10(sol_S6_mstarts_eps_1[,3]), col = "black", cex = 0.9, lwd = 3, lty = 2)
  lines(mytimes_S6_mstarts_eps_1,log10(sol_S6_mstarts_eps_1[,4]), col = "black", cex = 0.9, lwd = 3, lty = 3)
  lines(mytimes_S6_mstarts_eps_1,log10(sol_S6_mstarts_eps_1[,2]), col = "black", cex = 0.9, lwd = 3, lty = 4)
  
  abline(h = log10(Active_S6_ss_eps_1), col = c("grey"),lwd = 2, lty = 5)
  
### FIGURE B | Steady States | Epsilon = 1 #####################################
  
   segments(x0 = -10, x1 = 0,
           y0 = log10(sumP_S6_ss_eps_1), y1 = log10(sumP_S6_ss_eps_1),
           lwd = 3, lty = 2, col = "black")
  segments(x0 = -10, x1 = 0,
           y0 = log10(sumNP_S6_ss_eps_1), y1 = log10(sumNP_S6_ss_eps_1),
           lwd = 3, lty = 3, col = "black")
  segments(x0 = -10, x1 = 0,
           y0 = log10(fit_S6_ss_value_eps_1[1,]), y1 = log10(fit_S6_ss_value_eps_1[1,]),
           lwd = 3, lty = 4, col = "black")
  
  segments(x0 = -10, x1 = 0,
           y0 = log10(sol_S6_ss_eps_1[,5]), y1 = log10(sol_S6_ss_eps_1[,5]),
           lwd = 3, lty = 1, col = "black")
  segments(x0 = -10, x1 = 0,
           y0 = log10(sol_S6_ss_eps_1[,3]), y1 = log10(sol_S6_ss_eps_1[,3]),
           lwd = 3, lty = 2, col = "black")
  segments(x0 = -10, x1 = 0,
           y0 = log10(sol_S6_ss_eps_1[,4]), y1 = log10(sol_S6_ss_eps_1[,4]),
           lwd = 3, lty = 3, col = "black")
  segments(x0 = -10, x1 = 0,
           y0 = log10(sol_S6_ss_eps_1[,2]), y1 = log10(sol_S6_ss_eps_1[,2]),
           lwd = 3, lty = 4, col = "black")
  
  mtext(expression(bold("B")), padj=-0.4, adj=-0.15, cex=1.4, font=1)
  
  segments(x0 = 0, x1 = 0,
           y0 = -5, y1 = 0.8,
           lwd = 3, lty = 3, col = "olivedrab")
  
  segments(x0 = 10, x1 = 10,
           y0 = -5, y1 = 0.8,
           lwd = 3, lty = 3, col = "olivedrab")
  
  ### Title for Start Treatment & End Treatment
  title1 = "Start"
  title2 = "End"
  mtext(side=3, line=-2.4, adj=0.07, cex=0.75, font = 2, title1)
  mtext(side=3, line=-2.4, adj=0.17, cex=0.75, font = 2, title2)
  
################################################################################
### Treatment Stops | Epsilon = 1
################################################################################
  
### Needed to merge the plots[par(new = TRUE)] #################################
  
  par(new = TRUE)
  lines(mytimes_S6_mstops_eps_1,log10(sol_S6_mstops_eps_1[,5]), 
        xlab='', ylab='',
        ylim = c(-5,2), xlim = c(-10,100),
        pch=1, tck = 0.02,
        type= 'l', col = "black", cex = 0.9, lwd = 3, lty = 1, 
        yaxt="n", xaxt="n", bty = "l")
  par(new = TRUE)
  lines(mytimes_S6_mstops_eps_1,log10(sol_S6_mstops_eps_1[,3]), col = "black", cex = 0.9, lwd = 3, lty = 2)
  par(new = TRUE)
  lines(mytimes_S6_mstops_eps_1,log10(sol_S6_mstops_eps_1[,4]), col = "black", cex = 0.9, lwd = 3, lty = 3)
  par(new = TRUE)
  lines(mytimes_S6_mstops_eps_1,log10(sol_S6_mstops_eps_1[,2]), col = "black", cex = 0.9, lwd = 3, lty = 4)
  
  # X-axis
  axis(1, las=1, at=seq(-10, 100, by=10), tck = 0.02, cex.axis=1.1, font.axis=1, 
       labels = TRUE)
  
  # Y-axis
  axis(2, las=2, at = seq(-5,2), tck = 0.02, cex.axis=1.1, font.axis=2,
       labels = expression(10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2))
  
  S <- expression("Susceptible")
  I <- expression("Progressive")
  L <- expression("Non-progressive")
  A <- expression("Active")
  
  legend(50, 1.6,
         legend = as.expression(c(S,L,I,A)),
         col = c("black", "black", "black", "black"),
         cex = 0.85, pch=c(-1,-1,-1,-1), lwd = c(2,2,2,2),
         lty=c(4,3,2,1), bty = "n", border = FALSE)
  
  legend(-15, 2.05, 
         expression(paste(beta*T==0.01, ",", ~~gamma==0.05, ",", ~~epsilon==1)),
         text.col = "black", text.font= 2, cex = 1, 
         bty = "n")
  
################################################################################
### FIGURE 3D | Linear Scale
################################################################################   
  
  plot(mytimes_S6_mstarts_eps_1,norm_active_S6_mstarts_eps_1, 
       xlab='', ylab='',
       ylim = c(0,2), xlim = c(-10,100),
       type = "l", lty = 1, col = "blue", cex = 0.9, lwd = 3,
       yaxt="n", xaxt="n", bty = "l")
  
  abline(h = 1, col = c("grey"),lwd = 2, lty = 5)
  
### Treatment Stops | Epsilon = 1 ##############################################
  
### Needed to merge the plots[par(new = TRUE)] #################################
  
  par(new = TRUE)
  lines(auc_mstops_dataset_S6_condition_1_eps_1$mytimes_S6_mstops_eps_1,
        auc_mstops_dataset_S6_condition_1_eps_1$norm_active_S6_mstops_eps_1, 
        col = "blue", cex = 0.9, lwd = 3, lty = 1)
  par(new = TRUE)
  lines(auc_mstops_dataset_S6_condition_2_eps_1$mytimes_S6_mstops_eps_1,
        auc_mstops_dataset_S6_condition_2_eps_1$norm_active_S6_mstops_eps_1,  
        col = "red", cex = 0.9, lwd = 3, lty = 1)  
  mtext("Years after treating Mtb infection", side=1, line=1.9, cex=1.2)
  mtext("Scaled TB cases", side=2, line=2.8, cex=1.2)
  
  mtext(expression(bold("D")), padj=-0.75, adj=-0.15, cex=1.4, font=1)
  
  segments(x0 = -10, x1 = 0,
           y0 = 1, y1 = 1,
           cex = 0.9, lwd = 3, lty = 1, col = "black")

  segments(x0 = 0, x1 = 0,
           y0 = 0, y1 = 1.6,
           lwd = 3, lty = 3, col = "olivedrab")
  
  segments(x0 = 10, x1 = 10,
           y0 = 0, y1 = 1.6,
           lwd = 3, lty = 3, col = "olivedrab")
  
  ### Title for Start Treatment & End Treatment
  
  title1 = "Start"
  title2 = "End"
  mtext(side=3, line=-2.4, adj=0.07, cex=0.75, font = 2, title1)
  mtext(side=3, line=-2.4, adj=0.172, cex=0.75, font = 2, title2)
  
  # X-axis
  
  axis(1, las=1, at=seq(-10, 100, by=10), tck = 0.02, cex.axis=1.1, font.axis=1, 
       labels = TRUE)
  
  # Y-axis
  axis(2, las=1, at = c(0, 0.25, 0.50, 0.75, 1, 1.25, 1.5, 1.75, 2), tck = 0.02,
       cex.axis=1, font.axis=1.2,
       labels = expression(0, 0.25, 0.50, 0.75, 1, 1.25, 1.5, 1.75, 2))
  
  ### Legend ###
  
  legend(-15, 2, 
         expression(paste(beta*T==0.01, ",", ~~gamma==0.05, ",", ~~epsilon==1)),
         text.col = "black", text.font= 2, cex = 1, 
         bty = "n")

  #### Text labeling Steady States ####
  text(x = -3, y = 1.04, label = "0.09%",
       col = "black",   # Color of the text
       font = 2,      # Bold face
       cex = 0.75)
  
  #### Text labeling Treatment Starts ####
  text(x = 18, y = 0.75, label = "0.07%",
       col = "blue",   # Color of the text
       font = 2,      # Bold face
       cex = 0.75)
  
  #### Text labeling Treatment Stops ####
  text(x = 20, y = 1.12, label = "0.1%",
       col = "red",   # Color of the text
       font = 2,      # Bold face
       cex = 0.75)
}

################################################################################
### Repeat the same for one stage model for HIGH Infection and Recovery rates
### Parametric values are given in Supplemental Table S1
################################################################################


################################################################################
### FIGURE 4
################################################################################

{
 par(mfrow = c(1,1), mar = c(3,4.1,1.8,1), mgp=c(2.0,0.7,0), 
      xaxs = "i", yaxs = "i", cex.lab = 1.2, cex.axis = 1.2)
  
################################################################################
### FIGURE 4B | Epsilon = 1 (PEAK)
################################################################################ 
  
  plot(mytimes_S6_mstarts_eps_1,norm_active_S6_mstarts_eps_1,
       xlab='', ylab='',
       ylim = c(0,2), xlim = c(-10,100),
       type = "l", lty = 5, col = "blue", cex = 0.9, lwd = 3,
       yaxt="n", xaxt="n", bty = "l")
  
  par(new = TRUE)
  lines(auc_mstops_dataset_S6_condition_1_eps_1$mytimes_S6_mstops_eps_1,
        auc_mstops_dataset_S6_condition_1_eps_1$norm_active_S6_mstops_eps_1, col = "blue", cex = 0.9, lwd = 3, lty = 5)
  par(new = TRUE)
  lines(auc_mstops_dataset_S6_condition_2_eps_1$mytimes_S6_mstops_eps_1,
        auc_mstops_dataset_S6_condition_2_eps_1$norm_active_S6_mstops_eps_1, col = "red", cex = 0.9, lwd = 3, lty = 5)

  ### Repeat again the same with epsilon=0.95/0.5/0 ###
  
  
  mtext("Years after treating Mtb infection", side=1, line=1.9, cex=1.4)
  mtext("Scaled TB cases", side=2, line=2.8, cex=1.4)
  
  mtext(expression(bold("A")), padj=-0.6, adj=-0.1, cex=1.6, font=1)
  
  abline(h = 1, col = c("grey"),lwd = 2, lty = 5)
  
  segments(x0 = -10, x1 = 0,
           y0 = norm_active_S6_ss_eps_1, y1 = norm_active_S6_ss_eps_1,
           cex = 0.9, lwd = 3, lty = 1, col = "black")
  
  segments(x0 = 0, x1 = 0,
           y0 = 0, y1 = 1.6,
           lwd = 3, lty = 3, col = "olivedrab")
  
  segments(x0 = 10, x1 = 10,
           y0 = 0, y1 = 1.6,
           lwd = 3, lty = 3, col = "olivedrab")
  
  ### Title for Start Treatment & End Treatment
  
  title1 = "Start"
  title2 = "End"
  mtext(side=3, line=-5, adj=0.07, cex=0.85, font = 2, title1)
  mtext(side=3, line=-5, adj=0.172, cex=0.85, font = 2, title2)
  
  # X-axis
  
  axis(1, las=1, at=seq(-10, 100, by=10), tck = 0.02, cex.axis=1.1, font.axis=1, 
       labels = TRUE)
  
  # Y-axis
  axis(2, las=1, at = c(0, 0.25, 0.50, 0.75, 1, 1.25, 1.5, 1.75, 2), tck = 0.02,
       cex.axis=1, font.axis=1.2,
       labels = expression(0, 0.25, 0.50, 0.75, 1, 1.25, 1.5, 1.75, 2))

  legend(-15, 2, 
         expression(paste(beta*T==0.01, ",", ~~gamma==0.05)),
         text.col = "black", text.font= 2, cex = 1.2, 
         bty = "n")
  
  E1 <- expression(paste(epsilon==1, ~("AOC"==2)))
  
  ### Repeat again the same with epsilon=0.95/0.5/0 ###
  
  legend(15, 0.85,
         legend = as.expression(c(E1)),
         col = c("blue"), cex = 1,
         pch=c(-1), lwd = c(2), lty=c(5),
         bty = "n", border = FALSE)
  
  E11 <- expression(paste(epsilon==1, ~("AUC"==1.6)))
  
  ### Repeat again the same with epsilon=0.95/0.5/0 ###
  
  legend(40, 1.2,
         legend = as.expression(c(E11)),
         col = c("red"), cex = 1,
         pch=c(-1), lwd = c(2), lty=c(5),
         bty = "n", border = FALSE)
}


