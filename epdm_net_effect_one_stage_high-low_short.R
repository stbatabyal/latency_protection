# R codes for Net Effect of the Treatment Active TB Individuals with one stage: Saikat Batabyal #

library(deSolve)
library(rootSolve)
library(bayestestR)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(scales)

################################################################################
### Model Frame (ODE)
################################################################################

mtb_epdm <- function(time,state,parameters) {
  with(as.list(c(state,parameters)),{
    dS.dt <- birth - mu*S - beta*S + gamma*L + delta_L*L + delta_L*I #Susceptible
    dI.dt <- eta*beta*S + beta*(1-epsilon)*L - theta*I - beta*I - lambda*I - mu*I - delta_L*I #Progressive_Stage_1
    dL.dt <- (1-eta)*beta*S + lambda*I - beta*(1-epsilon)*L - gamma*L - mu*L - delta_L*L #Non-progressive_Stage_1
    dA.dt <- theta*I + beta*I - mu*A - delta_R*A #Final Progressive Stage
    return(list(c(dS.dt, dI.dt, dL.dt, dA.dt)))
  })
}


################################################################################
### Figure 5D. Fully Protective(EPSILON=1) with LOW Infection and Recovery rates
################################################################################

################################################################################
### Prior to the Treatment 
################################################################################

### STEADY STATES ##############################################################

set.seed(601)
nboot <- 100
sd0 <- 0.2

vec_set_6A <- replicate(nboot, {
  list(birth = 0.013,
       mu = 0.013,
       beta = exp(rnorm(1, log(0.01), sd0)),
       gamma = exp(rnorm(1, log(0.05), sd0)),
       lambda = exp(rnorm(1, log(1), sd0)),
       eta = 0.05,
       theta = 2,
       epsilon = 1,
       delta_R = 0.3)
})

init1_set_6A  <- c(S = 0.7, I = 0.0006, L = 0.12, A = 0.004) #initial conditions for odes
time1_set_6A  <- seq(0,100,by=0.01) #time period

ss1_set_6A <- vector(ncol(vec_set_6A), mode="list")

for (k in 1:100){#range of values for beta
  ss1_set_6A[[k]] <- steady(y=init1_set_6A,times=0,func=mtb_epdm,
                            parms=c(vec_set_6A[,k],
                                    delta_L = 0))
}

names(ss1_set_6A) <- c()  ## to get "parameter" value incorporated in results
sol_ss1_set_6A <- dplyr::bind_rows(lapply(ss1_set_6A,as.data.frame),.id="")

### SOLUTION ###################################################################

res1_set_6A <- vector(ncol(vec_set_6A),mode="list")

for (k in 1:100){ #range of values for beta
  res1_set_6A[[k]] <- ode(y=ss1_set_6A[[k]][["y"]],times=time1_set_6A,func=mtb_epdm,
                          parms=c(vec_set_6A[,k],
                                  delta_L = 0))
}

names(res1_set_6A) <- c()  ## to get "parameter" value incorporated in results
sol1_set_6A <- dplyr::bind_rows(lapply(res1_set_6A,as.data.frame),.id="")

### Vector Array - Active Cases ################################################

mylist1_set_6A <- list() #create an empty list

for (k in 1:100) {
  vec1_set_6A <- res1_set_6A[[k]][,5] #pre-allocate a numeric vector
  mylist1_set_6A[[k]] <- vec1_set_6A #put all vectors in the list
}
df1_set_6A <- do.call("rbind", mylist1_set_6A) #combine all vectors into a matrix

res1_active_set_6A <- as.data.frame(t(df1_set_6A))

### Normalised Active Cases ####################################################

mylist2_set_6A <- list() #create an empty list

for (k in 1:100) {
  vec2_set_6A <- res1_set_6A[[k]][,5]/ss1_set_6A[[k]][["y"]][["A"]] #preallocate a numeric vector
  mylist2_set_6A[[k]] <- vec2_set_6A #put all vectors in the list
}
df2_set_6A <- do.call("rbind",mylist2_set_6A) #combine all vectors into a matrix

norm1_active_set_6A <- as.data.frame(t(df2_set_6A))

################################################################################
### During Treatment
################################################################################

### STEADY STATES ##############################################################

time2_set_6A <- seq(0,10,0.01) #data taken from years 1-10

ss2_set_6A <- vector(ncol(vec_set_6A),mode="list")

for (k in 1:100){ #range of values for "parameter"
  ss2_set_6A[[k]] <- steady(y=ss1_set_6A[[k]][["y"]],times=c(0,10),func=mtb_epdm,
                            parms=c(vec_set_6A[,k],
                                    delta_L = 1.5), method = "runsteady")
}

names(ss2_set_6A) <- c()  ## to get "parameter" value incorporated in results
sol_ss2_set_6A <- dplyr::bind_rows(lapply(ss2_set_6A,as.data.frame),.id="")

### SOLUTION ###################################################################

res2_set_6A <- vector(ncol(vec_set_6A),mode="list")

for (k in 1:100){ #range of values for "parameter"
  res2_set_6A[[k]] <- ode(y=ss1_set_6A[[k]][["y"]],times=time2_set_6A,func=mtb_epdm,
                          parms=c(vec_set_6A[,k],
                                  delta_L = 1.5))
}

names(res2_set_6A) <- c()  ## to get "parameter" value incorporated in results
sol2_set_6A <- dplyr::bind_rows(lapply(res2_set_6A,as.data.frame),.id="")

### Vector Array - Active Cases ################################################

mylist3_set_6A <- list() #create an empty list

for (k in 1:100) {
  vec3_set_6A <- res2_set_6A[[k]][,5] #pre-allocate a numeric vector
  mylist3_set_6A[[k]] <- vec3_set_6A #put all vectors in the list
}
df3_set_6A <- do.call("rbind",mylist3_set_6A) #combine all vectors into a matrix

res2_active_set_6A <- as.data.frame(t(df3_set_6A))

### Normalised Active Cases ####################################################

mylist4_set_6A <- list() #create an empty list

for (k in 1:100) {
  vec4_set_6A <- res2_set_6A[[k]][,5]/ss1_set_6A[[k]][["y"]][["A"]] #pre-allocate a numeric vector
  mylist4_set_6A[[k]] <- vec4_set_6A #put all vectors in the list
}
df4_set_6A <- do.call("rbind",mylist4_set_6A) #combine all vectors into a matrix

norm2_active_set_6A <- as.data.frame(t(df4_set_6A))

################################################################################
### After stopping the treatment
################################################################################

### SOLUTION ###################################################################

time3_set_6A <- seq(10,100,0.01) #data taken from years 10-100

res3_set_6A <- vector(ncol(vec_set_6A),mode="list")

for (k in 1:100){ #range of values for "parameter"
  res3_set_6A[[k]] <- ode(y=ss2_set_6A[[k]][["y"]],times=time3_set_6A,func=mtb_epdm,
                          parms=c(vec_set_6A[,k],
                                  delta_L = 0))
}

names(res3_set_6A) <- c()  ## to get "parameter" value incorporated in results
sol3_set_6A <- dplyr::bind_rows(lapply(res3_set_6A,as.data.frame),.id="")

### Vector Array - active cases w.r.t "parameter" ##############################

mylist5_set_6A <- list() #create an empty list

for (i in 1:100) {
  vec5_set_6A <- res3_set_6A[[k]][,5] #pre-allocate a numeric vector
  mylist5_set_6A[[k]] <- vec5_set_6A #put all vectors in the list
}
df5_set_6A <- do.call("rbind",mylist5_set_6A) #combine all vectors into a matrix

res3_active_set_6A <- as.data.frame(t(df5_set_6A))

### Normalised Active Cases ####################################################

mylist6_set_6A <- list() #create an empty list

for (k in 1:100) {
  vec6_set_6A <- res3_set_6A[[k]][,5]/ss1_set_6A[[k]][["y"]][["A"]] #pre-allocate a numeric vector
  mylist6_set_6A[[k]] <- vec6_set_6A #put all vectors in the list
}
df6_set_6A <- do.call("rbind",mylist6_set_6A) #combine all vectors into a matrix

norm3_active_set_6A <- as.data.frame(t(df6_set_6A))

################################################################################
### Calculate: Area Under the Curve (AUC) & Area Over the Curve (AOC) ##########
################################################################################

### DECREASE: Part-1, AUC (during the treatment) ###

auc1_list_set_6A <- list() #create an empty list

for (k in 1:100) {
  auc1_set_6A <- area_under_curve(time2_set_6A,norm2_active_set_6A[,k], method = "spline") #pre-allocate a numeric vector
  auc1_list_set_6A[[k]] <- auc1_set_6A #put all vectors in the list
}
auc_1_set_6A <- do.call("rbind",auc1_list_set_6A) #combine all vectors into a matrix

dat_auc_1_set_6A <- as.data.frame(auc_1_set_6A)

### Finding AOC (benefited blue area during treatment) #########################

aoc1_list_set_6A <- list() #create an empty list

for (k in 1:100) {
  aoc1_set_6A <- 10 - dat_auc_1_set_6A[k,1]  #pre-allocate a numeric vector
  aoc1_list_set_6A[[k]] <- aoc1_set_6A #put all vectors in the list
}
aoc_1_set_6A <- do.call("rbind",aoc1_list_set_6A) #combine all vectors into a matrix

dat_aoc_1_set_6A <- as.data.frame(aoc_1_set_6A)

### DECREASE: Part-2, Finding AUC ###

norm_mstops_set_6A <-cbind.data.frame(norm3_active_set_6A,time3_set_6A)
attach(norm_mstops_set_6A)

### INCREASE: Part-1, Finding AOC (due to stopping the treatment) ###

auc2_list_set_6A <- list() #create an empty list

for (k in 1:100) {
  condition_1_set_6A <- norm3_active_set_6A[,k] < 1 
  norm_mstops_condition_1_set_6A <- norm_mstops_set_6A[1:(min(which(condition_1_set_6A == FALSE))-1),]
  auc2_set_6A <- area_under_curve(norm_mstops_condition_1_set_6A$time3_set_6A,
                                  norm_mstops_condition_1_set_6A[,k], method = "spline") #pre-allocate a numeric vector
  auc2_list_set_6A[[k]] <- auc2_set_6A #put all vectors in the list
}
auc_2_set_6A <- do.call("rbind",auc2_list_set_6A) #combine all vectors into a matrix

dat_auc_2_set_6A <- as.data.frame(auc_2_set_6A)

### TIME with respect to each AUC ##############################################

time_stop_list1_set_6A <- list() #create an empty list

for (k in 1:100) {
  condition_1_set_6A <- norm3_active_set_6A[,k] < 1 
  norm_mstops_condition_1_set_6A <- norm_mstops_set_6A[1:(min(which(condition_1_set_6A == FALSE))-1),]
  time_stop_1_set_6A <- max(norm_mstops_condition_1_set_6A$time3_set_6A)
  time_stop_list1_set_6A[[k]] <- time_stop_1_set_6A #put all vectors in the list
}
time_stop_T1_set_6A <- do.call("rbind",time_stop_list1_set_6A) #combine all vectors into a matrix

dat_time_stop_T1_set_6A <- as.data.frame(time_stop_T1_set_6A)

### Finding AOC (benefited blue area after stopping treatment) #################

aoc2_list_set_6A <- list() #create an empty list

for (k in 1:100) {
  aoc2_set_6A <- (dat_time_stop_T1_set_6A[k,1]-10) - dat_auc_2_set_6A[k,1]  #pre-allocate a numeric vector
  aoc2_list_set_6A[[k]] <- aoc2_set_6A #put all vectors in the list
}
aoc_2_set_6A <- do.call("rbind",aoc2_list_set_6A) #combine all vectors into a matrix

dat_aoc_2_set_6A <- as.data.frame(aoc_2_set_6A)

### INCREASE: Part-2, Finding AOC (for the peak, red curve) ####################

auc3_list_set_6A <- list() #create an empty list

for (k in 1:100) {
  condition_2_set_6A <- norm3_active_set_6A[,k] > 1
  norm_mstops_condition_2_set_6A <- norm_mstops_set_6A[(min(which(condition_2_set_6A == TRUE))):nrow(norm_mstops_set_6A),]
  
  auc3_set_6A <- area_under_curve(norm_mstops_condition_2_set_6A$time3_set_6A,
                                  norm_mstops_condition_2_set_6A[,k], method = "spline") #pre-allocate a numeric vector
  auc3_list_set_6A[[k]] <- auc3_set_6A #put all vectors in the list
}
auc_3_set_6A <- do.call("rbind",auc3_list_set_6A,) #combine all vectors into a matrix

dat_auc_3_set_6A <- as.data.frame(auc_3_set_6A)

### TIME with respect to each AUC ##############################################

time_stop_list2_set_6A <- list() #create an empty list

for (k in 1:100) {
  condition_2_set_6A <- norm3_active_set_6A[,k] > 1
  norm_mstops_condition_2_set_6A <- norm_mstops_set_6A[(min(which(condition_2_set_6A == TRUE))):nrow(norm_mstops_set_6A),]
  time_stop_2_set_6A <- min(norm_mstops_condition_2_set_6A$time3_set_6A)
  time_stop_list2_set_6A[[k]] <- time_stop_2_set_6A #put all vectors in the list
}

time_stop_T2_set_6A <- do.call("rbind",time_stop_list2_set_6A) #combine all vectors into a matrix

dat_time_stop_T2_set_6A <- as.data.frame(time_stop_T2_set_6A)

### Finding Area Under the Curve ###############################################

aoc3_list_set_6A <- list() #create an empty list

for (k in 1:100) {
  aoc3_set_6A <- dat_auc_3_set_6A[k,1] - (100-dat_time_stop_T2_set_6A[k,1])   #pre-allocate a numeric vector
  aoc3_list_set_6A[[k]] <- aoc3_set_6A #put all vectors in the list
}
aoc_3_set_6A <- do.call("rbind",aoc3_list_set_6A) #combine all vectors into a matrix

dat_aoc_3_set_6A <- as.data.frame(aoc_3_set_6A)

################################################################################
### Calculate Total Area
################################################################################

auc_decrease_set_6A <- -(dat_aoc_1_set_6A + dat_aoc_2_set_6A)
auc_increase_set_6A <- dat_aoc_3_set_6A

per_auc_decrease_set_6A <- auc_decrease_set_6A*10
per_auc_increase_set_6A <- auc_increase_set_6A*10

total_auc_set_6A <- auc_decrease_set_6A + auc_increase_set_6A
total_auc_set_6A
per_total_auc_set_6A <- total_auc_set_6A*10

################################################################################
### Re-scaling with percentage of Total Area
################################################################################

x1_set_6A <- per_total_auc_set_6A$V1
x2_set_6A <- per_auc_decrease_set_6A$V1
x3_set_6A <- per_auc_increase_set_6A$V1

################################################################################
##### SCALING INTO LOG_SCALE in GG-PLOT ########################################
################################################################################

### NET ###
DF1_net_stage_one_low_eps_1 <- data.frame(X=x1_set_6A)
Den_net_stage_one_low_eps_1 <- density(DF1_net_stage_one_low_eps_1$X)
DF2_net_stage_one_low_eps_1 <- data.frame(x = Den_net_stage_one_low_eps_1$x, 
                                           y = Den_net_stage_one_low_eps_1$y)
### AOC ###
DF1_aoc_stage_one_low_eps_1 <- data.frame(X=x2_set_6A)
Den_aoc_stage_one_low_eps_1 <- density(DF1_aoc_stage_one_low_eps_1$X)
DF2_aoc_stage_one_low_eps_1 <- data.frame(x = Den_aoc_stage_one_low_eps_1$x, 
                                           y = Den_aoc_stage_one_low_eps_1$y)

### AUC ###
DF1_auc_stage_one_low_eps_1 <- data.frame(X=x3_set_6A)
Den_auc_stage_one_low_eps_1 <- density(DF1_auc_stage_one_low_eps_1$X)
DF2_auc_stage_one_low_eps_1 <- data.frame(x = Den_auc_stage_one_low_eps_1$x, 
                                           y = Den_auc_stage_one_low_eps_1$y)

# plotmath expression wrapped in expression()

P_low_eps_1 <- "Mtb~infection~is~fully~protective~'('~epsilon==1~','~beta*T==0.01~','~gamma==0.05~')'"
L1 <- "Reduction~during~MTBI~treatment~program~(AOC)"
L2 <- "Increase~after~program~stop~(AUC)"
L3 <- "Net~effect~','~Delta==AUC-AOC"

# Plot  
plot4 <-ggplot()+
  geom_line(aes(x = x, y = y), data = DF2_net_stage_one_low_eps_1, color='black', size=1)+
  geom_line(aes(x = x, y = y), data = DF2_aoc_stage_one_low_eps_1, color='blue', linetype="dashed", size=1)+
  geom_line(aes(x = x, y = y), data = DF2_auc_stage_one_low_eps_1, color='red', linetype="dotted", size=1)+
  geom_segment(aes(x = 0, y = 10^-3, xend = 0, yend = 10^0), color='grey', linetype="dashed", size=1)+
  # set axes titles and labels
  scale_x_continuous(name = "Percentage change in TB cases due to program") +
  scale_y_log10(name = "Fraction", breaks = trans_breaks("log10", function(x) 10^x, 4),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(10^-3, 10^0), expand = c(0, 0))+
  coord_cartesian(xlim = c(-200, 50), ylim = c(10^-3, 10^0), expand = FALSE)+ # Set axis limits
  # Legends
  geom_segment(aes(x = -195, y = 0.5, xend = -175, yend = 0.5), color='blue', linetype="dashed", size=1)+
  geom_segment(aes(x = -195, y = 0.32, xend = -175, yend = 0.32), color='red', linetype="dotted", size=1)+ 
  geom_segment(aes(x = -195, y = 0.2, xend = -175, yend = 0.2), color='black', linetype="solid", size=1)+  
  annotate("text", -135, 0.8, label = P_low_eps_1, size = 9.4/.pt, parse = TRUE)+
  annotate("text", -120, 0.5, label = L1, size = 9.4/.pt, parse = TRUE)+
  annotate("text", -134, 0.32, label = L2, size = 9.4/.pt, parse = TRUE)+
  annotate("text", -140, 0.2, label = L3, size = 9.4/.pt, parse = TRUE)+
  # Theme adjustments:
  theme_classic() + # eliminate default background 
  theme(panel.grid.major = element_blank(), # eliminate major grids
        panel.grid.minor = element_blank(), # eliminate minor grids
        # set font family for all text within the plot ("serif" should work as "Times New Roman")
        # note that this can be overridden with other adjustment functions below
        text = element_text(family="Helvetica"),
        # adjust X-axis title
        axis.title.x = element_text(size = 13),
        # adjust X-axis labels; also adjust their position using margin (acts like a bounding box)
        # using margin was needed because of the inwards placement of ticks
        axis.text.x = element_text(size = 11, face = "bold", margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),
        # adjust Y-axis title
        axis.title.y = element_text(size = 14),
        # adjust Y-axis labels
        axis.text.y = element_text(size = 11, face = "bold", margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
        # length of tick marks - negative sign places ticks inwards
        axis.ticks.length = unit(-1.4, "mm"),
        # width of tick marks in mm
        axis.ticks = element_line(size = 0.5)
  ) 

plot4

# Repeat the same process for epsilon=0 for low transmission 
# Repeat the same process for epsilon=1 for high transmission
# Repeat the same process for epsilon=0 for high transmission


