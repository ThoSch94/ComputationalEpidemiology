library(deSolve)
library(reshape2)
library(ggplot2)
library(xkcd)
library(ggpubr)

sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- - beta * S * I - beta_treatment * S * T
    dI <- beta * S * I - gamma * I  - I * treatment_rate + beta_treatment * S * T
    dR <-                 gamma * I  + T * gamma_treated
    dT <- I * treatment_rate - T * gamma_treated
    
    return(list(c(dS, dI, dR, dT)))
  })
}


calculate_and_plot = function(sir, parameters, title, init){
  
  ## Time frame
  times      <- seq(0, 100, by = 0.01)
  ## Solve using ode (General Solver for Ordinary Differential Equations)
  out <- ode(y = init, times = times, func = sir, parms = parameters)
  
  ## change to data frame
  out <- as.data.frame(out)
  
  
  dat = melt(out,id.vars = c("time"))
  
  
  names(dat) = c("time","Group","value")
  
  
  plot = 
    ggplot(data = dat,aes(x = time,y = value/(max(value)))) +
    geom_line(aes(col = Group,group = Group),size=1.2) +
    theme_minimal() +
    labs(title=paste(paste("SIR Model", title)))+
    xlab ("Time")  +
    ylab("Proportion of Population")+
    theme_classic() + 
    theme(text = element_text(size = 20)) +
    ylim(0,1)+ 
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "black")) 
  return(plot)
}


### Set parameters
## Proportion in each compartment
init<- c(S = 1-1e-3, I = 1e-3, R = 0.0, T = 0)
parameter_matrix= matrix(nrow = 3, ncol = 5)
## beta: infection parameter; gamma: recovery parameter
colnames(parameter_matrix) = c("beta", "gamma", "treatment_rate", "gamma_treated", "beta_treatment")
parameter_matrix[1,] = c(0.15, 0.003, 0.2, 0.003, 0.15)
parameter_matrix[2,] = c(0.15, 0.003, 0.2, 0.03, 0.15)
parameter_matrix[3,] = c(0.15, 0.003, 0.2, 0.1, 0.15)

plot1 = calculate_and_plot(sir, parameter_matrix[1,], "gamme_t = 0.1, beta_t = 2", init)
plot2 = calculate_and_plot(sir, parameter_matrix[2,], "gamme_t = 0.3, beta_t = 2", init)
plot3 = calculate_and_plot(sir, parameter_matrix[3,], "gamme_t = 0.6, beta_t = 2", init)
arranged_plot = ggarrange(plot1, plot2, plot3, ncol = 3, labels="AUTO")
ggsave("part1_plot.png", arranged_plot, width = 30, height = 6, units = "in", dpi = 300)


# part 2
init<- c(S = 1-1e-6, I = 1e-6, R = 0.0, T = 0)
parameter_matrix= matrix(nrow = 9, ncol = 5)
## beta: infection parameter; gamma: recovery parameter
colnames(parameter_matrix) = c("beta", "gamma", "treatment_rate", "gamma_treated", "beta_treatment")
parameter_matrix[1,] = c(2, 0.1, 0.2, 0.1, 2)
parameter_matrix[2,] = c(2, 0.1, 0.2, 0.1, 1.3)
parameter_matrix[3,] = c(2, 0.1, 0.2, 0.1, 0.6)
parameter_matrix[4,] = c(2, 0.1, 0.2, 0.2, 2)
parameter_matrix[5,] = c(2, 0.1, 0.2, 0.3, 1.3)
parameter_matrix[6,] = c(2, 0.1, 0.2, 0.2, 0.6)
parameter_matrix[7,] = c(2, 0.1, 0.2, 0.6, 2)
parameter_matrix[8,] = c(2, 0.1, 0.2, 0.6, 1.3)
parameter_matrix[9,] = c(2, 0.1, 0.2, 0.6, 0.6)
plot1 = calculate_and_plot(sir, parameter_matrix[1,], "gamme_t = 0.1, beta_t = 2", init)
plot2 = calculate_and_plot(sir, parameter_matrix[2,], "gamme_t = 0.1, beta_t = 1.3", init)
plot3 = calculate_and_plot(sir, parameter_matrix[3,], "gamme_t = 0.1, beta_t = 0.6", init)
plot4 = calculate_and_plot(sir, parameter_matrix[4,], "gamme_t = 0.3, beta_t = 2", init)
plot5 = calculate_and_plot(sir, parameter_matrix[5,], "gamme_t = 0.3, beta_t = 1.3", init)
plot6 = calculate_and_plot(sir, parameter_matrix[6,], "gamme_t = 0.3, beta_t = 0.6", init)
plot7 = calculate_and_plot(sir, parameter_matrix[7,], "gamme_t = 0.6, beta_t = 2", init)
plot8 = calculate_and_plot(sir, parameter_matrix[8,], "gamme_t = 0.6, beta_t = 1.3", init)
plot9 = calculate_and_plot(sir, parameter_matrix[9,], "gamme_t = 0.6, beta_t = 0.6", init)
arranged_plot = ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol = 3, nrow = 3,labels="AUTO")
ggsave("part2_plot.png", arranged_plot, width = 30, height = 18, units = "in", dpi = 300)

