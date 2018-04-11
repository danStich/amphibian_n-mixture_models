# Front-end ---------------------------------------------------------------
  library(R2jags)
  library(plyr)
  library(lubridate)
  library(plotrix)
  library(MASS)

# MCMC Settings
  ni = 250000
  nb = 50000
  nt = 50
  nc = 3
  
# ------------------------------------------------------------------------------
# 2016 WF data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vPools_2016.csv')

# Look at the first few lines of data
  head(herps)

# Create a new variable that indicates whether the pond was in the lower or the
# upper site
  for(i in 1:nrow(herps)){
    if(herps$pond[i] <= 8) {herps$site[i]=0
    } else {
      herps$site[i]=1
    }
  }

# Fix the date formatting so we can use it in R if needed
  herps$date = as.POSIXct(as.character(herps$date), format='%m/%d/%Y')

# Define day of year based on the new date column
  herps$day = yday(herps$date)

# Define week number based on the date column
  herps$week = as.numeric(as.factor(week(herps$date)))
  #herps = herps[herps$week==2 | herps$week==3, ]
  #herps$week = as.numeric(as.factor((herps$week)))

# Form initial values for wf abundance
  Nst = ddply(herps, c("pond", "week"), summarize, count=sum(wf))#[,2]+1

# Make a pond x week matrix
  Nst2 =  xtabs(count~pond+week, data=Nst)
  dimnames(Nst2)=NULL

# Get weekdays from data
  herps$weekday = as.numeric(as.factor(weekdays(herps$date)))

# Make a matrix of counts by weekday and week
  Nst3 = xtabs(wf~pond + weekday + week, data = herps)
  dimnames(Nst3) = NULL

# Get pond as a numerically represented factor
  herps$pond = as.numeric(as.factor(herps$pond))

# 2016 WF mixture model specification ----------------------------
# Define the base model
  modelString = "
    model{
      # Likelihood
        for(i in 1:nponds){
          for(t in 1:nweeks){
            N[i,t] ~ dpois(lambda[i, t])
          }
        }

      # Random effect of site and day on abundance
        for(z in 1:nponds){
          for(i in 1:7){
            for(t in 1:nweeks){
              C[z,i,t] ~ dbin(p[pond[z], week[t]], N[pond[z], week[t]])
            }
          }
        }

      # Priors
        for(i in 1:nponds){
          #meanP[i] ~ dbeta(1, 1)

          for(t in 1:nweeks){
            log(lambda[i,t]) <- alpha[i, t]
            alpha[i,t] ~ dnorm(0, .1)T(-10, 10)
            #p[i,t] <- meanP[i]
            p[i,t] ~ dbeta(1, 5)
          }
        }
    }
  "

# Write the model to a file
  writeLines(modelString, "wfmodel.txt")

# 2016 WF model settings and calibration ------------------------------------------
# Bundle the data for JAGS
  herps.data = list(
    nponds = length(unique(herps$pond)),
    nweeks = length(unique(herps$week)),
    C = Nst3,
    pond = unique(herps$pond),
    week = sort(unique(herps$week))
  )

# Provide initial values for the Gibbs sampler
  inits = function(){
    list(
      N = Nst2,
      #meanP = rbeta(15, 1, 1),
      alpha = structure(rnorm(60, 0, 1), .Dim=c(15,4)),
      p = structure(rbeta(60, 1, 5), .Dim=c(15,4))
    )
  }

# Parameters to monitor during estimation
  params = c(
    #'meanP',
    'N',
    'p'
  )

# # MCMC Settings
#   ni = 25000
#   nb = 5000
#   nt = 50
#   nc = 3

# Run the model
  wf.model16 = jags(data = herps.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "wfmodel.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(wf.model16, digits = 2)

# 2016 WF Results --------------------------------------------------------------
# Get the results
  wfres16 = data.frame(wf.model16$BUGSoutput$sims.list)
  
# Get abundance estimates for each week
  wk1 = wfres16[ , 1:15]
  wk2 = wfres16[ , 16:30]
  wk3 = wfres16[ , 31:45]
  wk4 = wfres16[ , 46:60]

# Now, get the week for max abundance of species in each of the 15 pools and put
# the posterior for that week into a list object.
  check = c()
  for(i in 1:15){
    check[i] =
      which(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))==
            max(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))))
  }

  big = list(wk1, wk2, wk3, wk4)

  wfres16 = vector(mode='list', length=15)
  for(i in 1:15){
    wfres16[[i]] = big[[check[i]]][,i]
  }
  wfres16 = do.call(cbind, wfres16)
  
# Make boxplot of abundance by pool in upper and lower sites
  par(mar=c(5,6,1,1))
  boxplot(wfres16, ylim=c(0, 700), outline = FALSE,
          names=c(as.character(seq(1,8,1)),
            c("9", "10", "11", "18", "19", "20", "21")
          ),
          cex.axis = 1.15,
          yaxt='n',
          col='gray'
  )
  axis(2, las=2, cex.axis=1.15)
  mtext('Estimated abundance (N)', side=2, cex=1.25, line=5)
  mtext('Vernal pool', side=1, cex=1.25, line = 4)
  abline(v=8.5, col='black', lwd=2, lty=2)
  text(x = 1, y=700, "Lower pools")
  text(x = 14.5, y=700, "Upper pools")

# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      wflowers16 = apply(wfres16[ , 1:8], 1, mean)
      mean(wflowers16)
      quantile(wflowers16, probs = c(0.05, 0.95))

    # Mean abundance in lower pools
      wfuppers16 = apply(wfres16[ , 9:15], 1, mean)
      mean(wfuppers16)
      quantile(wfuppers16, probs = c(0.05, 0.95))  
 
  # Make the plots
    differences = wfuppers16-wflowers16
    par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,4,1,1))
    hist(differences, col='gray87', border='gray60', yaxt='n', xaxt='n',
         main='', xlab='',cex.lab=1.25)
    axis(1, pos=0, cex=1.15)
    axis(2, pos=-200, las=2, cex=1.15)
    mtext(expression(paste("Mean difference (N"["upper"], "- N"["lower"],")")),
          side = 1, cex=1.25, line = 4)
    abline(v=mean(differences), col='black', lwd=3, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='black', lwd=2, lty=2)
    text(x=600, y=700, expression(paste(mu, "= " , "203")), adj=0)
    text(x=600, y=650, expression(paste("95% CRI = 100-567")), adj=0)


# 2016 Jeff data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vPools_2016.csv')

# Look at the first few lines of data
  head(herps)

# Create a new variable that indicates whether the pond was in the lower or the
# upper site
  for(i in 1:nrow(herps)){
    if(herps$pond[i] <= 8) {herps$site[i]=0
    } else {
      herps$site[i]=1
    }
  }

# Fix the date formatting so we can use it in R if needed
  herps$date = as.POSIXct(as.character(herps$date), format='%m/%d/%Y')

# Define day of year based on the new date column
  herps$day = yday(herps$date)

# Define week number based on the date column
  herps$week = as.numeric(as.factor(week(herps$date)))
  #herps = herps[herps$week==2 | herps$week==3, ]
  #herps$week = as.numeric(as.factor((herps$week)))

# Form initial values for wf abundance
  Nst = ddply(herps, c("pond", "week"), summarize, count=sum(jeff))#[,2]+1

# Make a pond x week matrix
  Nst2 =  xtabs(count~pond+week, data=Nst)
  dimnames(Nst2)=NULL

# Get weekdays from data
  herps$weekday = as.numeric(as.factor(weekdays(herps$date)))

# Make a matrix of counts by weekday and week
  Nst3 = xtabs(jeff~pond + weekday + week, data = herps)
  dimnames(Nst3) = NULL

# Get pond as a numerically represented factor
  herps$pond = as.numeric(as.factor(herps$pond))

# 2016 Jeff model specification ----------------------------
# Define the base model
  modelString = "
    model{
      # Likelihood
        for(i in 1:nponds){
          for(t in 1:nweeks){
            N[i,t] ~ dpois(lambda[i, t])T(0,2000)
          }
        }

      # Random effect of site and day on abundance
        for(z in 1:nponds){
          for(i in 1:7){
            for(t in 1:nweeks){
              C[z,i,t] ~ dbin(p[pond[z], week[t]], N[pond[z], week[t]])
            }
          }
        }

      # Priors
        for(i in 1:nponds){
          #meanP[i] ~ dbeta(1, 1)

          for(t in 1:nweeks){
            log(lambda[i,t]) <- alpha[i, t]
            alpha[i,t] ~ dnorm(0, .1)T(-10, 10)
            #p[i,t] <- meanP[i]
            p[i,t] ~ dbeta(1, 5)
          }
        }
    }
  "

# Write the model to a file
  writeLines(modelString, "jeffmodel.txt")

# 2016 Jeff model settings and calibration ------------------------------------------
# Bundle the data for JAGS
  herps.data = list(
    nponds = length(unique(herps$pond)),
    nweeks = length(unique(herps$week)),
    C = Nst3,
    pond = unique(herps$pond),
    week = sort(unique(herps$week))
  )

# Provide initial values for the Gibbs sampler
  inits = function(){
    list(
      N = Nst2,
      alpha = structure(rnorm(60, 0, 1), .Dim=c(15,4)),
      p = structure(rbeta(60, 1, 5), .Dim=c(15,4))
    )
  }

# Parameters to monitor during estimation
  params = c(
    'N',
    'p'
  )

# # MCMC Settings
#   ni = 25000
#   nb = 5000
#   nt = 50
#   nc = 3

# Run the model
  jeff.model16 = jags(data = herps.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "jeffmodel.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(jeff.model16, digits = 2)

# 2016 Jeff results ------------------------------------------------------------
# Collect results from model object
  jres = data.frame(jeff.model16$BUGSoutput$sims.list)
# Get abundance estimates for each week
  wk1 = jres[ , 1:15]
  wk2 = jres[ ,16:30]
  wk3 = jres[ ,31:45]
  wk4 = jres[ ,46:60]

# Now, get the week for max abundance of species in each of the 15 pools and put
# the posterior for that week into a list object.
  check = c()
  for(i in 1:15){
    check[i] =
      which(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))==
            max(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))))
  }

  big = list(wk1, wk2, wk3, wk4)

  jeffres16 = vector(mode='list', length=15)
  for(i in 1:15){
    jeffres16[[i]] = big[[check[i]]][,i]
  }
  jeffres16 = do.call(cbind, jeffres16)
  
# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      jefflowers16 = apply(jeffres16[ , 1:8], 1, mean)
      mean(jefflowers16)
      quantile(jefflowers16, probs = c(0.05, 0.95))
    # Mean abundance in lower pools
      jeffuppers16 = apply(jeffres16[ , 9:15], 1, mean)
      mean(jeffuppers16)
      quantile(jeffuppers16, probs = c(0.05, 0.95))  

# Test for differences in abundance based on presence/absence of streams
  streams = apply(jeffres16[ , c(1,2,3,6,7)], 1, mean)
  nostream = apply(jeffres16[ , c(4,5,8)], 1, mean)
  quantile(streams-nostream, probs=c(0.025, 0.975))
# Test for differences in abundance based on hydroperiod
  full = apply(jeffres16[ , c(1,3,6,7)], 1, mean)
  empty = apply(jeffres16[ , c(2,4,5,8)], 1, mean)
  quantile(full-empty, probs=c(0.025, 0.975))        
    
# 2016 Spot data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vPools_2016.csv')

# Look at the first few lines of data
  head(herps)

# Create a new variable that indicates whether the pond was in the lower or the
# upper site
  for(i in 1:nrow(herps)){
    if(herps$pond[i] <= 8) {herps$site[i]=0
    } else {
      herps$site[i]=1
    }
  }

# Fix the date formatting so we can use it in R if needed
  herps$date = as.POSIXct(as.character(herps$date), format='%m/%d/%Y')

# Define day of year based on the new date column
  herps$day = yday(herps$date)

# Define week number based on the date column
  herps$week = as.numeric(as.factor(week(herps$date)))
  #herps = herps[herps$week==2 | herps$week==3, ]
  #herps$week = as.numeric(as.factor((herps$week)))

# Form initial values for wf abundance
  Nst = ddply(herps, c("pond", "week"), summarize, count=sum(spot))#[,2]+1

# Make a pond x week matrix
  Nst2 =  xtabs(count~pond+week, data=Nst)
  dimnames(Nst2)=NULL

# Get weekdays from data
  herps$weekday = as.numeric(as.factor(weekdays(herps$date)))

# Make a matrix of counts by weekday and week
  Nst3 = xtabs(spot~pond + weekday + week, data = herps)
  dimnames(Nst3) = NULL

# Get pond as a numerically represented factor
  herps$pond = as.numeric(as.factor(herps$pond))

# 2016 Spot model specification ----------------------------
# Define the base model
  modelString = "
    model{
      # Likelihood
        for(i in 1:nponds){
          for(t in 1:nweeks){
            N[i,t] ~ dpois(lambda[i, t])
          }
        }

      # Random effect of site and day on abundance
        for(z in 1:nponds){
          for(i in 1:7){
            for(t in 1:nweeks){
              C[z,i,t] ~ dbin(p[pond[z], week[t]], N[pond[z], week[t]])
            }
          }
        }

      # Priors
        for(i in 1:nponds){
          #meanP[i] ~ dbeta(1, 1)

          for(t in 1:nweeks){
            log(lambda[i,t]) <- alpha[i, t]
            alpha[i,t] ~ dnorm(0, .1)T(-10, 10)
            #p[i,t] <- meanP[i]
            p[i,t] ~ dbeta(1, 5)
          }
        }
    }
  "

# Write the model to a file
  writeLines(modelString, "spotmodel.txt")

# 2016 Spot model settings and calibration ------------------------------------------
# Bundle the data for JAGS
  herps.data = list(
    nponds = length(unique(herps$pond)),
    nweeks = length(unique(herps$week)),
    C = Nst3,
    pond = unique(herps$pond),
    week = sort(unique(herps$week))
  )

# Provide initial values for the Gibbs sampler
  inits = function(){
    list(
      N = Nst2,
      #meanP = rbeta(15, 1, 1),
      alpha = structure(rnorm(60, 0, 1), .Dim=c(15,4)),
      p = structure(rbeta(60, 1, 5), .Dim=c(15,4))
    )
  }

# Parameters to monitor during estimation
  params = c(
    #'meanP',
    'N',
    'p'
  )

# # MCMC Settings
#   ni = 25000
#   nb = 5000
#   nt = 50
#   nc = 3

# Run the model
  spot.model16 = jags(data = herps.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "spotmodel.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(spot.model16, digits = 2)


# 2016 Spot results ------------------------------------------------------------
# Collect results from model object
  sres = data.frame(spot.model16$BUGSoutput$sims.list)
# Get abundance estimates for each week
  wk1 = sres[ , 1:15]
  wk2 = sres[ ,16:30]
  wk3 = sres[ ,31:45]
  wk4 = sres[ ,46:60]

# Now, get the week for max abundance of species in each of the 15 pools and put
# the posterior for that week into a list object.
  check = c()
  for(i in 1:15){
    check[i] =
      which(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))==
            max(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))))
  }

  big = list(wk1, wk2, wk3, wk4)

  spotres16 = vector(mode='list', length=15)
  for(i in 1:15){
    spotres16[[i]] = big[[check[i]]][,i]
  }
  spotres16 = do.call(cbind, spotres16)

# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      spotlowers16 = apply(spotres16[ , 1:8], 1, mean)
      mean(spotlowers16)
      quantile(spotlowers16, probs = c(0.025, 0.975))
    # Mean abundance in lower pools
      spotuppers16 = apply(spotres16[ , 9:15], 1, mean)
      mean(spotuppers16)
      quantile(spotuppers16, probs = c(0.025, 0.975))  

# Test for differences in abundance based on presence/absence of streams
  streams = apply(spotres16[ , c(1,2,3,6,7)], 1, mean)
  nostream = apply(spotres16[ , c(4,5,8)], 1, mean)
  quantile(streams-nostream, probs=c(0.025, 0.975))

# Test for differences in abundance based on hydroperiod
  full = apply(spotres16[ , c(1,3,6,7)], 1, mean)
  empty = apply(spotres16[ , c(2,4,5,8)], 1, mean)
  quantile(full-empty, probs=c(0.025, 0.975))         
        
# 2016 ALL SPP FIGURES --------------------------------------------------------- 
# Box plots of abundance at each pool
  # WF
    par(mfrow=c(3,1), oma=c(1,2.5,1,1))
    boxplot(wfres16, ylim=c(0, 700), outline = FALSE,
            names=c(as.character(seq(1,8,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    abline(v=8.5, col='black', lwd=2, lty=2)
    text(x=11, y=650, "Upper pools", cex=1.25)
    text(x=4, y=650, "Lower pools", cex=1.25)  
    text(x=.5, y=650, '(a)', cex=1.5)
  # Jeffs
    boxplot(jeffres16, ylim=c(0, 125), outline = FALSE,
            names=c(as.character(seq(1,8,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    mtext('Estimated abundance (N)', side=2, cex=1.25, line=5)
    abline(v=8.5, col='black', lwd=2, lty=2)
    text(x=.5, y=110, '(b)', cex=1.5)
  # Spotted salamanders
    boxplot(spotres16, ylim=c(0, 500), outline = FALSE,
            names=c(as.character(seq(1,8,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    mtext('Vernal pool', side=1, cex=1.25, line = 4)
    abline(v=8.5, col='black', lwd=2, lty=2)
    text(x=.5, y=450, '(c)', cex=1.5)   
    
# Histograms of differences between upper and lower sites   
  par(mfrow=c(3,1), oma=c(1,1,1,1))
  # Woodfrogs
    differences = wfuppers16-wflowers16
    hist(differences, col='gray87', border='gray60', yaxt='n', xaxt='n',
         main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-100, 400), breaks=25)
    axis(1, pos=0, cex=1.15)
    abline(v=mean(differences), col='black', lwd=3, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='black', lwd=2, lty=2)
    text(x=300, y=800, expression(paste(mu, "= " , "77")), adj=0, cex=1.5)
    text(x=300, y=700, expression(paste("95% CRI = 21-234")), adj=0, cex=1.5)
    text(-100, 800, '(a)', cex=1.5)#; axis(2)
  # Jefferson salamanders
    differences = jeffuppers16-jefflowers16
    hist(differences, col='gray87', border='gray60', yaxt='n', xaxt='n',
         main='', xlab='', ylab = '', xlim=c(-100, 400), breaks=25)
    axis(1, pos=0, cex=1.15)
    mtext(side=2, 'Frequency', cex=1.25, line=3)
    abline(v=mean(differences), col='black', lwd=3, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='black', lwd=2, lty=2)
    text(x=300, y=1500, expression(paste(mu, "= " , "203")), adj=0, cex=1.5)
    text(x=300, y=1350, expression(paste("95% CRI = - 6-96")), adj=0, cex=1.5)    
    text(-100, 1500, '(b)', cex=1.5)#; axis(2)
  # Spotted salamander
    differences = spotuppers16-spotlowers16
    hist(differences, col='gray87', border='gray60', yaxt='n', xaxt='n',
         main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-100, 400), breaks=25)
    axis(1, pos=0, cex=1.15)
    mtext(expression(paste("Mean difference (N"["upper"], "- N"["lower"],")")),
          side = 1, cex=1.25, line = 4)
    abline(v=mean(differences), col='black', lwd=3, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='black', lwd=2, lty=2)
    text(x=300, y=700, expression(paste(mu, "= " , "76")), adj=0, cex=1.5)
    text(x=300, y=600, expression(paste("95% CRI = 9-224")), adj=0, cex=1.5)    
    text(-100, 700, '(c)', cex=1.5)#; axis(2)

# ------------------------------------------------------------------------------
# 2017 WF data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vPools_2017.csv')

# Look at the first few lines of data
  head(herps)

# Create a new variable that indicates whether the pond was in the lower or the
# upper site
  for(i in 1:nrow(herps)){
    if(herps$pond[i] <= 8) {herps$site[i]=0
    } else {
      herps$site[i]=1
    }
  }

# Fix the date formatting so we can use it in R if needed
  herps$date = as.POSIXct(as.character(herps$date), format='%m/%d/%Y')

# Define day of year based on the new date column
  herps$day = yday(herps$date)

# Define week number based on the date column
  herps$week = as.numeric(as.factor(week(herps$date)))
  #herps = herps[herps$week==2 | herps$week==3, ]
  #herps$week = as.numeric(as.factor((herps$week)))

# Form initial values for wf abundance
  Nst = ddply(herps, c("pond", "week"), summarize, count=sum(wf))#[,2]+1

# Make a pond x week matrix
  Nst2 =  xtabs(count~pond+week, data=Nst)
  dimnames(Nst2)=NULL

# Get weekdays from data
  herps$weekday = as.numeric(as.factor(weekdays(herps$date)))

# Make a matrix of counts by weekday and week
  Nst3 = xtabs(wf~pond + weekday + week, data = herps)
  dimnames(Nst3) = NULL
  #Nst3[Nst3==0] = NA

# Get pond as a numerically represented factor
  herps$pond = as.numeric(as.factor(herps$pond))

# 2017 WF mixture model specification ----------------------------
# Define the base model
  modelString = "
    model{
      # Likelihood
        for(i in 1:nponds){
          for(t in 1:nweeks){
            N[i,t] ~ dpois(lambda[i, t])
          }
        }

      # Random effect of site and day on abundance
        for(z in 1:nponds){
          for(i in 1:7){
            for(t in 1:nweeks){
              C[z,i,t] ~ dbin(p[pond[z], week[t]], N[pond[z], week[t]])
            }
          }
        }

      # Priors
        for(i in 1:nponds){
          #meanP[i] ~ dbeta(1, 1)

          for(t in 1:nweeks){
            log(lambda[i,t]) <- alpha[i, t]
            alpha[i,t] ~ dnorm(0, 0.1)T(-10, 10)
            #p[i,t] <- meanP[i]
            p[i,t] ~ dbeta(1, 1)
          }
        }
    }
  "

# Write the model to a file
  writeLines(modelString, "wfmodel17.txt")

# 2017 WF model settings and calibration ------------------------------------------
# Bundle the data for JAGS
  herps.data = list(
    nponds = length(unique(herps$pond)),
    nweeks = length(unique(herps$week)),
    C = Nst3,
    pond = unique(herps$pond),
    week = sort(unique(herps$week))
  )

# Provide initial values for the Gibbs sampler
  inits = function(){
    list(
      N = Nst2,
      #meanP = rbeta(15, 1, 1),
      alpha = structure(
              rnorm(length(unique(herps$week))*length(unique(herps$pond)), 0, 1),
              .Dim=c(length(unique(herps$pond)),length(unique(herps$week)))),
      p = structure(
          rbeta(length(unique(herps$week))*length(unique(herps$pond)), 1, 5),
         .Dim=c(length(unique(herps$pond)),length(unique(herps$week))))
    )
  }

# Parameters to monitor during estimation
  params = c(
    #'meanP',
    'N',
    'p'
  )

# # MCMC Settings
#   ni = 25000
#   nb = 5000
#   nt = 50
#   nc = 3

# Run the model
  wf.model17 = jags(data = herps.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "wfmodel17.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(wf.model17, digits = 2)

# 2017 Results WF --------------------------------------------------------------
# Get the results
  wfres17 = data.frame(wf.model17$BUGSoutput$sims.list)
  
# Get abundance estimates for each week
  wk1 = wfres17[ , 1:14]
  wk2 = wfres17[ , 15:28]
  wk3 = wfres17[ , 29:43]
  wk4 = wfres17[ , 43:56]
  wk5 = wfres17[ , 57:70]

# Now, get the week for max abundance of species in each of the 15 pools and put
# the posterior for that week into a list object.
  check = c()
  for(i in 1:14){
    check[i] =
      which(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]), mean(wk5[,i]))==
            max(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]), mean(wk5[,i]))))
  }

  big = list(wk1, wk2, wk3, wk4, wk5)

  wfres17 = vector(mode='list', length=14)
  for(i in 1:14){
    wfres17[[i]] = big[[check[i]]][,i]
  }
  wfres17 = do.call(cbind, wfres17)
  
# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    wflowers = apply(wfres17[ , 1:7], 1, mean)
    mean(wflowers)
    quantile(wflowers, probs = c(0.05, 0.95))
  # Upper pools
    wfuppers = apply(wfres17[ , 8:14], 1, mean)
    mean(wfuppers)
    quantile(wfuppers, probs = c(0.05, 0.95))  
 
# Test for differences in abundance based on presence/absence of streams
  streams = apply(wfres17[ , c(1,2,3,6,7)], 1, mean)
  nostream = apply(wfres17[ , c(4,5)], 1, mean)
  quantile(streams-nostream, probs=c(0.05, 0.95))

# Test for differences in abundance based on hydroperiod
  full = apply(wfres17[ , c(1,3,6,7)], 1, mean)
  empty = apply(wfres17[ , c(2,4,5)], 1, mean)
  quantile(full-empty, probs=c(0.05, 0.95))    

# 2017 Jeff data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vPools_2017.csv')

# Look at the first few lines of data
  head(herps)

# Create a new variable that indicates whether the pond was in the lower or the
# upper site
  for(i in 1:nrow(herps)){
    if(herps$pond[i] <= 8) {herps$site[i]=0
    } else {
      herps$site[i]=1
    }
  }

# Fix the date formatting so we can use it in R if needed
  herps$date = as.POSIXct(as.character(herps$date), format='%m/%d/%Y')

# Define day of year based on the new date column
  herps$day = yday(herps$date)

# Define week number based on the date column
  herps$week = as.numeric(as.factor(week(herps$date)))
  #herps = herps[herps$week==2 | herps$week==3, ]
  #herps$week = as.numeric(as.factor((herps$week)))

# Form initial values for wf abundance
  Nst = ddply(herps, c("pond", "week"), summarize, count=sum(jeff))

# Make a pond x week matrix
  Nst2 =  xtabs(count~pond+week, data=Nst)
  dimnames(Nst2)=NULL

# Get weekdays from data
  herps$weekday = as.numeric(as.factor(weekdays(herps$date)))

# Make a matrix of counts by weekday and week
  Nst3 = xtabs(jeff~pond + weekday + week, data = herps)
  dimnames(Nst3) = NULL

# Get pond as a numerically represented factor
  herps$pond = as.numeric(as.factor(herps$pond))

# 2017 Jeff model specification ----------------------------
# Define the base model
  modelString = "
    model{
      # Likelihood
        for(i in 1:nponds){
          for(t in 1:nweeks){
            N[i,t] ~ dpois(lambda[i, t])
          }
        }

      # Random effect of site and day on abundance
        for(z in 1:nponds){
          for(i in 1:7){
            for(t in 1:nweeks){
              C[z,i,t] ~ dbin(p[pond[z], week[t]], N[pond[z], week[t]])
            }
          }
        }

      # Priors
        for(i in 1:nponds){
          #meanP[i] ~ dbeta(1, 1)

          for(t in 1:nweeks){
            log(lambda[i,t]) <- alpha[i, t]
            alpha[i,t] ~ dnorm(0, 0.1)T(-10, 10)
            #p[i,t] <- meanP[i]
            p[i,t] ~ dbeta(1, 1)
          }
        }
    }
  "

# Write the model to a file
  writeLines(modelString, "jeffmodel17.txt")

# 2017 Jeff model settings and calibration ------------------------------------------
# Bundle the data for JAGS
  herps.data = list(
    nponds = length(unique(herps$pond)),
    nweeks = length(unique(herps$week)),
    C = Nst3,
    pond = unique(herps$pond),
    week = sort(unique(herps$week))
  )

# Provide initial values for the Gibbs sampler
  inits = function(){
    list(
      N = Nst2,
      #meanP = rbeta(15, 1, 1),
      alpha = structure(
              rnorm(length(unique(herps$week))*length(unique(herps$pond)), 0, 1),
              .Dim=c(length(unique(herps$pond)),length(unique(herps$week)))),
      p = structure(
          rbeta(length(unique(herps$week))*length(unique(herps$pond)), 1, 5),
         .Dim=c(length(unique(herps$pond)),length(unique(herps$week))))
    )
  }

# Parameters to monitor during estimation
  params = c(
    #'meanP',
    'N',
    'p'
  )

# # MCMC Settings
#   ni = 25000
#   nb = 5000
#   nt = 50
#   nc = 3

# Run the model
  jeff.model17 = jags(data = herps.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "jeffmodel17.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(jeff.model17, digits = 2)

# 2017 Jeff results ------------------------------------------------------------
# Collect results from model object
  jeffres17 = data.frame(jeff.model17$BUGSoutput$sims.list)
  
# Get abundance estimates for each week
  wk1 = jeffres17[ , 1:14]
  wk2 = jeffres17[ , 15:28]
  wk3 = jeffres17[ , 29:43]
  wk4 = jeffres17[ , 43:56]
  wk5 = jeffres17[ , 57:70]

# Now, get the week for max abundance of species in each of the 15 pools and put
# the posterior for that week into a list object.
  check = c()
  for(i in 1:14){
    check[i] =
      which(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]), mean(wk5[,i]))==
            max(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]), mean(wk5[,i]))))
  }

  big = list(wk1, wk2, wk3, wk4, wk5)

  jeffres17 = vector(mode='list', length=14)
  for(i in 1:14){
    jeffres17[[i]] = big[[check[i]]][,i]
  }
  jeffres17 = do.call(cbind, jeffres17)
  
# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      jefflowers = apply(jeffres17[ , 1:7], 1, mean)
      mean(jefflowers)
      quantile(jefflowers, probs = c(0.025, 0.975))
    # Mean abundance in lower pools
      jeffuppers = apply(jeffres17[ , 8:14], 1, mean)
      mean(jeffuppers)
      quantile(jeffuppers, probs = c(0.025, 0.975))  

# Test for differences in abundance based on presence/absence of streams
  streams = apply(jeffres17[ , c(1,2,3,6,7)], 1, mean)
  nostream = apply(jeffres17[ , c(4,5)], 1, mean)
  quantile(streams-nostream, probs=c(0.025, 0.975))

# Test for differences in abundance based on hydroperiod
  full = apply(jeffres17[ , c(1,3,6,7)], 1, mean)
  empty = apply(jeffres17[ , c(2,4,5)], 1, mean)
  quantile(full-empty, probs=c(0.025, 0.975))    

# 2017 Spot data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vPools_2017.csv')

# Look at the first few lines of data
  head(herps)

# Create a new variable that indicates whether the pond was in the lower or the
# upper site
  for(i in 1:nrow(herps)){
    if(herps$pond[i] <= 8) {herps$site[i]=0
    } else {
      herps$site[i]=1
    }
  }

# Fix the date formatting so we can use it in R if needed
  herps$date = as.POSIXct(as.character(herps$date), format='%m/%d/%Y')

# Define day of year based on the new date column
  herps$day = yday(herps$date)

# Define week number based on the date column
  herps$week = as.numeric(as.factor(week(herps$date)))
  #herps = herps[herps$week==2 | herps$week==3, ]
  #herps$week = as.numeric(as.factor((herps$week)))

# Form initial values for wf abundance
  Nst = ddply(herps, c("pond", "week"), summarize, count=sum(spot))

# Make a pond x week matrix
  Nst2 =  xtabs(count~pond+week, data=Nst)
  dimnames(Nst2)=NULL

# Get weekdays from data
  herps$weekday = as.numeric(as.factor(weekdays(herps$date)))

# Make a matrix of counts by weekday and week
  Nst3 = xtabs(spot~pond + weekday + week, data = herps)
  dimnames(Nst3) = NULL

# Get pond as a numerically represented factor
  herps$pond = as.numeric(as.factor(herps$pond))

# 2017 Spot model specification ----------------------------
# Define the base model
  modelString = "
    model{
      # Likelihood
        for(i in 1:nponds){
          for(t in 1:nweeks){
            N[i,t] ~ dpois(lambda[i, t])T(0,2000)
          }
        }

      # Random effect of site and day on abundance
        for(z in 1:nponds){
          for(i in 1:7){
            for(t in 1:nweeks){
              C[z,i,t] ~ dbin(p[pond[z], week[t]], N[pond[z], week[t]])
            }
          }
        }

      # Priors
        for(i in 1:nponds){
          #meanP[i] ~ dbeta(1, 1)

          for(t in 1:nweeks){
            log(lambda[i,t]) <- alpha[i, t]
            alpha[i,t] ~ dnorm(0, 0.1)T(-7, 7)
            #p[i,t] <- meanP[i]
            p[i,t] ~ dbeta(1, 1)
          }
        }
    }
  "

# Write the model to a file
  writeLines(modelString, "spotmodel17.txt")

# 2017 Spot model settings and calibration ------------------------------------------
# Bundle the data for JAGS
  herps.data = list(
    nponds = length(unique(herps$pond)),
    nweeks = length(unique(herps$week)),
    C = Nst3,
    pond = unique(herps$pond),
    week = sort(unique(herps$week))
  )

# Provide initial values for the Gibbs sampler
  inits = function(){
    list(
      N = Nst2,
      #meanP = rbeta(15, 1, 1),
      alpha = structure(
              rnorm(length(unique(herps$week))*length(unique(herps$pond)), 0, 1),
              .Dim=c(length(unique(herps$pond)),length(unique(herps$week)))),
      p = structure(
          rbeta(length(unique(herps$week))*length(unique(herps$pond)), 1, 5),
         .Dim=c(length(unique(herps$pond)),length(unique(herps$week))))
    )
  }

# Parameters to monitor during estimation
  params = c(
    #'meanP',
    'N',
    'p'
  )

# # MCMC Settings
#   ni = 25000
#   nb = 5000
#   nt = 50
#   nc = 3

# Run the model
  spot.model17 = jags(data = herps.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "spotmodel17.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(spot.model17, digits = 2)


# 2017 Spot results ------------------------------------------------------------
  spotres17 = data.frame(spot.model17$BUGSoutput$sims.list)
  
# Get abundance estimates for each week
  wk1 = spotres17[ , 1:14]
  wk2 = spotres17[ , 15:28]
  wk3 = spotres17[ , 29:43]
  wk4 = spotres17[ , 43:56]
  wk5 = spotres17[ , 57:70]

# Now, get the week for max abundance of species in each of the 15 pools and put
# the posterior for that week into a list object.
  check = c()
  for(i in 1:14){
    check[i] =
      which(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]), mean(wk5[,i]))==
            max(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]), mean(wk5[,i]))))
  }

  big = list(wk1, wk2, wk3, wk4, wk5)

  spotres17 = vector(mode='list', length=14)
  for(i in 1:14){
    spotres17[[i]] = big[[check[i]]][,i]
  }
  spotres17 = do.call(cbind, spotres17)
  
# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      spotlowers = apply(spotres17[ , 1:7], 1, mean)
      mean(spotlowers)
      quantile(spotlowers, probs = c(0.025, 0.975))
    # Mean abundance in lower pools
      spotuppers = apply(spotres17[ , 8:14], 1, mean)
      mean(spotuppers)
      quantile(spotuppers, probs = c(0.025, 0.975))  

# Test for differences in abundance based on presence/absence of streams
  streams = apply(spotres17[ , c(1,2,3,6,7)], 1, mean)
  nostream = apply(spotres17[ , c(4,5)], 1, mean)
  quantile(streams-nostream, probs=c(0.025, 0.975))

# Test for differences in abundance based on hydroperiod
  full = apply(spotres17[ , c(1,3,6,7)], 1, mean)
  empty = apply(spotres17[ , c(2,4,5)], 1, mean)
  quantile(full-empty, probs=c(0.025, 0.975))   
      
# 2017 ALL SPP FIGURES --------------------------------------------------------- 
# Box plots of abundance at each pool
  # WF
    par(mfrow=c(3,1), oma=c(3.5,5,2,1), mar=c(2,2,1,1))
    boxplot(wfres17, ylim=c(0, 1500), outline = FALSE,
            names=c(as.character(seq(1,7,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    abline(v=7.5, col='black', lwd=2, lty=2)
    text(x=11, y=1350, "Upper pools", cex=1.25)
    text(x=4, y=1350, "Lower pools", cex=1.25)  
    text(x=.5, y=1350, '(d)', cex=1.5)
  # Jeffs
    boxplot(jeffres17, ylim=c(0, 200), outline = FALSE,
            names=c(as.character(seq(1,7,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    mtext('Estimated abundance (N)', side=2, cex=1.25, line=5)
    abline(v=7.5, col='black', lwd=2, lty=2)
    text(x=.5, y=180, '(e)', cex=1.5)
  # Spotted salamanders
    boxplot(spotres17, ylim=c(0, 1500), outline = FALSE,
            names=c(as.character(seq(1,7,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    mtext('Vernal pool', side=1, cex=1.25, line = 4)
    abline(v=7.5, col='black', lwd=2, lty=2)
    text(x=.5, y=1350, '(f)', cex=1.5)   
    
# Histograms of differences between upper and lower sites   
  par(mfrow=c(3,1), oma=c(3,3,1,1), mar=c(2,1,1,1))
  # Woodfrogs
    differences = wfuppers-wflowers
    hist(differences, col='gray87', border='gray60', xaxt='n', yaxt='n',
         main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-500, 1500), breaks=25)
    axis(1, pos=0, cex=1.15)
    abline(v=mean(differences), col='black', lwd=3, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='black', lwd=2, lty=2)
    #text(x=300, y=800, expression(paste(mu, "= " , "77")), adj=0, cex=1.5)
    #text(x=300, y=700, expression(paste("95% CRI = 21-234")), adj=0, cex=1.5)
    text(-180, 500, '(a)', cex=1.5)#; axis(2)
  # Jefferson salamanders
    differences = jeffuppers-jefflowers
    hist(differences, col='gray87', border='gray60', xaxt='n', yaxt='n',
         main='', xlab='', ylab = '', xlim=c(-500, 1500), breaks=25)
    axis(1, pos=0, cex=1.15)
    mtext(side=2, 'Frequency', cex=1.25, line=2)
    abline(v=mean(differences), col='black', lwd=3, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='black', lwd=2, lty=2)
    #text(x=300, y=1500, expression(paste(mu, "= " , "203")), adj=0, cex=1.5)
    #text(x=300, y=1350, expression(paste("95% CRI = - 6-96")), adj=0, cex=1.5)    
    text(-180, 600, '(b)', cex=1.5)#; axis(2)
  # Spotted salamander
    differences = spotuppers-spotlowers
    hist(differences, col='gray87', border='gray60', xaxt='n', yaxt='n',
         main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-500, 1500), breaks=45)
    axis(1, pos=0, cex=1.15)
    mtext(expression(paste("Mean difference (N"["upper"], "- N"["lower"],")")),
          side = 1, cex=1.25, line = 4)
    abline(v=mean(differences), col='black', lwd=3, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='black', lwd=2, lty=2)
    #text(x=300, y=700, expression(paste(mu, "= " , "76")), adj=0, cex=1.5)
    #text(x=300, y=600, expression(paste("95% CRI = 9-224")), adj=0, cex=1.5)    
    text(-180, 150, '(c)', cex=1.5)#; axis(2)

# ------------------------------------------------------------------------- 
# Multi-year figures ------------------------------------------------------
# Box plots of abundance at each pool in each year
  # Wood frogs
    # 2016
      par(mfrow=c(3,2), oma=c(3,4,1,1), mar=c(2,2,1,2))
      boxplot(wfres16, ylim=c(0, 700), outline = FALSE,
              names=c(as.character(seq(1,8,1)),
                c("9", "10", "11", "18", "19", "20", "21")
              ),
              cex.axis = 1.15,
              yaxt='n',
              col='gray'
      )
      axis(2, las=2, cex.axis=1.15)
      abline(v=8.5, col='black', lwd=2, lty=2)
      text(x=11, y=650, "Upper pools", cex=1.25)
      text(x=4, y=650, "Lower pools", cex=1.25)  
      text(x=.5, y=650, '(a)', cex=1.5)
      mtext(side=3, '2016')
    # 2017  
      boxplot(wfres17, ylim=c(0, 1500), outline = FALSE,
              names=c(as.character(seq(1,7,1)),
                c("9", "10", "11", "18", "19", "20", "21")
              ),
              cex.axis = 1.15,
              yaxt='n',
              col='gray'
      )
      axis(2, las=2, cex.axis=1.15)
      abline(v=7.5, col='black', lwd=2, lty=2)
      text(x=11, y=1350, "Upper pools", cex=1.25)
      text(x=4, y=1350, "Lower pools", cex=1.25)  
      text(x=.5, y=1350, '(d)', cex=1.5)
      mtext(side=3, '2017')
  # Jeffs
    # 2016
      boxplot(jeffres16, ylim=c(0, 125), outline = FALSE,
              names=c(as.character(seq(1,8,1)),
                c("9", "10", "11", "18", "19", "20", "21")
              ),
              cex.axis = 1.15,
              yaxt='n',
              col='gray'
      )
      axis(2, las=2, cex.axis=1.15)
      mtext('Estimated abundance (N)', side=2, cex=1.25, line=3)
      abline(v=8.5, col='black', lwd=2, lty=2)
      text(x=.5, y=110, '(b)', cex=1.5)
    # 2017
      boxplot(jeffres17, ylim=c(0, 200), outline = FALSE,
              names=c(as.character(seq(1,7,1)),
                c("9", "10", "11", "18", "19", "20", "21")
              ),
              cex.axis = 1.15,
              yaxt='n',
              col='gray'
      )
      axis(2, las=2, cex.axis=1.15)
      abline(v=7.5, col='black', lwd=2, lty=2)
      text(x=.5, y=180, '(e)', cex=1.5)
  # Spotted salamanders
    # 2016
      boxplot(spotres16, ylim=c(0, 500), outline = FALSE,
              names=c(as.character(seq(1,8,1)),
                c("9", "10", "11", "18", "19", "20", "21")
              ),
              cex.axis = 1.15,
              yaxt='n',
              col='gray'
      )
      axis(2, las=2, cex.axis=1.15)
      abline(v=8.5, col='black', lwd=2, lty=2)
      text(x=.5, y=450, '(c)', cex=1.5)       
    # 2017
      boxplot(spotres17, ylim=c(0, 1500), outline = FALSE,
              names=c(as.character(seq(1,7,1)),
                c("9", "10", "11", "18", "19", "20", "21")
              ),
              cex.axis = 1.15,
              yaxt='n',
              col='gray'
      )
      axis(2, las=2, cex.axis=1.15)
      mtext('Vernal pool', side=1, cex=1.25, line = 4, adj=-.25)
      abline(v=7.5, col='black', lwd=2, lty=2)
      text(x=.5, y=1350, '(f)', cex=1.5)   
  
# Histograms of differences between upper and lower sites   
  # Set graphical parameters
    par(mfrow=c(3,2), oma=c(3,4,3,1), mar=c(2,1,1,1))
  # Woodfrogs
    # 2016
      differences = wfuppers16-wflowers16
      hist(differences, col='gray87', border='gray60', yaxt='n', xaxt='n',
           main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-100, 800))#, breaks=25)
      axis(1, pos=0, cex=1.15)
      abline(v=mean(differences), col='black', lwd=2, lty=1)
      abline(v=quantile(differences, probs=c(0.05, 0.95)),
             col='black', lwd=1, lty=2)  
      text(700, 4500, '(a)', cex=1.5)#; axis(2)
      mtext(side=3, '2016', line=2)
    # 2017
      differences = wfuppers-wflowers
      hist(differences, col='gray87', border='gray60', xaxt='n', yaxt='n',
           main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-100, 800))#, breaks=25)
      axis(1, pos=0, cex=1.15)
      abline(v=mean(differences), col='black', lwd=2, lty=1)
      abline(v=quantile(differences, probs=c(0.05, 0.95)),
             col='black', lwd=1, lty=2)  
      text(700, 5300, '(d)', cex=1.5)#; axis(2)
      mtext(side=3, '2017', line=2)      
  # Jefferson salamanders
    # 2016
      differences = jeffuppers16-jefflowers16
      hist(differences, col='gray87', border='gray60', yaxt='n', xaxt='n',
           main='', xlab='', ylab = '', xlim=c(-100, 800))#, breaks=25)
      axis(1, pos=0, cex=1.15)
      mtext(side=2, 'Frequency', cex=1.25, line=2)
      abline(v=mean(differences), col='black', lwd=2, lty=1)
      abline(v=quantile(differences, probs=c(0.05, 0.95)),
             col='black', lwd=1, lty=2)  
      text(700, 9000, '(b)', cex=1.5)#; axis(2)
    # 2017    
      differences = jeffuppers-jefflowers
      hist(differences, col='gray87', border='gray60', xaxt='n', yaxt='n',
           main='', xlab='', ylab = '', xlim=c(-100, 800))#, breaks=25)
      axis(1, pos=0, cex=1.15)
      abline(v=mean(differences), col='black', lwd=2, lty=1)
      abline(v=quantile(differences, probs=c(0.05, 0.95)),
             col='black', lwd=1, lty=2)  
      text(700, 9000, '(e)', cex=1.5)#; axis(2)
  # Spotted salamander
    # 2016
      differences = spotuppers16-spotlowers16
      hist(differences, col='gray87', border='gray60', yaxt='n', xaxt='n',
           main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-100, 800), breaks=25)
      axis(1, pos=0, cex=1.15)
      abline(v=mean(differences), col='black', lwd=2, lty=1)
      abline(v=quantile(differences, probs=c(0.05, 0.95)),
             col='black', lwd=1, lty=2)  
      text(700, 4500, '(c)', cex=1.5)#; axis(2)
    # 2017    
      differences = spotuppers-spotlowers
      hist(differences, col='gray87', border='gray60', xaxt='n', yaxt='n',
           main='', xlab='',cex.lab=1.25, ylab='', xlim=c(-100, 800))#, breaks=45)
      axis(1, pos=0, cex=1.15)
      mtext(expression(paste("Mean difference (N"["upper"], "- N"["lower"],")")),
            side = 1, cex=1.25, line = 4, adj=-1.25)
      abline(v=mean(differences), col='black', lwd=2, lty=1)
      abline(v=quantile(differences, probs=c(0.05, 0.95)),
             col='black', lwd=1, lty=2)  
      text(700, 2700, '(f)', cex=1.5)#; axis(2)    
      
# Multi-year derived quantities -----------------------------------------------
# Wood frogs
  # Test for differences in abundance based on presence/absence of streams
    wfstreams16 = apply(wfres16[ , c(1,2,3,6,7)], 1, mean)
    wfnostream16 = apply(wfres16[ , c(4,5,8)], 1, mean)
    wfdiffstream16 = wfstreams16-wfnostream16
    mean(wfdiffstream16)
    quantile(wfdiffstream16, probs=c(0.05, 0.95))
  # Test for differences in abundance based on hydroperiod
    wffull16 = apply(wfres16[ , c(1,3,6,7)], 1, mean)
    wfempty16 = apply(wfres16[ , c(2,4,5,8)], 1, mean)
    wfdiffhydro16 = wffull16-wfempty16
    mean(wfdiffhydro16)
    quantile(wfdiffhydro16, probs=c(0.05, 0.95))
  # Test for differences in abundance based on presence/absence of streams
    wfstreams17 = apply(wfres17[ , c(1,2,3,6,7)], 1, mean)
    wfnostream17 = apply(wfres17[ , c(4,5,8)], 1, mean)
    wfdiffstream17 = wfstreams17-wfnostream17
    mean(wfdiffstream17)
    quantile(wfdiffstream17, probs=c(0.05, 0.95))
  # Test for differences in abundance based on hydroperiod
    wffull17 = apply(wfres17[ , c(1,3,6,7)], 1, mean)
    wfempty17 = apply(wfres17[ , c(2,4,5,8)], 1, mean)
    wfdiffhydro17 = wffull17-wfempty17
    mean(wfdiffhydro17)
    quantile(wfdiffhydro17, probs=c(0.05, 0.95))
    
  # Calculate averages across years
    mean((wfdiffhydro16+wfdiffhydro17)/2)  
    quantile((wfdiffhydro16+wfdiffhydro17)/2, c(0.05, 0.95))
    
    mean((wfdiffstream16+wfdiffstream17)/2) 
    quantile((wfdiffstream16+wfdiffstream17)/2, c(0.05, 0.95)) 
  
# Jefferson salamanders  
  # Test for differences between sites  
  # Calculate mean and 90% Credible intervals in lower and upper pools overall
    # 2016
      # Mean abundance in lower pools
        jefflowers16 = apply(jeffres16[ , 1:8], 1, mean)
        mean(jefflowers16)
        quantile(jefflowers16, probs = c(0.05, 0.95))
      # Mean abundance in upper pools
        jeffuppers16 = apply(jeffres16[ , 9:15], 1, mean)
        mean(jeffuppers16)
        quantile(jeffuppers16, probs = c(0.05, 0.95)) 
      # Differences
        jeffdiff16 = jeffuppers16-jefflowers16
        mean(jeffdiff16)   
        quantile(jeffdiff16, probs = c(0.05, 0.95))        
         
    # 2017    
      # Mean abundance in lower pools
        jefflowers17 = apply(jeffres17[ , 1:7], 1, mean)
        mean(jefflowers17)
        quantile(jefflowers17, probs = c(0.05, 0.95))
      # Mean abundance in upper pools
        jeffuppers17 = apply(jeffres17[ , 8:14], 1, mean)
        mean(jeffuppers17)
        quantile(jeffuppers17, probs = c(0.05, 0.95))  
      # Differences
        jeffdiff17 = jeffuppers17-jefflowers17
        mean(jeffdiff17)   
        quantile(jeffdiff17, probs = c(0.05, 0.95))         
    # Across years
      mean(c(jeffdiff16, jeffdiff17))
      quantile(c(jeffdiff16, jeffdiff17), probs = c(0.05, 0.95))
        
  # Test for differences in abundance based on presence/absence of streams
    jeffstreams16 = apply(jeffres16[ , c(1,2,3,6,7)], 1, mean)
    jeffnostream16 = apply(jeffres16[ , c(4,5,8)], 1, mean)
    jeffdiffstream16 = jeffstreams16-jeffnostream16
    mean(jeffdiffstream16)
    quantile(jeffdiffstream16, probs=c(0.05, 0.95))
  # Test for differences in abundance based on hydroperiod
    jefffull16 = apply(jeffres16[ , c(1,3,6,7)], 1, mean)
    jeffempty16 = apply(jeffres16[ , c(2,4,5,8)], 1, mean)
    jeffdiffhydro16 = jefffull16-jeffempty16
    mean(jeffdiffhydro16)
    quantile(jeffdiffhydro16, probs=c(0.05, 0.95))
  # Test for differences in abundance based on presence/absence of streams
    jeffstreams17 = apply(jeffres17[ , c(1,2,3,6,7)], 1, mean)
    jeffnostream17 = apply(jeffres17[ , c(4,5,8)], 1, mean)
    jeffdiffstream17 = jeffstreams17-jeffnostream17
    mean(jeffdiffstream17)
    quantile(jeffdiffstream17, probs=c(0.05, 0.95))
  # Test for differences in abundance based on hydroperiod
    jefffull17 = apply(jeffres17[ , c(1,3,6,7)], 1, mean)
    jeffempty17 = apply(jeffres17[ , c(2,4,5,8)], 1, mean)
    jeffdiffhydro17 = jefffull17-jeffempty17
    mean(jeffdiffhydro17)
    quantile(jeffdiffhydro17, probs=c(0.05, 0.95))
    
  # Calculate averages across years
    mean((jeffdiffhydro16+jeffdiffhydro17)/2)  
    quantile((jeffdiffhydro16+jeffdiffhydro17)/2, c(0.05, 0.95))
    
    mean((jeffdiffstream16+jeffdiffstream17)/2) 
    quantile((jeffdiffstream16+jeffdiffstream17)/2, c(0.05, 0.95)) 
 
# Spotted salamanders  
  # Test for differences between sites  
  # Calculate mean and 90% Credible intervals in lower and upper pools overall
    # 2016
      # Mean abundance in lower pools
        spotlowers16 = apply(spotres16[ , 1:8], 1, mean)
        mean(spotlowers16)
        quantile(spotlowers16, probs = c(0.05, 0.95))
      # Mean abundance in upper pools
        spotuppers16 = apply(spotres16[ , 9:15], 1, mean)
        mean(spotuppers16)
        quantile(spotuppers16, probs = c(0.05, 0.95)) 
      # Differences
        spotdiff16 = spotuppers16-spotlowers16
        mean(spotdiff16)   
        quantile(spotdiff16, probs = c(0.05, 0.95))        
         
    # 2017    
      # Mean abundance in lower pools
        spotlowers17 = apply(spotres17[ , 1:7], 1, mean)
        mean(spotlowers17)
        quantile(spotlowers17, probs = c(0.05, 0.95))
      # Mean abundance in upper pools
        spotuppers17 = apply(spotres17[ , 8:14], 1, mean)
        mean(spotuppers17)
        quantile(spotuppers17, probs = c(0.05, 0.95))  
      # Differences
        spotdiff17 = spotuppers17-spotlowers17
        mean(spotdiff17)   
        quantile(spotdiff17, probs = c(0.05, 0.95))         

    # Across years
      mean(c(spotdiff16, spotdiff17))
      quantile(c(spotdiff16, spotdiff17), probs = c(0.05, 0.95))        
            
  # Test for differences in abundance based on presence/absence of streams
    spotstreams16 = apply(spotres16[ , c(1,2,3,6,7)], 1, mean)
    spotnostream16 = apply(spotres16[ , c(4,5,8)], 1, mean)
    spotdiffstream16 = spotstreams16-spotnostream16
    mean(spotdiffstream16)
    quantile(spotdiffstream16, probs=c(0.05, 0.95))
  # Test for differences in abundance based on hydroperiod
    spotfull16 = apply(spotres16[ , c(1,3,6,7)], 1, mean)
    spotempty16 = apply(spotres16[ , c(2,4,5,8)], 1, mean)
    spotdiffhydro16 = spotfull16-spotempty16
    mean(spotdiffhydro16)
    quantile(spotdiffhydro16, probs=c(0.05, 0.95))
  # Test for differences in abundance based on presence/absence of streams
    spotstreams17 = apply(spotres17[ , c(1,2,3,6,7)], 1, mean)
    spotnostream17 = apply(spotres17[ , c(4,5,8)], 1, mean)
    spotdiffstream17 = spotstreams17-spotnostream17
    mean(spotdiffstream17)
    quantile(spotdiffstream17, probs=c(0.05, 0.95))
  # Test for differences in abundance based on hydroperiod
    spotfull17 = apply(spotres17[ , c(1,3,6,7)], 1, mean)
    spotempty17 = apply(spotres17[ , c(2,4,5,8)], 1, mean)
    spotdiffhydro17 = spotfull17-spotempty17
    mean(spotdiffhydro17)
    quantile(spotdiffhydro17, probs=c(0.05, 0.95))
    
  # Calculate averages across years
    mean((spotdiffhydro16+spotdiffhydro17)/2)  
    quantile((spotdiffhydro16+spotdiffhydro17)/2, c(0.05, 0.95))
    
    mean((spotdiffstream16+spotdiffstream17)/2) 
    quantile((spotdiffstream16+spotdiffstream17)/2, c(0.05, 0.95))            