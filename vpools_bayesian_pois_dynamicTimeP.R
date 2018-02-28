# Front-end ---------------------------------------------------------------
  library(R2jags)
  library(plyr)
  library(lubridate)
  library(plotrix)
  library(MASS)

# WF data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vpools.csv')

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

# WF mixture model specification ----------------------------
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
  writeLines(modelString, "wfmodel.txt")

# WF model settings and calibration ------------------------------------------
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

# MCMC Settings
  ni = 3100
  nb = 100
  nt = 3
  nc = 3

# Run the model
  wf.model = jags(data = herps.data,
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
  print(wf.model, digits = 2)

# Results WF --------------------------------------------------------------
# Get the results
  wfres = data.frame(wf.model$BUGSoutput$sims.list)
  
# Get abundance estimates for each week
  wk1 = wfres[ , 1:15]
  wk2 = wfres[ , 16:30]
  wk3 = wfres[ , 31:45]
  wk4 = wfres[ , 46:60]

# Now, get the week for max abundance of species in each of the 15 pools and put
# the posterior for that week into a list object.
  check = c()
  for(i in 1:15){
    check[i] =
      which(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))==
            max(c(mean(wk1[,i]), mean(wk2[,i]), mean(wk3[,i]), mean(wk4[,i]))))
  }

  big = list(wk1, wk2, wk3, wk4)

  wfres = vector(mode='list', length=15)
  for(i in 1:15){
    wfres[[i]] = big[[check[i]]][,i]
  }
  wfres = do.call(cbind, wfres)
  
# Make boxplot of abundance by pool in upper and lower sites
  par(mar=c(5,6,1,1))
  boxplot(wfres, ylim=c(0, 700), outline = FALSE,
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
      wflowers = apply(wfres[ , 1:8], 1, mean)
      mean(wflowers)
      quantile(wflowers, probs = c(0.025, 0.975))

    # Mean abundance in lower pools
      wfuppers = apply(wfres[ , 9:15], 1, mean)
      mean(wfuppers)
      quantile(wfuppers, probs = c(0.025, 0.975))  
 
  # Make the plots
    differences = wfuppers-wflowers
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

# Test for differences in abundance based on presence/absence of streams
streams = apply(wfres[ , c(1,2,3,6,7)], 1, mean)
nostream = apply(wfres[ , c(4,5,8)], 1, mean)
quantile(streams-nostream, probs=c(0.025, 0.975))

# Test for differences in abundance based on hydroperiod
full = apply(wfres[ , c(1,3,6,7)], 1, mean)
empty = apply(wfres[ , c(2,4,5,8)], 1, mean)
quantile(full-empty, probs=c(0.025, 0.975))    

# Jeff data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vpools.csv')

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

# Jeff model specification ----------------------------
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

# Jeff model settings and calibration ------------------------------------------
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

# MCMC Settings
  ni = 31000
  nb = 1000
  nt = 30
  nc = 3

# Run the model
  jeff.model = jags(data = herps.data,
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
  print(jeff.model, digits = 2)

# Jeff results ------------------------------------------------------------
# Collect results from model object
  jres = data.frame(jeff.model$BUGSoutput$sims.list)
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

  jeffres = vector(mode='list', length=15)
  for(i in 1:15){
    jeffres[[i]] = big[[check[i]]][,i]
  }
  jeffres = do.call(cbind, jeffres)
  
# Make boxplot of abundance by pool in upper and lower sites
  par(mar=c(5,6,1,1))
  boxplot(jeffres, ylim=c(0, 125), outline = FALSE,
          names=c(as.character(seq(1,8,1)),
            c("9", "10", "11", "18", "19", "20", "21")
          ),
          cex.axis = 1.15,
          yaxt='n',
          ylim = c(0, 2000),
          col='gray'
  )
  axis(2, las=2, cex.axis=1.15)
  mtext('Estimated abundance (N)', side=2, cex=1.25, line=5)
  mtext('Vernal pool', side=1, cex=1.25, line = 4)
  abline(v=8.5, col='black', lwd=2, lty=2)
  text(x = 1, y=120, "Lower pools")
  text(x = 14.5, y=120, "Upper pools")

# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      jefflowers = apply(jeffres[ , 1:8], 1, mean)
      mean(jefflowers)
      quantile(jefflowers, probs = c(0.025, 0.975))

    # Mean abundance in lower pools
      jeffuppers = apply(jeffres[ , 9:15], 1, mean)
      mean(jeffuppers)
      quantile(jeffuppers, probs = c(0.025, 0.975))  
 
  # Make the plots
    differences = jeffuppers-jefflowers
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

# Test for differences in abundance based on presence/absence of streams
streams = apply(jeffres[ , c(1,2,3,6,7)], 1, mean)
nostream = apply(jeffres[ , c(4,5,8)], 1, mean)
quantile(streams-nostream, probs=c(0.025, 0.975))

# Test for differences in abundance based on hydroperiod
full = apply(jeffres[ , c(1,3,6,7)], 1, mean)
empty = apply(jeffres[ , c(2,4,5,8)], 1, mean)
quantile(full-empty, probs=c(0.025, 0.975))        
    
# Spot data manipulation -------------------------------------------------------
# Read in the data file
  herps = read.csv('vpools.csv')

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

# Spot model specification ----------------------------
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
  writeLines(modelString, "spotmodel.txt")

# Spot model settings and calibration ------------------------------------------
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

# MCMC Settings
  ni = 31000
  nb = 1000
  nt = 30
  nc = 3

# Run the model
  spot.model = jags(data = herps.data,
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
  print(spot.model, digits = 2)


# Spot results ------------------------------------------------------------
# Collect results from model object
  sres = data.frame(spot.model$BUGSoutput$sims.list)
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

  spotres = vector(mode='list', length=15)
  for(i in 1:15){
    spotres[[i]] = big[[check[i]]][,i]
  }
  spotres = do.call(cbind, spotres)
# Make boxplot of abundance by pool in upper and lower sites
  par(mar=c(5,6,1,1))
  boxplot(spotres, ylim=c(0, 500), outline = FALSE,
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
  text(x = 1, y=500, "Lower pools")
  text(x = 14.5, y=500, "Upper pools")

# Calculate the differences in abunance at upper and lower pools
  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      spotlowers = apply(spotres[ , 1:8], 1, mean)
      mean(spotlowers)
      quantile(spotlowers, probs = c(0.025, 0.975))

    # Mean abundance in lower pools
      spotuppers = apply(spotres[ , 9:15], 1, mean)
      mean(spotuppers)
      quantile(spotuppers, probs = c(0.025, 0.975))  
 
  # Make the plots
    differences = spotuppers-spotlowers
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

# Test for differences in abundance based on presence/absence of streams
streams = apply(spotres[ , c(1,2,3,6,7)], 1, mean)
nostream = apply(spotres[ , c(4,5,8)], 1, mean)
quantile(streams-nostream, probs=c(0.025, 0.975))

# Test for differences in abundance based on hydroperiod
full = apply(spotres[ , c(1,3,6,7)], 1, mean)
empty = apply(spotres[ , c(2,4,5,8)], 1, mean)
quantile(full-empty, probs=c(0.025, 0.975))         
      
# ALL SPP FIGURES --------------------------------------------------------- 
# Box plots of abundance at each pool
  # WF
    par(mfrow=c(3,1), oma=c(1,2.5,1,1))
    boxplot(wfres, ylim=c(0, 700), outline = FALSE,
            names=c(as.character(seq(1,8,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    abline(v=8.5, col='black', lwd=2, lty=2)
    mtext(side=3, "Upper pools", adj=.9)
    mtext(side=3, "Lower pools", adj=0.1)  
    text(x=.5, y=650, '(a)', cex=1.5)
  # Jeffs
    boxplot(jeffres, ylim=c(0, 125), outline = FALSE,
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
    boxplot(spotres, ylim=c(0, 500), outline = FALSE,
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
    differences = wfuppers-wflowers
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
    differences = jeffuppers-jefflowers
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
    differences = spotuppers-spotlowers
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
