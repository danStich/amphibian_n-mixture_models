# Front-end ---------------------------------------------------------------
  library(R2jags)
  library(plyr)
  library(lubridate)
  library(plotrix)

# Data manipulation -------------------------------------------------------
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

# Drop data from low-observation period
# herps = subset(herps, date>='2016-03-25' & date <='2016-04-08')

# Define day of year based on the new date column
  herps$day = yday(herps$date)

# Bayesian binomial mixture model specification ----------------------------
# Define the base model
  modelString = "
    model {

      # Likelihood
        for(i in 1:nponds){
          N[i] ~ dnegbin(p.N[i], r[i])
        }

      # Random effect of site and day on abundance
        for(i in 1:n){
          C[i] ~ dbin(p[pond[i]], N[pond[i]])
        }

      # Priors
        for(i in 1:nponds){
          p.N[i] ~ dbeta(1, 1)
          p[i] ~ dbeta(1, 1)
          log(r[i]) <- lr[i]
        }

        for(i in 1:nponds){
          lr[i] ~ dnorm(0, 0.0001)
        }
        

    }
  "

# Write the model to a file
  writeLines(modelString, "model.txt")

# Model settings and calibration wf ----------------------------------------
# Bundle the data for JAGS
  herps.data = list(
    n = nrow(herps),
    nponds = length(unique(herps$pond)),
    C = herps[,3],
    pond = rep(seq(1,15,1), 25)
  )

# Provide initial values for the Gibbs sampler
  # Form initial values for wf abundance
  Nst = ddply(herps, "pond", summarize, count=sum(wf))[,2]+1
  inits = function(){
    list(
      N = Nst,
      p = rbeta(15, 1, 1),
      p.N = rbeta(15, 1, 1),
      lr = rnorm(15, 1, 10)
    )
  }

# Parameters to monitor during estimation
  params = c(
    'N',
    'p',
    'p.N',
    'r'
  )

# MCMC Settings
  ni = 3300000
  nb = 300000
  nt = 3000
  nc = 3

# Run the model
  herps.model = jags(data = herps.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "model.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(herps.model, digits = 2)

# Get the results
  res = data.frame(herps.model$BUGSoutput$sims.list)

# Quick look at some convergence diagnostics
  # plot(res$N.13[1:(nrow(res)/3)], type='l', col='blue')
  # lines(res$N.13[((nrow(res)/3)+1):((nrow(res)/3)*2)], type='l', col='red')
  # lines(res$N.13[((nrow(res)/3)*2+1):((nrow(res)/3)*3)], type='l', col='green')

# Results wf-----------------------------------------------------------------
# Abundance
  # First, a table of the results that can be put in Excel
    results.table = herps.model$BUGSoutput$summary[1:15, c(1,2,3,7)]
    results.table = apply(results.table, c(1,2), round, digits=2)
    write.table(results.table, 'populationEstimates.csv', quote=FALSE,
                sep=",")

  # Boxplots by pool
    par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,6,1,1))
    boxplot(res[ , c(1:15)],
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
    abline(v=8.5, col='blue', lwd=2, lty=2)
    text(x = 1, y=2000, "Lower pools")
    text(x = 14.5, y=2000, "Upper pools")

  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      lowers = apply(res[ , 1:8], 1, mean)
      mean(lowers)
      quantile(lowers, probs = c(0.025, 0.975))

    # Mean abundance in lower pools
      uppers = apply(res[ , 9:16], 1, mean)
      mean(uppers)
      quantile(uppers, probs = c(0.025, 0.975))

  # # Plot the posterior estimates for population abundance at upper and lower
  # # pools as density curves
  #  # Set up plotting window
  # 	 par(mfrow=c(2,1), oma=c(5,6,.5,1))
  # 	 par(mar=c(1,1,1.5,1))
  #
  # 	 # Lower pools
  #   	 plot(density(lowers), col='red', xlim=c(0,1000), main='',
  #   	      yaxt = 'n', lwd=2)
  #   	 axis(2, las=2)
  #      text(x=800, y=.19, expression(paste(mu, "= " , "5")), adj=0)
  #      text(x=800, y=.17, expression(paste("95% CRI = 1-26")), adj=0)
  #   	 mtext("Posterior density", side=2, line = 5, cex=1.5, adj=-2)
  #
  # 	 # Upper pools
  #   	 plot(density(uppers), col='blue', xlim=c(0,1000), main='',
  #   	      yaxt = 'n', lwd=2)
  #   	 axis(2, las=2)
  #      text(x=800, y=.007, expression(paste(mu, "= " , "207")), adj=0)
  #      text(x=800, y=.006, expression(paste("95% CRI = 104-600")), adj=0)
  #   	 mtext("Mean estimated abundance (N)", side=1, line=4, cex=1.25)

  # # Plot the posterior estimates for population abundance at upper and lower
  # # pools as histograms
  #  # Set up plotting window
  # 	 par(mfrow=c(2,1), oma=c(5,6,.5,1))
  # 	 par(mar=c(1,1,1.5,1))
  #
  # 	 # Lower pools
  #   	 hist(lowers, col='red', xlim=c(0,1000), main='',
  #   	      yaxt = 'n', ylim=c(0,600))
  #   	 axis(2, las=2)
  #   	 legend(x=800, y=600, legend=c('Lower pools', 'Upper pools'),
  #   	        fill=c('red', 'blue'), bty='n')
  #      text(x=800, y=200, expression(paste(mu, "= " , "5")), adj=0)
  #      text(x=800, y=100, expression(paste("95% CRI = 1-26")), adj=0)
  #
  # 	 # Upper pools
  #   	 hist(uppers, col='blue', xlim=c(0,1000), main='',
  #   	      yaxt = 'n', ylim=c(0,600))
  #   	 axis(2, las=2)
  #      text(x=800, y=200, expression(paste(mu, "= " , "207")), adj=0)
  #      text(x=800, y=100, expression(paste("95% CRI = 104-600")), adj=0)
  #   	 mtext("Frequency", side=2, line = 5, cex=1.5, adj=2)
  #   	 mtext("Mean estimated abundance (N)", side=1, line=4, cex=1.5)

  # Calculate the differences in abunance at upper and lower pools
    differences = uppers-lowers
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

# # Detection probabilities
#   # Boxplots by pool
#     par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,6,1,1))
#     boxplot(res[ , c(17:31)],
#             names=c(as.character(seq(1,8,1)),
#               c("9", "10", "11", "18", "19", "20", "21")
#             ),
#             cex.axis = 1.15,
#             yaxt='n',
#             ylim = c(0, 1),
#             col='gray'
#     )
#     axis(2, las=2, cex.axis=1.15)
#     mtext('Detection probability (p)', side=2, cex=1.25, line=5)
#     mtext('Vernal pool', side=1, cex=1.25, line = 4)
#     abline(v=8.5, col='blue', lwd=2, lty=2)
#     text(x = 1, y=1, "Lower pools")
#     text(x = 14.5, y=1, "Upper pools")

# # Histogram of differences in detection probabilities between upper and lower
# # pools
#   # Mean for lower
#     lowerP = apply(res[ , 17:20, 22:23], 1, mean)
#     mean(lowerP)
#   # Mean for upper pools
#     upperP = apply(res[ , 26:31], 1, mean)
#     mean(upperP)
#
#   # Histogram
#   differences = upperP-lowerP
#   par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,6,1,1))
#   hist(differences, col='gray60', border='gray40', yaxt='n', xaxt='n',
#        main='', xlab='',cex.lab=1.25)
#   axis(1, pos=0, cex=1.15)
#   axis(2, las=2, pos=-0.12, cex=1.15)
#   mtext(expression(paste("Mean difference (p"["upper"], "- p"["lower"],")")),
#         side = 1, cex=1.25, line = 4)
#   abline(v=mean(differences), col='blue', lwd=2, lty=1)
#   abline(v=quantile(differences, probs=c(0.025, 0.975)),
#          col='red', lwd=2, lty=2)

# Test for differences in abundance based on presence/absence of streams
streams = apply(herps.model$BUGSoutput$sims.list$N[ , c(1,2,3,6,7)], 1, mean)
nostream = apply(herps.model$BUGSoutput$sims.list$N[ , c(4,5,8)], 1, mean)
quantile(streams-nostream, probs=c(0.025, 0.975))

# Test for differences in abundance based on hydroperiod
full = apply(herps.model$BUGSoutput$sims.list$N[ , c(1,3,6,7)], 1, mean)
empty = apply(herps.model$BUGSoutput$sims.list$N[ , c(2,4,5,8)], 1, mean)
quantile(full-empty, probs=c(0.025, 0.975))


# Model settings and calibration jeff ----------------------------------------
# Bundle the data for JAGS
  jeff.data = list(
    n = nrow(herps),
    nponds = length(unique(herps$pond)),
    C = herps[,4],
    pond = rep(seq(1,15,1), 25)
  )

# Provide initial values for the Gibbs sampler
  # Form initial values for wf abundance
  Nst = ddply(herps, "pond", summarize, count=sum(jeff))[,2]+1
  inits = function(){
    list(
      N = Nst,
      p = rbeta(15, 1, 1),
      p.N = rbeta(15, 1, 1),
      r = rgamma(15, 0.01, 0.01)
    )
  }

# Parameters to monitor during estimation
  params = c(
    'N',
    'p',
    'p.N',
    'r'
  )

# MCMC Settings
  ni = 15000
  nb = 5000
  nt = 50
  nc = 3

# Run the model
  jeff.model = jags(data = jeff.data,
                     inits = inits,
                     parameters.to.save = params,
                     model.file = "model.txt",
                     n.iter = ni,
                     n.burnin = nb,
                     n.thin = nt,
                     n.chains = nc,
                     working.directory = getwd()
                     )

# Print the model results
  print(jeff.model, digits = 2)

# Get the results
  res = data.frame(jeff.model$BUGSoutput$sims.list)
  #hist(res$N.1)

# Quick look at some convergence diagnostics
  plot(res$N.13[1:(nrow(res)/3)], type='l', col='blue')
  lines(res$N.13[((nrow(res)/3)+1):((nrow(res)/3)*2)], type='l', col='red')
  lines(res$N.13[((nrow(res)/3)*2+1):((nrow(res)/3)*3)], type='l', col='green')

# Results jeff-----------------------------------------------------------------
# Abundance
  # First, a table of the results that can be put in Excel
    results.table = jeff.model$BUGSoutput$summary[1:15, c(1,2,3,7)]
    results.table = apply(results.table, c(1,2), round, digits=2)
    write.table(results.table, 'populationEstimates.csv', quote=FALSE,
                sep=",")

  # Boxplots by pool
    par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,6,1,1))
    boxplot(res[ , c(1:15)],
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
    abline(v=8.5, col='blue', lwd=2, lty=2)
    text(x = 1, y=2000, "Lower pools")
    text(x = 14.5, y=2000, "Upper pools")

  # Calculate mean and 95% Credible intervals in lower and upper pools overall
    # Mean abundance in lower pools
      lowers = apply(res[ , 1:8], 1, mean)
      mean(lowers)
      quantile(lowers, probs = c(0.025, 0.975))

    # Mean abundance in lower pools
      uppers = apply(res[ , 9:16], 1, mean)
      mean(uppers)
      quantile(uppers, probs = c(0.025, 0.975))

  # Plot the posterior estimates for population abundance at upper and lower
  # pools as density curves
   # Set up plotting window
  	 par(mfrow=c(2,1), oma=c(5,6,.5,1))
  	 par(mar=c(1,1,1.5,1))

  	 # Lower pools
    	 plot(density(lowers), col='red', xlim=c(0,1000), main='',
    	      yaxt = 'n', lwd=2)
    	 axis(2, las=2)
       text(x=800, y=.19, expression(paste(mu, "= " , "5")), adj=0)
       text(x=800, y=.17, expression(paste("95% CRI = 1-26")), adj=0)
    	 mtext("Posterior density", side=2, line = 5, cex=1.5, adj=-2)

  	 # Upper pools
    	 plot(density(uppers), col='blue', xlim=c(0,1000), main='',
    	      yaxt = 'n', lwd=2)
    	 axis(2, las=2)
       text(x=800, y=.007, expression(paste(mu, "= " , "207")), adj=0)
       text(x=800, y=.006, expression(paste("95% CRI = 104-600")), adj=0)
    	 mtext("Mean estimated abundance (N)", side=1, line=4, cex=1.25)

  # Plot the posterior estimates for population abundance at upper and lower
  # pools as histograms
   # Set up plotting window
  	 par(mfrow=c(2,1), oma=c(5,6,.5,1))
  	 par(mar=c(1,1,1.5,1))

  	 # Lower pools
    	 hist(lowers, col='red', xlim=c(0,1000), main='',
    	      yaxt = 'n', ylim=c(0,600))
    	 axis(2, las=2)
    	 legend(x=800, y=600, legend=c('Lower pools', 'Upper pools'),
    	        fill=c('red', 'blue'), bty='n')
       text(x=800, y=200, expression(paste(mu, "= " , "5")), adj=0)
       text(x=800, y=100, expression(paste("95% CRI = 1-26")), adj=0)

  	 # Upper pools
    	 hist(uppers, col='blue', xlim=c(0,1000), main='',
    	      yaxt = 'n', ylim=c(0,600))
    	 axis(2, las=2)
       text(x=800, y=200, expression(paste(mu, "= " , "207")), adj=0)
       text(x=800, y=100, expression(paste("95% CRI = 104-600")), adj=0)
    	 mtext("Frequency", side=2, line = 5, cex=1.5, adj=2)
    	 mtext("Mean estimated abundance (N)", side=1, line=4, cex=1.5)

  # Calculate the differences in abunance at upper and lower pools
    differences = uppers-lowers
    par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,6,1,1))
    hist(differences, col='gray60', border='gray40', yaxt='n', xaxt='n',
         main='', xlab='',cex.lab=1.25)
    axis(1, pos=0, cex=1.15)
    axis(2, pos=0, las=2, cex=1.15)
    mtext(expression(paste("Mean difference (N"["upper"], "- N"["lower"],")")),
          side = 1, cex=1.25, line = 4)
    abline(v=mean(differences), col='blue', lwd=2, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='red', lwd=2, lty=2)
    text(x=1000, y=400, expression(paste(mu, "= " , "202")), adj=0)
    text(x=1000, y=375, expression(paste("95% CRI = 98-590")), adj=0)

# Detection probabilities
  # Boxplots by pool
    par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,6,1,1))
    boxplot(res[ , c(17:31)],
            names=c(as.character(seq(1,8,1)),
              c("9", "10", "11", "18", "19", "20", "21")
            ),
            cex.axis = 1.15,
            yaxt='n',
            ylim = c(0, 1),
            col='gray'
    )
    axis(2, las=2, cex.axis=1.15)
    mtext('Detection probability (p)', side=2, cex=1.25, line=5)
    mtext('Vernal pool', side=1, cex=1.25, line = 4)
    abline(v=8.5, col='blue', lwd=2, lty=2)
    text(x = 1, y=1, "Lower pools")
    text(x = 14.5, y=1, "Upper pools")

  # Histogram of differences in detection probabilities between upper and lower
  # pools
    # Mean for lower
      lowerP = apply(res[ , 17:20, 22:23], 1, mean)
      mean(lowerP)
    # Mean for upper pools
      upperP = apply(res[ , 26:31], 1, mean)
      mean(upperP)

    # Histogram
    differences = upperP-lowerP
    par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(5,6,1,1))
    hist(differences, col='gray60', border='gray40', yaxt='n', xaxt='n',
         main='', xlab='',cex.lab=1.25)
    axis(1, pos=0, cex=1.15)
    axis(2, las=2, pos=-0.12, cex=1.15)
    mtext(expression(paste("Mean difference (p"["upper"], "- p"["lower"],")")),
          side = 1, cex=1.25, line = 4)
    abline(v=mean(differences), col='blue', lwd=2, lty=1)
    abline(v=quantile(differences, probs=c(0.025, 0.975)),
           col='red', lwd=2, lty=2)
