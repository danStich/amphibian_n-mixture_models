# Front-end ---------------------------------------------------------------
  library(R2jags)
  library(plyr)
  library(lubridate)
  library(reshape)
  library(reshape2)

# Data manipulation -------------------------------------------------------
# Daily count data
  # Read in the herps data
    herps = read.csv('vpools_2016.csv')
    # herps2 = read.csv('vPools_2017.csv')
    # herps2 = herps2[1:6]
    # herps3 = rbind(herps, herps2)
  # Convert the dataframe from wide format to long format
    vp = melt(herps, id.vars=1:2)
    
  # Change the 'values' column to 'count'
    names(vp)[3:4] = c('species', 'count')
    vp = vp[vp$species!='newt', ]
    
  # Create a column for occupancy
    vp$occ = 0
    vp$occ[vp$count>0]=1
  
  # Convert pond to a numerically represented factor for looping in JAGS
    vp$pond = as.numeric(as.factor(vp$pond))
    
  # Create a new variable that indicates whether the pond was in the lower or 
  # upper site and fill it in.
    for(i in 1:nrow(vp)){
      if(vp$pond[i] <= 8) {vp$site[i]=0
      } else {
        vp$site[i]=1
      }
    }
  
  # Fix the date formatting so we can use it in R if needed
    vp$date = as.POSIXct(as.character(vp$date), format='%m/%d/%Y')
  
  # Define day of year based on the new date column
    vp$day = yday(vp$date)
    
  # Create a column for year
    vp$year = year(vp$date)
      
  # Create reps within each site based on day and year
    # Split the dataframe by site. This makes a list of dataframes split out by
    # sampling location
      temp.rep = split(vp, list(vp$year, vp$pond))
    # Create a new column for each site that contains numeric ordering of
    # sample reps
      for(i in 1:length(temp.rep)){
        temp.rep[[i]]$rep = as.numeric(as.factor(temp.rep[[i]]$date))
      }
    # Put all of the dfs in the list back into a single dataframe by row binding
      vp = do.call('rbind', temp.rep)
    
  # Define week number based on the date column
    vp$week = as.numeric(as.factor(week(vp$date)))
    #vp = vp[vp$week==2 | vp$week==3, ]
    vp$week = as.numeric(as.factor((vp$week)))
  
  # Get weekdays from data
    vp$weekday = as.numeric(as.factor(weekdays(vp$date)))
    
# Habitat data  
  # Read in the habitat data    
    habitat = read.csv('vpoolsHabitat.csv')
    
  # Have a look at the first few rows of data
    head(habitat)
    
  # Get pond as a numerically represented factor
    habitat$pool = as.numeric(as.factor(habitat$pool))  
    names(habitat)[1] = 'pond'
    
  # Merge the habitat data with the trap counts
    vhabs = merge(vp, habitat)
    
  # Order the data as a check
    vhabs = vhabs[with(vhabs, order(date, species, pond)), ]
    #head(vhabs)  

# Create data for occupancy model(s)
# Get the total counts of observations for each species
  total.count = tapply(vp$occ, vp$species, sum)

# Find the number of unique species
  n.species = unique(vp$species)
  n.species

# Define a variable 'n' for the number of unique sampling dates
  n = length(n.species)
 
# Find the number of unique sampling locations
  n.ponds = unique(vp$pond)
  
# Number of sampling reps
  nreps = length(unique(vp$rep))
  
# J is the number of sampled ponds
  J = length(unique(vp$site))
  J

# The detection/non-detection data is reshaped into a three dimensional
# array X where the first dimension, j, is the point; the second
# dimension, k, is the rep; and the last dimension, i, is the length class.
  junk.melt=melt(vhabs, id.var=c("species", "site", "pond"), measure.var="occ")
  X = cast(junk.melt, site ~ pond ~ species, fun.aggregate = max,
           fill=NA, add.missing = TRUE)

# K is a vector of length J indicating the number of reps at each point j
  K = rep(ncol(X[,,1]), J)

# nsites is the number of sites
  nsites=length(unique(vp$site))  
  
# Basic model -------------------------------------------------------------
# Write the model code
  modelstring = "
    model {
  
    for (i in 1:n) {
    
      # Create priors for rep i
      u[i] ~ dnorm(0, 0.001) # Prior for occupancy on logit scale
      v[i] ~ dnorm(0, 0.001) # Prior for detection on logit scale
    
      # Derived estimates- invert link function for estimates on real scale
      occ[i] = exp(u[i])/(1+exp(u[i]))
      det[i] = exp(v[i])/(1+exp(v[i]))
    
      # Create a loop to estimate the Z matrix (true occurrence for species i
      # at point j.
      for (j in 1:J) {
  
        logit(psi[j,i]) <- u[i] # Can put gear in here
        mu.psi[j,i] <- psi[j,i]
        Z[j,i] ~ dbern(mu.psi[j,i])
  
        # Create a loop to estimate detection for species i at point k
        # during sampling period k.
        for (k in 1:K[j]) {
          logit(p[j,k,i]) <-  v[i] # Can put gear in here
          mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
          X[j,k,i] ~ dbern(mu.p[j,k,i])
        } # Close k
  
      } # Close j
    } # Close i
  
    }
  "
  
# Write the model to a file  
  writeLines(modelstring, 'vpOcc.txt')

# Specify data, initial values, and model settings
  # Load all the data
  vp.data = list(n=n, J=J, K=K, X=X)
  
  # Specify the parameters to be monitored
  vp.params = c('occ', 'det')

  # Specify initial values 
  vp.inits = function() {
    list(
      u=rnorm(n, 0, 1),
      v=rnorm(n, 0, 1),
      Z = structure(rep(1, J*n),.Dim=c(J,n))
    )
  }
  
  # Specify MCMC settings
    n.iter = 1100
    n.burnin = 100
    n.chains = 3
    n.thin = 10
    
# Run the model and call the results ?fit?
vpOcc = jags(data=vp.data, inits=vp.inits, 
             parameters.to.save=vp.params, model.file = "vpOcc.txt",
             n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin,
             working.directory = getwd())
  
# Print the model results
print(vpOcc, digits=3)

# Date-specific estimates -------------------------------------------------
# The detection/non-detection data is reshaped into a three dimensional
# array X where the first dimension, j, is the point; the second
# dimension, k, is the rep; and the last dimension, i, is the length class.
  junk.melt=melt(vhabs, id.var=c("species", "pond", "rep",'site'), measure.var="occ")
  X = cast(junk.melt, site ~ pond ~ species ~ rep, fun.aggregate = max,
           fill=NA, add.missing = TRUE)
  
# Write the model code
  modelstring = "
    model {
    for (g in 1:n){
      v[g] ~ dnorm(0, 0.001) # Prior for detection on logit scale
      det[g] <- exp(v[g])/(1+exp(v[g]))
    }

    for (i in 1:nreps) {
    
      # Create a loop to estimate the Z matrix (true occurrence for length i
      # at point j.
      for (g in 1:n){  

        # Create priors for length species i
          u[i,g] ~ dnorm(0, 0.001) # Prior for occupancy on logit scale

        # Derived estimates- invert link function for estimates on real scale
          occ[i,g] <- exp(u[i,g])/(1+exp(u[i,g]))

        for (j in 1:J) {

          # Create a loop to estimate detection for species i at point k
          # during sampling period k.
      
          for (k in 1:K[j]) {
  
          # Occupancy model
            logit(psi[j,k,g,i]) <- u[i,g] # Can put gear in here
            mu.psi[j,k,g,i] <- psi[j,k,g,i]
            Z[j,k,g,i] ~ dbern(mu.psi[j,k,g,i])
  
          # Detection model
            logit(p[j,k,g,i]) <-  v[g]
            mu.p[j,k,g,i] <- p[j,k,g,i]*Z[j,k,g,i]
            X[j,k,g,i] ~ dbern(mu.p[j,k,g,i])
          
          } # Close k
        } # Close j
      } # Close g
    } # Close i
  
    }
  "
  
# Write the model to a file  
  writeLines(modelstring, 'vpOccInd.txt')
  
# Specify data, initial values, and model settings
  # Load all the data
    vp.data = list(n=n, J=J, K=K, X=X,
                 nreps=nreps)
  
  # Specify the parameters to be monitored
    vp.params = c('occ', 'det')

  # Specify initial values 
    vp.inits = function() {
      list(
        u=structure(rnorm(n*nreps,0,1), dim=c(nreps,n)),
        v=rnorm(n, 0, 1),#structure(rnorm(n*nreps,0,1), dim=c(nreps,n)),
        Z = structure(rep(1,J*K[1]*nreps*n), dim = c(J,K[1],n,nreps))
      )
    }
  
  # Specify MCMC settings
    n.iter = 1100
    n.burnin = 100
    n.chains = 3
    n.thin = 3
    
# Run the model
vpOccInd = jags(data=vp.data, inits=vp.inits, 
             parameters.to.save=vp.params, model.file = "vpOccInd.txt",
             n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin,
             working.directory = getwd())
  
# Print the model results
print(vpOccInd, digits=3)

# Make a plot of the results

