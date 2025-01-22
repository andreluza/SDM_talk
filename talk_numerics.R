
# -------------------------------------------------------------------------------
# Tasks:
# 1 - estimate the number of occupied sites (\psi) using different models (GAM, GLM, SOM)
# 2 - estimate occupancy probability per site (\psi_i) using different models (SRE, GAM, GLM, SOM)
# 3 - estimate the effect of covariates on \psi_i using different models (SRE, GAM, GLM, SOM)

# To do so we will simulate a data set which is subjected to imperfect detection, as usual in any ecological data set. 
# When needed (for GLMs, GAM) we will aggregate the occupancy data to single-visit data
# These simulations are based on the books of Kéry & Schaub 2012 & Kéry & Royle 2016 AHM I


# -------------------------------------------------------------------------------

# If you do not have the packages installed in your R, run for example for one package:
install.packages("raster")
install.packages("unmarked")
install.packages("mgcv")
install.packages("jagsUI")

# Load necessary libraries
rm(list=ls()) # clear the workspace
library(raster)
library(unmarked)
library(mgcv)
library(jagsUI)

####___________ Task 1 _______________####

#### 1) Detection-non-detection data simulation (without covariates) ####

# Choose sample sizes and prepare observed data array y_ij
set.seed(24) # Set the seed so we all get same data set and results
M <- 100 # Number of sites
J <- 3 # Number of presence/absence measurements (n surveys)
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
psi <- 0.8 # Probability of occupancy/proportion of occupied sites
p <- 0.5 # Probability of detection in occupied sites

# ----- Ecological model
# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi) # R has no Bernoulli, so Binomial with size=1 do the trick

# ----- Observation model
# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = z, prob = p) # detection conditional on true occurrence
}

# ----- 2) Which model can recover the true values of \psi and p?
# Try the GLM
# apply (y,1,max) will aggregate data to a single-visit data. Any 1 will result in a site with detection
m1 <- glm (y ~1 ,
           data = data.frame (y = apply (y,1,max)),
           family="binomial")
summary(m1)
plogis(m1$coefficients) # backtransform (inverse logit function), it's the same as 1/(1+exp(beta0))

# Try the GAM
# apply (y,1,max) will aggregate data so that any 1 will result in a site with detection
m1_gam <- gam (y ~1 ,
           data = data.frame (y = apply (y,1,max)),
           family="binomial")
summary(m1_gam)
plogis(m1_gam$coefficients) # backtransform (inverse logit function)

# Try the hierarchical site-occupancy model
# hierarchical model that enables to estimate psi and p

umf <- unmarkedFrameOccu(y = y) # Create unmarked data frame
summary(umf) # Summarize data frame
(fm1 <- occu(~1 ~1, data = umf)) # Fit model

# Get estimates on probability scale
backTransform(fm1, "state")
backTransform(fm1, "det")

# table of comparison  (sum(hat(psi)))
data.frame (True =  psi, 
            GLM = plogis(m1$coefficients), # backtransform (inverse logit function)
            GAM = plogis(m1_gam$coefficients) , # backtransform (inverse logit function)
            SOM = backTransform(fm1, "state")@estimate) # backtransform (inverse logit function)


# What's the conclusion here?


####___________ Task 2 _______________####


####  Simulate more complex data influence by site/observation covariates ####

# same M and J as before
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Create a covariate called vegHt
vegHt <- array(runif(M * J, -1, 1), dim = c(M, J))

# Create a covariate called edgeDist
edgeDist <- sort(runif(M, -1, 1)) # Create a gradient

# Define spatial grid (for visualization)
n <- sqrt(M) # Define grid dimensions (assuming M is a perfect square)
coords <- expand.grid(x = 1:n, y = 1:n) # create all combinations of coordinates
coords$x <- coords$x/max(coords$x) # standardize coordinates
coords$y <- coords$y/max(coords$y) # standardize coordinates

# Create raster layers for predictors
edgeDist_raster <- rasterFromXYZ(cbind(coords, edgeDist))
vegHt_raster <- rasterFromXYZ(cbind(coords, apply(vegHt, 1, mean)))

# Plot predictors to see the environment in the simulated 
par(mfrow=c(2,3))
plot(edgeDist_raster, main = "Edge distance")
plot(vegHt_raster, main = "Average vegHt")
plot(0)
plot(rasterFromXYZ(cbind(coords, vegHt[,1])),main = "Average vegHt,j=1")
plot(rasterFromXYZ(cbind(coords, vegHt[,2])),main = "Average vegHt,j=2")
plot(rasterFromXYZ(cbind(coords, vegHt[,3])),main = "Average vegHt,j=3")

# Choose true parameter values for occupancy model and compute occupancy
beta0 <- 0 # Logit-scale intercept
beta1 <- 2 # Logit-scale slope for edgeDist (positive effect)
psi <- plogis(beta0 + beta1 * edgeDist) # True occupancy probability
plogis(beta0) # inverse logit, the same as 1/(1+exp(beta0))

# plot the true relationships
par(mfrow=c(2,3),mar=c(4,4,4,4))
plot(edgeDist, psi, 
     ylim = c(0,1), type = "p", pch = 19,ylab=expression(psi[i])) # Plot psi relationship

# Rasterize occupancy probability (psi)
psi_raster <- rasterFromXYZ(cbind(coords, psi))
plot(psi_raster, main = expression(paste("True Occupancy Probability ",psi[i],""))) # plot the true occupancy

# ----- Ecological model
# Now visit each site and observe presence/absence perfectly
z <- rbinom(M, 1, psi) # True presence/absence

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2 # Logit-scale intercept
alpha1 <- -3 # Logit-scale slope for vegHt - detection probability decreases with vegHt
p <- plogis(alpha0 + alpha1 * vegHt) # Detection probability
mean(p)
plot(p ~ vegHt, ylim = c(0,1),pch=19,ylab=expression(p[ij])) # Look at relationship

# ----- Observation model
# Take J measurements at each site
for(j in 1:J) {
  y[,j] <- rbinom(M, z, p[,j]) # sampling conditional on true occurrence
}

# Rasterize observation-level covariates (varying over I sites and J surveys)
y_raster_1 <- rasterFromXYZ(cbind(coords, y[,1]))
y_raster_2 <- rasterFromXYZ(cbind(coords, y[,2]))
y_raster_3 <- rasterFromXYZ(cbind(coords, y[,3]))
plot(y_raster_1, main = expression(paste("Occupancy data ",y[i],"")))
plot(y_raster_2, main = expression(paste("Occupancy data ",y[i],"")))
plot(y_raster_3, main = expression(paste("Occupancy data ",y[i],"")))


# ----- Which model can recover the true values of \psi_i?

# --- 1. Try the Species Range Envelope Model ---
# Presence-only model based on threshold of environmental conditions better suited for the species

# Define environmental thresholds for SRE
edgeDist_threshold <- range(edgeDist)   * c(0.2, 0.8) # Arbitrary thresholds
vegHt_threshold <- range(vegHt)   * c(0.2, 0.8)  # Arbitrary thresholds

# Apply SRE logic
sre_prediction <- ifelse(edgeDist > edgeDist_threshold[1] & edgeDist < edgeDist_threshold[2] &
                           apply(vegHt, 1, mean) > vegHt_threshold[1] &
                           apply(vegHt, 1, mean) < vegHt_threshold[2], 1, 0)

# Create raster for SRE predictions
sre_raster <- rasterFromXYZ(cbind(coords, sre_prediction))
par(mfrow=c(4,3),mar=c(5,4,3,6))
plot(sre_raster, main = "SRE Predicted Distribution")
plot(psi_raster, main = "True Distribution")
plot(psi_raster - sre_raster, main = "Difference")

# --- 2. Try the GLM Model ---
glm_model <- glm (y ~vegHt+edgeDist ,
                  data = data.frame (y = apply (y,1,max),
                                     vegHt = rowMeans(vegHt),
                                     edgeDist = edgeDist),family="binomial")
summary(glm_model)
plogis(glm_model$coefficients[1])

# predictions GLM
predictions_glm <- plogis(glm_model$coefficients [1] + glm_model$coefficients[2] * rowMeans(vegHt) +   glm_model$coefficients [3] * edgeDist)

# Create and save raster for GLM predictions
glm_raster <- rasterFromXYZ(cbind(coords, predictions_glm))
plot(glm_raster, main = "GLM Predicted Distribution")
plot(psi_raster, main = "True Distribution")
plot(psi_raster - glm_raster, main = "Difference")

# --- 3. Try the GAM Model ---
gam_model <- gam (y ~s(vegHt)+s(edgeDist) ,
                  data = data.frame (y = apply (y,1,max),
                                     vegHt = rowMeans(vegHt),
                                     edgeDist = edgeDist),family="binomial")
summary(gam_model)
plogis(gam_model$coefficients[1])

# predictions GAM
predictions_gam <- plogis(predict(gam_model))

# Create and save raster for GAM predictions
gam_raster <- rasterFromXYZ(cbind(coords,predictions_gam))
plot(gam_raster, main = "GAM Predicted Distribution")
plot(psi_raster, main = "True Distribution")
plot(psi_raster - gam_raster, main = "Difference")

# --- 4. Try the Hierarchical Occupancy Model ---
# Format data and summarize
umf <- unmarkedFrameOccu(y = y, # Pres/Abs measurements
                         siteCovs = data.frame(edgeDist = edgeDist), # site-specific covs.
                         obsCovs = list(vegHt = vegHt)) # obs-specific covs.
summary(umf)

# Fit model and extract estimates
# Detection covariates follow first tilde, then occupancy covariates
summary(fm2 <- occu(~vegHt ~edgeDist, data=umf))
predictions_SOM <- predict (fm2,type = "state")$Predicted

# Create and save raster for SOM predictions
SOM_raster <- rasterFromXYZ(cbind(coords, predictions_SOM))
plot(SOM_raster, main = "SOM Predicted Distribution")
plot(psi_raster, main = "True Distribution")
plot(psi_raster - SOM_raster, main = "Difference")

# table of comparison of all models (sum(hat(psi_i)))
data.frame (True =  sum(values(psi_raster))/M, 
            SRE = sum(values(sre_raster))/M,
            GLM = sum(values(glm_raster))/M,
            GAM = sum(values(gam_raster))/M,
            SOM = sum(values(SOM_raster))/M)

# What's the conclusion here?


# ----- Comparison of covariate effects -  regression coefficients
plot.gam(gam_model,se = T,pages = 1,add=T) # see estimates (linear scale)

par(mfrow=c(3,2),mar=c(5,4,4,6))
plot(p ~ vegHt, main="True", ylim = c(0,1),pch=19,ylab=expression(p[ij])) # Look at relationship
plot(edgeDist, psi,
     ylim = c(0,1), type = "p", pch = 19,ylab=expression(psi[i])) # Plot psi relationship

# GLM
plot(seq(from = -1, to = 1,length.out=M)  , 
     plogis(glm_model$coefficients[1] + 
            glm_model$coefficients[2] * seq(from = -1, to = 1,length.out=M) +
            glm_model$coefficients[3] * 0),
     ylim = c(0,1),
     main="GLM", 
     ylab=expression(paste (psi[i],".", p)),
     xlab="")

plot(seq(from = -1, to = 1,length.out=M)  , 
     plogis(glm_model$coefficients[1] + 
            glm_model$coefficients[2] * 0+
            glm_model$coefficients[3] *  seq(from = min(edgeDist), to = max(edgeDist),length.out=M) ),
     ylim = c(0,1),
     ylab=expression(paste (psi[i],".", p)),
     xlab="")

# SOM
# detection
plot(seq(from = -1, to = 1,length.out=M)  , 
     plogis(fm2@estimates@estimates$det@estimates[1]+ 
              fm2@estimates@estimates$det@estimates[2] * seq(from = -1, to = 1,length.out=M)) ,
     ylim = c(0,1),
     ylab=expression(paste (p[ij])),
     main="SOM", 
     xlab="")

# state
plot(seq(from = -1, to = 1,length.out=M)  , 
     plogis(fm2@estimates@estimates$state@estimates[1]+ 
              fm2@estimates@estimates$state@estimates[2] * seq(from = -1, to = 1,length.out=M)) ,
     ylim = c(0,1),
     ylab=expression(paste (psi[i])),
     xlab="")


# What's the conclusion here?

# end
