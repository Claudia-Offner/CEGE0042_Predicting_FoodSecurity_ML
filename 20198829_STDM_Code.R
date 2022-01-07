# PREDICTING FOOD INSECURITY WITH RANDOM FORESTS
# Module: CEGE0042
# Student no.: 20198829
# Date: 26 April 2021

# This code was developed to analyze and predict food insecurity in Chad
# whilst accounting for space and time autocorrelations.

# To run this code, set the work directory to folder containing the provided data.
setwd("C:/Users/offne/Downloads/20198829_STDM_Data/")

#### import libraries ####
library(rgdal)
library(sp)
library(sf)
library(zoo)
library(reshape)
library(ggplot2)
library(lattice)
library(tmap)
library(spdep)
library(pdp)
library(dplyr) 
library(knitr)
library(ranger)
library(ggcorrplot) 
library(caret)

#### 1. READ DATA ####
df <- readOGR(dsn="chad_data.shp")
# set projection
df <- spTransform(df, CRS("+proj=sinu +ellps=WGS84 +datum=WGS84 +no_defs")) 



#### 2. CLEAN DATA ####

# Remove all predictors except fews_ipc
df <- subset (df, select = -c(country, fews_ha, fws_prj_n, fws_prj_n_, fws_prj_m, fws_prj_m_))
# Remove anom variables to optimize independence
df <- subset (df, select = -c(ndvi_nm, rain_nm, et_anom))

# Set column names
colnames(df@data)[1] <- "admin_code"
colnames(df@data)[2] <- "admin_name"
colnames(df@data)[5] <- "year_month"
colnames(df@data)[8] <- "fews_ipc"
colnames(df@data)[9] <- "ndvi_mean"
colnames(df@data)[10] <- "rain_mean"
colnames(df@data)[11] <- "et_mean"
colnames(df@data)[12] <- "acled_count"
colnames(df@data)[13] <- "acled_fatalities"
colnames(df@data)[14] <- "p_staple_food"
colnames(df@data)[16] <- "ruggedness_mean"
colnames(df@data)[17] <- "area_a"
colnames(df@data)[18] <- "pasture_pct"
colnames(df@data)[19] <- "cropland_pct"

# Remove unnecessary years
df <- df[df$year != 2007,]
df <- df[df$year != 2008,]
df <- df[df$year != 2009,]
df <- df[df$year != 2020,]

# Extend NA's where applicable
na_extension <- na.fill(df$fews_ipc, "extend")
df$fews_ipc <- na_extension
na_extension <- na.fill(df$et_mean, "extend")
df$et_mean <- na_extension

# Round Target Feature Classifications
df$fews_ipc <- round(df$fews_ipc)

# Replace strange characters in Admin Names
df$admin_name <- gsub('©', 'C', df$admin_name)

# order alphabetically
df <- df[order(df$admin_name), ]



#### ** Matrix Transformations ** ####

# Transform target variable to matrices for spatial temporal exploration 
# https://www.r-statistics.com/2012/01/aggregation-and-restructuring-data-from-r-in-action/

# Transpose few_ipc into time series
df_ms <- cast(df@data, admin_name~year_month, value = "fews_ipc") # in months
df_ys <-cast(df@data, admin_name~year, mean, value = "fews_ipc") # in years
s <- subset(df@data, select = c("admin_name", "admin_code", "centx", "centy"))
s <- s[!duplicated(s), ] # remove duplicates
df_ms <- merge(s, df_ms, by='admin_name') 
df_ys <- merge(s, df_ys, by='admin_name') 

# Transform into matrix
df_m_matrix <- data.matrix(df_ms[,5:ncol(df_ms)])
df_y_matrix <- data.matrix(df_ys[,5:ncol(df_ys)])
rownames(df_m_matrix) <- df_ms[,1]
rownames(df_y_matrix) <- df_ys[,1]



#### 3. NON SPATIAL-TEMPORAL ANALYSIS ####

# Look at data
colnames(df@data) 
head(df@data)
summary(df)
for (i in colnames(df@data)){
  print(paste(i, sd(df[[i]])))
}

# Examine predictive feature (few_ipc)
par(fig=c(0,1,0,0.8), new=TRUE)
barplot(prop.table(table(df$fews_ipc)), ylim=c(0,1), xlab = "Fews_ipc score", ylab = "Frecuency (%)", col="#6699FF")
par(fig=c(0,1,0.56,1), new=TRUE)
boxplot(df$fews_ipc,
        col = "#FF5050",
        border = "black",
        horizontal = TRUE,
        notch = TRUE,
        axes = FALSE,
        main = "Enhanced Fews_ipc Barplot")

# Examine correlation between features and target variable
ggcorrplot(cor(df@data[,8:19]), legend.title = "Correlation", colors = c("#ff5050", "#FFFFFF", "#6699FF"))+
  labs(title = "Correlation Matrix of Model Parameters")



#### 3. SPATIAL ANALYSIS ####

# Spatial distribution over annual average
bound <- readOGR(dsn="chad.shp")
bound <- spTransform(bound, CRS("+proj=sinu +ellps=WGS84 +datum=WGS84 +no_defs"))
df2 <- merge(bound, df_ys, by='admin_code') 
df2 <- subset(df2, select=-c(centx, centy, admin_name))
d <- melt(df2@data,
            id.vars = "admin_code",
            variable.name = "Year",
            value.name = "Value")
df2@data <-d
df2$level <- "Low"
df2$level[which(d$value >= 2)] <- "High"

tm_shape(df2) +
  tm_polygons("level", title="Food Insecurity Level", palette=c("#ff5050", "#ffe5e5")) +
  tm_compass(position=c("left","top"))+
  tm_scale_bar()+
  tm_facets("variable")

# QUANTIFY SPATIAL DEPENDENCIES

# Re-order admins
data <- data.frame(admin_code = df@data$admin_code, admin_name = df@data$admin_name)
data <- data[!duplicated(data), ]
bound <- merge(bound, data, by='admin_code') 
bound <- bound[order(bound$admin_name), ]
rownames(bound@data) <- bound@data$admin_name

# Check Spatial Dependency
# Moran's Global I
W <- nb2listw(poly2nb(bound))
spaceMean <- rowMeans(df_m_matrix)
moran.test(x=spaceMean, listw=W)

# Moran's Local I Clusters
moranCluster <- function(shape, W, var, alpha=0.05, p.adjust.method="bonferroni")
{
  Ii <- localmoran(var, W, p.adjust.method=p.adjust.method)
  shape$Ii <- Ii[,"Ii"]
  shape$Iip <- Ii[,"Pr(z > 0)"]
  shape$sig <- shape$Iip<alpha
  # Scale the data to obtain low and high values
  shape$scaled <- scale(var) # high low values at location i
  shape$lag_scaled <- lag.listw(W, shape$scaled) # high low values at neighbours j
  shape$lag_cat <- factor(ifelse(shape$scaled>0 & shape$lag_scaled>0, "High-High",
                                 ifelse(shape$scaled>0 & shape$lag_scaled<0, "High-Low",
                                        ifelse(shape$scaled<0 & shape$lag_scaled<0, "Low-Low",
                                               ifelse(shape$scaled<0 & shape$lag_scaled<0, "Low-High", "Equivalent")))))
  shape$sig_cluster <- as.character(shape$lag_cat)
  shape$sig_cluster[!shape$sig] <- "Non-significant"
  shape$sig_cluster <- as.factor(shape$sig_cluster)
  results <- data.frame(Ii=shape$Ii, pvalue=shape$Iip, type=shape$lag_cat, sig=shape$sig_cluster)
  
  return(list(results=results))
}

clusters <- moranCluster(bound, W=W, var=spaceMean)$results
bound$Ii_cluster <- clusters$sig
tm_shape(bound) + tm_polygons(col="Ii_cluster", palette=c("#ff5050", "#6699FF"), title = "Moran's Local I Clusters")
d <- subset(bound@data, select=-c(3))
bound@data <- d

# So there is a strong spatial dependency



#### 3. TEMPORAL ANALYSIS ####

# Look at the mean few_ipc across country over time (in years)
d<- as.data.frame(colMeans(df_y_matrix))
ggplot(data = d, aes(y=colMeans(df_y_matrix), x = seq(2010, 2019, by=1)))+
  geom_line(color = "#6699FF", size = 2)+
  scale_x_continuous(breaks=seq(2010, 2019, by=1))+
  labs(x = "Year", y = "Fews_ipc Score", title = "Yearly Time Series of Food Insecurity")

# Look at the mean few_ipc across country over time (in months)
d<- as.data.frame(colMeans(df_m_matrix))
a <- 1:120
ggplot(data = d, aes(y=colMeans(df_m_matrix), x = a))+
  geom_line(color = "#ff5050", size = 2)+
  scale_x_continuous(breaks=a[seq(1, length(a),12)], labels = c("2010_01", "2011_01", "2012_01", "2013_01", "2014_01","2015_01", "2016_01", "2017_01", "2018_01", "2019_01"))+
  labs(x = "Month", y = "Fews_ipc Score", title = "Monthly Time Series of Food Insecurity")

# Order heat map by admin food security means to examine effects
ggplot(df@data) +
  geom_tile(aes(x=year_month, y=reorder(admin_name, centy, mean, order=TRUE), fill = fews_ipc)) +
  scale_fill_gradient(name = "Classification",low = "blue", high = "red") +
  labs(x = "Month", y = "Admin")

# QUANTIFY TEMPORAL DEPENDENCIES

# Check Temporal Dependency
timeMean <- colMeans(df_m_matrix)

bacf <- acf(timeMean, lag.max=50, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=0.20, linetype='dashed', col = 'blue')+
  geom_hline(yintercept=-0.20, linetype='dashed', col = 'blue')+
  scale_y_continuous(breaks=seq(-1.0, 1.0, by=0.2))+
  labs(y="ACF", x="Lag", title="ACF")

bacf <- acf(diff(timeMean), lag.max=50, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=0.20, linetype='dashed', col = 'blue')+
  geom_hline(yintercept=-0.20, linetype='dashed', col = 'blue')+
  scale_y_continuous(breaks=seq(-1.0, 1.0, by=0.2))+
  labs(y="ACF", x="Lag", title="ACF Differenced")

bacf <- pacf(timeMean, lag.max=50, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=0.20, linetype='dashed', col = 'blue')+
  geom_hline(yintercept=-0.20, linetype='dashed', col = 'blue')+
  scale_y_continuous(breaks=seq(-1.0, 1.0, by=0.2))+
  labs(y="PACF", x="Lag", title="PACF")

bacf <- pacf(diff(timeMean), lag.max=50, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=0.20, linetype='dashed', col = 'blue')+
  geom_hline(yintercept=-0.20, linetype='dashed', col = 'blue')+
  scale_y_continuous(breaks=seq(-1.0, 1.0, by=0.2))+
  labs(y="PACF", x="Lag", title="PACF Differenced")

#### 4. MODEL PRE-PROCESSING ####

# 1. TIME LAGS
# Add lag of 6 for general time trends
# Add lag of 12 for seasonal time trends
# https://statisticsglobe.com/create-lagged-variable-by-group-in-r

df2 <- data.frame(df@data)
df2 <- df2 %>%
  group_by(admin_name) %>%
  dplyr::mutate(time1 = dplyr::lag(fews_ipc, n = 6, default = NA)) %>% # general time  lag of 6 by admin_code
  dplyr::mutate(time2 = dplyr::lag(fews_ipc, n = 12, default = NA)) %>% # seasonal time lag of 12 by admin_code
  as.data.frame()

# 2. SPACE LAG
## https://wlm.userweb.mwn.de/R/wlmRspma.htm

# A. Spatial weight matrix
# Queen contiguity matrix (we consider observations that share a vertex to be neighbors)
W <- nb2listw(poly2nb(bound), style="B") # binary weight order (nb2list)
sp_w_matrix <- listw2mat(W) #store as matrix
# Name cols rows
rownames(sp_w_matrix) <- bound@data$admin_name
colnames(sp_w_matrix) <- bound@data$admin_name

#B. Extract list of neighbors for each admin_name
nest_list <- vector("list", nrow(sp_w_matrix))
for(i in 1:nrow(sp_w_matrix)){
  row <- sp_w_matrix[i, ] == 1 # extract values equal to 1
  matched_col <- colnames(sp_w_matrix)[row == T] # extract col names of the values
  nest_list[[i]] <- matched_col # append them to list
}
names(nest_list) <- rownames(sp_w_matrix)

#C. For each list of neighbors (nest_list), at every time interval (df_m_matrix), take the neighbors average of fews_ipc 
nest_av <- vector("list", nrow(df_m_matrix))
names(nest_av) <- rownames(df_m_matrix)
for (list in nest_list) {
  
  timeseries <- c()
  
  for(month in 1:ncol(df_m_matrix)){
    
    num <- c()
    
    for(neighbor in list) {
      
      for(admin in 1:nrow(df_m_matrix)) {
        
        if(rownames(df_m_matrix)[admin] == neighbor) {
          
          num <- c(num, df_m_matrix[admin, month])
        }
      }
    }
    av <- mean(num)
    
    timeseries <- c(timeseries, av)
    
  }
  nest_av[list] <- list(timeseries)
}

# D. Store nested averages in final data frame
sp_w_var <- as.data.frame(nest_av)
rownames(sp_w_var) <- colnames(df_m_matrix)
sp_w_var <- t(sp_w_var)
sp_w_var <- melt(sp_w_var)
sp_w_var <- sp_w_var[order(sp_w_var$X1), ]
df2$spacelag <- sp_w_var$value

# Normalize features
process <- preProcess(as.data.frame(df2[,9:22]), method=c("range"))
df2[,9:22] <- predict(process, as.data.frame(df2[,9:22]))

# Set fews_ipc to binary classification - splitting on third quartile
y <- df2$fews_ipc
yClass <- y
yClass[which(y>=2)] <- 1 # high insecurity
yClass[which(y<2)] <- -1 # low insecurity
df2$fews_ipc <- as.factor(yClass)

# Ommit NA values
df2 <- df2 %>% na.omit()

# Remove unnecessary features
df2 <- subset(df2, select = -c(1:7))



#### 5. RANDOM FOREST CLASSIFICATION ####

# A. SPLIT TRAIN/TEST SETS 

# Random 70% split in train and test sets (time series is accounted for by lags)
set.seed(101)
n <- nrow(df2)
trainInd <- sort(sample(1:n, n*.7)) 
train <- df2[trainInd,]
test <- df2[-trainInd,]

# B. OPTIMISE MODEL HYPERPARAMETERS
# https://uc-r.github.io/random_forests 

# Hyperparameter grid search (112)
features <- train[,2:15]

hyper_grid <- expand.grid(
  mtry       = seq(1, 14, by = 1),
  node_size  = seq(3, 9, by = 1),
  sampe_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0
)

nrow(hyper_grid)

##!! Commented out for computational efficiency !!##
#for(i in 1:nrow(hyper_grid)) {

# train model
#  model <- ranger(
#    formula         = fews_ipc ~ ., 
#    data            = train, 
#    num.trees       = 500,
#    mtry            = hyper_grid$mtry[i],
#    min.node.size   = hyper_grid$node_size[i],
#    sample.fraction = hyper_grid$sampe_size[i],
#    seed            = 123
#  )

# add OOB error to grid
#  hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
#}

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

# Optimized hyper parameters
# mtry av = 10
# node_size av = 4
# sampe_size = 0.8

# C. FIT OPTIMISED MODEL
model <- ranger(
  formula   = fews_ipc ~ ., 
  data      = train, 
  num.trees = 500,
  mtry      = 10,
  min.node.size = 4,
  sample.fraction = 0.800,
  classification = TRUE,
  importance = 'impurity'
)

# D. RESULTS

# Evaluate model
1 - model$prediction.error # accuracy

# Plot Confusion matrix
cm <- model$confusion.matrix
ggplot(as.data.frame(cm), aes(y=reorder(true, Freq), x=predicted)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low="#6699FF", high="#ff5050")+
  labs(x="Predicted", y="True", fill="Frequency", title="Model Confusion Matrix")

# Plot Variable importance
imp <- as.data.frame(model$variable.importance)
colnames(imp) <- c("importance")
imp$variable <- rownames(imp)

ggplot(imp, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  labs(x = "", y="RMSE", title= "Variable Importance Scores")+
  guides(fill=F)+
  scale_fill_gradient(low="#ff5050", high="#6699FF")

# Make prediction
predicted <- predict(model, test)

# Calculate accuracy score
CM <- table(predicted[["predictions"]], test$fews_ipc)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy


