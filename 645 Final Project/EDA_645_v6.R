directory <- "c:/users/beaus/Desktop/STAT 645 Project"
setwd(directory)

require("aTSA")
require("splines")

sst_df <- read.csv("sea-surface-temp_fig-1_epa.csv", skip=7, header = FALSE)
ace_df <- read.csv("cyclones_fig-2_epa.csv", skip=7, header=FALSE)
tna_df <- read.table(text=paste0(head(readLines("https://www.esrl.noaa.gov/psd/data/correlation/tna.data"), -2), collapse="\n"), skip = 1, header = FALSE)

tna_df <- head(tna_df, -1)
tna_df <- tail(tna_df, -2)
tna_df$SST <- rowMeans(tna_df[,-1])
tail(tna_df)

kaggle_df <- read.csv("atlantic_kaggle.csv")
plot(kaggle_df$Date, kaggle_df$Maximum.Wind)
plot(kaggle_df$Date, kaggle_df$Minimum.Pressure, col = rgb(0,1,0, 0.2), ylim = c(850,1025))
kaggle_df$Year <- as.numeric(substr(kaggle_df$Date, start=1, stop=4))

Median.Pressure <- tapply(kaggle_df$Minimum.Pressure, kaggle_df$Year, FUN=median)
Maximum.Wind <- tapply(kaggle_df$Maximum.Wind, kaggle_df$Year, FUN=median)

Median.Pressure1 <- tail(Median.Pressure, -29)
Maximum.Wind1 <- tail(Maximum.Wind, -29)
data <- cbind(Median.Pressure1, Maximum.Wind1, sst_df)

sst_df2 <- tail(sst_df, -70)
#newData <- data.frame(sst_df2$V1, sst_df2$V4, ace_df$V2)
#colnames(newData) <- c("Year", "SST", "ACE")

df <- data.frame(tna_df$V1, tna_df$SST, ace_df$V2)
colnames(df) <- c("Year", "SST", "ACE")

######################

y = df$ACE
x = df$SST

# Exponential fit.
f <- function(x,a,b) {a * exp(b * x)}
fit <- nls(y ~ f(x,a,b), start = c(a=1, b=1)) 
co <- coef(fit)
plot(df$SST, df$ACE, main="ACE vs. North Atlantic SSTs", xlab="SST", ylab="ACE", pch=16)
legend("bottomright", legend=c("linear", "exponential", "polynomial"), pch=16, col = c("red", "green", "blue"))
curve(f(x, a=co[1], b=co[2]), add = TRUE, col="green", lwd=2) 

# Linear fit.
linModel <- lm(y ~ x, df)
abline(linModel, col="red", lwd=2)
summary(linModel)

# Polynomial fit.
f <- function(x,a,b,c,d) {(a*x^3) + (b*x^2) + (c*x) + d}
fit2 <- nls(y ~ f(x,a,b,c,d), start = c(a=1, b=1, c=1, d=1)) 
co <- coef(fit2)
curve(f(x, a=co[1], b=co[2], c=co[3], d=co[4]), add = TRUE, col="blue", lwd=2) 

# Correlation test.
cor.test(x, y)

##normaility assumption
qqnorm(linModel$residuals)
qqline(linModel$residuals, col="red", lwd=3)

##Equal variance
plot(linModel$fitted.values, linModel$residuals)
abline(h=0, col="red" )

summary(linModel)

#######################

# Do a stationarity test.
stationary.test(df$ACE)
plot(df$Year, df$ACE, type="l", xlab="Year", ylab="ACE", main="ACE Time Series")
acf(df$ACE, type="correlation", plot=TRUE, main="ACF of ACE")
ACE_diff <- diff(as.matrix(ace_df$V2))
ACE_diff <- c(NA, ACE_diff)
df$ACE_diff <- ACE_diff
head(df)

ACE_lag <- c(NA, ace_df$V2[-nrow(ace_df)])
df$ACE_lag <- ACE_lag
head(df)

# Analyze an lm model.
ACE_model <- lm(ACE~ACE_lag, df)
summary(ACE_model)
cor.test(df$ACE, df$ACE_lag)

# Analyze another.
diff_model <- lm(ACE~ACE_diff, df)
summary(diff_model)
cor.test(df$ACE, df$ACE_diff)

# Spline model.
splineModel <- lm(ACE~bs(Year, knots=c(1967, 1987, 1997)), data=df)
plot(splineModel)
summary(splineModel)

logACE <- log(df$ACE)
df$logACE <- logACE

#3 cutpoints at ages 25 ,50 ,60
fit<-lm(logACE ~ bs(Year,degree=1),data = df)
summary(fit)

#Plotting the Regression Line to the scatterplot   
#plot(df$Year,df$logACE,col="grey",xlab="Year",ylab="logACE")
#points(year.grid,predict(fit,newdata = list(Year=year.grid)),col="darkgreen",lwd=2,type="l")
#adding cutpoints
#abline(v=c(1960,1980,2000),lty=2,col="darkgreen")

# Do a least-squares fit. Log the ACE values and store them in df.
lsModel <- lm(logACE~Year, data=df)
plot(df$Year, df$logACE)
abline(lsModel)

summary(lsModel)

########################################

# MLE calculations.
fn <- function(theta) {
  sum ( 0.5*(df$ACE - theta[1])^2/theta[2] + 0.5* log(theta[2]) )
}

optim(theta <- c(0, 1), fn)

########################################

#step function analysis

df$YearsM10Plus <- df$Year - 1959
df$YearsM10Plus[df$YearsM10Plus < 0] <- 0
df$YearsM20Plus <- df$Year - 1969
df$YearsM20Plus[df$YearsM20Plus < 0] <- 0
df$YearsM30Plus <- df$Year - 1979
df$YearsM30Plus[df$YearsM30Plus < 0] <- 0
df$YearsM40Plus <- df$Year - 1989
df$YearsM40Plus[df$YearsM40Plus < 0] <- 0
df$YearsM50Plus <- df$Year - 1999
df$YearsM50Plus[df$YearsM50Plus < 0] <- 0
df$YearsM60Plus <- df$Year - 2009
df$YearsM60Plus[df$YearsM60Plus < 0] <- 0

##########################################

##choosing knots using variable selection

# Model 1.
splineModel <- lm(logACE~Year+YearsM10Plus+YearsM20Plus+YearsM30Plus+YearsM40Plus+YearsM50Plus+YearsM60Plus,
                  data=df)
summary(splineModel)
BIC(splineModel)
splineModel1 <- lm(logACE~Year+YearsM20Plus+YearsM30Plus+YearsM40Plus+YearsM50Plus+YearsM60Plus,
                  data=df)
BIC(splineModel1)
summary(splineModel1)

# Model 2.
splineModel2 <- lm(logACE~Year+YearsM30Plus+YearsM40Plus+YearsM50Plus+YearsM60Plus,
                   data=df)
BIC(splineModel2)
summary(splineModel2)

# Model 3. 
splineModel3 <- lm(logACE~Year+YearsM30Plus+YearsM40Plus+YearsM50Plus,
                   data=df)
BIC(splineModel3)
summary(splineModel3)

# Model 4.
splineModel3 <- lm(logACE~Year+YearsM40Plus+YearsM50Plus,
                   data=df)
BIC(splineModel3)
summary(splineModel3)

plot(splineModel3)
df$predict <- predict(splineModel3)
preds <- predict(splineModel3, se.fit = TRUE)
df$predict2 <- preds$se.fit
df$lwr <- preds$fit - 1.96*preds$se.fit
df$upr <- preds$fit + 1.96*preds$se.fit
summary(splineModel3)
BIC(splineModel3)

# xtable
xtable::xtable(summary(splineModel3))


##create X on p 359
model.matrix(splineModel3)

yearlims <- range(df$Year)
year.grid <- seq(from=yearlims[1], to=yearlims[2])
plot(df$Year, df$logACE, xlab="Year", ylab="log(ACE)", 
     main="Linear model with spline", pch=20)
points(year.grid, predict(splineModel3), type="l", lwd=2, col="blue")
abline(v=1989, lty=2)
abline(v=1999, lty=2)

plot(splineModel3, pch=19)

# Graph linear model with spline in ggplot.
polydata <- data.frame(c(df$Year, rev(df$Year)), c(df$lwr, rev(df$upr)), c(df$SST, df$SST))
names(polydata) <- c("x", "y", "SST")
ggplot(data=df, aes(x=Year, y=log(ACE), color=SST)) +
  geom_point() +
  geom_line(aes(x=Year, y=predict), lwd=2) +
  geom_polygon(data=polydata, aes(x=x, y=y), color=NA, fill="blue", alpha=0.2) +
  ggtitle("Linear Model with Splines and 95% Confidence Band") +
  theme(plot.title = element_text(hjust=0.5), text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(1950, 2010, 10)) +
  scale_y_continuous(breaks=seq(3, 6, 0.5))

# Plot a spline against log(ACE).
plot(df$Year, df$YearsM30Plus, type="l", xlab="Year", ylab=expression('(x-30)'['+']), main="Graphical depiction of the spline")

########################################

plot(predict(splineModel3, newdata = data.frame(Year=2017:2020, 
                                           YearsM30Plus=37:40)), type="l")

