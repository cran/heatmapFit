require(ROCR)
set.seed(459871)

# generate data set
x<-runif(100, min=0, max=10)
y<-ifelse(runif(100, min=0, max=1)<pnorm(0.15*x-2),1,0)
dat2<-as.data.frame(cbind(y,x))

# run a correctly specified model on the data
probit.model<-glm(y~x, family=binomial(link=probit), data=dat2)
probit.predict.pr<-predict(probit.model, type="response", newdata=data.frame(x=seq(from=0, to=10, by=0.1)))

# make a plot of the data generating process and data set
par(mfrow=c(1,2))
plot(y~x, main="Probit Model")
lines(probit.predict.pr~seq(from=0, to=10, by=0.1), lty=1)
x.true<-seq(from=0, to=10, by=0.1)
y.true<-pnorm(0.15*x.true-2)
lines(y.true~x.true, lty=2)
legend("right", pch=c(1,-1, -1), lty=c(0,1,2), legend=c("data", "predicted Pr(y=1)", "true Pr(y=1)"))
summary(probit.model)

# plot the ROC curve for that model
probit.predict<-predict(probit.model)
pred <- prediction(probit.predict, dat2$y)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, col=rainbow(10), main="ROC Curve")
abline(0, 1, lty=2)
legend("bottomright", lty=c(1,2), col=c(2,1), legend=c("ROC curve", "random expectation"))

# area under the ROC curve for the model
performance(pred, measure = "auc")

# determine the percent correctly predicted
probit.predict.pr<-predict(probit.model, type="response")
pcp.opt<-0
t.opt<-0
for(j in seq(from=0, to=1, by=0.05)){

	y.pred<-ifelse(probit.predict.pr>j, 1, 0)
	pcp<-sum(y.pred==y)/length(y)
	if(pcp>pcp.opt){pcp.opt<-pcp; t.opt<-j}	

}
cat("The optimum percent correctly predicted is ", pcp.opt, " with a threshold of ", t.opt, ".\n", sep="")

# determine the ePCP
epcp<-sum(ifelse(y==1, (1/(length(y)))*probit.predict.pr, (1/(length(y)))*(1-probit.predict.pr)))
cat("The expected percent correctly predicted is ", epcp, ".\n", sep="")

# create a heatmap fit plot and calculate the heatmap statistic
out<-heatmap.fit(y~x, fam=binomial(link=probit), dat=dat2, reps=1000, span.l="aicc")


