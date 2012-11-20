
\name{heatmapFit}
\alias{heatmapFit}
\alias{heatmap.fit}
\title{Heatmap Fit Statistic For Binary Dependent Variable Models}
\description{Generates a fit plot for diagnosing misspecification in GLM models, and calculates the related heatmap fit statistic (Esarey and Pierce, 2012).}
\usage{heatmap.fit(form, fam, dat, reps=1000, span.l="aicc",color=F)}
\arguments{
 \item{form}{Formula for the GLM model to be estimated.}
 \item{fam}{A family argument for the GLM model. A link should be specified 
	(e.g., fam = family( binomial( link = probit ))).}
 \item{dat}{Data set on which the GLM model is to be estimated.}
 \item{reps}{Number of repetitions (default = 1000).}
 \item{span.l}{Bandwidth for the nonparametric fit between y and y-hat (predicted from the binomial model). Defaults to "aicc", calculation of an AICc-minimizing bandwidth. Other options are "gcv", which minimizes the generalized cross-validation statistic, or a numerical bandwidth.}
 \item{color}{Whether the plot should be in color (= TRUE) or grayscale (the default, = FALSE).}
}
\value{
\item{heatmap.obs.p}{A vector of the bootstrapped p-values for each observation in the data set (where the p-value gives the one-tailed proportion of bootstrap replicates with a larger deviation between predicted and empirical probability at that point).}
}
\references{Esarey, Justin and Andrew Pierce (2012). "Assessing Fit Quality and Testing for Misspecification in Binary Dependent Variable Models." Political Analysis 20(4): 480-500.}
\details{This function plots the degree to which a GLM model's predicted probabilities are a good (in-sample) predictor of empirical probabilities. For example, if a model predicts that Pr(y=1)=k\%, about k\% of observations with this predicted probability should have y = 1. Lowess smoothing (with an automatically-selected optimum bandwidth) is used to estimate empirical probabilities in the data set and to overcome sparseness of the data. Systematic deviations are distinguished from sampling variation via bootstrapping of the distribution under the null that the model is an accurate predictor, with p-values indicating the one-tailed proportion of bootstrap samples that are less-extreme than the observed deviation. The plot shows GLM predicted probabilities on the x-axis and smoothed empirical probabilities on the y-axis, with a rug of points indicating the location of sample observations. The ideal fit is a 45-degree line. The shading of the plotted line indicates the degree to which fit deviations are larger than expected due to sampling variation.

	A summary statistic for fit (the "heatmap statistic") is also reported. This statistic is the proportion of the sample in a region with one-tailed p-value less than or equal to 10\%. Finding more than 20\% of the dataset with this p-value in this region is diagnostic of misspecification in the model.

	More details for the technique are given in Esarey and Pierce, "Assessing Fit Quality and Testing for Misspecification in Binary Dependent Variable Models," Political Analysis 2012.}
\examples{

set.seed(459871)

# generate data set
x<-runif(100, min=0, max=10)
y<-ifelse(runif(100, min=0, max=1)<pnorm(0.15*x-2),1,0)
dat2<-as.data.frame(cbind(y,x))

# create fit plot
out<-heatmap.fit(y~x, fam=binomial(link=probit), dat=dat2, reps=1000, span.l="aicc")

}
 
 