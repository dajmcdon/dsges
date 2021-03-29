## Process output
setwd('~/Dropbox/Daniel/DSGEsSuck/')
## Series flipping
load('output/Perm/permOut.Rdata')
load('data/SWDataDM.Rdata')
source('code/functions/modelsol.R')
source('code/functions/reparlogp.R')
source('code/functions/SimsCodesFort.R')
source('code/priorsetup.R')
require(QZ)
require(FKF)
perms = perm(7,7)
SWPerm = c(3,1,2,6,5,7,4)
truePerm = which.min(rowSums((perms-matrix(SWPerm,nrow=nrow(perms), ncol=7,byrow=TRUE))^2))
series = rownames(y[truePerm])
usualorder = output[truePerm,]
testerr = output[,8:14]
baseline = matrix(usualorder[8:14][perms],ncol=7) ## series in matching order
perc = log(testerr)-log(baseline)
unsort = t(apply(perms,1,function(x) match(SWPerm,x)))
percent.improvement = matrix(NA, nrow(output), 7)
test.error.sort = matrix(NA, nrow(output), 7)
for(i in 1:nrow(output)){
  percent.improvement[i,] = perc[i,unsort[i,]]
  test.error.sort[i,] = testerr[i,unsort[i,]]
}
require(xtable)
perc.errs = rowMeans(percent.improvement)
round(percent.improvement[which.min(perc.errs),],2)
perc.errs[which.min(perc.errs)]

## Figs
pdf(file='gfx/permutations/AvgPercentImprovement.pdf',width=12,height=8)
par(mar=c(5,5,3,1))
plot(sort(perc.errs), pch=19, main='Average Percent Improvement',ylab='mean log ratio',bty='n',xaxt='n',xlab=paste('sorted indices,', round(mean(perc.errs<0),2)*100,'% improved'),las=1)
abline(h=0,col=2)
dev.off()
pdf(file='gfx/permutations/ParamsBoxPlot.pdf',width=12,height=8)
par(mar=c(1,5,3,1))
boxplot(scale(output[,15:50],center=FALSE,scale=prior$sd),las=1,main='Parameter estimates',xaxt='n')
points(output[1,15:50]/prior$sd,col=2,pch=19,cex=.5)
dev.off()
pdf(file='gfx/permutations/ParamsBoxPlotDev.pdf',width=12,height=8)
par(mar=c(1,5,3,1))
boxplot(scale(output[,15:50],center=output[1,15:50],scale=output[1,15:50]), las=1,main='% deviation Parameter estimates from baseline',xaxt='n')
dev.off()
pdf(file='gfx/permutations/SeriesPercentImprovement.pdf',width=12,height=8)
par(mar=c(3,5,3,1))
boxplot(percent.improvement, main='Percent Improvement',ylab='log ratio',bty='n',names=series,las=1)
abline(h=0,col=2)
dev.off()


sc = apply(y,1,var)
scaled.errors = matrix(NA,nrow(output),7)
for(i in 1:nrow(output)) scaled.errors[i,] = test.error.sort[i,]/sc
sc.errs = rowMeans(scaled.errors)

print(xtable(scaled.errors[c(1,which.min(sc.errs)),]))
cat('Best scaled MSE',sc.errs[c(1,which.min(sc.errs))],'\n')

pdf(file='gfx/permutations/AvgScaledError.pdf',width=12,height=8)
par(mar=c(5,5,3,1))
plot(sort(log(sc.errs)), pch=19, main='Average Scaled MSE',ylab='log scale',bty='n',xaxt='n',xlab=paste('sorted indices,', round(mean(sc.errs<sc.errs[1]),2)*100,'% improved'),las=1)
abline(h=log(sc.errs[1]),col=2)
dev.off()
pdf(file='gfx/permutations/SeriesScaledError.pdf',width=12,height=8)
par(mar=c(3,5,3,1))
boxplot(log(scaled.errors), main='Scaled MSE',ylab='log scale',bty='n',names=series,las=1)
points(log(scaled.errors[1,]), pch=19, col=2,cex=1)
dev.off()


# llike = double(nrow(perms))
# for(ii in 1:nrow(perms)){
#   if(ii %% 50 == 0) print(ii)
#   llike[ii] = getlogLike(output[ii,15:50], y[perms[ii,],1:200], prior)
# }
# save(llike, file='output/Perm/llikeCalc.Rdata')
llikeML = double(nrow(perms))
for(ii in 1:nrow(perms)){
  if(ii %% 50 == 0) print(ii)
  llikeML[ii] = getlogLike(output[ii,15:50], y[perms[ii,],1:200], prior, ML=TRUE)
}
save(llikeML, file='output/Perm/llikeCalcML.Rdata')

llikeAll = double(nrow(perms))
for(ii in 1:nrow(perms)){
   if(ii %% 50 == 0) print(ii)
   llikeAll[ii] = getlogLike(output[ii,15:50], y[perms[ii,],1:231], prior)
 }
load('output/Perm/llikeCalc.Rdata')
llikeTest = llikeAll - llike 
save(llikeTest, file='output/Perm/llikeTestCalc.Rdata')



load('output/Perm/llikeCalc.Rdata')
pdf(file='gfx/permutations/TrainLogLike.pdf',width=12,height=8)
par(mar=c(5,5,3,1))
plot(sort(llike), main='Negative Penalized Log Likelihood',ylab='',
     xlab=paste('sorted indices,', round(mean(llike<llike[1]),2)*100,'% improved'),
     bty='n',las=1)
abline(h=llike[1],col=2)
dev.off()

load('output/Perm/llikeTestCalc.Rdata')
pdf(file='gfx/permutations/PredLogLike.pdf',width=12,height=8)
par(mar=c(5,5,3,1))
plot(sort(llikeTest), main='Negative Predictive Likelihood',ylab='',
     xlab=paste('sorted indices,', round(mean(llikeTest<llikeTest[1]),2)*100,'% improved'),
     bty='n',las=1)
abline(h=llikeTest[1],col=2)
dev.off()

top20sc = sort(sc.errs, index.return=TRUE)$ix[1:20]
top20perc = sort(perc.errs, index.return=TRUE)$ix[1:20]
top20llike = sort(llike, index.return=TRUE)$ix[1:20]

mat = cbind(matrix(rownames(y)[perms[top20sc,]],nrow=20),matrix(rownames(y)[perms[top20perc,]],nrow=20),matrix(rownames(y)[perms[top20llike,]],nrow=20))


require(xtable)
print(xtable(mat))

cat('Percent improve Scaled = ',round(mean(sc.errs<sc.errs[1]),2)*100,'\n')
cat('Percent improve Average % = ',round(mean(perc.errs<0),2)*100,'\n')
cat('Percent improve Llike = ',round(mean(llike<llike[1]),2)*100,'\n')

## Plotting forecasts
load('output/Perm/llikeCalc.Rdata')
load('output/Perm/llikeCalcML.Rdata')
predPerms = c(truePerm, which.min(llike), which.min(llikeML), which.min(llikeTest), 
              which.min(perc.errs), which.min(sc.errs))
errs = list()
preds = list()
for(i in 1:length(predPerms)){
  ypermed = y[perms[predPerms[i],],]
  filt = getFilterOutput(ypermed, output[predPerms[i],15:50])
  errs[[i]] = filt$vt[unsort[predPerms[i],],]
  preds[[i]] = y-errs[[i]]
}
for(i in 1:nrow(y)){
  #pdf(file=paste0('gfx/permutations/forecasts',rownames(y)[i],'.pdf'),width=12,height=8)
  #ylims=range(lapply(errs, function(x) range(x[i,201:231])))
  ylims=range(y[i,201:231], lapply(preds, function(x) range(x[i,201:231])))
  par(mar=c(5,3,1,1))
  xx = seq(2006.25,2013.75,by=.25)
  plot(xx, y[i,201:231], ty='l',col=1, las=1, bty='n',
       main='', ylab='', xlab=rownames(y)[i], lwd=1, ylim=ylims)
  cols=c(2,3,3,4,4,4)
  ltys = c(1,1,2,1,2,3)
  for(j in 1:length(preds)) lines(xx, preds[[j]][i,201:231], col=cols[j],
                                  lty=ltys[j], lwd=.5)
  #dev.off()
}
xtable(matrix(series[perms[predPerms,]],ncol=7))
testerrTab = scaled.errors[predPerms,]
xtable(cbind(testerrTab,sc.errs[predPerms],llike[predPerms]))
pdf(file='gfx/permutations/inSampLabor.pdf',width=12,height=8)
ylims=range(y[4,], preds[[1]][4,])
par(mar=c(5,3,1,1))
xx = seq(to=2013.75,by=.25,length=231)
plot(NA, NA, ty='l',col=1, las=1, bty='n',
     main='', ylab='', xlab='hours worked', lwd=1, ylim=ylims,xlim=range(xx))
rect(2006,-1000,2015,1000,col='lightgrey',border='transparent')
lines(xx, y[4,], col=1)
lines(xx, preds[[1]][4,], col=2)
fewdev.off()
