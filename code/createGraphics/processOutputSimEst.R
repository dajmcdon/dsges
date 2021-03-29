## Process output
setwd('~/Dropbox/Daniel/DSGEsSuck')
load('output/SimEst/summaryStats.Rdata')
attach(SummaryStats)
pdf(file='gfx/simulateEstimate/mspe.pdf',width=12,height=8)
boxplot(mspe~nestim, main='Average Mean Sq predicton error',
        las=1,xlab='n train (1000 test)')
lines(by(mspe,nestim,mean),col=1)
lines(by(msbe,nestim,mean),col=2)
lines(by(msbe,nestim,median),col=2,lty=2)
dev.off()
pdf(file='gfx/simulateEstimate/mspeSc.pdf',width=12,height=8)
boxplot(sc.mspe~nestim, main='Scaled Average Mean Sq predicton error',las=1,xlab='n train (1000 test)',ylim=range(sc.mspe[nestim==760]))
lines(by(sc.mspe,nestim,mean),col=1)
lines(by(sc.msbe,nestim,mean),col=2)
lines(by(sc.msbe,nestim,median),col=2,lty=2)
dev.off()
pdf(file='gfx/simulateEstimate/mste.pdf',width=12,height=8)
boxplot(mste~nestim, main='Average Mean Sq training error',las=1,xlab='n train')
dev.off()
pdf(file='gfx/simulateEstimate/msteSc.pdf',width=12,height=8)
boxplot(sc.mste~nestim, main='Scaled Average Mean Sq training error',las=1,xlab='n train')
dev.off()
pdf(file='gfx/simulateEstimate/parmerror.pdf',width=12,height=8)
boxplot(perror~nestim, main='Average L2 parameter error',las=1,xlab='n train',ylim=c(0,1))
dev.off()
llike.per.obs = llike/nestim
pdf(file='gfx/simulateEstimate/llikeObs.pdf',width=12,height=8)
boxplot(llike.per.obs~nestim, main='Penalized negative log likelihood per observation',las=1,xlab='n train')
dev.off()
pdf(file='gfx/simulateEstimate/logPred.pdf',width=12,height=8)
boxplot(-logPred/1000~nestim, main='Negative log predictive density/1000', las=1,xlab='n train (1000 predictions)')
dev.off()
detach(SummaryStats)

plot.series.errs <- function(x, base, save.file.name, mainlab, xlab1='n train (1000 test)',
                             col.base=2){
  for(i in 1:7){
    jj = i+5
    ser = names(x)[jj]
    pdf(file=paste0('gfx/simulateEstimate/',save.file.name,ser,'.pdf'), width=12,height=8)
    boxplot(x[[jj]]~x$nestim, main=paste(mainlab, ser), las=1,xlab=xlab1, 
            ylim=range(x[[jj]][x$nestim==760]))
    lines(by(x[[jj]],x$nestim,mean),col=1)
    lines(by(base[[jj]],base$nestim,mean),col=col.base)
    lines(by(base[[jj]],base$nestim,median),col=col.base,lty=2)
    dev.off()
  } 
}

load('output/SimEst/predErrs.Rdata')
load('output/SimEst/baselineErrs.Rdata')
plot.series.errs(pred.errs.sc,base.errs.sc,'predSC','Scaled MSPE:')
plot.series.errs(base.errs.sc,base.errs.sc,'baseSC','Scaled MSPE (true parameters):')

load('output/SimEst/parmEst.Rdata')
load('output/StdPvec.Rdata')
names(pvec) = c('stderrea',  'stderreb',  'stderreg',	'stderreqs',	'stderrem',	'stderrepinf',	'stderrew',	'crhoa',	'crhob',	'crhog',	'crhoqs',	'crhoms',	'crhopinf',	'crhow',	'cmap',	'cmaw',	'csadjcost',	'csigma',	'chabb',	'cprobw',	'csigl',	'cprobp',	'cindw',	'cindp',	'czcap',	'cfc',	'crpi',	'crr',	'cry',	'crdy',	'constepinf',	'constebeta',	'constelab',	'ctrend',	'cgy',	'calfa')
express <- function(char.expressions){
  return(parse(text=paste(char.expressions,collapse=";")))
}
parmGreek = c('sigma[a]','sigma[b]','sigma[g]','sigma[I]', 'sigma[r]', 'sigma[p]','sigma[w]',
              'rho[a]', 'rho[b]', 'rho[g]', 'rho[I]', 'rho[r]', 'rho[p]','rho[w]', 'mu[p]',
              'mu[w]', 'phi1', 'sigma[c]', 'h', 'xi[w]', 'sigma[l]', 'xi[p]', 'iota[w]',
              'iota[p]', 'Psi', 'Phi', 'r[pi]', 'rho', 'r[y]', 'r[Delta][y]', 'bar(pi)', 
              '100(beta^{-1} -1)', 'bar(l)', 'bar(gamma)', 'rho[ga]', 'alpha')
              

plot.parms <- function(x, pvec, xlab1 ='n train', col.base=2){
  for(ii in 1:36){
    jj = ii+5
    ylim = range(x[jj],pvec[ii])
    pdf(file=paste0('gfx/simulateEstimate/parmEstBox/',names(pvec)[ii],'.pdf'),
        width=12,height=8)
    boxplot(x[[jj]]~x$nestim, ylim=ylim,las=1, xlab=xlab1, main = express(parmGreek[ii]))
    abline(h=pvec[ii],col=col.base)
    dev.off()
  }
}

plot.parms(parmEst,pvec)

