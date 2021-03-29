## make table of parameters
setwd('~/Dropbox/DSGEs suck/')
load('output/StdPvec.Rdata')
load('data/SWDataDM.Rdata')
source('code/functions/reparlogp.R')
source('code/priorsetup.R')
names(pvec) = c('stderr.ea',  'stderr.eb',	'stderr.eg',	'stderr.eqs',	'stderr.em',	'stderr.epinf',	'stderr.ew',	'crhoa',	'crhob',	'crhog',	'crhoqs',	'crhoms',	'crhopinf',	'crhow',	'cmap',	'cmaw',	'csadjcost',	'csigma',	'chabb',	'cprobw',	'csigl',	'cprobp',	'cindw',	'cindp',	'czcap',	'cfc',	'crpi',	'crr',	'cry',	'crdy',	'constepinf',	'constebeta',	'constelab',	'ctrend',	'cgy',	'calfa')
df = prior[c(1,9,10,7,8)]
df$post = pvec
df = as.data.frame(df)
rownames(df)=names(pvec)

dynare.pvec = c(0.45871111106042,  0.24009239419021,	0.529613696317831,	0.4531916867673,	0.244603658393596,	0.140828947359983,	0.242517400032316,	0.957117060263353,	0.21856277590412,	0.976549211012774,	0.711665010893796,	0.148212585537272,	0.889559028271276,	0.966638241757864,	0.699074613590939,	0.835926051950719,	5.74892326430474,	1.38632937132328,	0.714246124427773,	0.701497680306459,	1.84221903167925,	0.651017699411878,	0.583114729890919,	0.241090935891256,	0.544365625097931,	1.60962669099295,	2.04837837701555,	0.811258325894124,	0.0889336233601046,	0.222121568801256,	0.789144390304714,	0.16813331075723,	0.536389578153887,	0.431084643746677,	0.519794486733075,	0.190957432132584)

sw.pvec = c(.45,.24,.52,.45,.24,.14,.24,.95,.18,.97,.71,.12,.9,.97,.74 ,.88, 5.48,1.39,.71,.73, 1.92,.65,.59,.22,.54,1.61,2.03,.81,.08,.22,.81,.16,-.1,.43, .52,.19)

df$swpost = sw.pvec
source('code/functions/modelsol.R')
source('code/functions/SimsCodesFort.R')

stopifnot(require(xtable))
stopifnot(require(QZ))
stopifnot(require(FKF))

xtable(df)
getlogLike(pvec,y,prior)
getlogLike(sw.pvec,y,prior)
(getlogLike(pvec,y,prior,TRUE)-getlogLike(sw.pvec,y,prior,TRUE))/getlogLike(sw.pvec,y,prior,TRUE)

