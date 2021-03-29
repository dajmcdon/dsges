
scales <- function(vec){
    return(vec[53:59])
}

getClustPred <- function(job, res, lnames){
    sc = scales(res)
    preds = res[9:15]/sc
    out = split(preds, 1:7)
    names(out) = lnames
    return(out)
}

getClustTrain <- function(job, res, lnames){
    scales = scales(res)
    train = res[2:8]/scales
    out = split(train, 1:7)
    names(out) = lnames
    return(out)
}

getClustBase <- function(job, res, lnames){
    sc = scales(res)
    base = res[60:66]/sc
    out = split(base, 1:7)
    names(out) = lnames
    return(out)
}

getClustParms <- function(job, res){
    parms = res[16:51]
    out = split(parms, 1:36)
    return(out)
}

getClustOutput <- function(job, res, parms){
    aa = mspe(res)
    bb = mste(res)
    gg = msbe(res)
    sc = scales(res)
    cc = scaled.mspe(res, sc)
    dd = scaled.mste(res, sc)
    ii = scaled.msbe(res, sc)
    ee = parmerror(res, parms)
    ff = llike(res)
    hh = logpred(res)
    return(list(mspe=aa, mste=bb, msbe=gg, sc.mspe=cc, sc.mste=dd, sc.msbe=ii,
                perror=ee, llike=ff, logPred=hh))
}

scaled.mspe <- function(vec, scales){
  return(mean(vec[9:15]/scales))
}

scaled.mste <- function(vec, scales){
  return(mean(vec[2:8]/scales))
}

scaled.msbe <- function(vec, scales){
    return(mean(vec[60:66]/scales))
}

mspe <- function(vec){
  return(mean(vec[9:15]))
}

mste <- function(vec){
  return(mean(vec[2:8]))
}

msbe <- function(vec){
    return(mean(vec[60:66]))
}

parmerror <- function(vec, parms){
  return(mean((vec[16:51]-parms)^2))
}

llike <- function(vec){
  return(vec[1])
}

logpred <- function(vec){
    return(vec[52])
}


