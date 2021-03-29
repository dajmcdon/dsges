qz.sims <- function(a=array(1,dim=1,1),b=array(1,dim=1,1)) {
  ## R does not load all of lapack in its own lapack.so library.  The routine
  ## needed here, zgges, therefore has to be loaded directly.
  ## if(!is.loaded("zgges")){
  ##   filename <- "/usr/lib/liblapack.dylib"
  ##   ##filename <- "/usr/lib/liblapack.so" #not working with atlas-compiled R.
  ##   ##if (file.access(filename)) 
  ##   ##  filename <- "/usr/lib/liblapackgf-3.so"
  ##   dyn.load(filename, now=FALSE)
  ## }
  N<-dim(a)[1]
  SDIM<-as.integer(1);
  ALPHA<-vector("complex",N)
  BETA <-vector("complex",N)
  q <- array(data=0+0i,dim=c(N,N))
  z <- array(data=0+0i,dim=c(N,N))
  LWORK <- as.integer(33*N)
  WORK <- vector("complex",length=LWORK)
  RWORK <- vector("numeric",length=8*N)
  INFO <- 0
## -------------------------------
## Notes:  The reordering could be done inside zgges.f.  This would require writing a fortran module that contained both the
## routine to check the size or sign of the real part of the root and a separate program, that sets div in the module.  (All this to
## get arund the fact that zgges wants the comparison function to have just two arguments.)  Also, to implement this in f77, and thereby
## guarantee wider portability, one would have to use a common block and ENTRY, rather than a module and double,save.  
  ##-------------------------------
  ## Returns a,b,q,z s.t. q %*% t(Conj(q)) = z %*% t(Conj(z)) = I, a and b upper triangular,
  ## q %*% a %*% t(Conj(z)) = A, and q %*% b %*% t(Conj(z))=B.  The diagonal of b should be
  ## non-negative and real (a normalization).  gev is a two-column array containing the generalized
  ## eigenvalues of the pair A,B.
  ##
  ## rc=0 indicates return without problems.
  ## rc=-18 indicates not enough work space
  ## was allocated.  The current code attempts to be overgenerous in allocating workspace
  ## based on empirical testing of the lapack code with random matrices.  If rc=-18 ever now
  ## occurs in practice, the code could be modified to first test for required space.
  ## Also, if the current version's space requirements are a problem, the code can surely
  ## be modified to use less workspace, with some possible loss in efficiency.
  ## rc=-j indicates an illegal argument type for the j'th argument to the lapack routine.  This should
  ## not occur.
  ## -------------------------------------
  ## The call to dyn.load works on  SuSE 9 and 10 systems with the distribution's installation
  ## of lapack.  If your lapack is at a different location, modify accordingly.
  ## -------------------------------------
  ## The first two V's indicate that we want q and z computed.  The quoted
  ## "N" indicates we want no sorting.  (This would require producing a
  ## different fortran routine for each possible choice of the "stake" around
  ## which we are sorting.) The first dum is where the sort criterion routine's
  ## name would go.  N is the 
  ## dimension of all the matrices.  SDIM would be the number of eigenvalues
  ## for which the sort criterion is true, but since we're not sorting it
  ## will come back 0, always.  ALPHA and BETA are the generalized eigenvalue
  ## vectors (the diagonals of the returned a and b).  WORK and RWORK
  ## are workspace, the former complex, the latter real.  The last dum is
  ## a placeholder for logical workspace needed only if a sort is done.
  ## INFO is 0 if all goes well, -j if the j'th argument has an incorrect
  ## form, j>0 if the algorithm failed, but ALPHA(k) and BETA(k) are accurate
  ## for k>j, N+1 if there was another sort of failure.
  ## ------------------------------
  out<-.Fortran("zgges","V","V","N","dum",N,as.complex(a),N,as.complex(b),N,
                SDIM,ALPHA,BETA,q,N,z,N,WORK,LWORK,RWORK,"dum",INFO)
  gev<-matrix(c(out[[11]],out[[12]]),nrow=N,ncol=2)
  ## if you run into problems with workspace, uncomment the line below and use the
  ## message it produces to modify the setting of LWORK above.
  ## cat("workspace needed:",out[[17]][1],"\n")
  return(list(a=matrix(out[[6]],nrow=N,ncol=N),
              b=matrix(out[[8]],nrow=N,ncol=N),
              q=matrix(out[[13]],nrow=N,ncol=N),
              z=matrix(out[[15]],nrow=N,ncol=N),
              gev=gev,
              rc=out[[21]]))
}
    



qzdiv <- function(stake, qzlist){
    a = qzlist$a
    b = qzlist$b
    q = qzlist$q
    z = qzlist$z
    gev = qzlist$gev
    root =  abs(cbind(diag(a), diag(b)))
    toosmall = root[,1]<1e-13
    root[,2] = root[,2]/root[,1]
    root[toosmall,2] = -1
    ## root = abs(cbind(diag(a), diag(b)))
    ## root[,1] = root[,1]-(root[,1]<1.e-13)*(root[,1]+root[,2])
    ## root[,2] = root[,2]/root[,1]
    bottom = !(root[,2] > stake | root[,2] < -.1)
    out = qz.ztgsen(a, b, q, z, bottom, ijob = 4L, want.Q = TRUE, want.Z = TRUE)
    qzlist$a = out$S
    qzlist$b = out$T
    qzlist$q = out$Q
    qzlist$z = out$Z
    qzlist$gev = cbind(out$ALPHA,out$BETA)
    return(qzlist)
}

gensys <- function(g0, g1, c0=matrix(0,dim(g0)[1],1), psi, pi, div=-1)
  {
    ##System given as
    ##        g0*y(t)=g1*y(t-1)+c0+psi*z(t)+pi*eta(t),
    ##with z an exogenous variable process and eta being endogenously determined
    ##one-step-ahead expectational errors.  Returned system is
    ##       y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) + loose*eta .
    ## If z(t) is i.i.d., the term involving fmat and fwt drops out.
    ## If the solution is unique (eu[2]==1) there is no "loose" term.  Otherwise
    ## loose characterizes the dimensions along which there is non-uniqueness.
    ## If div is omitted from argument list, a div>1 is calculated.
    ## eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu=[-2,-2] for coincident zeros.
    ## By Christopher A. Sims 2/24/2004, from earlier matlab code of same author.
    eu <- c(0,0)
    realsmall <- 1e-7
    fixdiv <- (div>0)
    n <- dim(g0)[1]
    nshock <- if (is.matrix(psi)) dim(psi)[2] else if (is.null(psi)) 0 else 1
    qzl <- qz.sims(g0,g1)
    zxz <- any((abs(diag(qzl$a))<realsmall) & (abs(diag(qzl$b))<realsmall))
    if (zxz) {
      ## "Coincident zeros.  Indeterminacy and/or nonexistence.\n"
      eu <- c(-2,-2)
      gev <- qzl$gev
      return(list(eu=eu,gev=gev))
    }
    zeroax <- abs(diag(qzl$a)) < realsmall
    unstabx <- abs(diag(qzl$a)) < (1-realsmall)*abs(diag(qzl$b)) # near unit roots don't count
    unstabx <- (! zeroax) & unstabx
    if (! fixdiv) {
      if (! any(unstabx)){
        div <- 1.01
      }else{
        div <- .5*(min(abs(diag(qzl$b)[unstabx]/diag(qzl$a)[unstabx]))+1)
      }
    }
    unstabx <- div*abs(diag(qzl$a))<= abs(diag(qzl$b))
    nunstab <- sum(unstabx)
    qzl <- qzdiv(div,qzl)
    qq <- t(Conj(qzl$q))                # to match matlab convention 
    gev <- qzl$gev
    ## note that this means that gev is not simply the diagonals of a nd b.  qzdiv
    ## changes the numbers on the diagonals (though not their ratios), but merely reorders
    ## the original gev.  
    if (nunstab==n){
      six <- NULL
      uix <- 1:n
    }else{
      if (nunstab==0){
        uix <- NULL
        six <- 1:n
      }else{
        uix <- (n-nunstab+1):n
        six <- 1:(n-nunstab)
      }
    }
    q1 <- qq[six,,drop=FALSE]
    q2 <- qq[uix,,drop=FALSE]
    z1 <- t(Conj(qzl$z[,six,drop=FALSE]))
    z2 <- t(Conj(qzl$z[,uix,drop=FALSE]))
    a2 <- qzl$a[uix,uix,drop=FALSE]
    b2 <- qzl$b[uix,uix,drop=FALSE]
    ## debug
    ## browser()
    etawt <- q2 %*% pi
    neta <- if (is.matrix(pi)) dim(pi)[2] else if (is.null(pi)) 0 else 1
    ndeta <- min(nunstab,neta)
    if(ndeta==0){
      ueta <- matrix(0,nunstab,0)
      deta <- vector("numeric",0)
      veta <- matrix(0,neta,0)
      bigev <- vector("logical",0)
    }else{
      sd <- svd(etawt)
      ueta <- sd$u; deta <- sd$d; veta <- sd$v
      bigev <- deta>realsmall
      ueta<-ueta[,bigev,drop=FALSE]
      veta<-veta[,bigev,drop=FALSE]
      deta<-deta[bigev]
    }
    eu[1] <- sum(bigev) >= nunstab
    ##----------------------------------------------------
    ## Note that existence and uniqueness are not just matters of comparing
    ## numbers of roots and numbers of endogenous errors.  These counts are
    ## reported below because usually they point to the source of the problem.
    ##------------------------------------------------------
    etawt1 <- q1 %*% pi
    ndeta1 <- min(n-nunstab,neta)
    if(ndeta1==0){
      ueta1 <- matrix(0,n-nunstab,0)
      deta1 <- vector("numeric",0)
      veta1 <- matrix(0,neta,0)
      bigev1 <- vector("logical",0)
    } else {
      sd <- svd(etawt1)
      ueta1<-sd$u
      deta1 <- sd$d
      veta1 <- sd$v
      bigev1 <- deta1 > realsmall
    }
    if (any(bigev1)) { #needed because empty dimensions are dropped after select
      ueta1 <- ueta1[,bigev1,drop=FALSE]
      veta1 <- veta1[,bigev1,drop=FALSE]
      deta1 <- deta1[bigev1]
      loose <- veta1-veta %*% t(Conj(veta)) %*% veta1
      svdl <- svd(loose)
      loose <- sum(abs(svdl$d)>realsmall*n)
      unq <- (loose==0)
    } else {
      ueta1 <- matrix(1,n-nunstab,0)
      veta1 <- matrix(1,neta,0)
      deta1 <- vector("complex",0)
      unq <- TRUE
    }
    if (unq) eu[2] <- 1
    ## } else
    ## {
    ##   cat("Indeterminacy.", loose, "loose endog errors.\n")
    ## }
    ## Note: if v is a vector of length n and m is an nxp matrix,
    ## v*m==diag(v)%*%m, m/v==solve(diag(v),m)==diag(v)\m (matlab notation)
    ##
    tmat <- cbind(diag(n-nunstab),-t(Conj((ueta %*% (t(Conj(veta))/deta)) %*% veta1 %*%
                                          (deta1 * t(Conj(ueta1)))))  )  
    G0<- rbind( tmat %*% qzl$a, cbind(matrix(0,nunstab,n-nunstab), diag(nunstab)))
    G1<- rbind(tmat %*% qzl$b, matrix(0,nunstab,n))
    ##----------------------
    ## G0 is always non-singular because by construction there are no zeros on
    ## the diagonal of a[1:(n-nunstab),1:(n-nunstab)], which forms G0's ul corner.
    ##-----------------------
    G0I <- solve(G0)
    G1 <- G0I%*%G1
    ##----------- uix can be empty, e.g. in indeterminate systems with no unstable roots ------------
    if(is.null(uix)){
      C <- G0I %*% tmat %*% qq %*% c0
      ## fmat <- matrix(0,0,0)
      ## fwt <- matrix(0, 0, nshock)
      impact <- G0I %*% tmat %*% qq %*% psi
  }else{
      C <- G0I %*% rbind(tmat%*% qq %*%c0,solve(qzl$a[uix,uix,drop=FALSE]-
                                                qzl$b[uix,uix,drop=FALSE],q2%*%c0) )
      impact <- G0I %*% rbind(tmat %*% qq %*% psi, matrix(0,nunstab, nshock))
      ## fmat <- solve(qzl$b[uix,uix,drop=FALSE], qzl$a[uix,uix,drop=FALSE])
      ## fwt <- -solve(qzl$b[uix,uix,drop=FALSE], q2 %*% psi)
    }
    ywt <- G0I[,uix,drop=FALSE]
    ##loose <- G0I %*% etawt1 %*% (diag(neta) - veta %*% t(Conj(veta)))
    loose <- G0I %*% qq %*% pi %*% (diag(neta) - veta %*% t(Conj(veta)))
    ## loose <- G0I %*% rbind(loose,matrix(0,nunstab,neta))  #(think the above is a mistaken remnant)
    ##-------------------- above are output for system in terms of z'y -------
    G1 <- Re(qzl$z %*% G1 %*% t(Conj(qzl$z)))
    C <- Re(qzl$z%*%C)
    impact <- Re(qzl$z%*%impact)
    ywt <- qzl$z%*%ywt
    loose <- Re(qzl$z %*% loose)
    vn <- dimnames(g0)[[2]]
    dimnames(G1) <- list(vn,vn)
    dimnames(C) <- list(vn,NULL)
    dimnames(impact)[[1]] <- vn
    dimnames(ywt)[[1]] <- vn
    dimnames(loose)[[1]] <- vn
    return(list(G1=G1,C=C,impact=impact,ywt=ywt,gev=gev,eu=eu,loose=loose)) # fmat=fmat,fwt=fwt,
  }

