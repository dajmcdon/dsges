gensolution <- function(pvec) {
  ## Function takes in a lengthy vector of parameters and returns matrices nearly ready for Kalman input
  stderr.ea = pvec[1]
  stderr.eb = pvec[2]
  stderr.eg = pvec[3]
  stderr.eqs = pvec[4]
  stderr.em = pvec[5]
  stderr.epinf = pvec[6]
  stderr.ew = pvec[7]
    ##
  crhoa = pvec[8]
  crhob = pvec[9]
  crhog = pvec[10]
    crhoqs = pvec[11]
    crhoms = pvec[12]
    crhopinf = pvec[13]
    crhow = pvec[14]
    cmap = pvec[15]
    cmaw = pvec[16]
    csadjcost = pvec[17]
    csigma = pvec[18]
    chabb = pvec[19]
    cprobw = pvec[20]
    csigl = pvec[21]
    cprobp = pvec[22]
    cindw = pvec[23]
    cindp = pvec[24]
    czcap = pvec[25]
    cfc = pvec[27]
    crpi = pvec[27]
    crr = pvec[28]
    cry = pvec[29]
    crdy = pvec[30]
    constepinf = pvec[31]
    constebeta = pvec[32]
    constelab = pvec[33]
    ctrend = pvec[34]
    cgy = pvec[35]
    calfa = pvec[36]
    
    nx = 52
    ny = 7
    ns = 7
    
    ## fixed parameters
    ctou = .025
    clandaw = 1.5
    cg = 0.18
    curvp = 10
    curvw = 10
    
    ## derived from steady state
    cbeta  =  (1 + constebeta / 100) ^ (-1)
    cgamma  =  ctrend / 100 + 1
    cpie    =  constepinf / 100 + 1
    
    clandap = cfc
    cbetabar = cbeta * cgamma ^ (-1 * csigma)
    cr = cpie / (cbeta * cgamma ^ (-1 * csigma))
    crk = (cbeta ^ (-1)) * (cgamma ^ csigma) - (1 - ctou)
    cw  =  (calfa ^ calfa * (1 - calfa) ^ (1 - calfa) / (clandap * crk ^
                                                             calfa)) ^ (1 / (1 - calfa))
    cikbar = (1 - (1 - ctou) / cgamma)
    cik = (1 - (1 - ctou) / cgamma) * cgamma
    clk = ((1 - calfa) / calfa) * (crk / cw)
    cky = cfc * (clk) ^ (calfa - 1)
    ciy = cik * cky
    ccy = 1 - cg - cik * cky
    crkky = crk * cky
    cwhlc = (1 / clandaw) * (1 - calfa) / calfa * crk * cky / ccy
    cwly = 1 - crk * cky
    
    
    conster = (cr - 1) * 100
    
    
    ##------------------------------------------
    ## Specify size of system of equations
    ##------------------------------------------
    nvar  =  nx          ## number of variables
    nshocks  =  7        ## number of shocks
    
    
    ##--------------------------------------------------
    ## Cast the system in Gensys form:
    ## G0*X(t)  =  G1*X(t-1) + Psi*e(t) + Pi*eta(t) + CC
    ##--------------------------------------------------
    
    ## Initializations
    G0  =  matrix(0, nvar, nvar)
    G1  =  matrix(0, nvar, nvar)
    CC  =  matrix(0, nvar, 1)
    Psi  =  matrix(0, nvar, nshocks)
    Pi  =  matrix(0, nvar, nvar)
    
    ## Indexing variables
    labobs  = 1
    robs  =  2
    pinfobs  = 3
    dy  =  4
    dc  =  5
    dinve  =  6
    dw  =  7
    ewma  =  8
    epinfma  = 9
    zcapf  = 10
    rkf  =  11
    kf  =  12
    pkf   =  13
    cf  =  14
    invef  =  15
    yf  =  16
    labf  =  17
    wf  =  18
    rrf  =  19
    mc  =  20
    zcap  =  21
    rk  =  22
    k  =  23
    pk  =  24
    c  =  25
    inve  =  26
    y  =  27
    lab  =  28
    pinf  =  29
    w  =  30
    r  =  31
    a  =  32
    b  =  33
    g  =  34
    qs  =  35
    ms  =  36
    spinf  =  37
    sw  =  38
    kpf  =  39
    kp   =  40
    
    ## dummies : for forward-looking variables##
    invef1  =  41
    rkf1    =  42
    pkf1    =  43
    cf1     =  44
    labf1   =  45
    inve1   =  46
    rk1     =  47
    pk1     =  48
    pinf1   =  49
    lab1    =  50
    c1      =  51
    w1      =  52
    
    ## shocks :
    ea     =  1
    eb     =  2
    eg     =  3
    eqs    =  4
    em     =  5
    epinf  =  6
    ew     =  7
    
    
    ## Define Log-Linearized Equations
    ##---------------------------------------------------------------
    ## For forward-looking variables, Replace E.t X.t+1 with
    ## y.t  =  E.t X.t+1 and
    ## add to the equation X.t  =  y.t-1 + eta.t (  =  X.t - E.t-1 X.t)
    ##---------------------------------------------------------------
    
    ##-------------------------------
    ## Flexible Economy
    ##-------------------------------
    
    ## 0*(1-calfa)*a + 1*a  =   calfa*rkf+(1-calfa)*(wf)
    G0[1, a]  =   0 * (1 - calfa) + 1
    G0[1, rkf]  =  -calfa
    G0[1, wf]  =  -(1 - calfa)
    
    ## zcapf  =   (1/(czcap/(1-czcap)))* rkf
    G0[2, zcapf]  =  1
    G0[2, rkf]  =  -(1 / (czcap / (1 - czcap)))
    
    ## rkf  =   (wf)+labf-kf
    G0[3, rkf]  =  1
    G0[3, wf]  =  -1
    G0[3, labf]  =  -1
    G0[3, kf]  =  1
    
    ## kf  =   kpf(-1)+zcapf
    G0[4, kf]  =  1
    G1[4, kpf]  =  1
    G0[4, zcapf]  =  -1
    
    
    ## invef  =  (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs
    G0[5, invef]  =  1
    G1[5, invef]  =  1 / (1 + cbetabar * cgamma)
    G0[5, invef1]  =  -(1 / (1 + cbetabar * cgamma)) * cbetabar * cgamma
    G0[5, pkf]  =  -(1 / (1 + cbetabar * cgamma)) * (1 / (cgamma ^ 2 * csadjcost))
    G0[5, qs]   =  -1
    
    ## invef1  =  invef(1)
    G0[6, invef]  =  1
    G1[6, invef1]  =  1
    Pi[6, invef]  =  1
    
    ## pkf  =  -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +
    ## ((1-ctou)/(crk+(1-ctou)))*pkf(1)
    G0[7, pkf]  =  1
    G0[7, rrf]  =  1
    G0[7, b]  =  0 - (1 / ((1 - chabb / cgamma) / (csigma * (1 + chabb /
                                                                 cgamma))))
    G0[7, rkf1]  =  -(crk / (crk + (1 - ctou)))
    G0[7, pkf1]  =  -((1 - ctou) / (crk + (1 - ctou)))
    
    ## rkf1  =  rkf[1)
    G0[8, rkf]  =  1
    G1[8, rkf1]  =  1
    Pi[8, rkf]  =  1
    
    ## pkf1  =  pkf[1]
    G0[9, pkf]  =  1
    G1[9, pkf1]  =  1
    Pi[9, pkf]  =  1
    
    ## cf  =  (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1)
    ## +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/
    ## (csigma*(1+chabb/cgamma))*(rrf+0*b) + b
    G0[10, cf]  =  1
    G1[10, cf]  =  (chabb / cgamma) / (1 + chabb / cgamma)
    G0[10, cf1]  =  -(1 / (1 + chabb / cgamma))
    G0[10, labf]  =  -((csigma - 1) * cwhlc / (csigma * (1 + chabb / cgamma)))
    G0[10, labf1]  =   ((csigma - 1) * cwhlc / (csigma * (1 + chabb / cgamma)))
    G0[10, rrf]  =  (1 - chabb / cgamma) / (csigma * (1 + chabb / cgamma))
    G0[10, b]  =  (1 - chabb / cgamma) / (csigma * (1 + chabb / cgamma)) *
        0 - 1
    
    ## cf1  =  cf[1]
    G0[11, cf]  =  1
    G1[11, cf1]  =  1
    Pi[11, cf]  =  1
    
    ## labff  =  labf[1]
    G0[12, labf]  =  1
    G1[12, labf1]  =  1
    Pi[12, labf]  =  1
    
    ## yf  =  ccy*cf+ciy*invef+g  +  crkky*zcapf
    G0[13, yf]  =  1
    G0[13, cf]  =  -ccy
    G0[13, invef]  =  -ciy
    G0[13, g]   =  -1
    G0[13, zcapf]  =  -crkky
    
    ## yf  =  cfc*( calfa*kf+(1-calfa)*labf +a )
    G0[14, yf]  =  1
    G0[14, kf]  =  -cfc * calfa
    G0[14, labf]  =  -cfc * (1 - calfa)
    G0[14, a]   =  -cfc
    
    ## wf  =  csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1)
    G0[15, wf]  =  1
    G0[15, labf]  =  -csigl
    G0[15, cf]  =  -(1 / (1 - chabb / cgamma))
    G1[15, cf]   =  -(chabb / cgamma) / (1 - chabb / cgamma)
    
    ## kpf  =   (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs
    G0[16, kpf]  =  1
    G1[16, kpf]  =  (1 - cikbar)
    G0[16, invef]  =  -(cikbar)
    G0[16, qs]   =  -(cikbar) * (cgamma ^ 2 * csadjcost)
    
    
    ##-------------------------------
    ## sticky price - wage economy
    ##-------------------------------
    
    ##  mc  =   calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a
    G0[17, mc]  =  1
    G0[17, rk]  =  -calfa
    G0[17, w]  =  -(1 - calfa)
    G0[17, a]   =  1 + 0 * (1 - calfa)
    
    ## zcap  =   (1/(czcap/(1-czcap)))* rk
    G0[18, zcap]  =  1
    G0[18, rk]  =  -(1 / (czcap / (1 - czcap)))
    
    ## rk  =   w+lab-k
    G0[19, rk]  =  1
    G0[19, w]  =  -1
    G0[19, lab]  =  -1
    G0[19, k]   =  1
    
    ## k  =   kp(-1)+zcap
    G0[20, k]  =  1
    G1[20, kp]  =  1
    G0[20, zcap]  =  -1
    
    ## inve  =  (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs
    G0[21, inve]  =  1
    G1[21, inve]  =  1 / (1 + cbetabar * cgamma)
    G0[21, inve1]  =  -(1 / (1 + cbetabar * cgamma)) * cbetabar * cgamma
    G0[21, pk]  =  -(1 / (1 + cbetabar * cgamma)) * (1 / (cgamma ^ 2 * csadjcost))
    G0[21, qs]   =  -1
    
    ##inve1  =  inve[1]
    G0[22, inve]  =  1
    G1[22, inve1]  =  1
    Pi[22, inve]  =  1
    
    ## pk  =  -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +
    ## (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1)
    G0[23, pk]  =  1
    G0[23, r]  =  1
    G0[23, pinf1]  =  -1
    G0[23, b]  =  0 - (1 / ((1 - chabb / cgamma) / (csigma * (1 + chabb /
                                                                  cgamma))))
    G0[23, rk1]  =  -(crk / (crk + (1 - ctou)))
    G0[23, pk1]  =  -((1 - ctou) / (crk + (1 - ctou)))
    
    ## rk1  =  rk[1]
    G0[24, rk]  =  1
    G1[24, rk1]  =  1
    Pi[24, rk]  =  1
    
    ## pk1  =  pk[1]
    G0[25, pk]  =  1
    G1[25, pk1]  =  1
    Pi[25, pk]  =  1
    
    ## pinf1  =  pinf[1]
    G0[26, pinf]  =  1
    G1[26, pinf1]  =  1
    Pi[26, pinf]  =  1
    
    ## c  =  (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1)
    ##     +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) -
    ## (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b
    G0[27, c]  =  1
    G1[27, c]  =  (chabb / cgamma) / (1 + chabb / cgamma)
    G0[27, c1]  =  -(1 / (1 + chabb / cgamma))
    G0[27, lab]  =  -((csigma - 1) * cwhlc / (csigma * (1 + chabb / cgamma)))
    G0[27, lab1]  =  ((csigma - 1) * cwhlc / (csigma * (1 + chabb / cgamma)))
    G0[27, r]  =  (1 - chabb / cgamma) / (csigma * (1 + chabb / cgamma))
    G0[27, pinf1]  =  -(1 - chabb / cgamma) / (csigma * (1 + chabb / cgamma))
    G0[27, b]  =  (1 - chabb / cgamma) / (csigma * (1 + chabb / cgamma)) *
        0 - 1
    
    ## lab1  =  lab[1]
    G0[28, lab]  =  1
    G1[28, lab1]  =  1
    Pi[28, lab]  =  1
    
    ## c1  =  c[1]
    G0[29, c]  =  1
    G1[29, c1]  =  1
    Pi[29, c]  =  1
    
    ## y  =  ccy*c+ciy*inve+g  +  1*crkky*zcap
    G0[30, y]  =  1
    G0[30, c]  =  -ccy
    G0[30, inve]  =  -ciy
    G0[30, g]   =  -1
    G0[30, zcap]  =  -1 * crkky
    
    ## y  =  cfc*( calfa*k+(1-calfa)*lab +a )
    G0[31, y]  =  1
    G0[31, k]  =  -cfc * calfa
    G0[31, lab]  =  -cfc * (1 - calfa)
    G0[31, a]   =  -cfc
    
    ## pinf  =   (1/(1+cbetabar*cgamma*cindp]] * ( cbetabar*cgamma*pinf(1] +cindp*pinf(-1]
    ##               +((1-cprobp]*(1-cbetabar*cgamma*cprobp]/cprobp]/((cfc-1]*curvp+1]*(mc]  ]  + spinf
    G0[32, pinf]  =  1
    G0[32, pinf1]  =  -(1 / (1 + cbetabar * cgamma * cindp)) * (cbetabar *
                                                                    cgamma)
    G1[32, pinf]  =  (1 / (1 + cbetabar * cgamma * cindp)) * cindp
    G0[32, mc]   =  -(1 / (1 + cbetabar * cgamma * cindp)) * ((1 - cprobp) *
                                                                  (1 - cbetabar * cgamma * cprobp) / cprobp) / ((cfc - 1) * curvp + 1)
    G0[32, spinf]  =  -1
    
    ##	      w  =   (1/(1+cbetabar*cgamma))*w(-1)
    ##               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
    ##               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
    ##               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
    ##               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
    ##               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
    ##               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w)
    ##               + 1*sw
    G0[33, w]  =  1 + (1 - cprobw) * (1 - cbetabar * cgamma * cprobw) /
        ((1 + cbetabar * cgamma) * cprobw) * (1 / ((clandaw - 1) * curvw + 1))
    G1[33, w]  =  (1 / (1 + cbetabar * cgamma))
    G0[33, w1]  =  -(cbetabar * cgamma / (1 + cbetabar * cgamma))
    G1[33, pinf]  =  (cindw / (1 + cbetabar * cgamma))
    G0[33, pinf]   =  (1 + cbetabar * cgamma * cindw) / (1 + cbetabar *
                                                             cgamma)
    G0[33, pinf1]  =  -(cbetabar * cgamma) / (1 + cbetabar * cgamma)
    G0[33, lab]  =  -(1 - cprobw) * (1 - cbetabar * cgamma * cprobw) / ((1 +
                                                                             cbetabar * cgamma) * cprobw) *
        (1 / ((clandaw - 1) * curvw + 1)) * csigl
    G0[33, c]  =  -(1 - cprobw) * (1 - cbetabar * cgamma * cprobw) / ((1 +
                                                                           cbetabar * cgamma) * cprobw) *
        (1 / ((clandaw - 1) * curvw + 1)) * (1 / (1 - chabb / cgamma))
    G1[33, c]  =  -(1 - cprobw) * (1 - cbetabar * cgamma * cprobw) / ((1 +
                                                                           cbetabar * cgamma) * cprobw) *
        (1 / ((clandaw - 1) * curvw + 1)) * ((chabb / cgamma) / (1 - chabb /
                                                                     cgamma))
    G0[33, sw]  =  -1
    
    ## w1  =  w[1]
    G0[34, w]  =  1
    G1[34, w1]  =  1
    Pi[34, w]  =  1
    
    ## r  =   crpi*(1-crr)*pinf + cry*(1-crr)*(y-yf) +crdy*(y-yf-y(-1)+yf(-1)) +crr*r(-1) +ms
    G0[35, r]  =  1
    G0[35, pinf]  =  -crpi * (1 - crr)
    G0[35, y]  =  -cry * (1 - crr) - crdy
    G0[35, yf]  =  cry * (1 - crr) + crdy
    G1[35, y]  =  -crdy
    G1[35, yf]  =  crdy
    G1[35, r]   =  crr
    G0[35, ms]  =  -1
    
    ## a  =  crhoa*a(-1)  + ea
    G0[36, a]  =  1
    G1[36, a]  =  crhoa
    Psi[36, ea]  =  stderr.ea
    
    ## b  =  crhob*b[-1] + eb
    G0[37, b]  =  1
    G1[37, b]  =  crhob
    Psi[37, eb]  =  stderr.eb
    
    ## g  =  crhog*(g[-1]) + eg + cgy*ea
    G0[38, g]  =  1
    G1[38, g]  =  crhog
    Psi[38, eg]  =  stderr.eg
    Psi[38, ea]  =  cgy * stderr.ea
    
    ## qs  =  crhoqs*qs[-1] + eqs
    G0[39, qs]  =  1
    G1[39, qs]  =  crhoqs
    Psi[39, eqs]  =  stderr.eqs
    
    ## ms  =  crhoms*ms[-1] + em
    G0[40, ms]  =  1
    G1[40, ms]  =  crhoms
    Psi[40, em]  =  stderr.em
    
    ## spinf  =  crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1)
    G0[41, spinf]  =  1
    G1[41, spinf]  =  crhopinf
    G0[41, epinfma]  =  -1
    G1[41, epinfma]  =  -cmap
    
    ## epinfma = epinf
    G0[42, epinfma]  =  1
    Psi[42, epinf]  =  stderr.epinf
    
    ## sw  =  crhow*sw[-1] + ewma - cmaw*ewma[-1]
    G0[43, sw]  =  1
    G1[43, sw]  =  crhow
    G0[43, ewma]  =  -1
    G1[43, ewma]  =  -cmaw
    
    ## ewma = ew
    G0[44, ewma]  =  1
    Psi[44, ew]  =  stderr.ew
    
    ## kp  =   [1-cikbar]*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs
    G0[45, kp]  =  1
    G1[45, kp]  =  1 - cikbar
    G0[45, inve]  =  -cikbar
    G0[45, qs]  =  -cikbar * cgamma ^ 2 * csadjcost
    
    ##-------------------------------
    ## measurement equations
    ##-------------------------------
    
    ##dy = y-y[-1]+ctrend
    G0[46, dy]  =  1
    G0[46, y]  =  -1
    G1[46, y]  =  -1
    CC[46]  =  ctrend
    
    ##dc = c-c[-1]+ctrend
    G0[47, dc]  =  1
    G0[47, c]  =  -1
    G1[47, c]  =  -1
    CC[47]  =  ctrend
    
    ##dinve = inve-inve[-1]+ctrend
    G0[48, dinve]  =  1
    G0[48, inve]  =  -1
    G1[48, inve]  =  -1
    CC[48]  =  ctrend
    
    ##dw = w-w[-1]+ctrend
    G0[49, dw]  =  1
    G0[49, w]  =  -1
    G1[49, w]  =  -1
    CC[49]  =  ctrend
    
    ##pinfobs  =  1*[pinf] + constepinf
    G0[50, pinfobs]  =  1
    G0[50, pinf]  =  -1
    CC[50]  =  constepinf
    
    ##robs  =     1*[r] + conster
    G0[51, robs]  =  1
    G0[51, r]  =  -1
    CC[51]  =  conster
    
    ##labobs  =  lab + constelab
    G0[52, labobs]  =  1
    G0[52, lab]  =  -1
    CC[52]  =  constelab
    
    ##----------------------------------------
    ## Create a solution space using gensys
    ##----------------------------------------
    if (anyNA(G0) ||
        anyNA(G1) ||
        any(is.infinite(G0)) ||
        any(is.infinite(G1)))
        return(list(eu = c(-2, -2)))
    system =  gensys(G0, G1, CC, Psi, Pi)
    G = system$G1
    C = system$C
    M = system$impact
    eu = system$eu
    
    ## C0[1]  =  CC[52]
    ## C0[2]  =  CC[51]
    ## C0[3]  =  CC[50]
    ## C0[4:7]  =  CC[49]
    
    simsout  =  list(
        G = G,
        C = C,
        M = M,
        eu = eu,
        nx = nx,
        ny = ny,
        ns = ns
    )
    return(simsout)
}

ss.model <- function(y, simsout, output = 'll', kill = 1e-10) {
    eu = simsout$eu
    if (output == 'll' && (eu[1] == -2 || eu[2] == -2)) return(-Inf)
    G = simsout$G
    C = simsout$C
    M = simsout$M
    nx = simsout$nx
    ny = simsout$ny
    ns = simsout$ns
    
    G[abs(G) < kill] = 0
    C[abs(C) < kill] = 0
    M[abs(M) < kill] = 0
    
    
    ## write the system in terms of
    ## yo_t+1 = c+ Z*x_t + G*v_t+1
    ## x_t+1 = d+ T*x_t   + H*w_t+1
    Zt = cbind(diag(ny), matrix(0, ny, nx - ny))   # it is not clear to me that the observations are in the first 7 (last 7?)
    dt = as.matrix(C)
    Tt = G
    HH = M %*% t(M)
    GG = diag(0, ny) # note: no observation errors?
    ct = matrix(0, ny, 1)
    
    
    a0 = double(nx)
    P0 = 10 * diag(nx) # as in the matlab code
    ## P0 = matrix(solve(diag(nx^2) - T %x% T) %*% c(HH), nx, nx)
    ## prior is the steady state (way too slow, c(HH) is wrong)
    out = fkf(a0, P0, dt, ct, Tt, Zt, HH, GG, y)
    if (output == 'll') {
        if (any(is.na(out$Ft)))
            out$logLik = -Inf
        return(out$logLik)
    }
    out$ssmats =  list(
        Zt = Zt,
        dt = dt,
        Tt = Tt,
        HH = HH,
        GG = GG,
        ct = ct,
        a0 = a0,
        P0 = P0
    )
    return(out)
}

getlogLike <-  function(parm, y, prior, ML = FALSE) {
    transformed = repar(parm, prior)
    lnprior = transformed$lnprior
    outparm = transformed$outparv
    
    simsout = gensolution(outparm)
    llike = ss.model(y, simsout) #regular log likelihood (would need to maximize
    ifelse(ML, return(-1 * llike), return(-1 * llike - 1 * lnprior))
}


estimate <- function(y, prior, nstarts = 10, trunc = 4, ML = FALSE, save.all = FALSE,
                     tol = 1e-6, track = 0, maxit = 100, badvals = 1e8) {
    init = t(as.matrix(start.vals(prior, nstarts-1)$unbounded))
    init = rbind(repar(prior$sw, prior, TRUE)$outparv, init)
    out1 = multistart(init, getlogLike, y=y, prior=prior,
                      method='Nelder-Mead',
                      control=list(trace=track, maxit=maxit, badval=badvals))
    out2 = multistart(init, getlogLike, y=y, prior=prior,
                      method='SANN',
                      control=list(trace=track, maxit=maxit, badval=badvals))
    minima = rbind.data.frame(out1, out2)
    parm = unlist(minima[which.min(minima$value), 1:ncol(init)])
    parm = repar(parm, prior)$outparv
    ## simsout = gensolution(parm)
    ## filt = ss.model(y, simsout, output='all')
    ## ytrunc = y[,-c(1:trunc)]
    ## pred = ytrunc - filt$vt[,-c(1:trunc)]
    ## train = rowMeans( filt$vt[,-c(1:trunc)]^2 )
    ## base = rowMeans( (y[,-c(1:trunc)] - rowMeans(y[,-c(1:trunc)]))^2 )
        
    out = list()
    out$parm = parm
    ## out$train = train
    ## out$baseline = base
    ## out$logpost = objective
    ## if(save.all){
    ##     out$pred = pred[,-c(1:trunc,n+1)]
    ##     out$y = ytrunc
    ##     out$ssout = filt
    ## }
    return(out)
}

perm <- function (n, r, v = 1:n) {
    ## gets matrix of all permutations
    if (r == 1)
        matrix(v, n, 1)
    else if (n == 1)
        matrix(v, 1, r)
    else{
        X <- NULL
        for (i in 1:n)
            X <- rbind(X, cbind(v[i], perm(n - 1, r - 1, v[-i])))
        X
    }
}

estim.pred.wrap <- function(v, y, prior, estim.obs = 1:ncol(y), pred.obs = NULL,
                            small = TRUE, save.preds = FALSE, trunc = 4, ...) {
    out = estim.pred(y[v, ], prior, estim.obs, pred.obs, small, save.preds, trunc, ...)
    return(out)
}

estim.pred <- function(y, prior, estim.obs = 1:ncol(y), pred.obs = NULL, small = TRUE,
                       save.preds = FALSE, trunc = 4, ...) {
    est = estimate(y[, estim.obs], prior, save.all = FALSE, trunc = trunc, ...)
        
    if (is.null(pred.obs)) return(est)
    ret = list()
    estim.obs = estim.obs[-c(1:trunc)]
    simsout = gensolution(est$parm)
    filt = ss.model(y, simsout, output = 'all')
    ytrunc = y[, -c(1:trunc)]
    pred = ytrunc - filt$vt[, -c(1:trunc)]
    ret$train.err = rowMeans(filt$vt[, estim.obs] ^ 2)
    ret$test.err = rowMeans(filt$vt[, pred.obs] ^ 2)
    if (small) return(c(ret$train.err, ret$test.err, est$parm))
    ## if(save.preds){
    ##     ret$base.train = rowMeans( (y[,estim.obs] - rowMeans(y[,estim.obs]))^2 )
    ##     ret$base.test = rowMeans( (y[,pred.obs] - rowMeans(y[,estim.obs]))^2 )
    ##     ret$pred = pred
    ##     ret$y = ytrunc
    ##     ret$parm = est$parm
    ##     ret$ssout = filt
    ##     ret$logpost = est$logpost
    ## }
    return(ret)
}

estim.pred.gen <- function(y, parm2start, prior, nestim = 100, npred = 1000, 
                           maxit = 1e4, trunc = 4) {
    init = repar(parm2start, prior, bound.2.unbound = TRUE)$outparv
    temp = optim(init, getlogLike, y = y[, 1:nestim], prior = prior,
            method = 'SANN', control = list(maxit = maxit))
    objective = temp$value
    parm = repar(temp$par, prior)$outparv
    estim.obs = 1:nestim
    estim.obs = estim.obs[-c(1:trunc)]
    simsout = gensolution(parm)
    filt = ss.model(y, simsout, output = 'all')
    filt2 = ss.model(y[, 1:nestim], simsout)
    filt3 = ss.model(y, gensolution(parm2start), output = 'all')
    logpred = filt$logLik - filt2
    ytrunc = y[, -c(1:trunc)]
    sc = apply(ytrunc, 1, var)
    train.err = rowMeans(filt$vt[, estim.obs] ^ 2)
    pred.obs = (nestim + 1):(nestim + npred)
    test.err = rowMeans(filt$vt[, pred.obs] ^ 2)
    base.err = rowMeans(filt3$vt[, pred.obs] ^ 2)
    return(c(objective, train.err, test.err, parm, logpred, sc, base.err))
}



generate <- function(nobs, pvec = priordraw(nsamp), nburn = 1000, kill = 1e-10) {
    simsout = gensolution(pvec)
    if (simsout$eu == -2 || simsout$eu == -2) stop('Parameter vector does not satisfy Sims.')
        
    G = simsout$G
    C = simsout$C
    M = simsout$M
    nx = simsout$nx
    ny = simsout$ny
    ns = simsout$ns
        
    G[abs(G) < kill] = 0
    C[abs(C) < kill] = 0
    M[abs(M) < kill] = 0
        
        
    ## write the system in terms of
    ## yo_t+1 = c+ Z*x_t + G*v_t+1
    ## x_t+1 = d+ T*x_t   + H*w_t+1
    Zt = cbind(diag(ny), matrix(0, ny, nx - ny))
    dt = as.matrix(C)
    Tt = G
    Ht = M
        
    ## Simulate
    y = matrix(NA, ny, nobs + nburn)
    x = matrix(0, nx, 1)
    for (i in 1:(nobs + nburn)) {
        w = matrix(rnorm(ns, 0, sd = pvec[1:7]), ns, 1)
        x = dt + Tt %*% x + Ht %*% w
        y[, i] = Zt %*% x
    }
    if (nburn > 0) y = y[, -c(1:nburn)]
    return(y)
}

getFilterOutput <- function(y, parm) {
    simsout = gensolution(parm)
    filt = ss.model(y, simsout, output = 'all')
    return(filt)
}
