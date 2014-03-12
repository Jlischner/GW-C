more off;

%# read input file:
%#-----------------
gwc_inp;

%# derived variables:
%#-------------------
ws    = dat(:,1);
ReSig = dat(:,2);
ImSig = dat(:,3);
dw = ws(2)-ws(1);
dt = ts(2)-ts(1);

if(do_interpolate == 1);
  wso = ws;
  ReSigo = ReSig;
  ImSigo = ImSig;
  
  dw = ( ws(2)-ws(1) )/Ninter;
  ws = [min(wso) : dwf : max(wso)]';
  ReSig = interp1(wso,ReSigo,ws);
  ImSig = interp1(wso,ImSigo,ws);
endif;

%# calculate GW spectral function:
%#--------------------------------
ReSig = ReSig-vxc-dE_hedin;
Sigma = ReSig + I*ImSig;
Agw = 1/pi * abs(ImSig) ./ ( (ws-elda-ReSig).^2 + ImSig.^2);

%# use sharper ImSig if necessary:
%# -------------------------------
ImSig1 = ImSig;
if(useCor == 1);
  ImSig = dat(:,4);
endif;
Gamma = abs( ImSig )/pi;

%# at which energy should the cumulant be evaluated
%# Ec = Eqp - Almbladh/Hedin recipe
%# Ec = elda - Aryasetiawan recipe
%# ------------------------------------------------
if( useEqp == 1);
  [a,b]=max(Agw( qp_thresh:length(Agw) ) );
  Eqp = ws( qp_thresh:length(Agw) )(b); %# note: problem when satellite peak higher than qp peak!
  Ec = Eqp;
else;
  Ec = elda;
endif;

if( min(ws) < Ec && Ec < max(ws))
  
  [a,IndEc] = min( abs(Ec-ws) );
  GammaEc   = Gamma(IndEc);
  DGammaEc  = (Gamma(IndEc+1) - Gamma(IndEc-1) )/(2*dw);
  dE   = ReSig(IndEc); %#dE=ReSig(E)-vxc
  etak = abs(ImSig1(IndEc)); %# etak = |Im Sigma(E)|

  DSigmaEc = ( Sigma(IndEc+1) - Sigma(IndEc-1) )/(2*dw);
  alphak  =  imag(DSigmaEc); %# alphak = Im dSigma/dw(E)
  gammak  = -real(DSigmaEc); %# gammak = -Re dSigma/dw(E)

else
  printf("error - bad freq. range!");
  break;
endif;

if(useEqp == 0);
  Eqp = Ec + dE; %# use on-shell self-energy correction
endif;

if(zeroalpha == 1);
  alphak = 0.; %# set asymmetry factor to zero
endif;

%# formula for cumulant function:
%# C(w) = ( Gamma(w) - Gamma(Ec) - w*dGamma/dw(Ec) )/(w-Ec)^2
%#----------------------------------------------------------------------
CS = Gamma - GammaEc - (ws-Ec)*DGammaEc ;
denom = ws-Ec;
CSpp = ( CS(IndEc+1)-2*CS(IndEc)+CS(IndEc-1) ) / dw^2; %# use analytical form at w=Ec
CS2 = CS ./ denom.^2;
CS2(IndEc) = CSpp/2;

%# now calculate full cumulant spectral function by fourier transform to real time:
%# compare Eq. (221) of Almbadh/Hedin for Cqp(t)
%# note that we have include the energy eigenvalue Eqp in Cqp (instead of just the shift dE)
%# -------------------------------------------------------------------------------------
Cqp    = -I*ts*Eqp - etak*abs(ts) + I*alphak*sign(ts) - gammak;
expCqp = exp(Cqp);

%# Csat(t) = int dw ImSig(w)/(w-Ec)^2 * e^{i(Ec-w)t} [Almbladh Eq. (219), Aryasetiawan PRL Eq. (8)]
%# -----------------------------------------------------------------------------------------------
Csat   = transpose( dw*sum( dmult(CS2, exp(-I*ws*ts')), 1) );
Csat .*= exp(I*Ec*ts);
expC   = exp(Csat) .* expCqp;

%# calculate GWC spectral function:
%# A(w) = 1/(2*pi) int dt e^{iwt} e^{ -iEc*t + C(t) }   [Aryasetiawan/Gunnarsson Eq. (178)]
%# ----------------------------------------------------------------------------------------
exptwN = exp(I*ts*wsc');
Ac = sum( dmult(expC, exptwN ),1); %# GW+C full spectral function
Ac = transpose(Ac)*dt/2/pi;

%# calculate effectiv GWC self energy:
G = sum( dmult( expC .* (ts<0) , exptwN ),1); %# GW+C full spectral function  
G = I*dt*transpose(G);
sigma_c = wsc - elda - 1./G;

%# hedin method: don't separate Cqp from Csat:
%# -------------------------------------------
beta  = interp1(ws-Ec,abs(ImSig),ws); %# need to shift ImSig! beta(w) = |ImSig(w+ek)|/pi
beta(isna(beta)) = 0.0;
beta /= pi;
wt = ws*ts';
f    = exp(-I*wt) + I*wt - ones(size(wt));
f    = dmult(1./ws.^2,f);

[a,b] = min(abs(ws));
f(b,:) = -ts.^2/2; %# use analytic result for w=0
Cret   = dw*sum( dmult(beta,f), 1); %# \int dw beta(w) * ( exp{-iwt} + iwt - 1 )/w^2
Eshift = -dw*real( sum(beta./(ws+I*eta) ) );

Acum = sum( dmult( exp(-I*(Eqp-Eshift)*ts' + Cret) , exptwN ), 1 );
Acum = transpose(Acum); 
Acum *= dt/(2*pi);

%# calculate first order cumulant contribution without working in time:
%# -------------------------------------------------------------------
Aqp  = (etak*cos(alphak)-(ws-Eqp)*sin(alphak) )./( (ws-Eqp).^2 + etak.^2);
Aqp *= exp(-gammak)/pi;
dA1s = dw*conv(Aqp,CS2);

wsi  = dw*[0:length(Aqp)+length(CS2)-2]'; %# dA1 lives on shifted grid wsi
wsi += 2*ws(1)-Ec;
dA1 = interp1(wsi,dA1s,wsc);

%# second order expansion without working in time:
%# -----------------------------------------------
CS2conv = dw*conv(CS2,CS2)/2; 
wsconv  = dw*[0:2*length(CS2)-2]';
wsconv += 2*ws(1)-2*Ec;

Aqp   = (etak*cos(alphak)-(wsconv-Eqp)*sin(alphak) )./( (wsconv-Eqp).^2 + etak.^2);
Aqp  *= exp(-gammak)/pi;
ddA   = dw*conv(Aqp, CS2conv);
wsdd  = dw*[0:2*length(wsconv)-2]';
wsdd += 2*wsconv(1); 
dA2   = interp1(wsdd,ddA,wsc);

%# quasiparticle part of GWC spectral function:
%# --------------------------------------------
Aqp  = (etak*cos(alphak)-(wsc-Eqp)*sin(alphak) )./( (wsc-Eqp).^2 + etak.^2);
Aqp *= exp(-gammak)/pi;

%# add 1st order and 2nd order satellite to Aqp:
%# ---------------------------------------------
A1 = Aqp + dA1;
A2 = A1  + dA2;

plot(wsc,Ac,'r-',ws,Agw,'b-',wsc,A2,'g-',wsc,Acum,'m-');
more on;

