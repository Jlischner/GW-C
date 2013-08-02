more off;

%# read input file:
%#-----------------
gwc_inp;

%# derived variables:
%#-------------------
ws = dat(:,1);
ReSig = dat(:,2);
ImSig = dat(:,3);

if(do_interpolate == 1);
  wso = ws;
  ReSigo = ReSig;
  ImSigo = ImSig;
  
  dwf = ( ws(2)-ws(1) )/Ninter;
  ws = [min(wso) : dwf : max(wso)]';
  ReSig = interp1(wso,ReSigo,ws);
  ImSig = interp1(wso,ImSigo,ws);
endif;

ReSig = ReSig-vxc-dE_hedin;
Sigma = ReSig + I*ImSig;
A = 1/pi * abs(ImSig) ./ ( (ws-elda-ReSig).^2 + ImSig.^2);
[a,b]=max(A( qp_thresh:length(A) ) );
Eqp = ws( qp_thresh:length(A) )(b); %# note: problem when satellite peak higher than qp peak!

ImSig1 = ImSig;
if(useCor == 1);
  ImSig = dat(:,4);
endif;
Gamma = abs( ImSig )/pi;
dw = ws(2)-ws(1);
dt = ts(2)-ts(1);

elda_save = elda;
if( useEqp == 1);
  elda = Eqp;
endif;

if( min(ws) < elda && elda < max(ws))
  
  [a,IndElda] = min( abs(elda-ws) );
  GammaElda   = Gamma(IndElda);
  DGammaElda  = (Gamma(IndElda+1) - Gamma(IndElda-1) )/(2*dw);
  dE = ReSig(IndElda); %#dE=ReSig(E)-vxc
  etak = abs(ImSig1(IndElda)); %# etak = |Im Sigma(E)|
  DSigmaElda = ( Sigma(IndElda+1) - Sigma(IndElda-1) )/(2*dw);
  
  alphak  =  imag(DSigmaElda); %# alphak = Im dSigma/dw(E)
  gammak  = -real(DSigmaElda); %# gammak = -Re dSigma/dw(E)

else
  printf("error - bad freq. range!");
  break;
endif;

if(zeroalpha == 1);
  alphak = 0.;
endif;

%# formula for cumulant function:
%# C(w) = ( Gamma(w+E)*Theta(mu-[w+E]) - Gamma(E) - w*dGamma/dw(E) )/w^2
%#----------------------------------------------------------------------
CS = Gamma.*( ws < eF )- GammaElda-(ws-elda)*DGammaElda ;
denom = ws-elda;
CS2 = CS ./ denom.^2;

%# doing interpolation to get smooth cumulant:
fitmin = -fitrange;
fitmax = +fitrange;
CS2A= CS2(IndElda + fitmin);
CS2B= CS2(IndElda + fitmax); 

for ii = fitmin:fitmax;
  
  wA = (fitmax-ii)/(fitmax-fitmin);
  wB = (ii-fitmin)/(fitmax-fitmin);
  CS2(IndElda + ii) = wA*CS2A + wB*CS2B;
endfor;

%# calculate first order cumulant contribution
Aqp2 = (etak*cos(alphak)-(ws-Eqp)*sin(alphak) )./( (ws-Eqp).^2 + etak.^2);
Aqp2 *= exp(-gammak)/pi;
dA1s = dw*conv(Aqp2,CS2);

%# dA1 lives on shifted grid wsi:
wsi  = dw*[0:length(Aqp2)+length(CS2)-2]';
wsi += 2*ws(1)-elda;
dA1 = interp1(wsi,dA1s,wsc);

%# second order expansion without working in time
CS2conv = dw*conv(CS2,CS2)/2; 
wsconv  = dw*[0:2*length(CS2)-2]';
wsconv += 2*ws(1)-2*elda;

Aqpconv = (etak*cos(alphak)-(wsconv-Eqp)*sin(alphak) )./( (wsconv-Eqp).^2 + etak.^2);
Aqpconv *= exp(-gammak)/pi;
ddA     = dw*conv(Aqpconv, CS2conv);
wsdd    = dw*[0:2*length(wsconv)-2]';
wsdd   += 2*wsconv(1); 
ddAI    = interp1(wsdd,ddA,wsc);

%# now calculate full cumulant spectral function by fourier transform to real time:
Cqp   = -I*ts*elda - etak*abs(ts) + I*alphak*sign(ts) - gammak;
expCqp= exp(Cqp);
exptw = exp(I*ts*ws');

Csat   = transpose( dw*sum( dmult(CS2, exp(-I*ws*ts')), 1) );
Csat .*= exp(I*elda*ts);
expC   = exp(Csat) .* expCqp;
exptwN = exp(I*ts*wsc');

Aqp  = sum( dmult(expCqp, exptwN ),1);
Aqp  = transpose(Aqp);
Aqp *= dt/2/pi;

Ac = sum( dmult(expC, exptwN ),1); %# GW+C full spectral function
Ac = transpose(Ac);
Ac *= dt/2/pi;

dA = dt* sum( dmult( expCqp .* Csat, exptwN), 1); %# first order satellite corr.
dA = transpose(dA);
dA /= 2*pi;

G = sum( dmult( expC .* (ts<0) , exptwN ),1); %# GW+C full spectral function  
G = transpose(G);
G *= I*dt;
%# self-energy from GW+C:
sigma_c = wsc - elda_save - 1./G;

%# 1st and 2nd order GW+C spectral functions:
A1 = Aqp + dA1;
A2 = A1  + ddAI;


elda = elda_save;

plot(wsc,Ac,'r-',ws,A,'b-',wsc,A2,'g-');
more on;

