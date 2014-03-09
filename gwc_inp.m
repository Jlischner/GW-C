%# sigma.dat contains self energy from a GW calculations
%# 1st col: freq grid (eV), 2nd col: ReSig (eV), 3rd col: ImSig (eV)
dat = load("sigma.dat");
elda = -1.3702; %# mean-field energy in eV
vxc  = -6.0844; %# exchange-correlation potential in eV;
eF   = 0.0; %# Fermi level in eV;
dE_hedin = 0; %# Hedin shift = Sigma(kF, elda_F) - Vxc(kF) in eV

ts  = [-60: 0.1 : 60]'; %# grid of times (plot ts versus expC as check)
wsc = [-20 :0.1: 20]'; %# freq grid for gwc spectral function (eV)
useEqp = 1; %# use Eqp in cumulant from Dyson eq.
do_interpolate = 0; %# interpolate sigma onto finer freq grid
Ninter = 10; %# fine grid spacing; used when do_interpolate = 1
zeroalpha = 0; %# set asymmetry to zero
fitrange  = 5; %# deal with w=0 in C(w)
qp_thresh = 1; %# helps finding the qp-peak if the satellite peak is higher
useCor = 0; %# use corrected ImSig to avoid negative spectral function