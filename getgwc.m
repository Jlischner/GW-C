function [wsc,Acum] = getgwc(dat,elda,vxc,dE_hedin,eta,ts,wsc,useEqp,do_interpolate,Ninter,zeroalpha,qp_thresh,useCor);

  %# make empty input file:
  fid = fopen("gwc_inp.m","w");
  fprintf(fid,"FakeInputFile=1;\n");
  fclose(fid);
  %# run GW+C code:
  gwc;
  
endfunction;