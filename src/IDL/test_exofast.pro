r = [0.1, 0.5, 0.9999, 1.0, 1.0001, 2.0, 10.0]
nr = 7 & nb = 10000 & u1 = 0.3 & u2 = 0.3
for i=0,nr-1 do begin &$
  b = (1+2*r[i])*(-1.+2.*(dindgen(nb)+0.5)/nb) &$
  exofast_occultquad,abs(b),u1,u2,r[i],muo1,mu0 &$
  plot,b,muo1 &$
endfor
