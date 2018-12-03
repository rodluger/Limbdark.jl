pro test_exofast
r = [0.01, 0.1, 0.5, 0.9999, 1.0, 1.0001, 2.0, 10.0, 100.0]
nr = 9 &  u1 = 0.3 & u2 = 0.3
epsilon = 1d-14 & delta = 1d-3

nb = 10000
for i=0,nr-1 do begin &$
;  b = (1+2*r[i])*(-1.+2.*(dindgen(nb)+0.5)/nb) &$
  if r[i] lt 1.0 then begin &$
    b= [linearspace(1e-15,epsilon,nb), linearspace(epsilon,delta,nb), linearspace(delta,r[i]-delta,nb),$
     -logarithmspace(delta,epsilon,nb) + r[i], linearspace(r[i]-epsilon,r[i]+epsilon,nb), logarithmspace(epsilon,delta,nb) + r[i],$
     linearspace(r[i]+delta,1-r[i]-delta,nb), -logarithmspace(delta,epsilon,nb) + (1-r[i]), linearspace(1-r[i]-epsilon,1-r[i]+epsilon,nb),$
     logarithmspace(epsilon,delta,nb) + (1-r[i]), linearspace(1-r[i]+delta,1+r[i]-delta,nb), - logarithmspace(delta,epsilon,nb) + (1+r[i]),$
     linearspace(1+r[i]-epsilon,1+r[i]-1e-15,nb)] &$
  endif else begin &$
    b= [linearspace(r[i]-1+1e-15,r[i]-1+epsilon,nb), logarithmspace(epsilon,delta,nb) + (r[i]-1), linearspace(r[i]-1+delta,r[i]-delta,nb),$
     -logarithmspace(delta,epsilon,nb) + r[i], linearspace(r[i]-epsilon,r[i]+epsilon,nb), logarithmspace(epsilon,delta,nb) + r[i], $
     linearspace(r[i]+delta,r[i]+1-delta,nb), -logarithmspace(delta,epsilon,nb) + (r[i]+1), linearspace(r[i]+1-epsilon,r[i]+1-1e-15,nb)] &$
  end &$
  TIC & exofast_occultquad,abs(b),u1,u2,r[i],muo1,mu0 & TOC &$
  TIC & exofast_occultquad_original,abs(b),u1,u2,r[i],muo1_original,mu0 & TOC &$
  off = max([r[i],1d0])*1d-2 &$
  plot,b,muo1,ys=1,xr=[min(b)-off,max(b)+off],xs=1 &$
  oplot,b,muo1_original,color=255,linestyle=2 &$
  print,'r: ',r[i],' max deviation/min(r^2,1): ',max(muo1-muo1_original)/min([r[i]^2,1d0]) &$
  c=get_kbrd(1)&$
  plot,b,muo1-muo1_original,ys=1,xr=[min(b)-off,max(b)+off],xs=1 &$
  c=get_kbrd(1)&$
endfor

return
end
