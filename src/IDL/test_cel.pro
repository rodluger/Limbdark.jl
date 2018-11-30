nk2 = 10
k2 = randomu(seed,nk2)
kc = sqrt(1d0-k2)
p = 2*randomu(seed,nk2)-1d0  & ps  = p & ps2 = p
a1 = 2*randomu(seed,nk2)-1d0 & a1s = a1 & a1s2 = a1
a2 = 2*randomu(seed,nk2)-1d0 & a2s = a2 
a3 = 2*randomu(seed,nk2)-1d0 & a3s = a3
b1 = 2*randomu(seed,nk2)-1d0 & b1s = b1 & b1s2 = b1
b2 = 2*randomu(seed,nk2)-1d0 & b2s = b2
b3 = 2*randomu(seed,nk2)-1d0 & b3s = b3
f1 = dblarr(nk2)
f2 = dblarr(nk2)
f3 = dblarr(nk2)
cel_bulirsch,k2,kc,ps,a1s,a2s,a3s,b1s,b2s,b3s,f1,f2,f3

; Now, compute these numerically:

nphi = 1000 & dphi = 0.5*!dpi/nphi
phi = dphi*(dindgen(nphi)+0.5)
c2 = cos(phi)^2 & s2 = 1-c2
for i=0,nk2-1 do begin & $
  if p[i] le 0.0 then begin & $
    p0 = sqrt((1d0-k2[i]-p[i])/(1d0-p[i])) & $
    a0 = (a1[i]-b1[i])/(1d0-p[i]) & $
    b0 = -(b1[i]-a1[i]*p[i])/(1d0-p[i])^2*k2[i]/p0+a0*p0 & $
    a1[i] = a0 & $
    b1[i] = b0*p0 & $
    p[i] = p0*p0 & $
  endif & $
  int1 = total(dphi*(a1[i]*c2+b1[i]*s2)/(c2+p[i]*s2)/sqrt(1d0-k2[i]*s2)) & $
  int2 = total(dphi*(a2[i]*c2+b2[i]*s2)/sqrt(1d0-k2[i]*s2)) & $
  int3 = total(dphi*(a3[i]*c2+b3[i]*s2)/sqrt(1d0-k2[i]*s2)) & $
  print,ps2[i],int1,f1[i],int1/f1[i],int2,f2[i],int2/f2[i],int3,f3[i],int3/f3[i] & $
endfor
