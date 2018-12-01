function logarithmspace,x1,x2,n &$
return,exp(alog(x1)+alog(x2/x1)*dindgen(n)/(n-1)) &$
end
