function sort3(a,b,c)
if c > b
  tmp = c; c = b; b = tmp
end
if c > a
  tmp = a; a = c; c = tmp
end
if b > a
  tmp = b; b = a; a = tmp
end
println(a," ",b," ",c)
return
end
