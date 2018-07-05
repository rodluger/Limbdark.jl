# Plots angles:
nb = 1001; nr = 1001
b=linspace(0.0,2,nr)
r=linspace(0.0,2,nr)
lambda = zeros(Float64,nr,nb)
phi = zeros(Float64,nr,nb)
for i=1:nr, j=1:nr
  if b[i] <= 1-r[j]
    phi[j,i]=pi/2; lambda[i,j]=pi/2
  elseif abs(1-r[j]) < b[i] && b[i] < 1+r[j]
    sphi = (1-r[j]^2-b[i]^2)/(2*b[i]*r[j])
    if abs(sphi) <= 1
      phi[j,i]=asin(sphi)
    else
      cphi = (1-b[i]+r[j])*(1+b[i]-r[j])*(b[i]+r[j]-1)*(1+b[i]+r[j])/(2*b[i]*r[j])
      println("r: ",r[j]," b: ",b[i]," sin(phi): ",sphi," cos(phi): ",cphi)
      phi[j,i] = acos(cphi)
    end
    slam = (1+(r[j]+b[i])*(b[i]-r[j]))/(2*b[i])
    if abs(slam) > 1
      println("r: ",r[j]," b: ",b[i]," sin(lam): ",slam)
    end
    lambda[j,i]=asin(slam)
  end
end
using PyPlot
clf()
img=imshow(phi,interpolation="nearest",origin="lower",extent=[0,2,0,2])
xlabel("Impact parameter")
ylabel("Radius ratio")
colorbar(img)
#read(STDIN,Char)
clf()
img=imshow(lambda,interpolation="nearest",origin="lower",extent=[0,2,0,2])
xlabel("Impact parameter")
ylabel("Radius ratio")
colorbar(img)
