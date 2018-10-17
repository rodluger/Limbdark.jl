function test_transit_poly_gradient_case(u_n,test_name)

r0 = [0.01,.1,0.5,0.99,1.0,1.01,2.0,10.,100.0]; n_u = length(u_n)
#r0 = [0.01,.1,0.5,0.999999,1.0,1.000001,2.0,10.,100.0]; n_u = length(u_n)
#r0 = [1.0]; n_u = length(u_n)
#r0_name =["0.01","0.1","0.5","0.999999","1","1.000001","2","10","100"]
r0_name =["0.01","0.1","0.5","0.99","1","1.01","2","10","100"]
#r0_name =["1"]
nb = 50

epsilon = 1e-12; delta = 1e-3
dfdrbu = zeros(n_u+2)
dfdrbu_big = zeros(BigFloat,n_u+2)
label_name=["r","b","u_1","u_2","u_3","u_4","u_5","u_6","u_7","u_8","u_9","u_10","u_11","u_12","u_13"]
floor = 1e-20
if ~skip_plots
  fig,axes = subplots(3,3)
end
fgrid = zeros(length(r0),nb*13)
@testset "$test_name" begin
for i=1:length(r0)
  global r=r0[i]
  if r < 1.0
    bgrid = abs.([linearspace(0.0,epsilon,nb); linearspace(epsilon,delta,nb); linearspace(delta,r-delta,nb);
     -logarithmspace(log10(delta),log10(epsilon),nb) .+ r; linearspace(r-epsilon,r+epsilon,nb); logarithmspace(log10(epsilon),log10(delta),nb) .+ r;
     linearspace(r+delta,1-r-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (1-r); linearspace(1-r-epsilon,1-r+epsilon,nb);
     logarithmspace(log10(epsilon),log10(delta),nb) .+ (1-r); linearspace(1-r+delta,1+r-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (1+r);linearspace(1+r-epsilon,1+r-1e-13,nb)])
     nticks = 14
     xticknames=[L"$10^{-13}$",L"$10^{-12}$",L"$10^{-3}$",L"$r-10^{-3}$",L"$r-10^{-12}$",L"$r+10^{-12}$",L"$r+10^{-3}$",
     L"$1-r-10^{-3}$",L"$1-r-10^{-12}$",L"$1-r+10^{-12}$",L"$1-r+10^{-3}$",L"$1+r-10^{-3}$",L"$1+r-10^{-12}$",L"$1+r-10^{-13}$"]
  else
    bgrid = [r-1+1e-13; logarithmspace(log10(epsilon),log10(delta),nb) .+ (r-1); linearspace(r-1+delta,r-delta,nb);
     -logarithmspace(log10(delta),log10(epsilon),nb) .+ r; linearspace(r-epsilon,r+epsilon,nb); logarithmspace(log10(epsilon),log10(delta),nb) .+ r;
     linearspace(r+delta,r+1-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (r+1); r+1-1e-13]
     nticks = 8
     xticknames=[L"$r-1+10^{-13}$",L"$r-1+10^{-3}$",L"$r-10^{-3}$",L"$r-10^{-12}$",L"$r+10^{-12}$",L"$r+10^{-3}$",
     L"$r+1-10^{-3}$",L"$r+1-10^{-13}$"]
  end
  nbgrid = length(bgrid)
  tp_grad_grid = zeros(nbgrid,n_u+2)
  tp_grad_array= zeros(n_u+2)
  tp_grid = zeros(nbgrid)
  tp_grid_big = zeros(nbgrid)
  tp_grad_grid_num = zeros(nbgrid,n_u+2)
  tp_grad_grid_ana = zeros(nbgrid,n_u+2)
  tp_grad_grid_big = zeros(nbgrid,n_u+2)
  for j=1:nbgrid
    tp=transit_poly(r,bgrid[j],u_n)
    tp_grid[j]=tp
    fgrid[i,j]=tp
    # Now compute with BigFloat finite difference (tests whether derivatives are coded properly):
    tp,tp_grad_array =  transit_poly_grad_num(r,bgrid[j],u_n,1e-18)
    if bgrid[j] == 0.0
      # Assign df/db exactly:
      tp_grad_array[2] =  0.0
    end
    tp_grad_grid_num[j,:] .=tp_grad_array
    tp_grid_big[j,:] .=tp
    # Compute analytic derivatives:
    tp = transit_poly!(r,bgrid[j],u_n,dfdrbu)
    tp_grad_grid_ana[j,:] .=dfdrbu
    # Compute derivatives in BigFloat (tests precision of derivatives):
    tp_big = transit_poly!(big(r),big(bgrid[j]),big.(u_n),dfdrbu_big)
    tp_grad_grid_big[j,:] .=dfdrbu_big
    @test isapprox(dfdrbu,tp_grad_array,atol=1e-6)
    if !isapprox(dfdrbu,tp_grad_array,atol=1e-6)
      println("r: ",r," bgrid: ",bgrid[j])
    end
    @test isapprox(dfdrbu,dfdrbu_big,atol=1e-8)
  end
  if ~skip_plots
# Now, make plots:
    ax = axes[i]
    y = abs.(asinh.(tp_grid) .-asinh.(tp_grid_big)); mask = y .<= floor; y[mask] .=floor
    ax[:semilogy](y,lw=1,label="flux")
    for n=1:n_u+2
      y = abs.(asinh.(tp_grad_grid_ana[:,n]) .-asinh.(tp_grad_grid_big[:,n])); mask = y .<= floor; y[mask] .=floor
      ax[:semilogy](y,lw=1,label=string(label_name[n]," vs. big"))
      y = abs.(asinh.(tp_grad_grid_ana[:,n]) .-asinh.(tp_grad_grid_num[:,n])); mask = y .<= floor; y[mask] .=floor
      ax[:semilogy](y,lw=1,label=string(label_name[n]," vs. num"))
    end
    if mod(i,3) == 0
      ax[:set_xlabel]("b values")
    end
    ax[:set_ylabel]("Derivative Error")
    ax[:axis]([0,nbgrid,floor,1])
    ax[:set_xticks](nb*linearspace(0,nticks-1,nticks))
    ax[:set_xticklabels](xticknames,rotation=45,fontsize=6)
    ax[:set_title](string("r = ",r0[i]),fontsize=6)
    ax[:legend](loc="upper right",fontsize=6)
  end
end
end
return fgrid
end
