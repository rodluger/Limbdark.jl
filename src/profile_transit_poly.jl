# Timing tests of transit_poly_g!:

include("transit_poly_struct.jl")
#trans = transit_init(0.1,0.5,[0.5,0.5],true)
#u_n = [0.5,0.25,0.125,0.0625,0.03125]
u_n = [0.5,0.25]
#trans = transit_init(0.1,0.5,u_n,true)
#n = 6
#u_n = 2.^(-linspace(1,n,n))
#trans = transit_init(0.1,0.5,[0.5,0.25,0.125,0.0625,0.03125],false)
trans = transit_init(0.1,0.5,u_n,true)
#Transit_Struct{Float64}(0.1, 0.5, [0.5, 0.5], 2, 3, [-0.25, 1.5, -0.125], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], false, [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0])

function transit_b!(trans::Transit_Struct{T},b::Array{T,1},flux::Array{T,1}) where {T <: Real}
#  flux = zeros(length(b))
#  gradient = zeros(length(b),2+trans.n)
  for i=1:length(b)
    trans.b = b[i]
#    flux[i]=transit_poly_g!(trans)
#    gradient[i,1:2]=trans.dfdrb
#    gradient[i,3:2+trans.n]=trans.dfdu
    if trans.grad
      flux[i]=transit_poly_g!(trans)
    else
      flux[i]=transit_poly_g(trans)
    end
  end
#  return flux,gradient
  return
end

nb = 10000000
b = abs.(linspace(-1.2,1.2,nb))
flux = zeros(nb)
@time transit_b!(trans,b,flux);
@time transit_b!(trans,b,flux);

#=
Profile.clear()
Profile.init(n=10^7,delay=0.01)
@profile transit_b(trans,b)
Profile.print()
=#
