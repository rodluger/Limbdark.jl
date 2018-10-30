include("transit_poly_struct.jl")
trans = transit_init(0.1,0.5,[0.5,0.5],false)
#trans = transit_init(0.1,0.5,[0.5,0.25,0.125,0.0625,0.03125],false)
#Transit_Struct{Float64}(0.1, 0.5, [0.5, 0.5], 2, 3, [-0.25, 1.5, -0.125], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], false, [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0 0.0; 0.0 0.0; 0.0 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0])

function transit_b(trans,b)
#  flux = zeros(length(b))
#  gradient = zeros(length(b),2+trans.n)
  for i=1:length(b)
    trans.b = b[i]
#    flux[i]=transit_poly!(trans)
#    gradient[i,:]=trans.dfdrbu
    transit_poly!(trans)
  end
#  return flux,gradient
  return
end
b = abs.(linspace(-1.2,1.2,10000000))

@time transit_b(trans,b);
@time transit_b(trans,b);

#=
Profile.clear()
Profile.init(n=10^7,delay=0.01)
@profile transit_b(trans,b)
Profile.print()
=#
