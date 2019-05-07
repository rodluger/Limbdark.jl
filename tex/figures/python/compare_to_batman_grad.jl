if VERSION >= v"0.7"
  using DelimitedFiles
end

include("../../../src/transit_poly_struct.jl")

function profile_transit_poly!(trans,flux,b)
    t1 = time_ns()
    for k=1:10
        for i=1:length(b)
          trans.b = b[i]
          flux[i] = transit_poly_g!(trans)
        end
    end
    elapsed = 1e-10 * (time_ns() - t1)
    return elapsed
end

b = readdlm("b.txt")
u = readdlm("u.txt")
u_n = zeros(length(u))
for i = 1:length(u)
    u_n[i] = u[i]
end
trans = transit_init(0.1, 0.0, u_n, true)
flux = zeros(length(b))
println(profile_transit_poly!(trans, flux, b))
