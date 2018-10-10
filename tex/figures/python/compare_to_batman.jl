using DelimitedFiles

include("../../../src/transit_poly_struct.jl")

function profile_transit_poly!(trans,flux,b)
    t1 = time_ns()
    for i=1:length(b)
      trans.b = b[i]
      flux[i] = transit_poly!(trans)
    end
    elapsed = time_ns() - t1
    return elapsed
end

b = readdlm("b.txt")
u_n = zeros(2)
u_n[1] = 0.4
u_n[2] = 0.26
trans = transit_init(0.1, 0.0, u_n, false)
flux = zeros(length(b))
println(1e-9 * profile_transit_poly!(trans, flux, b))
writedlm("flux.txt", flux)
