# Plot the various c_n components:

z = linspace(0,1,1000)

using PyPlot

clf()
plot(z,ones(z),label="n=0")

plot(z,z,label="n=1")

plot(z,4.*z.^2-2.,label="n=2")
plot(z,5.*z.^3-3.*z,label="n=3")
plot(z,6.*z.^4-4.*z.^2,label="n=4")
plot(z,7.*z.^5-5.*z.^3,label="n=5")
legend(loc="upper left")
savefig("plot_components.pdf", bbox_inches="tight")

