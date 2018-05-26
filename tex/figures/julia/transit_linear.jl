using PyPlot
x = linspace(0,2*pi,1000); y = sin.(3*x + 4*cos.(2*x))
plot(x, y, color="red", linewidth=2.0, linestyle="--")
savefig("transit_linear.pdf")
