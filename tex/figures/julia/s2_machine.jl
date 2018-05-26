using PyPlot
x = linspace(0,2*pi,1000); y = cos.(3*x + 4*tan.(2*x))
plot(x, y, color="green", linewidth=2.0, linestyle="--")
savefig("s2_machine.pdf")
