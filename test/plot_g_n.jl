import Limbdark

nu = 20
u=ones(nu)/nu

r=0.1
b=0.5

t = Limbdark.transit_init(r,b,u,false)

using PyPlot

plot(t.g_n)
xlabel(L"n")
ylabel(L"$g_n$")

savefig("g_n_vs_n.pdf",bbox_inches="tight")
