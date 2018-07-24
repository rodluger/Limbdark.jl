muc = 0.04

using PyPlot

mu = linspace(muc,1,1000)

plot(mu,(mu-muc)./sqrt.(1-2*muc*mu+muc^2),label=L"$\mu_\infty=(\mu-\mu_c)/\sqrt{1-2\mu_c\mu+\mu_c^2}$")
xlabel(L"$\mu$")
ylabel(L"$\mu_\infty$")
plot([0,1],[0,1],linestyle="--",label=L"$\mu_\infty=\mu$")
axis([0,1,0,1])
plot(mu,(mu-muc)/(1-muc),linestyle=":",label=L"$\mu_\infty=(\mu-\mu_c)/(1-\mu_c)$")
legend(loc="upper left")
