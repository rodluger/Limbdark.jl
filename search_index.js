var documenterSearchIndex = {"docs":
[{"location":"api.html#API","page":"API","title":"API","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"Modules = [Limbdark]\nOrder   = [:function, :type]","category":"page"},{"location":"api.html#Limbdark.compute_uniform!-Union{Tuple{Limbdark.Transit_Struct{T}}, Tuple{T}} where T<:Real","page":"API","title":"Limbdark.compute_uniform!","text":"compute_uniform(t)\n\nComputes the coefficient for the uniform disk case, s_0, given a TransitStruct instance t.\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.s2!-Union{Tuple{T}, Tuple{T,T,Vector{T}}} where T<:Real","page":"API","title":"Limbdark.s2!","text":"s2!(r,b,s2_grad)\n\nComputes the linear limb-darkening case, as well as the gradient, s2_grad=[ds_2/dr,ds_2/db] is a pre-allocated two-element array. Returns s2 and complete elliptic integrals, E(k) and  (E(m)-(1-m)K(m))/m, where m = k^2.\n\nExample\n\njulia> s2_grad=[0.0,0.0]\njulia> s2!(0.1,0.5,s2_grad)\n(2.067294367278038, 1.4726447391114554, 0.8111640472029077)\njulia> s2_grad\n2-element Array{Float64,1}:\n -0.53988  \n  0.0182916\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.s2-Union{Tuple{T}, Tuple{T,T}} where T<:Real","page":"API","title":"Limbdark.s2","text":"s2(r,b)\n\nCompute the integral of the linear limb-darkening case.\n\nExamples\n\njulia> s2(0.1,0.5)\n2.067294367278038\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.s2_ell-Union{Tuple{T}, Tuple{T,T}} where T<:Real","page":"API","title":"Limbdark.s2_ell","text":"s2_ell(r,b)\n\nCompute the integral of the linear limb-darkening case, s_2 from starry (which is S_1 in Agol & Luger) and also return two elliptic  integrals, E(k) and (E(m)-(1-m)K(m))/m, where m = k^2.\n\nExamples\n\njulia> s2_ell(0.1,0.5)\n(2.067294367278038, 1.4726447391114554, 0.8111640472029077)\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.sqarea_triangle-Union{Tuple{T}, Tuple{T,T,T}} where T<:Real","page":"API","title":"Limbdark.sqarea_triangle","text":"sqarea_triangle(a,b,c)\n\nFunction which computes sixteen times the square of the area of a triangle with sides a, b and c using Kahan method. How to compute (sixteen times the) area squared of triangle with  high precision (Goldberg 1991).\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.transit_init-Union{Tuple{T}, Tuple{T,T,Vector{T},Bool}} where T<:Real","page":"API","title":"Limbdark.transit_init","text":"transit_init(r,b,u_n,grad)\n\nInitialize and return a transit structure. This structure holds user-facing information about the transit as well as internal temporary variables used for computing the light curve. Pass an instance of this structure to transit_poly() to compute the actual flux.\n\nArguments\n\nr::Real: The radius of the occultor in units of the radius of the occulted body.\nb::Real: The (initial) impact parameter of the occultation.\nu_n::Array{Real,1}: The array of limb darkening coefficients.\ngrad::Bool: Compute the gradient of the flux as well?\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.transit_poly!-Union{Tuple{Limbdark.Transit_Struct{T}}, Tuple{T}} where T<:Real","page":"API","title":"Limbdark.transit_poly!","text":"transit_poly!(t)\n\nGiven a TransitStruct instance t, with Green's coefficients already initialized, computes a limb-darkened transit light curve with optional gradient with respect to r, b and u coefficients.\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.transit_poly!-Union{Tuple{T}, Tuple{T,T,Vector{T},Vector{T}}} where T<:Real","page":"API","title":"Limbdark.transit_poly!","text":"transit_poly!(r,b,u_n,dfdrbu)\n\nGiven a radius-ratio, impact parameter, vector of limb-darkening coefficients of size N, and pre-allocated vector for derivatives of size N+2, returns the flux (normalized to one for unocculted star) for a limb-darkened transit  light curve and its gradient with respect to the input parameters is returned.\n\nArguments\n\nr::Real: The radius of the occultor in units of the radius of the occulted body.\nb::Real: The (initial) impact parameter of the occultation.\nu_n::Array{Real,1}: The array of limb darkening coefficients.  Size N.\ndfdrbu::Array{Real,1}: Gradient of the flux with respect to all input parameters. Size is N+2\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.transit_poly-Union{Tuple{T}, Tuple{T,T,Vector{T}}} where T<:Real","page":"API","title":"Limbdark.transit_poly","text":"transit_poly(r,b,u_n)\n\nGiven a radius-ratio, impact parameter, vector of limb-darkening coefficients of size N, returns the flux (normalized to one for unocculted star) for a  limb-darkened transit.\n\nArguments\n\nr::Real: The radius of the occultor in units of the radius of the occulted body.\nb::Real: The (initial) impact parameter of the occultation.\nu_n::Array{Real,1}: The array of limb darkening coefficients.  Size N.\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.transit_poly_g!-Union{Tuple{Limbdark.Transit_Struct{T}}, Tuple{T}} where T<:Real","page":"API","title":"Limbdark.transit_poly_g!","text":"transit_poly_g!(t)\n\nGiven a TransitStruct instance t, with Green's coefficients already initialized, computes a limb-darkened transit light curve including gradient with respect to r, b and g coefficients.\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.transit_poly_g-Union{Tuple{Limbdark.Transit_Struct{T}}, Tuple{T}} where T<:Real","page":"API","title":"Limbdark.transit_poly_g","text":"transit_poly_g(t)\n\nGiven a TransitStruct instance t, computes a limb-darkened transit light  curve without gradient.\n\n\n\n\n\n","category":"method"},{"location":"api.html#Limbdark.Transit_Struct","page":"API","title":"Limbdark.Transit_Struct","text":"Transit_Struct\n\nStructure to hold arrays and other quantities for computing transit. \n\nMembers\n\nr::T: radius ratio\nb::T: impact parameter\nu_n::Array{T,1}: limb-darkening coefficients\nn::Int64: number of limb-darkening coefficients\nn_max::Int64: maximum number of terms in M_n\ng_n::Array{T,1}: Green's basis coefficients\nsn::Array{T,1}: Green's basis terms\nMn::Array{T,1}: integral over (k^2-sin^2(phi))^m2\nNn::Array{T,1}: integral over (k^2-sin^2(phi))^m2sin^2(phi)\ngrad::Bool: true for gradient; false for no gradient\ns2_grad::Array{T,1}: gradient of s2 with respect to (rb) - this is S_1 = s_2 from starry\ndsndr::Array{T,1}: derivatives of Green's basis with respect to r\ndsndb::Array{T,1}: derivatives of Green's basis with respect to b\ndgdu::Array{T,2}: derivatives g_n with respect to u_n\ndfdrb::Array{T,1}: derivative of flux with respect to r,b\ndfdg::Array{T,1}: derivative of flux with respect to g_n\ndfdu::Array{T,1}: derivative of flux with respect to u_n\njmax::Int64: maximum number of terms in series expansions of I_v and J_v\nMn_coeff::Array{T,3}: coefficients for series expansion of M_n\nNn_coeff::Array{T,2}: coefficients for series expansion of N_n\nninv::Array{T,1}: inverse of the integers n\nk2::T: k^2 = (1-(r-b)^2)(4br)\nk::T: k = sqrtk^2\nkc::T: k_c = sqrt1-k^2 (unless k  1, then it is sqrt1-1k^2)\nkck::T: k_c k\nkap0::T: kappa = sin^-1(k) (=kappa_0 in M&A)\npimkap1::T: pi-kappa_1\nsqarea::T: (1-(b-r)^2)((b+r)^2-1)\nkite_area2::T: 2A_kite\nEofk::T: E(k^2) is complete elliptic integral of first kind\nEm1mKdm::T: (E(m)-(1-m)K(m))m is complete elliptic integral with m=k^2\nonembmr2:: T: 1-(b-r)^2\nonembpr2:: T: 1-(b+r)^2\nonembmr2inv::T: 1(1-(b-r)^2)\nonemr2mb2::T: 1-r^2-b^2\nsqonembmr2::T: sqrt1-(b-r)^2\nfourbr::T: 4br\nsqbr::T: sqrtbr\nsqbrinv::T: 1sqrtbr\nfourbrinv::T: 1(4br)\nk2inv::T: 1k^2\nkc2::T: 1-k^2 (or 1-1k^2 if k  1)\nsqrt1mr2:: T: sqrt1-r^2\nden::T: 1(g_1 + 23 g_2)\nthird::T: 13\ntwothird:: T: 23\nsqr1mr::T: sqrtr(1-r) if r  1\n\n\n\n\n\n","category":"type"},{"location":"quickstart.html#Quick-start","page":"Quick start","title":"Quick start","text":"","category":"section"},{"location":"quickstart.html#A-simple-light-curve","page":"Quick start","title":"A simple light curve","text":"","category":"section"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"First, let's import Limbdark:","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"push!(LOAD_PATH, \"../../src/\") # hide\nimport Limbdark","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"Now let's define some quantities we'll use throughout these examples:","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"# Define an array of impact parameters\nnpts = 100\nb = zeros(npts)\nfor i = 1:npts\n    b[i] = 3.0 * ((i - 0.5) / float(npts) - 0.5)\nend\n\n# Define the occultor radius\nr = 0.1\n\n# Define the quadratic limb darkening coefficients\nu_n = [0.40, 0.26]\n\n# Define the struct to hold the transit info\ntrans = Limbdark.transit_init(r, 0.0, u_n, true)\nnothing # hide","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"To compute a light curve, we call the transit_poly!(trans) method:","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"using PyPlot # hide\nflux = zeros(npts)\nfor i = 1:npts\n    trans.b = abs(b[i])\n    flux[i] = Limbdark.transit_poly!(trans)\nend","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"Here's what the light curve looks like:","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"plot(b, flux, linewidth=2);\nxlabel(\"Impact parameter\"); # hide\nylabel(\"Flux\") # hide\nsavefig(\"example1.svg\"); # hide\nclose(); nothing # hide","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"(Image: )","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"We can also compute the derivatives with respect to r, b, and u_n:","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"dfdr = zeros(npts)\ndfdb = zeros(npts)\ndfdu1 = zeros(npts)\ndfdu2 = zeros(npts)\nfor i = 1:npts\n    trans.b = abs(b[i])\n    flux[i] = Limbdark.transit_poly!(trans)\n    dfdr[i] = trans.dfdrb[1]\n    dfdb[i] = trans.dfdrb[2]\n    dfdu1[i] = trans.dfdu[1]\n    dfdu2[i] = trans.dfdu[2]\nend","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"Here's what they all look like:","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"plot(b, dfdr, label=\"df/dr\")\nplot(b, dfdb, label=\"df/db\")\nplot(b, dfdu1, label=\"df/du1\")\nplot(b, dfdu2, label=\"df/du2\")\nxlabel(\"Impact parameter\"); # hide\nylabel(\"Derivative\") # hide\nlegend() # hide\nsavefig(\"example2.svg\") # hide\nclose(); nothing # hide","category":"page"},{"location":"quickstart.html","page":"Quick start","title":"Quick start","text":"(Image: )","category":"page"},{"location":"install.html#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"install.html","page":"Installation","title":"Installation","text":"Limbdark was developed with julia 0.7 but is fully compatible with julia 1.0. To install Limbdark, simply run the following from within julia:","category":"page"},{"location":"install.html","page":"Installation","title":"Installation","text":"import Pkg\nPkg.clone(\"https://github.com/rodluger/Limbdark.jl.git\")","category":"page"},{"location":"index.html#Limbdark.jl-Documentation","page":"Home","title":"Limbdark.jl Documentation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Fast, analytic, and numerically stable transit light curves for stars with arbitrary order limb darkening. Check out the latest draft of the paper here. The code is open source and available on GitHub.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"This documentation is under construction and is still being updated with examples, tutorials, and detailed API documentation. If you don't see what you're looking for, feel free to  open an issue.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Check out the Installation, Quick start, and API pages for more information.","category":"page"}]
}
