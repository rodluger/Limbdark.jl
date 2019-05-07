5/6/2019

This directory gives an initial exploration of fitting the
"power-2" limb-darkening model (Maxted 2018; Maxted & Gill
2019) as a function of the polynomial limb-darkening model
order.

The power-2 model is a power-law, I(mu) = 1-c_alpha(1-mu^alpha),
which has the advantage like quadratic limb-darkening of
only having two parameters, but it does a better job of
fitting stellar models near the limb.  A high-order polynomial
can do a decent job of approximating different power-laws,
but it degrades in precision near the limb.

In principle the polynomial law could be used as a stand-in
for the power-2 law, which has yet to be explored.
