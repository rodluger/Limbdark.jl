# This is performance and accuracy tests for the integrate_lightcurve function

using BenchmarkTools
using PyPlot
using PyCall

# Python modules used for comparison
exo = pyimport("exoplanet")
tt = pyimport("theano.tensor")
starry = pyimport("starry")

include("../../../src/integrate_lightcurve.jl")
include("../../../test/loglinspace.jl")

println("Running comparison...")

fig, axes = subplots(3,1, figsize=(8,12))
fig.subplots_adjust(hspace=.4)

topplot = axes[1]
middleplot = axes[2]
perfplot = axes[3]

function get_time_array(dataPoints::Real)
  #trans_dur = 2*sqrt((1+r)^2-b0^2)/v*24 # approx transit duration in hours
  return linearspace(0.0 - 0.5*get_duration_days(), 0.5*get_duration_days(), dataPoints)
end

function get_duration_days()
  return 2.0
end

function get_exposure_hours(dataPoints::Real)
  return 2/60 # TODO is it OK if this doesn't exactly cover the space of time of the duration? Overlap ok? Gaps ok?
  #return get_duration_days() * 24 / dataPoints # this would make it gapless without overlap always
end

function run_ALFM_test(r::Float64, b0::Float64, s::Float64, period::Float64, u_1::Float64, u_2::Float64, dataPoints::Int64, trials::Int64, graph::PyObject=nothing)

  time = get_time_array(dataPoints)

  expression, timeTaken, bytes, gctime, memallocs = @timed for k=1:trials

    # Setup all the parameters
    trans = transit_init(r, b0, [u_1, u_2], false)  # Initialize transit structure
    param = [0.0, s, b0]   # parameters of the transit: [t_0, v, b_0]
    t = Array{Float64,1}(undef,dataPoints)
    t .= time
    dt = get_exposure_hours(dataPoints)
    favg1 = Array{Float64,2}(undef, 7, dataPoints)
    tol = 1e-6
    maxdepth = 6
    neval_t = Array{Int64,1}(undef, dataPoints)
    depthmax = Array{Int64,1}(undef, dataPoints)

    integrate_lightcurve!(trans, param, t, dt, favg1, dataPoints, tol, maxdepth, neval_t, depthmax)

    if (graph != nothing && k==1)
      graph.plot(t, favg1[1,:], label="ALFM (n=$dataPoints)", linestyle = "--", lw=2)
    end
  end

  return timeTaken/trials
end

function run_exoplanet_test(orbit::PyObject, r::Ty, u_1::Ty, u_2::Ty, dataPoints::Int64, trials::Int64, graph::PyObject=nothing) where {Ty <: Real}

  t = get_time_array(dataPoints)
  texp = get_exposure_hours(dataPoints)

  expression, timeTaken, bytes, gctime, memallocs = @timed for k=1:trials

    lc = exo.StarryLightCurve([u_1, u_2]).get_light_curve(orbit=orbit, r=r, t=t, texp=texp, order=2).eval()

    if (graph != nothing && k==1)
      lc = lc.+1.0 # move relative flux to be on same scale as ALFM
      graph.plot(t, lc, label="exoplanet (n=$dataPoints)", linestyle = "-", lw=1)
    end
  end

  return timeTaken/trials
end

# Transit parameters
datapoints_vals = collect(logarithmspace(1.0, 6.0, 15))
period = 1000.0 # days
r = 0.1 # radius of planet in units of stellar radius TODO varying this produces unexpected results from limbdark?
b = 0.2 # impact parameter

# Quadratic limb darkening parameters
u_1 = 0.2; u_2 = 0.2

# Comparison configuration
runs = length(datapoints_vals)
trialsPerConfig = 10 # number of trials to run for a given parameter configuration

# Structures for holding results
num_datapoints = zeros(runs)
limbdark_time = zeros(runs)
exoplanet_time = zeros(runs)

# with these parameters, it will assume solar mass and solar radius
orbit = exo.orbits.KeplerianOrbit(period=period, b=b, ecc=0.0, omega=0.0)
#a = orbit.a.eval()[1] # in units of R_sun

# speed differences around this time are trivial (such as at 1.0)
v0 = orbit.get_planet_velocity(0.0)
s0 = sqrt((v0[1].eval()[1]^2.0) + (v0[2].eval()[1]^2.0) + (v0[3].eval()[1]^2.0))
#println("  Speed: $s0 solar_radii / day")

# which iteration to graph the light curves (higher is more data points)
runToGraph = 7

# Run the timing test
for i = 1:runs
  #period = period_vals[mod(i-1, length(period_vals)) + 1]
  dataPoints = convert(Int64, round(datapoints_vals[i]))

  graph_to_use::PyObject = nothing
  if (i == runToGraph)
    graph_to_use = topplot
  end

  num_datapoints[i] = dataPoints
  limbdark_time[i] = run_ALFM_test(r, b, s0, period, u_1, u_2, dataPoints, trialsPerConfig, graph_to_use)
  exoplanet_time[i] = run_exoplanet_test(orbit, r, u_1, u_2, dataPoints, trialsPerConfig, graph_to_use)
end

# extra test for different period, for debugging
r2 = 0.02
dataPoints = convert(Int64, round(datapoints_vals[runToGraph]))
run_ALFM_test(r2, b, s0, period, u_1, u_2, dataPoints, trialsPerConfig, middleplot)
run_exoplanet_test(orbit, r2, u_1, u_2, dataPoints, trialsPerConfig, middleplot)

# Prepare labels for the light curve comparisons (for accuracy debugging)
topplot.set_title("Light Curves (r=$r; b=$b; T=$period d)")
topplot.set_xlabel("Time [days]")
topplot.set_ylabel("Relative Flux")
topplot.legend(loc="lower right",fontsize=9)

middleplot.set_title("Light Curves (r=$r2; b=$b; T=$period d)")
middleplot.set_xlabel("Time [days]")
middleplot.set_ylabel("Relative Flux")
middleplot.legend(loc="lower right",fontsize=9)

# Prepare plot for comparing performance
perfplot.plot(num_datapoints, limbdark_time, label = "ALFM (2019)", linestyle = "-", lw=2, color="C0")
perfplot.plot(num_datapoints, exoplanet_time, label = "exoplanet", linestyle = "-", lw=2, color=:red)

perfplot.set_title("Time Integration Comparison (r=$r; b=$b)")
perfplot.set_xlabel("Data Points")
perfplot.set_ylabel("Time [s]")
perfplot.set_xscale(:log)
perfplot.set_yscale(:log)
perfplot.legend(loc="lower right",fontsize=9)
