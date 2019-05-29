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

topplot = axes[1]
middleplot = axes[2]
perfplot = axes[3]

# Duration is observed duration in hours
# Returns a time array in units of hours
function get_time_points(duration::Real)
  # Simulate transit data with a 2-minute cadence:
  tobs = 3 * duration  # Observation duration in hours
  nobs = convert(Int64,round(tobs * 30))  # Number of data points
  return linearspace(-1.5*duration, 1.5*duration, nobs)  # Time array in units of hours
end

function get_time_array(dataPoints::Real)
  #trans_dur = 2*sqrt((1+r)^2-b0^2)/v*24 # approx transit duration in hours
  return linearspace(-2.0, 2.0, dataPoints)
end


function run_ALFM_test(r::Float64, b0::Float64, aonr::Float64, period::Float64, u_1::Float64, u_2::Float64, dataPoints::Int64, trials::Int64, plot_curve::Bool)

  v = 2*pi*aonr/period # velocity
  println("ALFM velocity estimate: $v")

  time = get_time_array(dataPoints)
  nobs = length(time)

  expression, timeTaken, bytes, gctime, memallocs = @timed for k=1:trials

    # Setup all the parameters
    trans = transit_init(r, b0, [u_1, u_2], false)  # Initialize transit structure
    param = [0.0, v/24, b0]   # parameters of the transit: [t_0,v,b_0]  # TODO fixed speed is different than exoplanet
    t = Array{Float64,1}(undef,nobs)
    t .= time
    dt = 2/60  # Transit exposure time is 2 minutes, convert to hours.
    favg1 = Array{Float64,2}(undef, 7, nobs)
    tol = 1e-6
    maxdepth = 6
    neval_t = Array{Int64,1}(undef, nobs)
    depthmax = Array{Int64,1}(undef, nobs)

    integrate_lightcurve!(trans, param, t, dt, favg1, nobs, tol, maxdepth, neval_t, depthmax)

    if (plot_curve && k==1)
      topplot.plot(t, favg1[1,:], label="ALFM (n=$nobs)", linestyle = "-", lw=1)
    end
  end

  return timeTaken/trials
end

function run_exoplanet_test(orbit::PyObject, r::Ty, u_1::Ty, u_2::Ty, dataPoints::Int64, trials::Int64, plot_curve::Bool) where {Ty <: Real}

  t = get_time_array(dataPoints)
  nobs = length(t)

  expression, timeTaken, bytes, gctime, memallocs = @timed for k=1:trials

    lc = exo.StarryLightCurve([u_1, u_2]).get_light_curve(orbit=orbit, r=r, t=t, texp=2.0/60.0).eval()

    if (plot_curve && k==1)
      lc = lc.+1.0 # move relative flux to be on same scale as ALFM
      topplot.plot(t, lc, label="exoplanet (n=$nobs)", linestyle = "-", lw=1)
    end
  end

  return timeTaken/trials
end

# Transit parameters. Some vary, some are fixed

datapoints_vals = collect(logarithmspace(1.0, 4.5, 15)) # orbital period (~2.5 to 10000 days)
#period_vals = collect(logarithmspace(0.4, 3.0, 25)) # orbital period (~2.5 to 10000 days)
period = 1000.0 # days
r = 0.1 # radius of planet in units of stellar radius
b = 0.3 # impact parameter

# Quadratic limb darkening parameters
u_1 = 0.2; u_2 = 0.2

# Comparison configuration
runs = length(datapoints_vals)
trialsPerConfig = 3 # number of trials to run for a given parameter configuration

# Structures for holding results
num_datapoints = zeros(runs)
limbdark_time = zeros(runs)
exoplanet_time = zeros(runs)

# Run the timing test
for i = 1:runs
  #period = period_vals[mod(i-1, length(period_vals)) + 1]
  dataPoints = convert(Int64, round(datapoints_vals[i]))

  # with these parameters, it will assume solar mass and solar radius
  # TODO that probably isn't consistent with what we did above
  orbit = exo.orbits.KeplerianOrbit(period=period, b=b, ecc=0.0, omega=0.0)

  # TODO calculate velocity at t=0 and compare to the ALFM calculation. should be very close.
  a = orbit.a.eval()
  aonr = a / r

  v0 = orbit.get_planet_velocity(0.0)
  s0 = sqrt((v0[1].eval()[1]^2.0) + (v0[2].eval()[1]^2.0) + (v0[3].eval()[1]^2.0))
  println("velocity from exoplanet: $s0")
  # TODO make sure units are the same. Is there a discrepancy?


  plot_it = 3 <= i <= 3

  num_datapoints[i] = dataPoints
  limbdark_time[i] = run_ALFM_test(r, b, aonr, period, u_1, u_2, dataPoints, trialsPerConfig, plot_it)
  exoplanet_time[i] = run_exoplanet_test(orbit, r, u_1, u_2, dataPoints, trialsPerConfig, plot_it)

end

topplot.set_xlabel("Time")
topplot.set_ylabel("Relative Flux")
topplot.legend(loc="upper left",fontsize=10)

middleplot.set_xlabel("Time")
middleplot.set_ylabel("Relative Flux")
middleplot.legend(loc="upper left",fontsize=10)

# Prepare plot for comparing performance
perfplot.plot(num_datapoints, limbdark_time, label = "ALFM (2019)", linestyle = "-", lw=2, color="C0")
perfplot.plot(num_datapoints, exoplanet_time, label = "exoplanet", linestyle = "-", lw=2, color=:red)

perfplot.set_title("Time Integration Comparison (r=0.1; b=0.3)")
perfplot.set_xlabel("Data Points")
perfplot.set_ylabel("Time [s]")
perfplot.set_xscale(:log)
perfplot.set_yscale(:log)
perfplot.legend(loc="upper left",fontsize=10)
