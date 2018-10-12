module Random
  export seed!
  function seed!(n)
    return MersenneTwister(n)
  end
end
