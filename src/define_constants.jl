
if VERSION >= v"0.7"
  if !(eval(:(@isdefined third_float)))
    const third_float = 1.0/3.0
    const third_big= big(1)/big(3)
    const twothird_float = 2.0/3.0
    const twothird_big = big(2)/big(3)
  end
else
  if !(isdefined(:third_float))
    const third_float = 1.0/3.0
    const third_big= big(1)/big(3)
    const twothird_float = 2.0/3.0
    const twothird_big = big(2)/big(3)
  end
end
