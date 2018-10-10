s_2,s2_gradient= s2_grad_func(r,b)
s2_big,s2_grad_numeric = s2_grad_num(r,b)
s2_grad_ana = zeros(2)
diff1 = s2_grad_numeric-s2_gradient
s_2=s2!(r,b,s2_grad_ana)
diff2 = s2_grad_ana-s2_gradient
println("b : ",b," r: ",r," flux diff: ",s_2-s2_big," grad diff(num-auto): ",diff1," diff(ana-auto): ",diff2 )
println("gradient: ",s2_gradient," expected: ",-[2,-2./3.])
@test isapprox(s2_grad_numeric[1],s2_gradient[1],atol=1e-8)
@test isapprox(s2_grad_numeric[2],s2_gradient[2],atol=1e-8)
@test isapprox(s2_grad_ana[1],s2_gradient[1],atol=1e-8)
@test isapprox(s2_grad_ana[2],s2_gradient[2],atol=1e-8)
