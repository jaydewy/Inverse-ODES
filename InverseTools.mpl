InverseTools := module()
description "tools for solving inverse ODE problems"
option package;
export InverseODE, CompareInverseSoln;

InverseODE := proc(x, N, a, b, x_init := FAIL, t_init := 0)
# x is the target solution as a function of t
# N is the degree of the polynomial approximation desired
# (a,b) is the interval I on which to solve
# x_init is optional. It is the inital value x(0) = x0
# t_init is optional. It is such that x(t_init) = x_init
local Tx, t0, x0, dist, dD, pol, pol_temp, soln, eta, i, M;

pol := 0;
for i from 0 to N do
	pol_temp := pol;
	pol := pol_temp + (eta[i])*x^i
end do;

x0 := ifelse(evalb(x_init = FAIL), eta[N+1], x_init);

Tx := x0 + int(subs(t=s, pol), s=t_init..t); # Picard operator
dist := int((subs(t=s, x - Tx))^2, s=a..b); # Collage distance, function of parameters eta__i

M := ifelse(evalb(x_init = FAIL), N+1, N); # M depends on whether the problem has x0 constrained or variable
for i from 0 to M do
	dD[i] := diff(dist, eta[i]) = 0;
end do;

soln := solve([seq(dD[i], i=0..M)], {seq(eta[i], i=0..M)});

end proc;

CompareInverseSoln := proc(x, soln, a, b, x_init := FAIL)
# x is the target solution, as a function of t
# soln is the output of a call to InverseODE
# x_init is FAIL if the inverse problem was variable, otherwise is it the x(0) = x0 value
local g, g_temp, y, N, i;

N := ifelse(evalb(x_init = FAIL), numelems(soln) - 1, numelems(soln));
g := 0;
for i from 1 to N do
	g_temp := g;
	g := g_temp + rhs(soln[i])*y(t)^(i-1);
end do;
y := rhs(dsolve([diff(y(t),t) = g, y(0) = ifelse(evalb(x_init = FAIL), rhs(soln[N+1]), x_init)]));
evalf(sqrt(int((x - y)^2, t=a..b)));
end proc:

end module;
