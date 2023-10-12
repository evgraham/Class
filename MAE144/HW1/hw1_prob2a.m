% Solving Problem 2A of HW1
% Setting up polynomials of transfer function G(s)
a = RR_poly([1 0 -46 0 369 0 -324]); b = RR_poly([1 0 -29 0 100]); 
% Target denominator polynomial
f = RR_poly([1 20 154 576 1089 972 324]);
% Calling RR_diophantine to find the solution with smallest order for y(s)
[x, y] = RR_diophantine(a, b, f);
test=trim(a*x+b*y);
residual1=norm(f-test);
