% Solving Problem 2B of HW1
% Setting up polynomials of transfer function G(s)
a = RR_poly([1 0 -46 0 369 0 -324]); b = RR_poly([1 0 -29 0 100]); 
% Target denominator polynomial
f = RR_poly([1 20 154 576 1089 972 324]);
% Calling RR_diophantine to find the solution with smallest order for y(s)
[x, y] = RR_diophantine(a, b, f);
k=0; % Iteration counter
% Using a while loop to add a pole at s=-20 until the answer is proper

% Condition for semi-proper D(s): while y.n > x.n
while y.n >= x.n
    f = f.poly;
    % Adding an additional pole for every iteration
    f = RR_poly(conv(f, [1 20]));
    k=k+1;
    [x, y] = RR_diophantine(a, b, f);
end
test=trim(a*x+b*y);
residual1=norm(f-test);
