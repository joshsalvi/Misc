
function [u, v, z] = iterate(a, y0, n)
% Outputs everything needed for graphical analysis of the logistic map; the
% inputs are the parameter a for the logistic map, the starting value y0,
% and the number of iterations for the function to do
% Usage:: to interate 500 times with IC y0 = .32, and parameter .5, we type
% >> [u,v,z] = iterate(0.5, 0.32, 500); to plot the 400th to 450th
% points of the orbit, you?d do
% >> plot(z(400:450)). To do the "web plot", with the last 50 points, type
% >> plot(u(900:1000),v(900:1000))
z = zeros(n+1,1);
z(1) = y0;
for i = 2:n
z(i) = logistic(a,z(i-1));
end
u = zeros(2*n,1);u(1) = y0;
v = zeros(2*n,1);v(1) = 0;
for i=2:2*n
u(i) = z(ceil(i/2));
end
v = [0;u(3:2*n)];
u = u(1:2*n-1);
