function [z] = iterate2(a, y0, n)
% this one just computes the orbit, so we don?t need u and v of iterate.m
z = zeros(n+1,1);
z(1) = y0;
for i = 2:n
z(i) = logistic(a,z(i-1));
end
