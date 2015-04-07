function y=logistic(a,x);
% this takes as input the parameter a and the value of x and
% returns y = ax(1-x). The value of x can be a vector, in which
% case, the corresponding y vector is returned
% Usage: >> x = 0:.01:1; (100 equidistant points in [0, 1])
% >> y = logistic(3.2, x);
y=a*x.*(1-x);

