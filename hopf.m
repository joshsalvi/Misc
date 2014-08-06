function xdot=hopf(t,x,mu)
% Example: Hopf Bifurcation

global mu 

xdot(1)=x(2)+x(1)*(mu-(x(1)^2)-(x(2)^2));
xdot(2)=-x(1)+x(2)*(mu-(x(1)^2)-(x(2)^2));
xdot=xdot';



% "Complex and Chaotic Nonlinear Dynamics. 
% Advances in Economics and Finance, 
% Mathematics and Statistics" 
% T.Vialar, Springer 2009
% Copyright(c).
