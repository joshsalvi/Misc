%% Create a heart
clear all; close all;

%% *Two-Dimensional Heart*
% This will create a heart in two dimensions.

t = linspace(-500,500,20000);

x1 = 16*sin(t).^3;
y1 = 13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - 4*cos(t);

figure(1);
plot(x1,y1,'r*'); axis([-30 30 -20 20]);

%% *Two-Dimensional Heart, Oscillations*
% This will create a heart in two dimensions from a series of oscillations.

x2=[-2:.001:2];
y2=(sqrt(cos(x2)).*cos(200*x2)+sqrt(abs(x2))-0.7).*(4-x2.*x2).^0.01;

figure(2);
plot(x2,y2,'r');



%% *Three-Dimensional Heart*
% This will create a heart in three dimensions.

v = -1.5:.02:1.5;
[x3,y3,z3]=ndgrid(v,v,v);
w=(2*x3.^2+y3.^2+z3.^2-1).^3-(1/10)*x3.^2.*z3.^3-y3.^2.*z3.^3;

figure(3);
p=patch(isosurface(w,0));
set(p,'facecolor','r','edgecolor','none')
camlight; lighting phong
view(160,15)

