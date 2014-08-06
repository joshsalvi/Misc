%%% Reminder: Normal Forms of Supercritical Bifurcations.%%%%%%%%%%%
%%                                                                %%
%%  (1): dx/dt= mu - x^2 : Node-Saddle                            %%
%%  (2): dx/dt= mu x - x^2 : Transcritical (Stability Exchange)   %%
%%  (3): dx/dt= mu x - x ^3 : Pitchfork                           %%
%%  (4): dz/dt= (mu + i gamma)z - z|z|^2 (z=complex) : Hopf       %%
%%                                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example of a Nonlinear Dynamical System (in which a Hopf 
% bifurcation can occur for appropriate values of parameters).
% The dynamical system is defined as follows:
% (1): dx/dt=[µ-(x^2+y^2)]x-gamma y, 
% (2): dy/dt=gamma x+[µ-(x^2+y^2)]y.
% In a general way, it is said that: 
% . The stationary solution is z=0 (i.e. x=y=0). 
% . There exists another solution such that |z|^2 is independent from time,
% i.e. |z|^2=x^2+y^2=µ, this condition defines the equation of a circle (in
% the plane) whose radius is (µ)^1/2.
%
% The routine attempts to depict (in a discrete way) the orbit for -0.5<µ<3 
% while gamma varies from -1 to 1. 
%
% ATTENTION: Considered as an exercise the following question: 
% IS THIS DISCRETE SIMULATION ("Dynamic2.m") VALID OR NOT? 


clear all

figure(1)

disp(' ')
disp('PLEASE, ENLARGE FIGURE, THEN, PRESS ANY KEY TO CONTINUE,')
disp(' ')
disp('TO STOP SIMULATION, TYPE "CTRL-C" IN THE COMMAND WINDOW,')
disp(' ')

pause

s = 1000; 
k = 700;
x0 = rand(1,s); 
y0 = rand(1,s);
A = linspace(-0.5,3,s); 
x = zeros(k,s);
x(1,:) = x0;
y = zeros(k,s);
y(1,:) = y0;
for v=0:0.01:2.02;
    gamma=v-1;
    for s = 1:k;          
    x(s+1,:) = (A-((x(s,:).^2)+(y(s,:).^2))).*x(s,:)-gamma.*y(s,:);
    y(s+1,:) = gamma.*x(s,:)+(A-((x(s,:).^2)+(y(s,:).^2))).*y(s,:);       
    end     

    pause(0.00000001) 

Ha=ones(1000,1);
VX=Ha*A;
SX=VX(1:701,:);
NX=SX(467:701,:);
plot3(NX,x(floor(2.*end/3):end,:),y(floor(2.*end/3):end,:), '.','markersize',1);
view(-15,29)
axis([-0.5 3 -1.8 1.8 -1.8 1.8]);
title(['Orbit of the system: dx/dt=[\mu-(x^2+y^2)]x-\gammay, dy/dt=\gammax+[\mu-(x^2+y^2)]y with \gamma= ', num2str(gamma)]);
xlabel('\mu-abscissa: -0.5 \leq \mu \leq 3');
ylabel('X');
zlabel('Y');
grid;

%box;

end


% "Complex and Chaotic Nonlinear Dynamics. 
% Advances in Economics and Finance, 
% Mathematics and Statistics" 
% T.Vialar, Springer 2009.
% Copyright(c).