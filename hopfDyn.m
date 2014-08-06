% Example of Supercritical Hopf bifurcation
% From a stable node to a stable limit cycle shown in the Phase-plane.
% Phase-plane plots of the components [X(1),X(2)] from µ=-1 to µ=2.

%%% Reminder: Normal Forms of Supercritical Bifurcations.%%%%%%%%%%%
%%                                                                %%
%%  (1): dx/dt= µ - x^2 : Node-Saddle                             %%
%%  (2): dx/dt= µ x - x^2 : Transcritical (Stability Exchange)    %%
%%  (3): dx/dt= µ x - x ^3 : Pitchfork                            %%
%%  (4): dz/dt= (µ + i gamma)z - z|z|^2 with z=complex : Hopf     %%
%%                                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

disp('To stop simulation, Press "ctrl-c" in the command window')


global mu

for mu=-1:0.1:2;
    tfinal=100;
    
% First set of initial values 
x0a=[-1.5 -1.5];
x0b=[0 -1.5];
x0c=[1.5 -1.5];
x0d=[1.5 0];
x0e=[1.5 1.5];
x0f=[0 1.5];
x0g=[-1.5 1.5];
x0h=[-1.5 0];

% Initial values close to the origin (0,0).
x0ia=[-0.001 -0.001];
x0ib=[0 -0.001];
x0ic=[0.001 -0.001];
x0id=[0.001 0];
x0ie=[0.001 0.001];
x0if=[0 0.001];
x0ig=[-0.001 0.001];
x0ih=[-0.001 0];

[t xa]=ode23(@hopf,[0 tfinal],x0a);
[t xb]=ode23(@hopf,[0 tfinal],x0b);
[t xc]=ode23(@hopf,[0 tfinal],x0c);
[t xd]=ode23(@hopf,[0 tfinal],x0d);
[t xe]=ode23(@hopf,[0 tfinal],x0e);
[t xf]=ode23(@hopf,[0 tfinal],x0f);
[t xg]=ode23(@hopf,[0 tfinal],x0g);
[t xh]=ode23(@hopf,[0 tfinal],x0h);
[t xia]=ode23(@hopf,[0 tfinal],x0ia);
[t xib]=ode23(@hopf,[0 tfinal],x0ib);
[t xic]=ode23(@hopf,[0 tfinal],x0ic);
[t xid]=ode23(@hopf,[0 tfinal],x0id);
[t xie]=ode23(@hopf,[0 tfinal],x0ie);
[t xif]=ode23(@hopf,[0 tfinal],x0if);
[t xig]=ode23(@hopf,[0 tfinal],x0ig);
[t xih]=ode23(@hopf,[0 tfinal],x0ih);

pause(0.001)

plot(xa(:,1),xa(:,2),'b',xb(:,1),xb(:,2),'b',xc(:,1),xc(:,2),'b',xd(:,1),xd(:,2),'b',...
    xe(:,1),xe(:,2),'b',xf(:,1),xf(:,2),'b',xg(:,1),xg(:,2),'b',xh(:,1),xh(:,2),'b',...
    xia(:,1),xia(:,2),'r',xib(:,1),xib(:,2),'r',xic(:,1),xic(:,2),'r',xid(:,1),xid(:,2),'r',...
    xie(:,1),xie(:,2),'r',xif(:,1),xif(:,2),'r',xig(:,1),xig(:,2),'r',xih(:,1),xih(:,2),'r');
title(['Hopf bifurcation with -1 \leq \mu \leq 2; \mu: ', num2str(mu)]);
grid;
xlabel('X1');
ylabel('X2');

end



% "Complex and Chaotic Nonlinear Dynamics. 
% Advances in Economics and Finance, 
% Mathematics and Statistics" 
% T.Vialar, Springer 2009
% Copyright(c).

