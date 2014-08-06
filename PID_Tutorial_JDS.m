% What follows is a tutorial on PID control using MATLAB as a guide. Each
% section and figure includes comments to help the reader. Many of the
% functions require the control system toolbox.
% This tutorial was written by J.D. Salvi (2014), josh.salvi@gmail.com and
% adapted from various web sources.

%% PID in MATLAB
% One can first define a PID controller in MATLAB in two different ways.
% Method (1) defines the PID controller from its transfer function. Method
% (2) uses the pid() object to define a PID controller. If you do not have
% pid() installed, use the first method to accomplish the same task

% Method (1): PID controller from transfer function. Select the values for
% your PID controller below.
Kp = 0.5;   % proportional gain
Ki = 0.3;   % integral gain
Kd = 0.3;   % derivative gain

s = tf('s');    % define a continuous time transfer function
C = Kp + Ki/s + Kd*s   % transfer function of PID controller

% Method (2): PID controller with MATLAB's pid() function. Select the values 
% for your PID controller below.
Kp = 0.5;   % proportional gain
Ki = 0.3;   % integral gain
Kd = 0.3;   % derivative gain
C = pid(Kp,Ki,Kd)   % PID controller with PID function
tf(C)               % convert pid object to transfer function

% This section will output the continuous time transfer functions and
% parameters of the PID controller for each of the two methods above. The
% transfer function from Method (1) should yield the same answer as that
% calculated in Method (2).

%% Mass-Spring Damper Problem
% Here we wish to model the transfer function of a mass with added drag and
% an elastic element. The equation of motion of the system is: M*x_ddot +
% G*x_dot + K*x = F. Here we can model the transfer function and the
% subsequent step response for this system. We will then use the same
% system and add a PID controller to it.

% Below are the parameters of the system. Select your parameters
% accordingly.
M = 1.1;    % kg
G = 1;      % N·s/m
K = 0.1;    % N/m
F = 0.2;    % N

% We may then calculate the transfer function of this system from the
% Laplace transform of the system's equation of motion. 
s = tf('s');
Po = 1/(M*s^2 + G*s + K);   % mass-spring-damper transfer function

figure(1);
step(Po)        % plot the step response of the system

% You should now see a plot with the open loop (PID off) step response of the
% mass-spring damper system. Use this image from figure 1 for comparison
% with those in the next steps to see how adding each of the PID controls
% affects the step response.

%% Add Proportional Control
% We now take the same mass-spring-damper and add proportional control to
% the system. Run the previous session so that you use the same values for
% the mass, damping, stiffness, and force for consistency. Play around with
% different values of proportional gain to see what the response looks
% like, comparing it with the open loop response and responses from other
% values of proportional gain.

% Set the value of your proportional gain
Kp = 50;    % proportional gain

Vc = pid(Kp)    % pid object with your value of proportional gain
Vo = feedback(Vc*Po,1)     % provide feedback using the pid controller

t = 0:0.01:10;      % define a time vector
figure(2);
step(Vo,t);         % step response

% What you should see is that by increasing the proportional gain, Kp, the
% rise time of the step response decreases, the overshoot increases, and
% the steady-state error decreases. 

%% Add Proportional-Derivative Control
% We now add bothe proportional and derivative gains, as above. This
% analysis plots the step responses as you adjust two different values, Kp
% and Kd. Be sure to run the previous instances first for consistency. Try
% holding Kp constant and see what happens when you adjust Kd.

% Set the gain values;
Kp = 50;    % proportional gain
Kd = 5;     % derivative gain

Vc = pid(Kp,0,Kd)   % pid object with Kp and Kd
Vo = feedback(Vc*Po,1)   % provide feedback using the pid controller

t = 0:0.01:10;      % define a time vector
figure(3);
step(Vo,t);         % step response

% In this case, you should see that an increase in the derivative gain
% reduces the overshoot and the settling time for the step response.

%% Add Proportional-Integral Control
% Now we turn the derivative controll off and adjust the integral control.
% The same protocol from before follows. Again, try holding the
% proportional gain constant and adjust the integral gain to see what
% happens.

% Set the gain values;
Kp = 20;    % proportional gain
Ki = 60;    % integral gain

Vc = pid(Kp,Ki)   % pid object with Kp and Kd
Vo = feedback(Vc*Po,1)   % provide feedback using the pid controller

t = 0:0.01:10;      % define a time vector
figure(4);
step(Vo,t);         % step response

% As the integral gain increases, you should see that the integral gain
% reduces the rise time and increases the overshoot of the system. It may
% be difficult to find a good answer with a large value of proportional
% gain, as a result. Try making Kp small and then adjust Ki. A high
% integral gain should also greatly reduce the steady-state error.

%% Use a complete PID Controller
% Finally, we use all three gains and set up a complete PID controller. In
% this case, try to find parameters that eliminates overshoot and
% steady-state error with the fastest rise time possible.

% Set the gain values;
Kp = 20;    % proportional gain
Ki = 70;    % integral gain
Kd = 60;    % derivative gain

Vc = pid(Kp,Ki,Kd)   % pid object with Kp and Kd
Vo = feedback(Vc*Po,1)   % provide feedback using the pid controller

t = 0:0.01:10;      % define a time vector
figure(5);
step(Vo,t);         % step response

% Did you see any strange behavior for certain values? How can these be
% troubleshooted?

%% How to Design a PID Controller

% (1) Start with an open loop response to see what you need to adjust.
% (2) Add Kp to reduce rise time.
% (3) Add Kd to reduce the overshoot.
% (4) Add Ki to reduce steady-state error.
% (5) Adjust Kp, Ki, and Kd to achieve an optimal response that does not
% self-oscillate.


%% Automated PID Tuning
% MATLAB has a built-in function that allows you to adjust the parameters
% in a PID controller to achieve optimal parameters in a GUI interface.
% Open it in this dialog, click on "Show Parameters" to see what you are
% adjusting, and play with different parameters. This GUI is set up to take
% the mass-spring-damper system and optimize its step response.

close all;

pidtool(Po,Vc);         % CHOICE (1) - pick one

% pidtool() takes as its first argument the system you are controlling and
% in its second parameter the PID controller. You can also use ONLY a
% proportional controller, for example, using the following command:

pidtool(Po,'p');        % CHOICE (2) - pick one

%% PID Tuning Without a GUI
% You can also define a PID controller like the one above from the command
% line instead of the GUI interface above. Here, you can define a bandwith
% in rad/s and a phase margin in deg. This code will then output the PID
% controller with associated values for Kp, Ki, and Kd for these two
% parameters.

% Define your parameters
bwpid = 20;     % bandwidth, rad/s
pmpid = 90;     % phase margin, deg

% Set the options for your PID controller
opts = pidtuneOptions('CrossoverFrequency',bwpid,'PhaseMargin',pmpid);
% Find the optimal settings for your controller using pidtune()
[Vc, info] = pidtune(Po, 'pid', opts)

% This creates an output with values for Kp, Ki, and Kd. Additionally, it
% outputs your bandwidth and phase margin. Finally, info.stable will tell
% you whether or not the PID controller will be stable. Which parameters
% are optimal? Do they make sense from the work you did in the GUI? How
% about from the manually generated PID controller? Finally, which
% parameters lead to the system's loss of stability? You can use a code
% like this or any of the above to define a possible set of PID parameters
% for a given system. In the case of active hair-bundles, for example, you
% can use the mass-spring-damper system, but be aware that a nonlinearity
% is present in at least one of the terms (damping or stiffness), and there
% is an additional active force. If you are daring, try taking the Laplace
% transform of a mass-spring-damper system with a nonlinearity included,
% find the transfer function, and copy some of the code from this file.



