%% Generate two data sets: one with a single exponential and one with two
t = (1:0.01:100)';
N_o = 100;
Tau1 = 20;
Tau2 = 5;
A1 = 50;    % Magnitude of each component
A2 = 50;    % Magnitude of each component
y = N_o - (A1+A2)*exp(-t./Tau1) + randn(size(t));   % Single exponential
yy = N_o ...                                        % Double exponential
    - A1*exp(-t./Tau1) ...
    - A2*exp(-t./Tau2) ...
    + randn(size(t));
h1 = plot(t,y,'color',[0 0.2 0.8]);
hold on
h2 = plot(t,yy,'k');
axis([-10 100 0 110])
legend([h1 h2],'Single Exp.','Double Exp.','location','southeast')
xlabel('Time (ms)')
ylabel('Dependent Var.')
title('Model Data')

%% fit for single exponential
clc
[Taus,CI,R2vals,fitA,fitB,gofA,gofB,AIC,dAIC,BIC,dBIC,iterations,residuals,tel] = taucalc(t, y, 5);
fits = {'Single Exponential','Double Exponential'};
disp(['INPUT: Single exponential'])
disp(['Number of iterations: ' num2str(iterations)])
disp(['R^2 values: ' num2str(R2vals)]);disp(' ');
disp(['*** AIC ***']);disp(['***********']);
disp(['Single exp: ' num2str(AIC(1))]);
disp(['Double exp: ' num2str(AIC(2))]);
disp(' ');disp(['*** BIC ***']);disp(['***********']);
disp(['Single exp: ' num2str(BIC(1))]);
disp(['Double exp: ' num2str(BIC(2))]);
disp(' ');
disp(['Best model from AIC: ' fits{find(dAIC == 0)}])
disp(['Best model from BIC: ' fits{find(dBIC == 0)}])
%% fit for double exponential
[Taus,CI,R2vals,fitA,fitB,gofA,gofB,AIC,dAIC,BIC,dBIC,iterations,residuals,tel] = taucalc(t, yy, 5);
clc
fits = {'Single Exponential','Double Exponential'};
disp(['INPUT: Double exponential'])
disp(['Number of iterations: ' num2str(iterations)])
disp(['R^2 values: ' num2str(R2vals)]);disp(' ');
disp(['*** AIC ***']);disp(['***********']);
disp(['Single exp: ' num2str(AIC(1))]);
disp(['Double exp: ' num2str(AIC(2))]);
disp(' ');disp(['*** BIC ***']);disp(['***********']);
disp(['Single exp: ' num2str(BIC(1))]);
disp(['Double exp: ' num2str(BIC(2))]);
disp(' ');
disp(['Best model from AIC: ' fits{find(dAIC == 0)}])
disp(['Best model from BIC: ' fits{find(dBIC == 0)}])
