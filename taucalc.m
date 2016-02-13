function [Taus,CI,R2vals,fitA,fitB,gofA,gofB,AIC,dAIC,BIC,dBIC,iterations,residuals,tel] = taucalc(t, y, maxiter)
%
% Calculates exponential time constants obtained from single and double 
% exponential fits. Provides AIC value for each fit for model selection.
%
% INPUTS:
%   t:          Independent variable vector (ie. time)
%   y:          Dependent variable vector (ie. heat, position)
%   maxiter:    Maximum number of iterations to perform
%
% OUTPUTS:
%   Taus:       Exponential time constants (structure)
%   CI:         Confidence intervals (structure)
%   fitA:       Fitted model for fitA
%   fitB:       Fitted model for fitB
%   gofA:       Goodness of fit for fitA
%   gofB:       Goodness of fit for fitB
%   AIC:        Akaike's Information Criterion
%   residuals:  Residuals from fits (structure)
%   tel:        Time elapsed for fitting algorithm
%
% --------------------------
% Julien B. Azimzadeh
% jazimzadeh@rockefeller.edu
% --------------------------
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Initialize          %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ftA = fittype( 'a - b*exp(c*x)','independent','x','dependent','y' );                % First-order exponential
ftB = fittype( 'a - b*exp(c*x) - d*exp(e*x)','independent','x','dependent','y' );   % Second-order exponential
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       1st Order Exp           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rt1=0; rt2=0; rt3=0; clear outliers; warning off;
opts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares',...
        'Upper', [Inf Inf 0], 'Lower',[0 -Inf -Inf],'TolX',1e-10,'TolFun',1e-10,...
        'MaxFunEvals',4000, 'MaxIter',4000,'Robust','Off');
[fitA, gofA] = fit( t, y, ftA, opts);
Aconfint = confint(fitA);                           % Store confidence intervals
rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));    % do CI95 cross zero?
rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));    % do CI95 cross zero?
rt3 = sign(Aconfint(1,3)) + sign(Aconfint(2,3));    % do CI95 cross zero?
m = 0;

fdata = feval(fitA,t);                              % Evaluate fit
resdA = fdata - y;                                  % Calculate residuals

while rt1 == 0 || rt2 == 0 || rt3 == 0 || std(resdA) >= 4 || isnan(rt1) ||isnan(rt2) || isnan(rt3) 
    if m >= maxiter
        break;
    end
    
    m = m+1;
    fdata = feval(fitA,t);                          % Evaluate fit
    resdA = fdata - y;                              % Calculate residuals
    I = abs(resdA) > 3*std(resdA);                  % I = 1 for any residuals greater than 3 x SD
    outliers = excludedata(t,y,'indices',I);        % Find outliers
    opts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares',...
        'Upper', [Inf Inf 0], 'Lower',[0 -Inf -Inf],'TolX',1e-10,'TolFun',1e-10,...
        'MaxFunEvals',4000, 'MaxIter',4000,'Robust','Off','Exclude',outliers);
    [fitA, gofA] = fit( t, y, ftA, opts);           % Fit w/outliers removed
    Aconfint = confint(fitA);                       % Check 95% CI
    rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));% do CI95 cross zero?
    rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));
    rt3 = sign(Aconfint(1,3)) + sign(Aconfint(2,3));
end
iterA = m;                                          % Save number of iterations
TauA = -1/fitA.c;                                 % Time constant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       2nd Order Exp           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rt1=0; rt2=0; rt3=0; rt4=0; rt5=0; clear outliers; warning off; 
opts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares',...
        'Upper', [Inf Inf 0 Inf 0], 'Lower',[0 -Inf -Inf -Inf -Inf],'TolX',1e-10,'TolFun',1e-10,...
        'MaxFunEvals',4000, 'MaxIter',4000,'Robust','Off');
[fitB, gofB] = fit( t, y, ftB, opts);
Bconfint = confint(fitB);                           % Store confidence intervals
rt1 = sign(Bconfint(1,1)) + sign(Bconfint(2,1));    % do CI95 cross zero?
rt2 = sign(Bconfint(1,2)) + sign(Bconfint(2,2));    % do CI95 cross zero?
rt3 = sign(Bconfint(1,3)) + sign(Bconfint(2,3));    % do CI95 cross zero?
rt4 = sign(Bconfint(1,4)) + sign(Bconfint(2,4));    % do CI95 cross zero?
rt5 = sign(Bconfint(1,5)) + sign(Bconfint(2,5));    % do CI95 cross zero?
m = 0;

fdata = feval(fitB,t);                              % Evaluate fit
resdB = fdata - y;                                  % Calculate residuals

while rt1 == 0 || rt2 == 0 || rt3 == 0 || rt4 == 0 || rt5 == 0 || std(resdB) >= 4 || isnan(rt1) ||isnan(rt2) || isnan(rt3) || isnan(rt4) || isnan(rt5)
    if m >= maxiter
        break;
    end
    m = m+1;
    fdata = feval(fitB,t);                          % Evaluate fit
    resdB = fdata - y;                              % Calculate residuals
    I = abs(resdB) > 3*std(resdB);                  % I = 1 for any residuals greater than 3 x SD
    outliers = excludedata(t,y,'indices',I);        % Find outliers
    opts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares',...
        'Upper', [Inf Inf 0 Inf 0], 'Lower',[0 -Inf -Inf -Inf -Inf],'TolX',1e-10,'TolFun',1e-10,...
        'MaxFunEvals',4000, 'MaxIter',4000,'Robust','Off','Exclude',outliers);
    [fitB, gofB] = fit( t, y, ftB, opts);           % Fit w/outliers removed
    Bconfint = confint(fitB);                       % Check 95% CI
    rt1 = sign(Bconfint(1,1)) + sign(Bconfint(2,1));% do CI95 cross zero?
    rt2 = sign(Bconfint(1,2)) + sign(Bconfint(2,2));
    rt3 = sign(Bconfint(1,3)) + sign(Bconfint(2,3));
    rt4 = sign(Bconfint(1,4)) + sign(Bconfint(2,4));
    rt5 = sign(Bconfint(1,5)) + sign(Bconfint(2,5));
end
iterB = m;                                          % Save number of iterations
TauB = [-1/fitB.c -1/fitB.e];                       % Time constants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Outputs             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CI.A = confint(fitA);                               % 95% confidence intervals
CI.B = confint(fitB);                               % 95% confidence intervals
R2vals = [ gofA.rsquare gofB.rsquare]; % R-squared
iterations = [iterA iterB];                   % Iterations
tel = toc;                                          % Time elapsed
residuals.A = resdA;                                % Residuals from fitA
residuals.B = resdB;                                % Residuals from fitB
Taus.A = TauA;                                      % Time constants from fits
Taus.B = TauB;                                      % Time constants from fits
n = length(t);

% R2 = SSR/SST = SSR/(SSR+SSE)
% so, SSR = SSE * R2/(1-R2)
ssrA = gofA.sse * gofA.rsquare / (1 - gofA.rsquare);    % calculate residual sum of squares
ssrB = gofB.sse * gofB.rsquare / (1 - gofB.rsquare);

AIC(1,1) = log(1./n * (residuals.A'*residuals.A)) + (2.*(3+1))./n;       % AIC for Fit A (single exponential)
AIC(1,2) = log(1./n * (residuals.B'*residuals.B)) + (2.*(5+1))./n;       % AIC for Fit B (double exponential)
dAIC = AIC - min(AIC);  % change in AIC; dAIC = 0 is the best model
BIC(1,1) = log(1./n * (residuals.A'*residuals.A)) + (2*(3+1))./n * log(n);         % AIC for Fit A (single exponential)
BIC(1,2) = log(1./n * (residuals.B'*residuals.B)) + (2*(5+1))./n * log(n);         % AIC for Fit B (double exponential)
dBIC = BIC - min(BIC);


