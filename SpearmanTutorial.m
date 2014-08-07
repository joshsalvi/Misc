% Calculate Spearman's correlation - A Tutorial
% This script calculate's Spearman's correlation using MATLAB's built-in
% methods and also does it manually. Notes are provided for each step.

clear all; close all;

% INPUT a vector of stiffnesses
kvec(:,1) = [100 100 100 200 200 200 300 300 300];   % µN/m
% INPUT a vector of amplitudes
amplvec(:,1) = [24 26 20 19 22 10 9 8 10];           % nm
% The above vectors must have the same length!

% Create a scatter plot of the stiffnesses and vectors
figure(1);
subplot(2,1,1);
scatter(kvec,amplvec,'k.'); xlabel('Stiffness (µN/m)');ylabel('Amplitude (nm)');
xmin = min(kvec);xmax=max(kvec);ymin=min(amplvec);ymax=max(amplvec);
axis([xmin-50 xmax+50 ymin-5 ymax+5])


% Use MATLAB's built-in function to calculate Spearman's rho and the
% p-value associated with it. 
[rho,prho] = corr(kvec,amplvec,'type','Spearman');

% We now rank these parameters by hand
% First, create the ranked variables. This is done using a sorting method 
% in this script.
% To do so, sort the variables in order and assign each a ranking. For
% those that are duplicates in each vector, assign the mean of their
% possible ranks.


% Sort both the stiffnesses and the forces
kvecu = unique(kvec);   % Find the unique variables
[xs,z1] = sort(kvec);   % Sort the stiffness vector
[z1,z2] = sort(z1);     % Sort the indixes
r = (1:length(kvec))';
r=r(z2);

for i=1:length(kvecu)
    s=find(kvecu(i)==kvec);    
    r(s,1) = mean(r(s));    
end
kveci = r;

clear xs z1 z2 r s
amplvecu = unique(amplvec);
[xs,z1] = sort(amplvec);   % Sort the stiffness vector
[z1,z2] = sort(z1);     % Sort the indixes
r = (1:length(kvec))';
r=r(z2);

for i=1:length(amplvecu)
    s=find(amplvecu(i)==amplvec);    
    r(s,1) = mean(r(s));    
end
amplveci = r;


% Create a second plot of the ranked variables.
figure(1);
subplot(2,1,2);
plot(kveci,amplveci,'k.'); xlabel('Ranked Stiffness'); ylabel('Ranked Amplitude');
xmin = min(kveci);xmax=max(kveci);ymin=min(amplveci);ymax=max(amplveci);
axis([xmin-1 xmax+1 ymin-1 ymax+1])

% Calculate Spearman's rho

% In the presence of tied ranks simply find the Pearson's correlation
% coefficient of the ranked variables.
[rho2] = corr(kveci,amplveci,'type','Pearson');

% Calculate Pearson’s correlation coefficient
% for comparison
[rpears,ppears] = corr(kvec,amplvec,'type’,’Pearson’);

display(sprintf('%s %s’,’Pearsons r from corr() function: ',num2str(rpears)));
display(sprintf('%s %s','Spearmans rho from corr() function: ',num2str(rho)));
display(sprintf('%s %s','Spearmans rho from hand calculation: ',num2str(rho2)));
display(sprintf('%s %s','p-value from Pearsons correlation: ',num2str(rpears)));
%display(sprintf('%s %s','p-value from Spearmans correlation: ',num2str(prho)));
