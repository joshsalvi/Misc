% INPUT PARAMETERS
%S = Xd;                 % Raw time series (S(:,1) is time in ms, S(:,2) is in nm) 
Fs = 10000;             % Sampling frequency
wsd = 501;              % Window size to remove drift
wss = 151;              % Window size for smoothing
N = 20;                 % Into how many bins would you like to divide total time trace for smoothing?
out = 0.2;              % Criterion for finding outliers, defined by mean(S)+out*std(S)

% Drift removal
ws2 = floor(wsd/2);
Savg = smooth(S(:,2),wsd,'sgolay',1);      % 1 = Moving Average, Anything greater is a higher degree spline
Sc = S(:,2) - Savg;

% Calculate the PSD, using a taper of choice
w = length(S(:,2))/N;   % Window size for smoothing (see above)
overlap = 0;            % Overlap of windows
NFFT = 2*Fs - 1;

choice = menu('Taper Choice','None (Boxcar)','Hann','Hamming','Blackman');

switch choice
    case 1
        h = spectrum.welch('Rectangular',w,overlap);
        Hpsd=psd(h,S(:,2),'NFFT',NFFT,'Fs',Fs,'ConfLevel',0.95);
        figure(1); loglog(Hpsd.Frequencies,Hpsd.Data);title('No Window');
        
    case 2
        h = spectrum.welch('Hann',w,overlap);
        Hpsd=psd(h,S(:,2),'NFFT',NFFT,'Fs',Fs,'ConfLevel',0.95);
        figure(2); loglog(Hpsd.Frequencies,Hpsd.Data);title('Hann Window');
        
    case 3
        h = spectrum.welch('Hamming',w,overlap);
        Hpsd=psd(h,S(:,2),'NFFT',NFFT,'Fs',Fs,'ConfLevel',0.95);
        figure(3); loglog(Hpsd.Frequencies,Hpsd.Data);title('Hamming Window');
        
    case 4
        h = spectrum.welch('Blackman',w,overlap);
        Hpsd=psd(h,S(:,2),'NFFT',NFFT,'Fs',Fs,'ConfLevel',0.95);
        figure(4); loglog(Hpsd.Frequencies,Hpsd.Data);title('Blackman Window');
        
end


% Fit the PSD to a Lorentzian

Pxx(:,2) = Hpsd.Data;
Pxx(:,1) = Hpsd.Frequencies;
wss2 = floor(wss/2);
PxxAvg = smooth(Pxx(:,2),'sgolay',1);
Pstd = zeros(1,length(Pxx));

for i = 1:wss2
    Pstd(i) = std(Pxx(1:wss2+1,2));
end
for i = wss2+1:length(Pxx)-(wss2+1)
    Pstd(i) = std(Pxx(i-wss2:i+wss2+1,2));
end
for i = length(Pxx)-(wss2+1):length(Pxx)
    Pstd(i) = std(Pxx(length(Pxx)-(wss2+1):length(Pxx),2));
end

outliers = abs(Pxx(:,2) - PxxAvg) > out*PxxAvg;

PxxOut(:,1) = Hpsd.Frequencies;
PxxOut(:,2) = Hpsd.Data;
PxxOut(outliers,:) = [];                % Removes outliers

subi = find(1<PxxOut(:,1) & PxxOut(:,1)<2*10^3);
x = PxxOut(subi,1);
y = PxxOut(subi,2);
figure(5);
loglog(x,y,'b');
hold on

figure(6);
plot(S(:,1),S(:,2),'k');

s = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[1 1]);
f = fittype('a/(f0sq+x^2)','options',s);
[c,gof,output] = fit(x,y,f)
plot(c,'g')
kB = 1.3806503*10^-23;
T = 295;
PSDunits = 10^-18; %m^2/Hz
gamma = kB*T/(pi^2*c.a*PSDunits)
k = sqrt(c.f0sq)*2*pi*gamma

