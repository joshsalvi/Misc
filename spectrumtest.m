delta_t = 10^-3;
total_t = 1;
%Subtract delta_t from the end to get an even number of data points.
time = 0:delta_t:total_t-delta_t;
Fs2 = 1/delta_t; %Sampling frequency

%Sinusoidal signal
NP = 10;%Number of periods
T = total_t/NP; %Period
w = 2*pi/T;
%X = cos(w.*time);%Without noise
X = 1*cos(w.*time) + 0.5.*randn(length(time),1)';%With noise

NFFT = 2^3*2^nextpow2(length(stim_center));

f = (Fs/2)*linspace(0,1,NFFT/2+1);
effective_delta_f = Fs/NFFT;
actual_delta_f = Fs/numel(time);

h_stim=psd(spectrum.periodogram('Rectangular'),X,'Fs',Fs2,'NFFT',NFFT);

sqrt(h_stim.Data(findnearest(h_stim.Frequencies,NP))*2)

L = length(X);
hfft  = fft(X,NFFT)/L;


