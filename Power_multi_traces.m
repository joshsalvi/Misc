function c = Power_multi_traces(trace_a,num_win,dt,paddingtimes)

%   This procedure calculates the averaged single-sided power spectral density with proper normalization, so
%   that if the signal is given in volts, the value is truely the power (with unit resistance). The output unit is (input trace unit)^2/Hz.
%   It uses multiple tapers to window the data, and thereafter performs
%   a padding. The procedure is inspired by the procedure from the chronux.org software 
%	package.
% 
%   The bandwidth of the spectrum 2W is determined by the relation:
%   2W =(2K-1)/T where N is the number of spectral estimates (num_win) and T 
%   is total length of the trace in time. So for instance, 5 windows and 2 sec trace gives
%   a bandwidth of 4.5 Hz.
%
%   The averaging is taken over both the trials (2nd dimension of trace_a) and
%   over multiple tapers (estimates of spectrum for a single trial).
%   These different estimates of the spectrum is also use to
%   calculate standard error of the estimates by boot strap technique.
% 
%	INPUT:
%  	"trace_a" is the array or matrix containing the trace that is the basis for
%	spectral estimation. The each row is a trial and each column is a data point sample
% 	at a period of "dt".  trace_a can be one-dimensional in the special case of a single
%	recording trace. Multiple traces give a better estimate but stationarity is required.
%	"num_win" is the number of tapers - 5-10 is usually good. 
%	"paddingtimes" gives you an artificial increase in spectral resolution- 5-10 times 
%	is usually good. This value slows down the estimation. 
%  
%   OUTPUT:
%   c(:,1) = frequency base
%   c(:,2) = averaged power spectral density (single sided) Unit is [(input trace unit)^2/Hz].
%   c(:,3) = standard error as a function of frequency. 2*SE provides a 95%
%   confidence limit.
%



if size(trace_a,2) > size(trace_a,1)  %Assuming more data points than traces
    trace_a=trace_a';
end

npp=size(trace_a,1);
np_traces=size(trace_a,2);
windows=dpss(npp,num_win/2);  % Defining window functions, using discrete prolate spherical sequences (Percival and Walden 1998).
fs=1./(dt);fmax=1./(2*dt)   ;% Maximum frequency with unit is Hz

for j=1:np_traces
    trace_b(:,j)=trace_a(:,j)-mean(trace_a(:,j)); % Removing off-set
end

for j=1:np_traces
    % Loop for spectral averaging over windows
    for i=1:num_win
        aa=fft(padding(trace_b(:,j).*windows(:,i),paddingtimes));
        bb=conj(aa);                            %Complex conjugation
        u=sum(windows(:,i).*windows(:,i)) ;     %/npp;  % the mean square of the window function (power from window).
    %    norm_a=fs*u*npp*paddingtimes/2  ;      %Normalization.  fs = sampling freq and u is the window,
        norm_a=u*npp*paddingtimes/2  ;          %Normalization.  fs = sampling freq and u is the window,                                                      %and the factor of 2 is for single sided power.
        pow_a(:,i+(j-1)*num_win)= aa.*bb/norm_a;              %np is length, which is increased when padding,  
    end
end

newnpp=size(pow_a,1);               % Length of padded power spectrum 
df=fmax./(round((newnpp/2)-1));     % Frequency resolution.
n_estimates=np_traces*num_win;

% Total mean power spectral density (power per frequency).
px=sum(pow_a,2)/(size(pow_a,2)*df);

% Delete one array spectrum:

baseset=1:n_estimates;

for k=1:n_estimates    
    setN=find(k~=baseset);
    pow_deleteone(:,k)=sum(pow_a(:,setN),2)/(n_estimates-1);
end

%  Calulating the standard error of the power spectral estimate (boot strap):

for ii=1:newnpp
    sigma_px(ii)=sqrt(((n_estimates-1)/n_estimates)*sum( (pow_deleteone(ii,:) - px(ii)).^2,2));
end

pxx=px(1:round(newnpp/2))   ;    % Disregarding the negative frequencies.
freqa=0:round((newnpp/2)-1);     % The frequency vector
sigma_px2=sigma_px(1:round(newnpp/2));

freq=freqa*df;  % Defining frequency trace with proper resolution, df.
c(:,1)=freq;
c(:,2)=pxx;
c(:,3)=sigma_px2';


function out=padding(xx,ntimes)
% this function add zeros around the trace xx, so that is is ntimes the
% original length. The mean of the trace is also removed.
% Rune W. Berg 2008

xx=xx-mean(xx); np=length(xx); xxp(1:np)=xx;
if ntimes > 0
    xxp(np+1:ntimes*np)=0;
end
out=xxp;
