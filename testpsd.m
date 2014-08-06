
Fs=10000;
t=linspace(0,1,Fs);
A=5;
f=30;
sigma=30;

y=A*sin(2*pi*f*t)+sigma*rand(1,length(t))+2*A*sin(2*pi*0.5*f*t)+0.5*A*sin(2*pi*3*f*t);

[q1(:,2),q1(:,1)]=pwelch(y,[],[],[],Fs);
q2=Power_multi_traces(y,1,1/Fs,1);
q3=Power_multi_traces(y,5,1/Fs,1);
q4=Power_multi_traces(y,10,1/Fs,1);
[q5(:,2),q5(:,1)] = spect(y,1,[],'as',Fs);

h=psd(spectrum.periodogram('Flat Top'),y,'Fs',Fs);
q6(:,1)=h.Frequencies;
q6(:,2)=h.Data;


figure(1);
plot(t,y);

figure(2);
loglog(q1(:,1),q1(:,2)); title('pwelch');

figure(3);
loglog(q2(:,1),q2(:,2)); title('MTM, T=1');

figure(4);
loglog(q3(:,1),q3(:,2)); title('MTM, T=5');

figure(5);
loglog(q4(:,1),q4(:,2)); title('MTM, T=10');

figure(6);
loglog(q5(:,1),q5(:,2)); title('asd');

figure(7);
loglog(q6(:,1),q6(:,2)); title('flat top');


p=find(q1(:,2)==max(q1(:,2)));
Xd_freq_pwelch=q1(p,1)
Xd_amp_pwelch=sqrt(q1(p,2)*(q1(p+1,1)-q1(p,1)))

p=find(q2(:,2)==max(q2(:,2)));
Xd_freq_MTM0=q2(p,1)
Xd_amp_MTM0=sqrt(q2(p,2)*(q1(p+1,1)-q1(p,1)))

p=find(q3(:,2)==max(q3(:,2)));
Xd_freq_MTM5=q3(p,1)
Xd_amp_MTM5=sqrt(q3(p,2)*(q3(p+1,1)-q3(p,1)))

p=find(q4(:,2)==max(q4(:,2)));
Xd_freq_MTM10=q4(p,1)
Xd_amp_MTM10=sqrt(q4(p,2)*(q4(p+1,1)-(q4(p+1,1,1)-q4(p,1))))

p=find(q5(10:length(q5(:,2)),2)==max(q5(10:length(q5(:,2)),2)));
Xd_freq_asd=q5(p+9,1)
Xd_amp_asd=q5(p+9,2)

p=find(q6(10:length(q6(:,2)),2)==max(q6(10:length(q6(:,2)),2)));
Xd_freq_FT=q6(p+9,1)
Xd_amp_FT=sqrt(q6(p+9,2).*(q6(p+9,1)))


%% Alpha correlation

phaselag = 27;

Fs=10000;
t=linspace(0,1,Fs);
A=5;
f=30;
sigma=30;

y = A*sin(2*pi*f*t) + 2*A*sin(2*pi*0.5*f*t) + 1.7*A*sin(2*pi*2*f*t) + sigma*rand(1,length(t));
y2 = A*sin(2*pi*f*t + phaselag*rand(1,length(t))) + 2*A*sin(2*pi*0.5*f*t + phaselag*rand(1,length(t))) + 1.7*A*sin(2*pi*2*f*t) + sigma*rand(1,length(t));

acalcorr = xcorr(y,y2,'unbiased');
z = smooth(acalcorr,length(acalcorr)/5);
acalcorr = acalcorr - z';

h=psd(spectrum.periodogram('Flat Top'),acalcorr,'Fs',Fs);
q(:,1)=h.Frequencies;
q(:,2)=h.Data;

a_corr = xcorr(y,'unbiased');
z = smooth(a_corr,length(a_corr)/5);
a_corr = a_corr - z';

h2=psd(spectrum.periodogram('Flat Top'),a_corr,'Fs',Fs);
q2(:,1)=h2.Frequencies;
q2(:,2)=h2.Data;

%p = find(q(10:length(q(:,2)),2)==max(q(10:length(q(:,2)),2)));
p = find(round(q(:,1))==2*f);
p2 = find(round(q2(:,1))==2*f);

Xd_freq_corr=q(p,1)
Xd_amp_corr=sqrt(q(p,2).*(q(p,1)))

Xd_freq_auto=q2(p,1)
Xd_amp_auto=sqrt(q2(p,2).*(q2(p,1)))

Xd_amp_auto_mean = mean(Xd_amp_auto);
Xd_amp_corr_mean = mean(Xd_amp_corr);

amp_ratio = Xd_amp_corr_mean/Xd_amp_auto_mean
