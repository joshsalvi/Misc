warning off

%clear all; close all;
close all;
%t = linspace(0,10,1e3);
%t = Timems;
t = linspace(0,length(CLASSIFICATION)*1e-4,length(CLASSIFICATION));

f = 10;
%y = sin(f*pi*t);
%y = XdInput1Displacementnm11;
%y = rand(1,length(t));
y = CLASSIFICATION(:,6);


figure(1);
plot(t,y);

yc = xcorr(y,'coeff');

ych = hilbert(yc);

yche = abs(ych);


index = find(yche(length(yche)/2:length(yche))<max(yche(length(yche)/2:length(yche)))/exp(1));

index1 = t(index(1));

[yf yfr] = fit(t',yche(length(yche)/2:length(yche)),'exp1');
figure(2);
text(max(t)/2,0.9*max(yche),'a = ');text(max(t)/2+0.07*max(t),0.9*max(yche),num2str(yf.a)); hold on;
text(max(t)/2,0.85*max(yche),'b = ');text(max(t)/2+0.07*max(t),0.85*max(yche),num2str(yf.b));
text(max(t)/2,0.8*max(yche),'tau = ');text(max(t)/2+0.07*max(t),0.8*max(yche),num2str(-1/yf.b));
text(max(t)/2,0.75*max(yche),'R2 = ');text(max(t)/2+0.07*max(t),0.75*max(yche),num2str(yfr.rsquare));
text(max(t)/2,0.7*max(yche),'tau2 = ');text(max(t)/2+0.07*max(t),0.7*max(yche),num2str(index1));
plot(t,yc(length(yc)/2:length(yc)),'k'); hold on;
plot(t,yche(length(yche)/2:length(yche)),'r');
plot(yf,'b');
