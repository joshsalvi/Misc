% bandpass filter

d = fdesign.bandpass('N,Fp1,Fp2,Ap', 30, (3090-20), (3090+20), .5, Fs);
Hd = design(d, 'cheby1');
ylow=filter(Hd,VarName1);

