%CT_saggital=CT_saggital(:,:,2:131);

D = smooth3(MRI_saggital(:,:,:));

X_range = 0:10:512;
Y_range = 0:10:512;
Z_range = 0:0.1:23;
[Xq,Yq,Zq]=meshgrid(X_range,Y_range,Z_range);
N= 2;

Dq = interp3(D,Xq,Yq,Zq);

  figure();
hq = vol3d('cdata',D,'texture','3D');
view(3);  
  axis tight;  daspect([1 1 .4])
  alphamap('rampup');
  alphamap(.1 .* alphamap);
  
  [f,v]=isosurface(Dq,400);
  vertface2obj(v,f,'jon-mri4.obj');
