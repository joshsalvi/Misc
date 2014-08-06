% Takes a pixelated square and rotates it randomly many times. Demonstrates
% that this results in something like a tent filter.

im=zeros(101);
im(27:75,27:75)=1;
im2=im;
for k=1:1000         
    im2=imrotate(im2,rand(1)*90,'bilinear','crop');
    imagesc(im2);
    drawnow;
end