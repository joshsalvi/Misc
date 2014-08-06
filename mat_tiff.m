MRI_transverse=uint16(MRI_transverse);
outputFileName = '/Users/jonathanfisher/Dropbox/Kickstarter/Imaging Data/MRI Images/Jons MRI/MRI_transverse.tif'
for K=1:length(MRI_transverse(1, 1, :))
imwrite(MRI_transverse(:, :, K), outputFileName, 'WriteMode',
'append','compression','none');
end