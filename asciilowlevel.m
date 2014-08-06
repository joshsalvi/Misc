%This m-file demonstrates how to read and write a simple ascii file 
%in MATLAB using the low level file I\O routines

%write the data file
fid = fopen('square_mat.txt', 'wt');  %open the file
fprintf(fid, '%s\n', 'This is a square matrix');   %write the header
fprintf(fid, '%i\t%i\t%i\n', [1 2 3; 4 5 6; 7 8 9]');   %write data to file
fclose(fid);   %close file
disp('Write of file is successful!... press spacebar to continue')

pause

%read the file back in
fid=fopen('square_mat.txt', 'rt');   %open the file
title = fgetl(fid);  %read in the header
[data,count]=fscanf(fid, '%i',[3,3]);   %read in data
fclose(fid);   %close the file
