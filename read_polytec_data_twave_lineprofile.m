% reduce_twave_2D_to_1D.m is an mfile for reading polytec data generated
% by a 2D scanning interferometer looking at an arc segment of a turn of the
% gerbil's cochlea.  It then reduces the arc data into a 1D line. The
% purpose is to create 1D data that can be used to generate a gain function
% (c.f. C. Shera, J. Acoust. Soc. Am. Volume 122, Issue 5, pp. 2738-2758).

% NOTE: This program requires the following mfiles in the same folder as this mfile:
% GetPointData.m
% Get0dBReference.m
% GetXYZCoordinates.m
% circfit.m

%EXAMPLE
%[pointdata1,img1,I1,timeaxis1]= reduce_twave_2D_to_1D('C:\ProgramFiles\MATLAB\R2008b\3khz\3khz', 'FFT', 'Vib','Displacement','Real & Imag.', 1.55e6,1);
% INPUTS:
% basename: text string describing file path, excluding filetype suffix
% domainname: name of the domain, e.g. 'FFT' or 'Time'
% channelname: name of the channel, e.g. 'Vib' or 'Ref1' or 'Vib &
% Ref1'
% signalname: name of the signal, e.g. 'Velocity' or 'Displacement'
% displayname: name of the display, e.g. 'Real' or 'Magnitude' or
% 'Samples'
% pix_per_m = pixels per meter. Polytec scanning software will save the
% scanning point locations in units of meters based on the inputed
% magnification. Even if the actual dimensions are completely different. To
% co-locate the scan points with the pixels in the background image
% pixels per meter must be known. That is simply one of the image
% dimensions (pix) divided by what the scanning software says is the max
% dislacement in m.
% firstfile: (1 = yes, 0 = no). This determines whether the spectrum window
% is loaded from file (if NO), or if you are given the opportunity to
% select by hand (if YES).
% newsettings: if you want to analyze the current file using the same
% angular geometry (arc, delta angle definitions) as the previous
% experiment, set this to 0. Otherwise should be 1.
% medfilt2d = (1 = yes, 0 = no). Do you want to perform 2D median filtering
% smoothing = (1 = yes, 0 = no). Do you want to perform 2D spatial
% smoothing? Uses sliding window averaging within nearest neighbor regions.
% of the data?

%Sample syntax:
%>> [pointdata1,img,I1,timeaxis1] =
%reduce_twave_2D_to_1D('D:\MATLAB_PATCH_1,5khz\1,5khz0',
%'FFT','Vib','Displacement', 'Real & Imag.', 1.55e6);
%
%J.A.N.F. 2011

function [pointdata, data1D,I] = read_polytec_data_twave_lineprofile('C:\Users\Josh\Documents\Lab\Hudspeth Lab\UV Interferometry\20110729\PassiveChinchilla1.svd', 'FFT', 'Vib', 'Displacement', 'Real & Imag.', 1.55e6,1, newsettings,medfilt2d,smoothing);

%PART 1: PULL ALL THE DATA FROM THE SVD FILE

%the polytec function GetPointData2.m reads the .svd file format and
%retrieves data for all points at all frequencies
%x = list of frequencies (1-by-N array, Hz)
%y = a M-by-2N array of complex Fourier amplitudes for each frequency and
%additionally for each point (there are M scan points). There are 2N FFT
%data points per scan place because the list of complex numbers is stored
%as a list of alternating real and imaginary numbers.
[x,y,usd] = GetPointData2([basename '.svd'], domainname, channelname, signalname, displayname, 0);

numpoints = size(y,1);

%now put complex values together as complex numbers (real and imag
%components are stored in neighboring data pairs in the polytec format)
if sum(findstr(displayname,'Imag'))~=0
    for j = 1:size(x,2)
        iy(:,j) = y(:,2*j);
        ry(:,j) = y(:,2*(j-1)+1);
        cy(:,j) = complex(ry(:,j),iy(:,j));
    end
else
    cy = y;
end

%we now begin adding information to the "master cluster" called pointdata.
%This cluster will eventually consist of:
%pointdata.info             = basic acquisition settings
%pointdata.location         = XY cartesion locations of scan points
%pointdata.FFT_data         = complex FFT amplitudes for all freqs for all points
%pointdata.spectrum         = averaged FFT spectrum (over all scanpoints)
%pointdata.spectrumwindow   = range of frequencies user has selected to use
%pointdata.medfiltFFT_data  = median filtered FFT data.
%

pointdata.info = usd;

%the polytec function GetXYZCoordinates.m gets the pointdata locations:
XYZ = GetXYZCoordinates([basename '.svd'], 0);

pointdata.location(:,:) = [XYZ(:,1) XYZ(:,2)];

%We assign the FFT data to a portion of this cluster
for i = 1:size(y,1)
    %pointdata.FFT_data has dimensions of (M,N,2), where
    % M is number of scan points
    % N is number of FFT lines
    % the 3rd dimension toggles between frequency (1) and FFT coef (2)
    pointdata.FFT_data(i,:,:) = [x',cy(i,:)'];
end
%pointdata.spectrum is the average FFT spectrum, averaged over all
%scan points. This is used for selecting spectral regions of interest.
pointdata.spectrum(:,:) = sum(abs(pointdata.FFT_data),1)./numpoints;


%PART 2: LET THE USER CHOOSE WHAT FREQUENCY TO VISUALIZE

if firstfile == 1
    
    for i = 1:2
        plot(pointdata.spectrum(:,1),pointdata.spectrum(:,2),'r')
        datacursormode on
        dcm_obj = datacursormode(gcf);
        user_entry = input('press return');
        c_info = getCursorInfo(dcm_obj);
        pointdata.spectrumwindow(i) = c_info.DataIndex;
    end
    
    spectwindow = pointdata.spectrum(pointdata.spectrumwindow(1):pointdata.spectrumwindow(2),1);
    
    %define the point in the FFT plot with the maximum magnitude
    max_spect_index = find(pointdata.spectrum(:,2)==max(pointdata.spectrum(pointdata.spectrumwindow(1):pointdata.spectrumwindow(2),2)));
    
    %which corresponds to a center frequency
    center_frequency = pointdata.spectrum(max_spect_index,1)
    
    
    
    range = pointdata.spectrumwindow;
    save('spect.mat','spectwindow','range','max_spect_index','center_frequency');
    
else
    spect = load('spect.mat','spectwindow','range','max_spect_index','center_frequency');
    spectwindow = spect.spectwindow;
    max_spect_index = spect.max_spect_index;
    pointdata.spectrumwindow = spect.range;
    center_frequency = spect.center_frequency;
end


%PART 3: SPATIAL FILTERING OF 2D SCANNING DATA

% First determine average point-to-point distance

nearests = zeros(1,numpoints);
for j=1:numpoints
    
    %Here we look through the list of points to see which is the closest to
    %our current point of interest
    distance = zeros(1,numpoints-1);
    for i = 1:numpoints
        if i==j
        else
            distance(i) = sqrt((pointdata.location(i,1)-pointdata.location(j,1))^2+(pointdata.location(i,2)-pointdata.location(j,2))^2);
        end
    end
    
    distance(find(distance==0))=NaN;
    
    nearests(j) = min(distance);
   % pause
end

%We multiply the nearest neighbor distance by 1.7 rather than sqrt(2)
%because the diagonal distance to nearest neighbors is frequently skewed
%due to non-cartesian scan point geometries.
neighbor_proximity = 1.7*(sum(nearests)/numpoints);

%First perform 2D median filtering, if that is selected in the input.

if medfilt2d ==1
    
    %allocate space for median filtered FFT_data
    pointdata.medfiltFFT_data = zeros(size(pointdata.FFT_data));
    
    %set the first dimension of the data equal to the same "timeaxis"
    pointdata.medfiltFFT_data(:,:,1) = pointdata.FFT_data(:,:,1);
    
    for k = 1:size(pointdata.location,1)
        nearest_neighbor = k;
        for p = 1:size(pointdata.location,1)
            
            if p==k
                
            else
                %a criterion for whether a point is a "neighbor"
                neighbor_val = sqrt((pointdata.location(k,1)-pointdata.location(p,1))^2+(pointdata.location(k,2)-pointdata.location(p,2))^2)<=neighbor_proximity;
                
                if neighbor_val ==1
                    nearest_neighbor = cat(1,nearest_neighbor,p);
                else
                end
                
            end
            
            %size(nearest_neighbor)
        end
        
        for g = 1:size(nearest_neighbor,1)
            if g==1
                %holdingpen will be the set of pointdata.FFT_data points
                %that are "nearest neighbors" of the central one you are
                %looking at (the k'th)
                holdingpen = zeros(size(nearest_neighbor,1),size(pointdata.FFT_data,2));
                holdingpen(g,:)= pointdata.FFT_data(nearest_neighbor(g),:,2);
            else
                holdingpen(g,:)= pointdata.FFT_data(nearest_neighbor(g),:,2);
            end
        end
        
        
        %set the FFT values of the central point to the median of the
        %points within the nearest neighbor region.
        
        pointdata.medfiltFFT_data(k,:,2) = median(holdingpen,1);
%         
%                 %some programming diagnostics
%                 for q = 1:size(nearest_neighbor,1)
%                     if q==1
%                         scatter(pointdata.location(nearest_neighbor(q),1),pointdata.location(nearest_neighbor(q),2),'r');
%                         hold on
%                     else
%                     end
%                     scatter(pointdata.location(nearest_neighbor(q),1),pointdata.location(nearest_neighbor(q),2),'r');
%         
%                 end
%         
%         
%         
%                 scatter(pointdata.location(k,1),pointdata.location(k,2),'b');
%                 % pointdata.location(k,:)
%                 set(gca,'xlim',[min(pointdata.location(:,1)) max(pointdata.location(:,1))]);
%                 set(gca,'ylim',[min(pointdata.location(:,2)) max(pointdata.location(:,2))]);
%                 hold off
% %         
        %some programming diagnostics
        %             size(nearest_neighbor)
        %             size(holdingpen)
        %             size(median(holdingpen,1))
        %             size(pointdata.medfiltFFT_data)
        
        
  %             pause
        
        clear nearest_neighbor
        
    end
    pointdata.FFT_data = pointdata.medfiltFFT_data;
end

%Now perform 2D central value averaging, if that is selected in the input. This
%replaces a central value with the average of the pixels in its nearest
%neighbor region.

if smoothing ==1
    
    %allocate space for median filtered FFT_data
    pointdata.smoothFFT_data = zeros(size(pointdata.FFT_data));
    
    %set the first dimension of the data equal to the same "timeaxis"
    pointdata.smoothFFT_data(:,:,1) = pointdata.FFT_data(:,:,1);
    
    for k = 1:size(pointdata.location,1)
        nearest_neighbor = k;
        for p = 1:size(pointdata.location,1)
            
            if p==k
                
            else
                %a criterion for whether a point is a "neighbor"
                neighbor_val = sqrt((pointdata.location(k,1)-pointdata.location(p,1))^2+(pointdata.location(k,2)-pointdata.location(p,2))^2)<=neighbor_proximity;
                
                if neighbor_val ==1
                    nearest_neighbor = cat(1,nearest_neighbor,p);
                else
                end
                
            end
            
            %size(nearest_neighbor)
        end
        
        for g = 1:size(nearest_neighbor,1)
            if g==1
                %holdingpen will be the set of pointdata.FFT_data points
                %that are "nearest neighbors" of the central one you are
                %looking at (the k'th)
                holdingpen = zeros(size(nearest_neighbor,1),size(pointdata.FFT_data,2));
                holdingpen(g,:)= pointdata.FFT_data(nearest_neighbor(g),:,2);
            else
                holdingpen(g,:)= pointdata.FFT_data(nearest_neighbor(g),:,2);
            end
        end
        
        %set the FFT values of the central point to the median of the
        %points within the nearest neighbor region.
        
        pointdata.smoothFFT_data(k,:,2) = mean(holdingpen,1);
        
        clear nearest_neighbor
        
    end
    pointdata.FFT_data = pointdata.smoothFFT_data;
end


close

FFTpoints_to_use = pointdata.FFT_data;


% PART 4: GENERATE INTERPOLATED 2D DATA, DISPLAY THAT

%convert point locations from units of meters to "pixels"
pointdata.location(:,:) = pix_per_m*pointdata.location(:,:);
pointdata.location(:,:) = round(pointdata.location(:,:));

%to match up with the image coordinates (1st quadrant coords, i.e. x and y
%positive), we translate the origin from the middle of the image (polytec's
%default origin position) to the image origin by adding (xdim/2, -ydim/2)
%to the (X,Y) location pairs. Note that Matlab's ordering of image coords
%is (Y,X).

xdim = 1344;
ydim = 1039;

for i = 1:size(pointdata.location,1)
    pointdata.location(i,:) = pointdata.location(i,:)+[round(xdim/2),round(-ydim/2)];
end

pointdata.location = abs(pointdata.location);

%now we generate a meshgrid on a finer scale than the data

%define the meshgrid delta
mdelta = 15 %pixels
pix_xdim = max(pointdata.location(:,1))-min(pointdata.location(:,1));
pix_ydim = max(pointdata.location(:,2))-min(pointdata.location(:,2));


tx = min(pointdata.location(:,1)):mdelta:max(pointdata.location(:,1));
ty = min(pointdata.location(:,2)):mdelta:max(pointdata.location(:,2));
[XI,YI] = meshgrid(tx,ty);

%construct a dataset that has only the center (modulation) frequency
%information. 

timeaxis = 0:2*pi/(100*center_frequency):2*pi/center_frequency;
%final_FFTpoints_to_use = zeros(numpoints,1);
timedomain_data = zeros(numpoints,size(timeaxis,2));

for i = 1:numpoints
    final_FFTpoints_to_use(i) = FFTpoints_to_use(i,max_spect_index,2);
    for j = 1:size(timeaxis,2)
        timedomain_data(i,j) = real(final_FFTpoints_to_use(i)*exp(-complex(0,1)*center_frequency*timeaxis(j)));
    end
end

ZI = griddata(pointdata.location(:,1),pointdata.location(:,2),final_FFTpoints_to_use(:),XI,YI,'cubic');

% size(pointdata.location)
% size(final_FFTpoints_to_use)
% pause
%

%relevantdata = real(final_FFTpoints_to_use.*exp(-complex(0,1)*center_frequency*timeaxis(i)));

A(:,:) = [min(min(timedomain_data)),max(max(timedomain_data))];
figure
for i = 1:size(timeaxis,2)
    %        Tri=delaunay(pointdata.location(:,1),pointdata.location(:,2));
    %        scatter3(pointdata.location(:,1),pointdata.location(:,2),timedomain_data(:,i),'b','markerface','r');
    
    %        trisurf(Tri,pointdata.location(:,1),pointdata.location(:,2),timedomain_data(:,i));
    %        alpha(0.5)
    %ZI = griddata(pointdata.location(:,1),pointdata.location(:,2),timedomain_data(:,i),XI,YI,'cubic');
    
    %surf(XI,YI,ZI) %,'FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',ZI)
    surf(XI,YI,real(ZI.*exp(-complex(0,1)*center_frequency*timeaxis(i))))
    %     size(XI)
    %     size(YI)
    %     size(ZI)
    %     scatter3(XI(1,:),YI(1,:),ZI(1,:));\
    %     plot3(x,y,z,'o'), hold off
    view(-45,60)
    set(gca,'CLim',A);
    zlim(A);
    I(i) = getframe;
end

%in a separate window plot the absolute magnitude of the 2D wave
B(:,:) = [min(abs(FFTpoints_to_use(:,max_spect_index,2))),max(abs(FFTpoints_to_use(:,max_spect_index,2)))];
figure
Z2 = griddata(pointdata.location(:,1),pointdata.location(:,2),abs(FFTpoints_to_use(:,max_spect_index,2)),XI,YI,'cubic');
surf(XI,YI,Z2)
view(70,60)
set(gca,'CLim',B);
zlim(B);
title('Absolute Magnitude')



% PART 5: GENERATE LINE PROFILE

%plot the coordinate points
figure, scatter(pointdata.location(:,1),pointdata.location(:,2));
axis equal
hold on

%If we have chosen to, we analyse the current data using saved geometry
%settings
if newsettings==0
    geometry = load('geometry.mat','xc','yc','R','A','radius_outer','radius_inner','theta1','theta2','dtheta1','dtheta2','delta_theta','numpoints1D');
    xc = geometry.xc;
    yc = geometry.yc;
    R = geometry.R;
    A = geometry.A;
    %     radius_outer = geometry.radius_outer;
    %     radius_inner = geometry.radius_inner;
    %     theta1 = geometry.theta1;
    %     theta2 = geometry.theta2;
    %     dtheta1 = geometry.dtheta1;
    %     dtheta2 = geometry.dtheta2;
    %     delta_theta = geometry.delta_theta;
    numpoints1D = geometry.numpoints1D;
else
    %here we define the ar
    title('choose the points that correspond to an arc');
    D = ginput;
    
    %use circfit to come up with a circular fit. The function circfit.m must be
    %in the matlab search paths.
    [xc, yc, R,A] = circfit(D(:,1),D(:,2))
    
end

%
% %subtract R from all the points
% for u = 1:size(pointdata.location,1)
%     pointdata.location1(u,:) = pointdata.location(u,:)-[xc,yc];
% end

 figure
  axis equal
 hold on
%  for i = 1:(size(ZI,1)*size(ZI,2))
%      scatter(XI(i),YI(i),'m');
%  end

%pointdata.location = pointdata.location1;

%create a scale to plot the circle on (on the same axes as the
%scatterplotted pointdata)
ts = min(min(pointdata.location(:,1))):max(max(pointdata.location(:,1)));

%plot the arc
plot(ts,sqrt(R^2-(ts-xc).^2)+yc,'g');

% define the angular resolution with which you will step through the
% linedata. This will be determined by your interpolation (from meshgrid).

delta_theta = 2*asin(mdelta/(2*R))
endangle = atan2((sqrt(R^2-(min(ts-xc)).^2)),(min(ts)-xc))
startangle = atan2((sqrt(R^2-(max(ts-xc)).^2)),(max(ts)-xc))


numpoints1D = round((endangle-startangle)/delta_theta);
data1D = zeros(1,numpoints1D);
% Here, for each delta_theta corresponding to one data point on the 1D data
% line, we search for the nearest datapoint from the interpolated data

%Step through the 1D data points
for datapoint1D = 1:numpoints1D
    theta = startangle+(datapoint1D-1)*delta_theta;
    
    pointlocation = [R*cos(theta)+xc,R*sin(theta)+yc];
    
    scatter(R*cos(theta)+xc,R*sin(theta)+yc,'r');
    
    %Here we look through the list of points to see which is the closest to
    %our current point of interest
    distance = zeros(1,size(ZI,1)*size(ZI,2));
    for i = 1:size(ZI,1)*size(ZI,2)
            distance(i) = sqrt((XI(i)-pointlocation(1))^2+(YI(i)-pointlocation(2))^2);       
    end
    
    distance(find(distance==0))=NaN;
% 
% XI(1,:)
% YI(:,1)
% size(ZI)
%find(distance==min(distance))
    data1D(datapoint1D) = ZI(find(distance==min(distance)));
 %   scatter(XI(find(distance==min(distance))),YI(find(distance==min(distance))),'r');
 %   pause
end

% save the geometry settings to use for other data sets
save('geometry.mat','xc','yc','R','A','numpoints1D');

%close all

% Visualize line profile
data1D = fliplr(data1D)
hold off
for i = 1:size(timeaxis,2)
    plot(real(data1D.*exp(-complex(0,1)*center_frequency*timeaxis(i))));

    if sum(findstr(signalname,'Velocity'))~=0
        ylim([-max(abs(data1D)) max(abs(data1D))])
    else
        ylim([-6e-7 6e-7])
    end
    %pause(0.005);
    I(i) = getframe(gcf);
end