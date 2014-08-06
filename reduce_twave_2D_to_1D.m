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

%Sample syntax:
%>> [pointdata1,img,I1,timeaxis1] =
%reduce_twave_2D_to_1D('D:\MATLAB_PATCH_1,5khz\1,5khz0',
%'FFT','Vib','Displacement', 'Real & Imag.', 1.55e6);


function [pointdata, A] = reduce_twave_2D_to_1D(basename, domainname, channelname, signalname, displayname,pix_per_m,firstfile, newsettings)

[x,y,usd] = GetPointData2([basename '.svd'], domainname, channelname, signalname, displayname, 0);

numpoints = size(y,1);

%now put complex values together as complex numbers (real and imag
%components are stored in neighboring data pairs in the polytec format)
if sum(findstr(displayname,'Imag'))~=0
    for j = 1:size(x,2)
        
        iy(:,j) = y(:,2*j);
        ry(:,j) = y(:,2*(j-1)+1);
        cy(:,j) = complex(ry(:,j),iy(:,j));
        % size(cy)
    end
else
    cy = y;
end

if sum(findstr(domainname,'FFT'))~=0
    %We assign the FFT data to a portion of this cluster and also take the inverse FFT to make the Time Domain data (TD)
    
    for i = 1:size(y,1)
        pointdata.FFT_data(i,:,:) = [x',cy(i,:)'];
    end
    pointdata.spectrum(:,:) = sum(abs(pointdata.FFT_data),1)./numpoints;
else
    pointdata.FFT_data = zeros(size(y,1),size(x,2),size(y,2));
    pointdata.FFT_data = [0];
    pointdata.spectrum(:,:) = [0];
end

pointdata.info = usd;

%get the pointdata locations:
XYZ = GetXYZCoordinates([basename '.svd'], 0);

pointdata.location(:,:) = [XYZ(:,1) XYZ(:,2)];


if sum(findstr(domainname,'FFT'))~=0
    %frequency domain case
    
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
        center_frequency = pointdata.spectrum(max_spect_index,1);
        
        range = pointdata.spectrumwindow;
        save('spect.mat','spectwindow','range','max_spect_index','center_frequency');
        
    else
        spect = load('spect.mat','spectwindow','range','max_spect_index','center_frequency');
        spectwindow = spect.spectwindow;
        max_spect_index = spect.max_spect_index;
        pointdata.spectrumwindow = spect.range;
        center_frequency = spect.center_frequency;
    end
    
    av_freq = sum(spectwindow)/size(spectwindow,1)
    %timeaxis = 0:2*pi/(100*av_freq):2*pi/av_freq;
     timeaxis = 0:2*pi/(100*av_freq):2*pi/av_freq;
    
    pointdata.TD_data = zeros(numpoints,size(timeaxis,2),2);
    for j = 1:numpoints
        %for j = 1
        for p = 1:size(spectwindow,1)
            %pointdata.FFT_data has three dimensions: (index, datapoint,
            %data).  Index is the number of the scanned point. datapoint is
            %the number of FFT lines. data is 2 pieces of info: the FFT
            %line frequency and the FFT coefficient.
            
            FFT_freq = pointdata.FFT_data(j,pointdata.spectrumwindow(1)+(p-1),1);
            FFT_coeff = pointdata.FFT_data(j,pointdata.spectrumwindow(1)+(p-1),2);
            
            pointdata.TD_data(j,:,1) = timeaxis;
            if p == 1
                p
                j
                
                pointdata.TD_data(j,:,2)= real(FFT_coeff*exp(-complex(0,1)*FFT_freq*timeaxis));
            else
                %following the initial frequency sinusoid, all others
                %(in the window) are summed on top of that.
                pointdata.TD_data(j,:,2)= pointdata.TD_data(j,:,2)+real(FFT_coeff*exp(-complex(0,1)*FFT_freq*timeaxis));
            end
        end
    end
else
    for j = 1:numpoints
        %Time-domain case
        size(x')
        size(cy(j,:)')
        pointdata.TD_data(j,:,:)= [x',cy(j,:)'];
    end
    
    timeaxis = x;
end


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

%now we generate a meshgrid on a finer scale than the data, say every 5
%pixels
pix_xdim = max(pointdata.location(:,1))-min(pointdata.location(:,1));
pix_ydim = max(pointdata.location(:,2))-min(pointdata.location(:,2));

%Now fit the data to a circle 

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
    radius_outer = geometry.radius_outer;
    radius_inner = geometry.radius_inner;
    theta1 = geometry.theta1;
    theta2 = geometry.theta2;
    dtheta1 = geometry.dtheta1;
    dtheta2 = geometry.dtheta2;
    delta_theta = geometry.delta_theta;
    numpoints1D = geometry.numpoints1D;
else
    
    %choose the points that correspond to an arc
    D = ginput;
    
    %use circfit to come up with a circular fit. The function circfit.m must be
    %in the matlab search paths.
    [xc, yc, R,A] = circfit(D(:,1),D(:,2))

end


%subtract R from all the points
for u = 1:size(pointdata.location,1)
    pointdata.location1(u,:) = pointdata.location(u,:)-[xc,yc];
end


figure, scatter(pointdata.location1(:,1),pointdata.location1(:,2));
axis equal
hold on

pointdata.location = pointdata.location1;

%create a scale to plot the circle on (on the same axes as the
%scatterplotted pointdata)
ts = min(min(pointdata.location1(:,1))):max(max(pointdata.location1(:,1)));

% plot(ts,sqrt(R^2-(ts).^2));


%Again, if using previous geometry settings, load those:
if newsettings ~=0
    %choose the points correspond to the outer, then inner radii
    title('choose the points correspond to the outer, then inner radii')
    D = ginput;
    radius_outer = sqrt((D(1,1))^2+(D(1,2))^2)
    radius_inner = sqrt((D(2,1))^2+(D(2,2))^2)
    plot(ts,sqrt(radius_outer^2-(ts).^2),'r');
    plot(ts,sqrt(radius_inner^2-(ts).^2),'g');

    %now select the beginning and ending of the arc; these points should be
    %selected in order of increasing angle (i.e. counter-clockwise)
    title('choose the points correspond to the beginning and ending of the arc segment')
    D = ginput;
    theta1 = atan2(D(1,2),D(1,1))
    theta2 = atan2(D(2,2),D(2,1))

    %now determine the angular resolution by finding the angle between
    %two neighboring closely packed data points.
    title('choose two consecutive data points along the rings arc')
    D = ginput;
    dtheta1 = atan2(D(1,2),D(1,1));
    dtheta2 = atan2(D(2,2),D(2,1));

    %dtheta now determines the total number of points on the 1D line, that is,
    %numpoints1D = abs(theta1-theta2)/delta_theta. I've arbitrarily set the delta
    %to correspond to the distance between two selected neighboring
    %datapoints.
    delta_theta = abs(dtheta1-dtheta2);
    numpoints1D = abs(theta1-theta2)/delta_theta;
else
end

% Here, for each delta_theta corresponding to one data point on the 1D data
% reduction, we search for included data points within the area on the ring
% sector.

%Here we step through the 1D data points
for datapoint1D = 1:numpoints1D
    theta = theta1+(datapoint1D-1)*delta_theta;
    counter = 1;
    % %     plot(-2000:0,tan(theta)*(-2000:0));
    %      plot(-2000:0,tan(theta+delta_theta)*(-2000:0));
    %      pause
    %Here we check the list of points to see if they are included in that
    %ring sector; points must be within the inner and outer radii,
    %and must be at an angle between theta and theta+delta_theta
    for n = 1:numpoints
        
        if sqrt(pointdata.location(n,1)^2+pointdata.location(n,2)^2)>=radius_inner && sqrt(pointdata.location(n,1)^2+pointdata.location(n,2)^2)<=radius_outer && atan2(pointdata.location(n,2),pointdata.location(n,1))<(theta+delta_theta) && atan2(pointdata.location(n,2),pointdata.location(n,1))>(theta)
            dArea_points(counter).FFT_data = pointdata.FFT_data(n,:,:);
            dArea_points(counter).index = n;
    %        scatter(pointdata.location(n,1),pointdata.location(n,2),'r');
            counter = counter+1;
        else
        end
    end
    
    if counter ~=1
        
        for p = 1:counter-1
            
 %           scatter(pointdata.location(dArea_points(p).index,1),pointdata.location(dArea_points(p).index,2),'r');
            %               size(pointdata.location)
            if p==1
                data1D(datapoint1D).FFT_data = dArea_points(p).FFT_data;
            else
                data1D(datapoint1D).FFT_data = data1D(datapoint1D).FFT_data+dArea_points(p).FFT_data;
            end
        end
        data1D(datapoint1D).FFT_data = data1D(datapoint1D).FFT_data/counter;
     %   pause
    else
        data1D(datapoint1D).FFT_data = zeros(1,size(pointdata.spectrum(:,2)),2);
    end
end

% save the geometry settings to use for other data sets
save('geometry.mat','xc','yc','R','A','radius_outer','radius_inner','theta1','theta2','dtheta1','dtheta2','delta_theta','numpoints1D');

close all

%now we pick off just the FFT coefficients and put them in a 1D variable
%which gets outputed by the function
A = zeros(1,size(data1D,2));
for i = 1:size(data1D,2)
    A(i)=data1D(i).FFT_data(1,max_spect_index,2);
end

A = boxcar(A,5);

%this following code can be used if you want to visualize the 1D data
%evolving in time
figure
for i = 1:size(timeaxis,2)
    plot(real(A.*exp(-complex(0,1)*center_frequency*timeaxis(i))));
    ylim([-max(abs(A)) max(abs(A))])
    pause(0.05);
%    I(i) = getframe(gcf);
end