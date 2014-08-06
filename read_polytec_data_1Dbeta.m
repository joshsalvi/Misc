% read_polytec_data_1D.m is an mfile for reading polytec data generated
% by a 2D scanning interferometer looking at an arc segment of a turn of the
% gerbil's cochlea.  It then reduces the arc data into a 1D line. The
% purpose is to create 1D data that can be used to generate a gain function
% (c.f. C. Shera, J. Acoust. Soc. Am. Volume 122, Issue 5, pp. 2738-2758).

% NOTE: This program requires that Polytec-created mfiles GetPointData.m and Get0dBReference.m be
% in the same folder as this mfile. 
%EXAMPLE
%[pointdata1,img1,I1,timeaxis1]= read_polytec_data('C:\ProgramFiles\MATLAB\R2008b\3khz\3khz', 'FFT', 'Vib','Displacement','Real & Imag.', 1.55e6,1);
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

% OUTPUTS:
% pointdata: contains all locations, spectral info, indexes for all points;
% I: movie frames for the traveling wave animation
% A: the complex velocity (or displacement) coefficients for each point

%Sample syntax: 
%for FFT
%[pointdata,I, A]= read_polytec_data_1Dbeta('E:\20oct10\5-30kHz-0.2k', 'FFT', 'Vib','Velocity','Real & Imag.', 1.55e6,1,1);
%for Time
%[pointdata,img,I,timeaxis] = read_polytec_data('g1-1_plus_timedata','Time','Usr','FFT / Vib / Velocity', 'Samples', 1.55e6);
%J.A.N.Fisher 2010

function [pointdata,I, A] = read_polytec_data_1Dbeta(basename, domainname, channelname, signalname, displayname,pix_per_m,firstfile,newsettings)

%We first need to know the number of points and the location, in the field
%of view, of these points. This we get from the txt file version of the
%.svd file.
[index, x, y, z, magnitude] = textread([basename '.txt'],'%d %f %f %d %f',10000,'headerlines',11);


numpoints = size(x,1);

pointdata.index = index;
pointdata.location(:,:) = [x y];
pointdata.bandmag = magnitude;

    [x,y,usd] = GetPointData2([basename '.svd'], domainname, channelname, signalname, displayname, 0);

    %3/24/10 commenting out this case structure
%     
% if sum(findstr(domainname,'FFT'))==0
%     %that is, if the domain name IS EQUAL to the string 'Time'
%     [x,y,usd] = GetPointData1([basename '.svd'], domainname, channelname, signalname, displayname, 0,1);
% %the following 4 lines are done just because matlab can't put together a
% %movie with thousands of frames due to memory limitations
%     y = y(:,1:200);
%     x = x(1:200);
%     size(x)
%     size(y)
% else
% 
%     for i = 1:numpoints
%         if i == 1
%         [x,y,usd] = GetPointData2([basename '.svd'], domainname, channelname, signalname, displayname, pointdata.index(i));    
% %    [x,y,usd] = GetPointData([basename '.svd'], domainname, channelname, signalname, displayname, 0);
%         else
%         [x,y1,usd] = GetPointData2([basename '.svd'], domainname, channelname, signalname, displayname, pointdata.index(i));
%         y = cat(1,y,y1);
%         end
% 
% 
%     end
% end

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



%     %pointdata.FFT_data(i,:,:) = [x',sqrt(cy'.*conj(cy'))];
% %    TD_data = ifft(cy','symmetric');
% %    pointdata.TD_data(i,:,:) = [x',TD_data];
% 
%     %figure,plot(pointdata.spectrum(:,1),pointdata.spectrum(:,2));
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
            pointdata.spectrumwindow = spect.range;
            max_spect_index = spect.max_spect_index;
            center_frequency = spect.center_frequency;
        end

        av_freq = sum(spectwindow)/size(spectwindow,1)
        %timeaxis = 0:2*pi/(100*av_freq):2*pi/av_freq;
        timeaxis = 0:2*pi/(100*av_freq):2*pi/av_freq;
        size(timeaxis)
        size(spectwindow)

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
%pointdata.location(:,:) = -pointdata.location(:,:);

%Mirror the data over one axis if the handedness is LEFT.
% if findstr(hand,'L') == 1
%     pointdata.location(:,1) = -pointdata.location(:,1);
% else
% end

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

%pointdata.location = abs(pointdata.location);

%now we generate a meshgrid on a finer scale than the data, say every 5
%pixels
pix_xdim = max(pointdata.location(:,1))-min(pointdata.location(:,1));
pix_ydim = max(pointdata.location(:,2))-min(pointdata.location(:,2));

tx = min(pointdata.location(:,1)):20:max(pointdata.location(:,1)); 
ty = min(pointdata.location(:,2)):20:max(pointdata.location(:,2));
[XI,YI] = meshgrid(tx,ty);

if sum(findstr(domainname,'FFT'))~=0
    %Frequency domain case
    A(:,:) = [min(min(pointdata.TD_data(:,:,2))),max(max(pointdata.TD_data(:,:,2)))];
    figure
    for i = 1:size(timeaxis,2)
    %    Tri=delaunay(pointdata.location(:,1),pointdata.location(:,2));
    %    scatter3(pointdata.location(:,1),pointdata.location(:,2),pointdata.TD_data(:,i,2),'b','markerface','r');

    %    trisurf(Tri,pointdata.location(:,1),pointdata.location(:,2),pointdata.TD_data(:,i,2));
        ZI = griddata(pointdata.location(:,1),pointdata.location(:,2),pointdata.TD_data(:,i,2),XI,YI,'cubic');

%        surf(XI,YI,ZI)
        surf(XI,YI,ZI,'FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',ZI)
    %     size(XI)
    %     size(YI)
    %     size(ZI)
    %    scatter3(XI(1,:),YI(1,:),ZI(1,:));\
    %    plot3(x,y,z,'o'), hold off
        view(-20.5,60)
        set(gca,'CLim',A);
        zlim(A);
        I(i) = getframe;
    end
else
    %Time-domain case
%    A(:,:) = [min(min(pointdata.TD_data(:,:,2))),max(max(pointdata.TD_data(:,:,2)))];
    A(:,:) = [-5e-7,5e-7];
    figure
    for i = 1:size(x,2)
    %    Tri=delaunay(pointdata.location(:,1),pointdata.location(:,2));
    %    scatter3(pointdata.location(:,1),pointdata.location(:,2),pointdata.TD_data(:,i,2),'b','markerface','r');

    %    trisurf(Tri,pointdata.location(:,1),pointdata.location(:,2),pointdata.TD_data(:,i,2));
        ZI = griddata(pointdata.location(:,1),pointdata.location(:,2),pointdata.TD_data(:,i,2),XI,YI,'cubic');


        surf(XI,YI,ZI)
    %     size(XI)
    %     size(YI)
    %     size(ZI)
    %    scatter3(XI(1,:),YI(1,:),ZI(1,:));\
    %    plot3(x,y,z,'o'), hold off
        view(-20.5,60)
        set(gca,'CLim',A);
        zlim(A);
        I(i) = getframe;
    end
end

%3/24/2010 I copied and pasted the following code for fitting a circle to
%the data.

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

% %flip startpoint if LEFT handed data
% if findstr(hand,'L')==1
%     startpoint(2) = -startpoint(2);
% else
% end

% if findstr(hand,'L')==1
%     startpoint = startpoint-[xc,yc];
% else
%     startpoint = startpoint-[yc,xc];
% end

figure, scatter(pointdata.location1(:,1),pointdata.location1(:,2));
axis equal
hold on


pointdata.location = pointdata.location1;

%create a scale to plot the circle on (on the same axes as the
%scatterplotted pointdata)
ts = min(min(pointdata.location1(:,1))):max(max(pointdata.location1(:,1)));

plot(ts,sqrt(R^2-(ts).^2));

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
            %ring sector
            for n = 1:numpoints
               
                if sqrt(pointdata.location(n,1)^2+pointdata.location(n,2)^2)>=radius_inner && sqrt(pointdata.location(n,1)^2+pointdata.location(n,2)^2)<=radius_outer && atan2(pointdata.location(n,2),pointdata.location(n,1))<(theta+delta_theta) && atan2(pointdata.location(n,2),pointdata.location(n,1))>(theta) 
                    dArea_points(counter).FFT_data = pointdata.FFT_data(n,:,:);
                    dArea_points(counter).index = n;
                    scatter(pointdata.location(n,1),pointdata.location(n,2),'r');
                    counter = counter+1;
                else
                end
            end

        if counter ~=1

            for p = 1:counter-1

                scatter(pointdata.location(dArea_points(p).index,1),pointdata.location(dArea_points(p).index,2),'r');
 %               size(pointdata.location)
                if p==1
                    data1D(datapoint1D).FFT_data = dArea_points(p).FFT_data;
                else
                data1D(datapoint1D).FFT_data = data1D(datapoint1D).FFT_data+dArea_points(p).FFT_data;
                end
            end
                data1D(datapoint1D).FFT_data = data1D(datapoint1D).FFT_data/counter;
            pause(0.01);
        else
            data1D(datapoint1D).FFT_data = zeros(1,size(pointdata.spectrum(:,2)),2);
        end
end

% save the geometry settings to use for other data sets
save('geometry.mat','xc','yc','R','A','radius_outer','radius_inner','theta1','theta2','dtheta1','dtheta2','delta_theta','numpoints1D');

A = zeros(1,size(data1D,2));
for i = 1:size(data1D,2)
    A(i)=data1D(i).FFT_data(1,max_spect_index,2);
end

A = boxcar(A,3);

figure
for i = 1:size(timeaxis,2)
    plot(real(A.*exp(-complex(0,1)*center_frequency*timeaxis(i))));
    ylim([-10e-4 10e-4])
%    xlim([0 size(data1D,2)])
    pause(0.01);
    I(i) = getframe(gcf);
end