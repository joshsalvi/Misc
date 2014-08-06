% XYZ = GetXYZCoordinates(filename, point)
% ----------------------------------------
% You have to install
% "http://www.mathworks.com/support/solutions/data/1-1BFEJ.html?solution=1-1BFEJ" for this to work.
% 
% Gets XYZ coordinates of the scan points from a polytec file.
%
% This is only possible for files containing 3D geometry or that have
%   a distance to the object specified. Otherwise there will be an
%   error message.
%
% filename is the path of the .svd file
% 
% point is the (1-based) index of the point to get the coordinates from. If point is
%   0 the coordiantes of all points will be returned. XYZ will contain the data of
%   point i at row index i.
%
% returns the xyz coordinates. columns correspond to the geometry X, Y, Z
%   in meter, rows to the point index.
%
function XYZ = GetXYZCoordinates(filename, point)
%

file = actxserver('PolyFile.PolyFile');
try
    invoke(file, 'Open', filename);
    infos=file.infos;
    measpoints = invoke(infos, 'Item', 'MeasPoints');
 
    if (point == 0)
        XYZ=zeros(measpoints.count,3);
        for i=1:measpoints.count
            measpoint=get(measpoints, 'Item', int32(i));
            [X,Y,Z]=invoke(measpoint,'CoordXYZ');
            XYZ(i,:)=[X,Y,Z];                  
        end
    else
        measpoint=get(measpoints, 'Item', int32(point));
        [X,Y,Z]=invoke(measpoint,'CoordXYZ');
        XYZ=[X,Y,Z];
    end
    invoke(file, 'Close');
    delete(file);
catch
    if file.IsOpen == 1
        invoke(file, 'Close');
    end
    delete(file);
    rethrow(lasterror);
end