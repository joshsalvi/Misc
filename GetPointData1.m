% [x,y,usd] = GetPointData(filename, domainname, channelname, signalname, displayname,
%   point, frame)
% 
% Gets original or user defined data from a polytec file.
%
% filename is the path of the .pvd or .svd file
% domainname is the name of the domain, e.g. 'FFT' or 'Time'
% channelname is the name of the channel, e.g. 'Vib' or 'Ref1' or 'Vib &
%   Ref1'
% signalname is the name of the signal, e.g. 'Velocity' or 'Displacement'
% displayname is the name of the display, e.g. 'Real' or 'Magnitude' or
%   'Samples'
% point is the (1-based) index of the point to get data from. If point is
%   0 the data of all points will be returned. y will contain the data of
%   point i at row index i.
% frame is the frame number of the data. for data acquired in MultiFrame
%   mode, 0 is the averaged frame and 1-n are the other frames. For user
%   defined datasets the frame number is in the range 1-n where n is the
%   number of frames in the user defined dataset. For all other data,
%   use frame number 0.
%
% returns x, the x axis values of the data
% returns y, the data. colomns correspond to the x-axis, rows to the point
%   index. for point = 0 rows for points that have no data is set to zeros.
% returns usd, a struct describing the signal
%
function [x,y,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point, frame)
%
file = actxserver('PolyFile.PolyFile');
invoke(file, 'open', filename);
% ptcBuildDefault | ptcBuildUser
pointdomains = invoke(file, 'GetPointDomainsEx', 607);
pointdomain = get(pointdomains, 'Item', domainname);
channel = get(pointdomain.channels, 'Item', channelname);
signal = get(channel.signals, 'Item', signalname);
display = get(signal.displays, 'Item', displayname);
%
xaxis = invoke(pointdomain, 'GetXAxis', display);
yaxes = invoke(pointdomain, 'GetYAxes', display);
yaxis = get(yaxes, 'Item', 1);
x = xaxis.Min:(xaxis.Max - xaxis.Min)/(xaxis.MaxCount - 1):xaxis.Max;
%
if ~isempty(strfind(channel.caps, 'ptcChannelCapsUser'))
    % get the user signal description of the signal
    polymatlab = actxserver('PolyFile.PolyMatlab');
    usersignal = invoke(polymatlab, 'Cast', 'IUserSignal', signal);
    usd.Name = usersignal.UserSignalDesc.Name;
    usd.Complex = usersignal.UserSignalDesc.Complex;
    usd.PowerSignal = usersignal.UserSignalDesc.PowerSignal;
    usd.Is3D = usersignal.UserSignalDesc.Is3D;
    usd.DbReference = usersignal.UserSignalDesc.DbReference;
    usd.XName = usersignal.UserSignalDesc.XName;
    usd.XUnit = usersignal.UserSignalDesc.XUnit;
    usd.XMin = usersignal.UserSignalDesc.XMin;
    usd.XMax = usersignal.UserSignalDesc.XMax;
    usd.XCount = usersignal.UserSignalDesc.XCount;
    usd.YName = usersignal.UserSignalDesc.YName;
    usd.YUnit = usersignal.UserSignalDesc.YUnit;
    usd.YMin = usersignal.UserSignalDesc.YMin;
    usd.YMax = usersignal.UserSignalDesc.YMax;
    delete(polymatlab);
else
    % create the user signal description from the axis info
    usd.Name = signal.name;
    usd.Complex = 0;
    usd.PowerSignal =  ~isempty(strfind(signal.name, 'AP')) | ~isempty(strfind(signal.name, 'CP'));
    usd.Is3D = ~isempty(strfind(channel.caps, 'ptcChannelCapsVector'));
    usd.DbReference = Get0dBReference(file, signal);
    usd.XName = xaxis.name;
    usd.XUnit = xaxis.unit;
    usd.XMin = xaxis.min;
    usd.XMax = xaxis.max;
    usd.XCount = xaxis.maxcount;
    usd.YName = yaxis.name;
    usd.YUnit = yaxis.unit;
    usd.YMin = yaxis.min;
    usd.YMax = yaxis.max;
end    
%
datapoints = pointdomain.datapoints;
if (point == 0)
    % get data of all points
    y = zeros(datapoints.count, xaxis.MaxCount);
    for i=1:datapoints.count
        datapoint = get(datapoints, 'Item', i);
        if ~isempty(strfind(channel.caps, 'ptcChannelCapsUser'))        
            ytemp = double(invoke(datapoint, 'GetData', display, frame));
            if (length(ytemp) > 0)
                % not all points might have user defined data attached.
                % polytec file access returns and empty array in this case.
                y(i,1:length(ytemp)) = ytemp;
            end
        else
            try
                % ignore errors because of invalid data, the result row
                % will contain zeros in this case
                y(i,:) = double(invoke(datapoint, 'GetData', display, frame));
            catch
            end
        end
    end
else
    datapoint = get(datapoints, 'Item', point);
    y = double(invoke(datapoint, 'GetData', display, frame));
end
%
invoke(file, 'close');
delete(file);