% Gets original or user defined data from all the points from a Polytec file.
%
% filename is the path of the .pvd or .svd file
% domainname is the name of the domain, e.g. 'FFT' or 'Time'
% channelname is the name of the channel, e.g. 'Vib' or 'Ref1' or 'Vib &
% Ref1'
% signalname is the name of the signal, e.g. 'Velocity' or 'Displacement'
% displayname is the name of the display, e.g. 'Real' or 'Magnitude' or
% 'Samples'
% point is the index of the point to get data from. If point is
% 0 the data of all points will be returned. y will contain the data of
% point i at row index i.
%
% returns x, the x axis values of the data
% returns y, the data. colomns correspond to the x-axis, rows to the point
% index.
% returns usd, a struct describing the signal
%
function [x,y,usd] = GetPointData(filename, domainname, channelname, signalname, displayname, point)
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
    %usd.DbReference = Get0dBReference(file, signal);
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
datapoints.count;
if (point == 0)
    % get data of all points
    %y = zeros(datapoints.count, xaxis.MaxCount);
    if (strcmp(displayname,'Real & Imag.')==1);
        y=zeros(datapoints.count,xaxis.MaxCount*2);
    else
        y = zeros(datapoints.count, xaxis.MaxCount);
    end
        
    for i=1:datapoints.count
        datapoint = get(datapoints, 'Item', i);
        y(i,:) = double(invoke(datapoint, 'GetData', display, 0));        
    end
else
    datapoint = get(datapoints, 'Item', point);
    y = double(invoke(datapoint, 'GetData', display, 0));
end
%
invoke(file, 'close');
delete(file);
