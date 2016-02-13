function importVICdata(datapath)
% This script imports data, without spacefiles, LabVIEW's
% variableloadclamp.vi to MATLAB. You will need to input the full directory
% name containing all files with your data, including the logfile.
%
% importVLCdata(datapath)
%
% datapath : string of the path containing your data
% MSVLC : Mech Stim or VLC? (1=VLC, 2=Mech Stim)
%
% Example:
% importVLCdata('/Users/username/Data/2014-09-30.01/Ear 1/Cell 2/')
%
% In some cases, the logfile will have incorrect indices. This becomes
% obvious for cases where "Fs" (the sampling rate) and "pulse" are not
% their expected values. If this occurs, modify the script and run it
% again.
%
% jsalvi@rockefeller.edu
%


% Find all appropriate files
file = dir(sprintf('%s%s',datapath,'*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
a = length(file);       % number of sessions in the directory
%{
% Find all appropriate files
file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.csv'));    % find the logfile
a = length(file);       % number of sessions in the directory
%}
% Import logdata
logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
if isstruct(logdata) == 0
%    logdata.data = logdata;
end
if isempty('logdata.textdata(isnan(logdata.data(:,8))==0,3)')==0
    comments = logdata.textdata(isnan(logdata.data(:,8))==0,3); % import comments
end

Fs = logdata.data(:,11);


% Note that some logfiles will have formatting issues and the indices
% listed above may be off by ±10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.



%Number of traces
ntrace=logdata.data(isnan(logdata.data(:,8))==0,8);
numavg=logdata.data(isnan(logdata.data(:,8))==0,10);
% Find all raw files
for j = 1:a
    rawfiles(j) = isempty(findstr(file(j).name,'raw'));
end
nonraw = find(rawfiles==1);raw = find(rawfiles==0);
ntraceraw = ntrace(ntrace~=0);
numavgraw = numavg(numavg~=0);

% Import the data, some of this may be redundant
for i = 1:a
    data2 = importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
    try
    if size(data2.data,2) == 5
        vec1{i} = data2.data(:,1)./1e3;
        Xc{i} = data2.data(:,2);
        Xd{i} = data2.data(:,3);
        Vc{i} = data2.data(:,4);
        Ir{i} = data2.data(:,5);
    elseif size(data2.data,2) == 4
        vec1{i} = data2.data(:,1)./1e3;
        Xd{i} = data2.data(:,2);
        Vc{i} = data2.data(:,3);
        Ir{i} = data2.data(:,4);
    elseif size(data2.data,2) == 3
        vec1{i} = data2.data(:,1)./1e3;
        Vc{i} = data2.data(:,2);
        Ir{i} = data2.data(:,3);
    else
        Vc{i} = data2.data(:,1);
        Ir{i} = data2.data(:,2);
    end
    headers{i} = data2.textdata;    % check headers to ensure Vr/Ir are correct, will be switched in current clamp
    end
end
clear data2

% Generate time vector
dt = 1./Fs;
for j = 1:a
    try
        tvec{j} = 0:dt(j):(length(vec1{j})-1)/Fs(j);
    end
    try
        comments{j} = ['(' num2str(j) '): ' comments{j}];
    end
end

clear i j 

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');

end
