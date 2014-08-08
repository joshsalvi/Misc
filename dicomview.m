function dicomview(dirName)
%DICOMVIEW opens a series of dicom images for viewing
%   DICOMVIEW('dir') opens a series of dicom images with from single
%   directory 'dir'. It creates a GUI with a slider for imagenumber and
%   image information displayed underneath.
%
%   qq=DCM; qq(2:end+1)=qq; qq(1).iFrame=[]; qq(1).image=ones(size(qq(1).image))*-32768; qq(1).imageAdj=qq(1).image; qq(1).contours = {};
%

if ischar(dirName)
    % Get filenames
    dirCont = dir([dirName '*.DCM']); % *.dcm directory contents
    if isempty(dirCont), dirCont = dir([dirName '*.dcm']); end
    Nframes = numel(dirCont); % number of frames
%     Nframes = 10;
    FN = dirCont(1).name;
    suffL = ceil(log10(Nframes)); % suffix length
    baseFN = fullfile(dirName,FN(1:(end-suffL-4))); % base filename without numeral suffixed
    % Generate analysis struct
    S = AnStruct(baseFN,suffL,Nframes);
else
    S = dirName;
end

% Store S
local_memory(S);

% Initialize GUI
local_init;


function S = AnStruct(baseFN,suffL,Nframes)
% Create struct used for all further analysis and containing all data

% Read DICOM images
S = struct(...
    'iFrame',{1:Nframes},...
    'FN',{},...
    'info',{},...
    'image',{},...
    'imageAdj',{},...
    'thrMin',{},...
    'thrMax',{},...
    'contourVal',{},...
    'contours',{},...
    'zoomXbase',{},...
    'zoomYbase',{},...
    'zoomX',{},...
    'zoomY',{}); % preallocation
for iFrame = 1:Nframes
    suff = sprintf(['%0' int2str(suffL) 'd'],iFrame-1);
    curFN = [baseFN suff '.dcm'];
    %
    D.iFrame = iFrame;
    D.FN = curFN;
    D.info = dicominfo(curFN);
    D.image = dicomread(curFN);
    D.imageAdj = D.image;
    D.thrMin = -32768;
    D.thrMax = 32767;
    D.contourVal = [];
    D.contours = [];
    D.zoomXbase = [];
    D.zoomYbase = [];
    D.zoomX = [];
    D.zoomY = [];
    %
    S(iFrame) = D;
end


function H = local_init
% Initialize GUI

S = local_memory;

H.fig = figure(...
    'Position',[185 142 840 800],...
    'IntegerHandle','off',...
    'NumberTitle','off',...
    'Toolbar','figure');

H.axes = axes(...
    'Parent',H.fig,...
    'Units','pixels',...
    'Position',[20 60 720 720],...
    'Box','on',...
    'XTick',[],...
    'YTick',[],...
    'NextPlot','add');

H.frmSldr = uicontrol(H.fig,...
    'Style','slider',...
    'Position',[20 20 720 20],...
    'Min',1,...
    'Max',numel(S),...
    'SliderStep',[1 10]/(numel(S)-1),...
    'Value',1,...
    'TooltipString','Slide to change frame number');

H.frmNr = uicontrol(H.fig,...
    'Style','text',...
    'Position',[750 20 70 20],...
    'FontSize',12,...
    'HorizontalAlignment','left',...
    'String','Frame');

H.thrMinSldr = uicontrol(H.fig,...
    'Style','slider',...
    'Position',[760 60 20 250],...
    'Min',-32768,...
    'Max',32767,...
    'SliderStep',[1 10]/65535,...
    'Value',-32768,...
    'TooltipString','Adjust minimum threshold');

H.thrMinVal = uicontrol(H.fig,...
    'Style','text',...
    'Position',[751 320 38 25],...
    'FontSize',10,...
    'String',{'Min';'-32768'});

H.thrMaxSldr = uicontrol(H.fig,...
    'Style','slider',...
    'Position',[800 60 20 250],...
    'Min',-32768,...
    'Max',32767,...
    'SliderStep',[1 10]/65535,...
    'Value',32767,...
    'TooltipString','Adjust maximum threshold');

H.thrMaxVal = uicontrol(H.fig,...
    'Style','text',...
    'Position',[791 320 38 25],...
    'FontSize',10,...
    'String',{'Max';'32767'});

H.cntrSldr = uicontrol(H.fig,...
    'Style','slider',...
    'Position',[760 365 20 250],...
    'Min',-32768,...
    'Max',32767,...
    'SliderStep',[1 10]/65535,...
    'Value',0,...
    'TooltipString','Contour detection threshold');

H.cntrVal = uicontrol(H.fig,...
    'Style','text',...
    'Position',[790 490 38 25],...
    'FontSize',10,...
    'String',{'Contour';'0'});

H.output = uicontrol(H.fig,...
    'Style','pushbutton',...
    'Position',[760 760 60 20],...
    'TooltipString','Output struct to workspace',...
    'String','Output');

H.info = uicontrol(H.fig,...
    'Style','pushbutton',...
    'Position',[760 735 60 20],...
    'TooltipString','Display information in prompt',...
    'String','Info');

H.original = uicontrol(H.fig,...
    'Style','pushbutton',...
    'Position',[760 710 60 20],...
    'TooltipString','Show original image',...
    'String','Original');

H.lockzoom = uicontrol(H.fig,...
    'Style','pushbutton',...
    'Position',[760 685 60 20],...
    'TooltipString','Lock zoom level',...
    'String','Lock zoom');

H.resetzoom = uicontrol(H.fig,...
    'Style','pushbutton',...
    'Position',[760 660 60 20],...
    'TooltipString','Reset zoom level',...
    'String','Reset zoom');

H.apply2all = uicontrol(H.fig,...
    'Style','pushbutton',...
    'Position',[760 635 60 20],...
    'TooltipString','Apply threshold to all slides',...
    'String','Apply');

set(H.fig,'DeleteFcn',@(src,event) local_outputCB(H,src,event))
set(H.frmSldr,'Callback',@(src,event) local_frmSldrCB(H,src,event))
set(H.thrMinSldr,'Callback',@(src,event) local_thrMinSldrCB(H,src,event))
set(H.thrMaxSldr,'Callback',@(src,event) local_thrMaxSldrCB(H,src,event))
set(H.cntrSldr,'Callback',@(src,event) local_cntrSldrCB(H,src,event))
set(H.output,'Callback',@(src,event) local_outputCB(H,src,event))
set(H.info,'Callback',@(src,event) local_infoCB(H,src,event))
set(H.original,'Callback',@(src,event) local_originalCB(H,src,event))
set(H.lockzoom,'Callback',@(src,event) local_lockzoomCB(H,src,event))
set(H.resetzoom,'Callback',@(src,event) local_resetzoomCB(H,src,event))
set(H.apply2all,'Callback',@(src,event) local_apply2allCB(H,src,event))

% Single call to intialize with first frame
local_frmSldrCB(H,H.frmSldr,[])


function local_frmSldrCB(object,src,event) %#ok<*INUSD>
% Frame slider callback

% Retreive S
S = local_memory;

% Frame
frmNum = round(get(src,'Value')); % frame number
curS = S(frmNum); % current frame

% Display frame number
set(object.frmNr,'String',['Frame ' int2str(curS.iFrame)]);

% Display window name
set(object.fig,'Name',curS.FN)

% Clear axes
cla(object.axes)

% Display adjusted image
imshow(curS.imageAdj,'Parent',object.axes)
axis(object.axes,'image')

% Store base zoom level
if isempty(curS.zoomXbase)
    [S.zoomXbase] = deal(get(object.axes,'XLim'));
    [S.zoomYbase] = deal(get(object.axes,'YLim'));
end

% Set zoom
if ~isempty(curS.zoomX)
    set(object.axes,'XLim',curS.zoomX,'YLim',curS.zoomY)
end

% Plot contours if available
if ~isempty(curS.contours)
    cellfun(@(C) plot(C(1,:),C(2,:),'r','LineWidth',2),curS.contours)
end

% Store analysis
local_memory(S);


function local_thrMinSldrCB(object,src,event)
% Threshold minimum slider callback

% Retreive S
S = local_memory;

% Frame
frmNum = round(get(object.frmSldr,'Value')); % frame number
curS = S(frmNum); % current frame
curIm = curS.image; % current image

% Thresholds
thrMin = get(src,'Value');
thrMax = get(object.thrMaxSldr,'Value');

% Correct indicator
set(object.thrMinVal,'String',{'Min';int2str(thrMin)})

% Thresholding
curIm(curIm<=thrMin) = -32768;
curIm(curIm>=thrMax) = 32767;

% Adjust analysis struct
S(frmNum).thrMin = thrMin;
S(frmNum).imageAdj = curIm;
[S.zoomX] = deal(get(object.axes,'XLim'));
[S.zoomY] = deal(get(object.axes,'YLim'));

% Store analysis
local_memory(S);

% Display changes
local_frmSldrCB(object,object.frmSldr,[])


function local_thrMaxSldrCB(object,src,event)
% Threshold maximum slider callback

% Retreive S
S = local_memory;

% Frame
frmNum = round(get(object.frmSldr,'Value')); % frame number
curS = S(frmNum); % current frame
curIm = curS.image; % current image

% Thresholds
thrMin = get(object.thrMinSldr,'Value');
thrMax = get(src,'Value');

% Correct indicator
set(object.thrMaxVal,'String',{'Max';int2str(thrMax)})

% Thresholding
curIm(curIm<=thrMin) = -32768;
curIm(curIm>=thrMax) = 32767;

% Adjust analysis struct
S(frmNum).thrMax = thrMax;
S(frmNum).imageAdj = curIm;
[S.zoomX] = deal(get(object.axes,'XLim'));
[S.zoomY] = deal(get(object.axes,'YLim'));

% Store analysis
local_memory(S);

% Display changes
local_frmSldrCB(object,object.frmSldr,[])


function local_cntrSldrCB(object,src,event)
% Contour detection slider callback

% Retreive S
S = local_memory;

% Contour detection threshold
thresh = get(src,'Value');

% Frame
frmNum = round(get(object.frmSldr,'Value')); % frame number
curS = S(frmNum); % current frame
curImAdj = curS.imageAdj; % current adjusted image

% Detect contours
C = contourc(double(curImAdj),[1 1]*thresh);
iStart = find(C(1,:)==thresh);
conts = arrayfun(@(indx) C(:,indx+1:indx+C(2,indx)),iStart,'UniformOutput',false);

% Save contours
S(frmNum).contourVal = thresh;
S(frmNum).contours = conts;
local_memory(S);

% Correct indicator
set(object.cntrVal,'String',{'Contour';int2str(thresh)})

% Display changes
local_frmSldrCB(object,object.frmSldr,[])


function local_outputCB(object,src,event)
% Output button callback

% Retreive S
S = local_memory;

% Write S to workspace
assignin('base','DCM',S)


function local_infoCB(object,src,event)
% Info button callback

% Retreive S
S = local_memory;

% Frame
frmNum = round(get(object.frmSldr,'Value')); % frame number
curS = S(frmNum); % current frame

% Display
disp(curS.info)


function local_originalCB(object,src,event)
% Original button callback

% Retreive S
S = local_memory;

% Frame
frmNum = round(get(object.frmSldr,'Value')); % frame number
curS = S(frmNum); % current frame

% Display image
hfig = figure('Visible','off');
hax = axes('Parent',hfig,'XTick',[],'YTick',[]);
imshow(curS.image,'Parent',hax)
set(hfig,'Position',[1032 494 420 420],'Visible','on');


function local_lockzoomCB(object,src,event)
% Lock zoom button callback

% Retreive S
S = local_memory;

% Get zoom
[S.zoomX] = deal(get(object.axes,'XLim'));
[S.zoomY] = deal(get(object.axes,'YLim'));

% Store analysis
local_memory(S);

% Redraw
local_frmSldrCB(object,object.frmSldr,[])


function local_resetzoomCB(object,src,event)
% Reset zoom button callback

% Retreive S
S = local_memory;

% Frame
frmNum = round(get(object.frmSldr,'Value')); % frame number
curS = S(frmNum); % current frame

% Display image
set(object.axes,'XLim',curS.zoomXbase,'YLim',curS.zoomYbase);

% Adjust S
[S.zoomX] = deal(S(frmNum).zoomXbase);
[S.zoomY] = deal(S(frmNum).zoomYbase);

% Store analysis
local_memory(S);


function local_apply2allCB(object,src,event)
% Apply to all button callback

% Retreive S
S = local_memory;

% Thresholds
thrMin = get(object.thrMinSldr,'Value');
thrMax = get(object.thrMaxSldr,'Value');
thresh = get(object.cntrSldr,'Value');

% Frame
for iFrame = 1:max([S.iFrame])
    curS = S(iFrame); % current frame
    curIm = curS.image; % current image
    % Thresholding
    curIm(curIm<=thrMin) = -32768;
    curIm(curIm>=thrMax) = 32767;
    % Adjust analysis struct
    S(iFrame).thrMin = thrMin;
    S(iFrame).thrMax = thrMax;
    S(iFrame).imageAdj = curIm;
    % Detect contours
    curImAdj = curIm; % current adjusted image
    C = contourc(double(curImAdj),[1 1]*thresh);
    iStart = find(C(1,:)==thresh);
    conts = arrayfun(@(indx) C(:,indx+1:indx+C(2,indx)),iStart,'UniformOutput',false);
    S(iFrame).contourVal = thresh;
    S(iFrame).contours = conts;
end

% Store analysis
[S.zoomX] = deal(get(object.axes,'XLim'));
[S.zoomY] = deal(get(object.axes,'YLim'));
local_memory(S);

% Display changes
local_frmSldrCB(object,object.frmSldr,[])


function D = local_memory(D)
% Store analysis struct for editing

persistent S
if nargin<1 % get
    D = S;
else % set
    S = D;
end










