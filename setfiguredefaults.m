function setfiguredefaults(varargin)
% Set the figure defaults, with axes order for N axes
%
% setfiguredefaults(N,dockyn)
%
% where N is the number of axes for your colormap
% If N is not included, the default line color will be black.
%
% If dockyn == 1, the figures will be docked by default. Otherwise,
% figures will remain undocked.
%
% Joshua Salvi
% jsalvi@rockefeller.edu
%
%
% NOTE: Requires the 'ametrine' colormap
%
figure;
if nargin > 0
    try
        set(0,'defaultAxesColorOrder',ametrine(varargin{1}))
    catch
        disp('Invalid number of axes. Defaulting to zero.');
        set(0,'defaultAxesColorOrder',[0 0 0])
    end
else
    set(0,'defaultAxesColorOrder',[0 0 0])
end
set(0,'defaultFigureColormap',jet)
set(0,'DefaultTextFontName','Helvetica Neue')
set(0,'DefaultTextFontUnits','Points')
set(0,'defaultlinelinewidth',2)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontName','Helvetica Neue')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesGridLineStyle',':')
set(0,'defaultFigureColor','White')
set(0,'DefaultFigurePaperType','a4letter')
set(0,'defaultFigurePosition',[1200 100 700 700])
set(0,'defaultAxesColor',[1 1 1])
set(0,'defaultLineMarker','None')
if nargin == 2
    if varargin{2} == 1
        set(0,'DefaultFigureWindowStyle','docked')
    else
        set(0,'DefaultFigureWindowStyle','normal')
    end
end
close;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% ametrine() %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmap=ametrine(n,varargin)
p=inputParser;
p.addParamValue('gamma',1.8, @(x)x>0);
p.addParamValue('minColor','none');
p.addParamValue('maxColor','none');
p.addParamValue('invert',0, @(x)x==0 || x==1);

if nargin==1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n);
elseif nargin>1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n, varargin{:});
else
    p.addParamValue('n',256, @(x)x>0 && mod(x,1)==0);
    p.parse();
end
config = p.Results;
n=config.n;
cP(:,1) = [30  60  150]./255; k(1)=1;  %cyan at index 1
cP(:,2) = [180 90  155]./255; k(3)=17; %purple at index 17
cP(:,3) = [230 85  65 ]./255; k(4)=32; %redish at index 32
cP(:,4) = [220 220 0  ]./255; k(5)=64; %yellow at index 64
for i=1:3
    f{i}   = linspace(0,1,(k(i+1)-k(i)+1))';  % linear space between these controlpoints
    ind{i} = linspace(k(i),k(i+1),(k(i+1)-k(i)+1))';
end
cmap = interp1((1:4),cP',linspace(1,4,64)); % for non-iso points, a normal interpolation gives better results
cmap = abs(interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,n)));
if config.invert
    cmap = flipud(cmap);
end
if ischar(config.minColor)
    if ~strcmp(config.minColor,'none')
        switch config.minColor
            case 'white'
                cmap(1,:) = [1 1 1];
            case 'black'
                cmap(1,:) = [0 0 0];
            case 'lightgray'
                cmap(1,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(1,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(1,:) = config.minColor;
end
if ischar(config.maxColor)
    if ~strcmp(config.maxColor,'none')
        switch config.maxColor
            case 'white'
                cmap(end,:) = [1 1 1];
            case 'black'
                cmap(end,:) = [0 0 0];
            case 'lightgray'
                cmap(end,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(end,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(end,:) = config.maxColor;
end
end
