function setfiguredefaults(N)
% Set the figure defaults, with axes order for N axes
%
% setfiguredefaults(N,dockyn)
%
% where N is the number of axes for your colormap
% If N is not included, the default line color will be black.
%
%
% Joshua Salvi
% jsalvi@rockefeller.edu
%
%
% NOTE: Requires the 'ametrine' colormap
%

set(0,'defaultFigureColormap',jet)
if exist('N')==1
    set(0,'defaultAxesColorOrder',ametrine(N))
else
    set(0,'defaultAxesColorOrder',[0 0 0])
end
set(0,'DefaultTextFontName','Times')
set(0,'DefaultTextFontUnits','Points')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontName','Times')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesGridLineStyle',':')
set(0,'defaultFigureColor','White')
set(0,'DefaultFigurePaperType','a4letter')
set(0,'defaultFigurePosition',[1200 100 700 700])
set(0,'defaultAxesColor',[1 1 1])
set(0,'defaultLineMarker','None')
