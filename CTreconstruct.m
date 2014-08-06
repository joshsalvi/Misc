function [faces, verts] = CTreconstruct(S,thr,fact,crop,FN)
%CTRECONSTRUCT - reconstructs a 3D representation of a .dcm CT image series
%   CTRECONSTRUCT('dir',thr) opens a series of DICOM images from single
%   directory 'dir'. It creates a 3D representation of the scanned
%   structure by thresholding the source images at thr (typically XX-XX).
%
%   [F, V] = CTreconstruct(DCM,[],2,{'x' {1:12}},'filename');
%   [F, V] = CTreconstruct(DCM,12000,[],{'xz' {1:12 6:8}});
%

% Defaults
if nargin < 2 || isempty(thr), thr = S(1).contourVal; end
if nargin < 3 || isempty(fact), fact = 1; end
if nargin < 4 || isempty(crop), crop = 'none'; end
if nargin < 5 || isempty(FN), FN = 'test'; end

% Construct 3D matrix
S1 = S(1);
rows = double(S1.info.Rows);
cols = double(S1.info.Columns);
stacks = numel(S);
MM = reshape([S.imageAdj],rows,cols,stacks);
MM = double(MM);

% Scale factors
Xscale = S1.info.PixelSpacing(1)*1e3;
Yscale = S1.info.PixelSpacing(2)*1e3;
Zscale = S1.info.SliceThickness*1e3;

% Base specifications
% Indices
Xind = 0:rows-1;
Yind = 0:cols-1;
Zind = [S.iFrame]-1;

% Crop and seal matrix
cropdim = crop{1};
cropind = crop{2};
if ~iscell(crop{2})
    cropind = {cropind};
end
cropX = strfind(cropdim,'x');
cropY = strfind(cropdim,'y');
cropZ = strfind(cropdim,'z');
if ~isempty(cropX)
    % Rretrieve 
    cropindX = cropind{cropX};
    % Index change x
    Xind = cropindX-1;
end
if ~isempty(cropY)
    % Rretrieve 
    cropindY = cropind{cropY};
    % Index change y
    Yind = cropindY-1;
end
if ~isempty(cropZ)
    % Rretrieve 
    cropindZ = cropind{cropZ};
    % Index change z
    Zind = cropindZ-1;
end

% Adjust matrix
TT = MM(Yind+1,Xind+1,Zind+1);

% Add blanks on all sides
sz = size(TT);
MM = zeros(sz(1)+2,sz(2)+2,sz(3)+2);
MM(2:end-1,2:end-1,2:end-1) = TT;

% Set axes
Xax = [Xind(1)-1 Xind Xind(end)+1]*Xscale;
Yax = [Yind(1)-1 Yind Yind(end)+1]*Yscale;
Zax = [Zind(1)-1 Zind Zind(end)+1]*Zscale;

% Create grid
[XX, YY, ZZ] = meshgrid(Xax,Yax,Zax);

% Resample matrix if scale factor ~= 1
if fact ~= 1
    rowsnew = numel(Xax);
    colsnew = numel(Yax);
    stacksnew = numel(Zax);
    Xnew = linspace(min(Xax),max(Xax),round(rowsnew/fact));
    Ynew = linspace(min(Yax),max(Yax),round(colsnew/fact));
    Znew = linspace(min(Zax),max(Zax),round(stacksnew/fact));
    [XXnew, YYnew, ZZnew] = meshgrid(Xnew,Ynew,Znew);
    MM = interp3(XX,YY,ZZ,MM,XXnew,YYnew,ZZnew);
else
    XXnew = XX;
    YYnew = YY;
    ZZnew = ZZ;
end

% Find largest structure and remove others
MM(MM<thr) = 0; % threshold matrix to binary
MM(MM>=thr) = 1; % threshold matrix to binary
CC = bwconncomp(MM,18); % compute connectivity
numPixels = cellfun(@numel,CC.PixelIdxList);
[junk, indx] = max(numPixels); %#ok<ASGLU>
NN = zeros(size(MM));
NN(CC.PixelIdxList{indx}) = 1;

% Compute surfaces
[faces, verts] = isosurface(XXnew,YYnew,ZZnew,NN,.9);

%vertface2obj(verts,faces,'test.obj')
stlwrite([FN '.stl'],faces,verts)




