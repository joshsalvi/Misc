function [n,binsca]=histmulti6(y,varargin)
% For fast multidimensional histograms for Matlab 6. This is my hack based
% upon histmulti5 (made for matlab5). I have changed the interface to use
% cell arrays, done a bug fix or two.
%
% Credits follow as taken as cut and pasted from its comments
%
% Andrew Diamond. Use at your own risk.
%
%
% CREDITS FROM histmulti5:
%
% Rewritten for Matlab 5 by Hans Olsson, Hans.Olsson@dna.lth.se
% Original Matlab 4 routine written by Hans Olsson, Hans.Olsson@dna.lth.se.
% Inspired by hist2 written by Kirill K. Pankratov, kirill@plume.mit.edu
% Which in turn was based on routines written by Hans.Olsson@dna.lth.se
%
%
% Note that histmulti5 is available through the Mathworks file exchange.
% It has the option of using a mex DLL at its core to speed it up. I have
% removed that option from this version because I assume that it won't be
% forward compatible with Matlab 6x. I haven't checked for the mex souce
% code. Anyway, it seems to be pretty fast.
%
% Note - for fastest results use equally space bins within a dimension
%(explicit or implied)
%
% Hacked fron histmulti5 by Andrew Diamond with various bug fixes and
%removed
% optional core mex hist5hel because my hacks could break it. Also updated
% for matlab 6 with cell arrays.
%
% [n,binsca]=histmulti6(y)
% y - <# points> x <# dimensions> data array to be histogramed (i.e. <#
%dimensions> histogram)
% n - # dimensions dimensional histogram where each dimension will be
%binned using 10 equaly sized
% bins covering its min to maximum element (i.e. as contained in
%y(:,col)
% binsca - 1x n dimensional histogram cell array s.t. binsca{i} are the
%bins used for the ith dimension.
%
% [n,binsca]=histmulti6(y, bindata)
% bindata can be one of the following
% 1) vector where the bindata(i) is the number of bins that will be used
%to histogram the ith dimension. The bins
% will be created using equal spacing (fast) from the min to max value
%of each dimension. Note that this implies
% that the # of columns in y is the same as the length of bindata.
% 2) Backwards compatability - A nan paddded matrix where the ith column
%is gives the bin edges of the ith
% dimension s.t. the jth bin of dimension i is bounded by [bin(j,i),
%bin(j+1,i)) which implies that
% if the first k elements of bin(i,j) are not nan then there are k-1
%bins. Note that this imples
% that the # of columns in y is the same as the # of columns in
%bindata.
%
% [n,binsca]=histmulti6(y, <bin vector for dim 1>, ..., <bin vector for dim
%n>)
% <bin vector for dim i> - 1xNi - vector specifying the Ni bin edges for
%the ith dimension. This
% implies that y should have should have n columns.
%
% Sample Usage 1)
% G2 = randn(10000,2); G2(:,2) = G2(:,2) .* 3;
%
% Bin1 = -3:0.1:3; % bin edges for first dimension
% Bin2 = -6:0.1:6; % bin edges for second dimension
% [H2, binsca] = histmulti6(G2, Bin1, Bin2); % H2 is the hisogram binsca
%is the same as {Bin1, Bin2}.
%
% Sample Usage 2)
% [H2, binsca] = histmulti6(G2, [100,200]); % make histogram with 100
%equal bins for dim 1 and 200 bins for dim 2 from min to max points.
%
% Note that both these examples are fast because they explicitly or
%implictly specify equal bins within each dimension.

if (nargin<1)
  help histmulti
  return
end;
if (nargin<2)
  bins=10*ones(1,size(y,2));
end;
if(length(varargin) == 1)
    bins = varargin{1};
elseif(length(varargin) > 1)
    nbins = zeros(1,length(varargin));
    for ivarargin = 1:length(varargin)
        sizei = size(varargin{ivarargin});
        if(~any(sizei == 1) | ~any(sizei > 1))
            error('Bins must be vectors');
        end
        nbins(ivarargin) = length(varargin{ivarargin})-1;
    end
    bins = zeros(max(nbins)+1,length(nbins)) + nan;
    for col=1:size(bins,2)
        bins(1:nbins(col)+1,col) = varargin{col}(:);
    end
end

if (length(bins)==1)
  bins=bins*ones(1,size(y,2));
end;
if (size(bins,2)~=size(y,2))
  error('Wrong number of bins in histmul6');
end;

if (size(bins,1)==1) % Andy Diamond - rehash bins to explicit form than do
standard code.
    nbins = bins;
    bins = zeros(max(nbins)+1,length(nbins)) + nan;
    for col=1:size(bins,2)
        bins(1:nbins(col)+1,col) = linspace(min(y(:,col)), max(y(:,col)),
nbins(col)+1)';
    end
end
% Transform elements into bin numbers.
% if (size(bins,1)==1)
% % All bins have equal size.
% nbins=bins;
% mincol=min(y);
% maxcol=max(y);
% % Construct bin-width:
% binwidth=(maxcol-mincol)./bins;
% x=nan*zeros(max(bins),size(nbins,2));
% prodold=1;
% for col=1:size(bins,2)
% % Construct the bins:
% x(1:nbins(col),col)=(mincol(col)+(0:nbins(col)-1)*binwidth(col))';
% % Transform data into bins in this dimension.
% add=floor((y(:,col)-mincol(col))/binwidth(col));
% % Ensure that it is valid bin, we subtract in order to simplify
% % the next line.
% add=min(add,nbins(col)-1);
% % Make one bin-number:
% if (col==1) S=add; else S=S+add*prodold; end;
% prodold=prodold*nbins(col);
% end;
% else
  global binorder;
  x=bins;
  % nbins=sum(~isnan(bins));
  nbins=sum(~isnan(bins))-1; % Andy Diamond - bins as edges N bin edgles ==>
%N-1 bins.
  % Index used by tobinh:
  I=1:size(y,1);
  prodold=1;
  for col=1:size(bins,2);
    x0=x(1,col);
 % dx=diff(x(1:nbins(col)));
    dx=diff(x(1:nbins(col), col));

    Ind=find(y(:,col) < x0 | y(:,col) >= x(nbins(col)+1,col)); % Andy
%Diamond - up front elimination of out of bounds data

    if(max(dx) - min(dx) < 10 .* eps) % Andy Diamond - floating point round
%error can exceed eps. Using 10eps because its practical for a histogram.
 % if (abs(dx-dx(1))<= dx(1)*10 .* eps)
      % Same size in this bin, make it faster:
      % We must here handle data outside all bins:
      add=min(nbins(col), floor((y(:,col)-x0)/dx(1)));
 % Ind=find((add<0)|(add>=nbins(col)));
      add(Ind)=nan*ones(size(Ind));
      binorder=add;
    else
      y(Ind,col) = nan; % Andy Diamond - Bug fix - nan-out out of bounds
data so it doesn't possibly store it.
      % Initialize binorder
      binorder=nan*ones(size(y,1),1);
      % run tobinh which stores the bin in the global variable binorder
      tobinh(y(:,col),I,x(:,col),1:nbins(col)+1,nbins(col)+1);
    end;
    % Make one bin-number
    if (col==1), S=binorder; else,S=S+binorder*prodold; end;
    prodold=prodold*nbins(col);
  end;
  S=S(~isnan(S));
% end;

% Transform bin numbers into histograms.
S=S+1;
% Compute one-dimensional histogram for integers S in the range
% 1..prod(nbins)
% This is the core of the routine and for large datasets this should
% be optimized.

% We have a special mex-file, but if that routine is not available
% we use a one-liner with sparse.
global DO_HAVEHELP
global HAVE_WARNED
if (isempty(DO_HAVEHELP)|(DO_HAVEHELP<3))
  DO_HAVEHELP=exist('hist5hel');
  if ((DO_HAVEHELP<3)&(isempty(HAVE_WARNED))&(exist('hist5hel.c')))
    warning(['Please compile ',which('hist5hel.c'),', and: clear mex']);
    HAVE_WARNED=1;
  end;
end;
 n=full(sparse(S,1,1,prod(nbins),1));

 % commented out by Andy Diamond
% if (DO_HAVEHELP>=3)
% n=hist5hel(S,prod(nbins));
% else
% n=full(sparse(S,1,1,prod(nbins),1));
% end;


% Finally make the output into the right form, do nothing if one
% dimension because reshape requires two dimensions.
if (length(nbins)>1)
  n=reshape(n,nbins);
end;

if(nargout > 1)
    binsca = cell(1,length(nbins));
    for col=1:size(bins,2)
        binsca{col} = bins(1:nbins(col)+1,col);
    end
end



function tobinh(y,I,x,binnr,nbins)
% Helper to place elements into bins.
% Place the y in different bins.
% The result is in the global variable binorder.
%
% We use binary search for the bins.
% Number of operations should be length(x)*log(nbins)*nbins
global binorder;
if (isempty(y)|(nbins==0))
  % Do nothing.
elseif (nbins==1)
  if (binnr==1)
    sel=y>=x(1);
    I=I(sel);
  end;
  binorder(I)=(binnr-1)*ones(size(I));
else
  i=fix(nbins/2);
  middle=x(binnr(i+1))
  sel=y<middle;
  tobinh(y(sel),I(sel),x,binnr(1:i),i);
  sel=y>=middle;
  tobinh(y(sel),I(sel),x,binnr(i+1:nbins),nbins-i);
end;
