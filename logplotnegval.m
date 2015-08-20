function logplotnegval(varargin)
% 
%  
% This function permits the plotting of negative values on a logarithmic
% plot.
%
%  logplotnegval(x,y,choice)
%
%       x,y:  x and y vectors
%       choice: Use this to plot as (1) loglog, (2) semilogX, (3) semilogY
%
% Programmer:
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%

if numel(varargin) < 3
    input1 = 1;
else
    input1 = varargin{3};
end

x = varargin{1};
y = varargin{2};



if input1 == 1
    X = sign(x).*log10(abs(x));
    Y = sign(y).*log10(abs(y));
    plot(X,Y,'.')
    yl = get(gca,'ytick');
    xl = get(gca,'xtick');
    %set(gca,'xticklabel',sign(xl).*10.^abs(xl))
    %set(gca,'yticklabel',sign(yl).*10.^abs(yl))
elseif input1 == 2
    X = sign(x).*log10(abs(x));    
    plot(X,y,'.')
    xl = get(gca,'xtick');
    %set(gca,'xticklabel',sign(xl).*10.^abs(xl))
elseif input1 == 3
    Y = sign(y).*log10(abs(y));
    plot(x,Y,'.')
    yl = get(gca,'ytick');
    set(gca,'yticklabel',sign(yl).*10.^abs(yl))
else
    disp('Incorrect input selected.');
    return;
end

end
