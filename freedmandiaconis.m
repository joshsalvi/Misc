function nb = freedmandiaconis(x)
% Implementation of the Freedman-Diaconis Rule
%
%   nb = freedmandiaconis(x)
%
%   x : 1xN array of data
%   nb : number of bins according to this rule
%
%   jsalvi@rockefeller.edu

bw = 2*iqr(x)/length(x)^(1/3);
nb = ceil((max(x) - min(x))/bw);

if nb < 4
    nb = 4;
elseif isinf(nb) == 1
    nb = 4;
elseif isnan(nb) == 1
    nb = 4;
end
if nb > 1e3
    nb = 1e3;
end

end
