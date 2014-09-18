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

if nb < 4 || isinf(nb) == 1 || isnan(nb) == 1
    nb = 4;
end

end
