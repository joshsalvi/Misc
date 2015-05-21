function r = drchrnd(a,n)
%
% Sample the Dirichlet distribution
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%


p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

end
