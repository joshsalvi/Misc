function [Y] = anovize(X)
Y = zeros(size(X,1),(size(X,2)+1)*(size(X,2)/2));
Y(:,1:size(X,2)) = X.^2;
offset = size(X,2);
for ii = 2:size(X,2)
Y(:,offset + (1:size(X,2)+1-ii)) = 2*X(:,1:end-ii+1).*X(:,ii:end);
offset = offset + size(X,2) + 1 - ii;
end
end
