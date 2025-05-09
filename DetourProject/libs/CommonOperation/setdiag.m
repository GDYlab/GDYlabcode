function X = setdiag(X,diav)
% make the diagonal of matrix X all equals to diav
% % s = diag(X);
% % news = ones(size(s));
% % news = news*diav;
% % X = X - diag(s) + diag(news);
% old code doesn't work for nan value
[m,n] = size(X);
for i = 1:min([m,n])
    X(i,i) = diav;
end

end