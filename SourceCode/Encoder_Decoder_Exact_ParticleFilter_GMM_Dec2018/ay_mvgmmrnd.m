function y = ay_mvgmmrnd(mu,sigma,p,n)
% MVGMMRND Random vectors from a mixture of multivariate normals.
% MU is an M-by-D matrix of means for the M component normals
% SIGMA is a D-by-D-by-M array of covariance matrices for the
% M component normals.
% P is an M-by-1 vector of component mixing probabilities.
% N is the desired number of random vectors.

[M,d] = size(mu);

% randomly pick from the components
[dum,compon] = histc(rand(n,1), [0; cumsum(p(:))./sum(p)]);

% generate random vectors from the selected components with a
% "stacked" matrix multiply
for i = 1:M
     Rt(i,:,:) = chol(sigma(:,:,i)); % holds the transposed cholesky factors
end
Z = repmat(randn(n,d), [1,1,d]);
y = squeeze(sum(Z.*Rt(compon,:,:),2)) + mu(compon,:);

% another way to generate the random vectors
% y = zeros(n,d);
% for i = 1:M
% mbrs = find(compon == i);
% ni = length(mbrs);
% y(mbrs,:) = randn(ni,d) * chol(sigma(:,:,i)) + repmat(mu(i,:),ni,1);
% end