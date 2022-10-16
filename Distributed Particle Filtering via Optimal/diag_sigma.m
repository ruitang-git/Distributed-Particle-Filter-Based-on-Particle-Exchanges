function sigma_r = diag_sigma(sigma)
% description:
%
% reshape the diagonal covariance matrix(vector form) generating from function 
% "fitgmdist" to matrix form
%
sigma_r = zeros(size(sigma,2),size(sigma,2),size(sigma,3));
for i = 1:size(sigma,3)
    sigma_r(:,:,i) = diag(sigma(:,:,i));
end