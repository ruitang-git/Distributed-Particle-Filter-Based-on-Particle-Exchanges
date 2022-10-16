%% sample from GMM
function p = GMsample(a,Mu,Sigma,N)
% 
% Input:
% N: sample size
% Output:
% p: particles
%
alphabet = 1:length(a); % row vector
prob = a'; % row vector
k = randsrc(N,1,[alphabet; prob]); % particle generate from kth Gaussian
d = size(Mu,1); % particle dimension
p = zeros(d,N); % particles 
for n = 1:N
    mu = Mu(:,k(n));
    sigma = Sigma(:,d*(k(n)-1)+1:d*k(n));
    R = chol(sigma);
    p(:,n) = mu + R.'*randn(d,1);
end
    
    
    