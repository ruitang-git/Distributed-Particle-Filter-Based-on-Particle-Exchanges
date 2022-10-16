%% Calculate the probability in GM Model
function p = GMprob(A, Mu, Sigma, x)
%
% Input:
% x: point
% others: statics
%
% Output:
% p: probabillity
%

C = length(A); % number of components
p = 0; % prob. initialization
d = size(Mu,1); % state dimension

for c = 1:C
    a = A(c); mu = Mu(:,c); sigma = Sigma(:,(c-1)*d+1:c*d);
    p = p + a*mvnpdf(x', mu', sigma);
end
