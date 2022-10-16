function [alpha,l, sigma_n, bestEstimate, y_hat] = GP_optfree(X, y, x, xstar)
% Dependencies: calcGP.m k_GP.m, GPtutorialFcn.m, hypSample.m
% clear;clc;close all;
if nargin < 4
    xstar = [];
    y_hat = [];
    bestEstimate = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%% 1. Initialization %%%%%%%%%%%%%%%%%%%%%%%%%

l = 1; sigma_n = 0.001;
%%%%%%%%%%%%%%%%%%%%% 3. Plot the regression %%%%%%%%%%%%%%%%%%%%%%%%%

% a) Plot Figure 2 in GPtutorial.pdf
% x_grid = min(X(:,1)):0.02:max(X(:,1)); y_grid = min(X(:,2)):0.02:max(X(:,2));
% [X_grid,Y_grid] = meshgrid(x_grid,y_grid);
% cord_star = [];
% for i = 1:size(X_grid,1)
%     for j = 1:size(X_grid,2)
%         cord_star = [cord_star; [X_grid(i,j) Y_grid(i,j)]];
%     end
% end
% [cord_star, bestEstimate_cord, bounds, K, kstar, alpha, y_hat] = calcGP (@k_GP, X, y, l, sigma_n, cord_star, 1);
if nargin == 4
    [xstar, bestEstimate, bounds, K, kstar, alpha, y_hat] = calcGP (@k_GP, X, y, l, sigma_n, xstar, 1);
    debug = 1;
else
    [~, ~, ~, ~, ~, alpha, ~] = calcGP (@k_GP, X, y, l, sigma_n, [], 1);
end
% figure
% figure; plot3(X(:,1),X(:,2),y,'.');hold on;
% plot3(x(1),x(2),0,'kx','LineWidth',2,'MarkerSize',10); hold on;
% if nargin == 4
%     bestEstimate = bestEstimate/sum(bestEstimate);
%     plot3(xstar(:,1),xstar(:,2),bestEstimate,'.');hold off;
% end
% hold off;
% bestEstimate_cord = bestEstimate_cord/sum(bestEstimate_cord);
% surf(X_grid,Y_grid,(reshape(bestEstimate_cord,[size(X_grid,2), size(X_grid,1)]))'); hold off;

end



function [covar, d_l, d_sigman] = k_GP (x1, x2, theta)
l = theta(1);
sigman = theta(2);
%  Covariance and gradients
covar = exp(-(norm(x1-x2))^2/(2*l^2));
d_l = covar * (l^-3) * (norm(x1-x2))^2; % Differentiate (2.16) from Rasmussen and Williams (2006)
s = sum(double(x1==x2));
switch s
    case 4
        covar = covar + sigman^2;
        d_sigman = 2*sigman;
    otherwise
        d_sigman = 0;
end
end

function [xstar, bestEstimate, bounds, K, V, alpha, y_hat] = calcGP (k, X, y, l, sigma_n, xstar, graphing)

% a) Initializations

% meany = mean(y); y = y - meany;
meany = 0;

n = length(y); lx = size(xstar,1); 
K = zeros(n); kstar = zeros(n,1);


K = exp(-(dist(X.').^2)/(2*l^2)) + sigma_n^2 * eye(n);

fstarbar = zeros(lx,1); V = zeros(lx,1);

% b) One-off calculations
diags = max(1e3*eps, 0); % beef up the diagonal if sigma_n = 0
L = chol (K + diags*eye(n),'lower'); 
alpha = L'\(L\y);
alpha = max(alpha,0);


% alpha = max(L'\(L\y),0); % but alpha can be very small ?
logpyX = -y'*alpha/2 - sum(log(diag(K)))/2 - n*log(2*pi)/2; % Log marginal likelihood

% c) xstar loop
for q = 1:lx
    for i = 1:n
        kstar(i) = k (X(i,:),xstar(q,:),[l sigma_n]);
    end
    % Mean of prediction
    fstarbar(q) = kstar' * alpha;
    % Variance of prediction
    v = L\kstar;
    V(q) = k (xstar(q,:),xstar(q,:),[l sigma_n]) - v'*v;
end
bounds = [fstarbar+1.96*sqrt(V) fstarbar-1.96*sqrt(V)]+meany;
bestEstimate = fstarbar+meany;
reg = 3*sqrt(V);

% 
% bestEstimate = bestEstimate/sum(bestEstimate);
% reg = reg/sum(reg);

y_hat = bestEstimate+reg;

debug = 1;
end
