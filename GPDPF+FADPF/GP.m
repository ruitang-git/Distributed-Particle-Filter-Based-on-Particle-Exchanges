% GP 
function [alpha,l, sigma_n, y_hat] = GP(X, y, x, xstar)
% Dependencies: calcGP.m k_GP.m, GPtutorialFcn.m, hypSample.m
% clear;clc;close all;
if nargin < 4
    xstar = [];
    y_hat = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%% 1. Initialization %%%%%%%%%%%%%%%%%%%%%%%%%

% a) Parameters
% sigma_f = 1; % signal variance
% % b) Inputs (example)
% X = [-1.5 -1 -.75 -.4 -.25 0]'; % N * d 
% y = .55*[-3 -2 -.6 .4 1 1.6]'; % N * 1 
% x_grid = -2:0.5:2; y_grid = x_grid;
% [X_grid,Y_grid] = meshgrid(x_grid,y_grid);
% cord = [];
% for i = 1:length(x_grid)
%     for j = 1:length(y_grid)
%         cord = [cord; [X_grid(i,j) Y_grid(i,j)]];
%     end
% end
% z =  sin(cord(:,2)).^2 + cos(cord(:,1));
% X = cord; y = z + 0.1*randn(length(z),1);
% X = x_P(1:2,:)'; y  = P_w';
% load X; load y; 
% 
% figure; plot3(X(:,1),X(:,2),y,'.')
%%%%%%%%%%%%%%%%%%% 2. Fit the ML parameters %%%%%%%%%%%%%%%%%%%%%%%%%

% a) Basic initializations
% options = optimset('GradObj','on');
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
Nstarts = 1;
% b) Choose several varied starting points
l_samp = hypSample ([1e-2 1e-1], Nstarts); % for l
sn_samp = hypSample ([1e-2 1e-1], Nstarts); % for sigma_n
inits = log([l_samp sn_samp]);
lb = [log(1e-2) log(1e-2)]; % lower bound
ub = [log(1e-1) log(1e-1)]; % upper bound

l_samp = hypSample ([1e-1 1], Nstarts); % for l
sn_samp = hypSample ([1e-1 1], Nstarts); % for sigma_n
inits = log([l_samp sn_samp]);
lb = [log(1e-1) log(1e-1)]; % lower bound
ub = [log(1) log(1)]; % upper bound

% lb = []; % lower bound
% ub = []; % upper bound
% c) Find candidates
paramVec = [];
for randomStart = 1:Nstarts
%     tic;
    [params, fval] = fmincon(@(params) GPtutorialFcn(params,X,y), inits(randomStart,:),[],[],[],[],lb,ub,[], options);
%     toc;
    paramVec = [paramVec; fval inits(randomStart,:) params];
end
paramVec(:,2:end) = exp(paramVec(:,2:end));
paramVec(:,1) = paramVec(:,1)/max(abs(paramVec(:,1)));
% d) Select best candidate
paramVec = sortrows(paramVec);
params = paramVec(1,size(inits,2)+2:end);
l = params(1); sigma_n = params(2);
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

% axis tight, v = axis; axis([-1.7 v(2:4)])
% v = axis; line ([v(1) v(1)], [v(3) v(4)], 'Color', 'k')
% xlabel ('x'), ylabel('y')
% % b) Calculate y*
% xindex = 452; % determined by inspection
% newX = xstar(xindex); newBounds = bounds(xindex,:);
% newEst = bestEstimate(xindex); newSigma = (newBounds(1)-newEst)/1.96;
% disp([newEst newSigma.^2]) % Answers reported in Step 2 of page 3 of GPtutorial.pdf
% % c) Update Figure 2 with y*
% errorbar (newX,newEst,newSigma,newSigma,'Color','green')
% plot (newX,newEst,'b.','MarkerSize',15)
% errorbar (X,y,sigma_n*ones(size(y)),sigma_n*ones(size(y)),'Color','red', 'Linestyle','none')
% plot (xstar,bestEstimate,'k'), plot (X,y,'k.', 'MarkerSize',15)
% % d) Plot Figure 1
% figure, hold on
% errorbar (X,y,sigma_n*ones(size(y)),sigma_n*ones(size(y)),'Color','red','Linestyle','none') % plot error bars
% errorbar (newX,newEst,newSigma,newSigma,'Color','green')
% text(newX-.01,newEst-newSigma-0.4,'?')
% plot (X,y,'k.', 'MarkerSize',15), plot (newX,newEst,'b.', 'MarkerSize',15)
% xlabel ('x'), ylabel('y'), axis(v)
end

function x = hypSample (bounds, N)

xmin = bounds(1); xmax = bounds(2);
F = rand(N,1);
x = xmin.^(1-F) .* xmax.^F;
end

function [fcn, grd] = GPtutorialFcn (params,X,y)

% a) Initializations
eparams = exp(params); 
% Variables:
l1 = eparams(1);
sigma_n = eparams(2); % avoid sigma_n being to small which may lead to K matrix non positive definite

n = length(y);
y = y - mean(y); % mean equals zero;

% b) Variance and its derivative
K = zeros(n); dKdl = zeros(n); dKds = zeros(n);

% tic;
K = exp(-(dist(X.').^2)/(2*l1^2));
dKdl = K .* (l1^-3) .* (dist(X.').^2); % Differentiate (2.16) from Rasmussen and Williams (2006)
K = K + sigma_n^2 * eye(n);
dKds = dKds + 2*sigma_n * eye(n);
% toc;


% c) Calculations
L = chol (K,'lower');
alpha = L'\(L\y);
invK = L'\(L\eye(size(L,1)));
% invK = pinv(K);
% alpha = invK*y; % alpha from page 114, not page 19, of Rasmussen and Williams (2006)
% alpha = max(alpha,0);
% d) Log marginal likelihood and its gradient 
logpyX = -y'*alpha/2 - sum(log(diag(L))) - n*log(2*pi)/2;
dlogp_dl = l1*trace((alpha*alpha' - invK)*dKdl)/2;
dlogp_ds = sigma_n*trace((alpha*alpha' - invK)*dKds)/2;
% NB: Since l = exp(p1), dlogp/dp1 = dlogp/dl dl/dp1 = dlogp/dl exp(p1) = dlogp/dl l, where p1 = params(1)

% e) Format output
fcn = -logpyX; grd = [-dlogp_dl -dlogp_ds];
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
meany = mean(y); y = y - meany;
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
% 
y_hat = bestEstimate+3*sqrt(V);  % f = mu + 3*sigma; required in resampling stage
end
