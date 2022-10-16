%% WEIGHTED Gaussian mixture learning

%
% ite = 100;
% while ite
%

function  [a, mu, Sigma] = GM(x, w, C)
%
% close all;
% clc;
% clear;
% load ('particles.mat')
% x = x_P_update;
% w = ones(1000,1)/1000;
% C = 2;
%

% Input:
% C: total number of components (adaptive)
% x: particles
% w: weight | column vector
%
% Output:
% a: C*1
% mu: d*C
% Sigma: d*(d*C)

% Questions here 
% 1. How to find the starting values ?
% 2. How to choose stopping criterion ?

% own code
% M = size(x,2);
d = size(x,1); % state dimension
% 
% mu = zeros(d, C); 
% a = ones(C, 1) / C;
% % initialize | how to initialize ?
% for c = 1:C
%   samples = x(:,randsample(M, 1));
%   mu(:, c) = samples;
%   Sigma(:,d*(c-1)+1:c*d) = 0.1 * rand(1,1) * cov(x');
%   Sigma(:,d*(c-1)+1:c*d) = Sigma(:,d*(c-1)+1:c*d) .* eye(d); % assume uncorrelated
% end
% 
% p = zeros(M,C); % M*C matirx
% 
% % loglikelihood = [];
% % epsilon = 1;
% ite = 100;
% while ite
%     for m = 1:M
%         for c = 1:C
%             R = Sigma(:,d*(c-1)+1:d*c);
%             p(m,c) = a(c)*mvnpdf(x(:,m)', mu(:,c)', R);
%         end
%         % normalize
%         p = p./sum(p,2);
%     end
%     for c = 1:C
%         a(c) = w' * p(:,c);
%         mu(:,c) = 1/a(c)*x*(w.*p(:,c));
%         Sigma(:,(c-1)*d+1:c*d) = 1/a(c)*((x-mu(:,c)).*(sqrt(p(:,c).*w))')...
%             *((x-mu(:,c)).*(sqrt(p(:,c).*w))')';
%         Sigma(:,d*(c-1)+1:c*d) = Sigma(:,d*(c-1)+1:c*d).*eye(d) + 0.01*eye(d);
%     end
%     % normalize
%     a = a/sum(a);
%     
% %     % compute loglikelihood
% %     Laccum = 0;
% %     for m = 1:M
% %         laccum = 0;
% %         for c = 1:C
% %             laccum = laccum + a(c)* mvnpdf(x(:,m), mu(:, c), Sigma(:,d*(c-1)+1:c*d));
% %         end
% %         Laccum = Laccum + log(laccum);
% %     end
% % 
% %     loglikelihood = [loglikelihood; Laccum];
% %     ite = ite + 1;
% %     
% %     % stopping condition
% %     % case when log-likelihood stop improving by more than eps
% %     if numel(loglikelihood) > 1
% %         if abs(loglikelihood(end)-loglikelihood(end-1)) <= epsilon
% %             break;
% %         end
% %     end
%     
%     ite = ite - 1;
% end


% build-in function
GMModel = fitgmdist(x',C,'RegularizationValue',0.1);
a = GMModel.ComponentProportion';
mu = GMModel.mu';
Sigma = reshape(GMModel.Sigma,[d,d*C]);

% plot
% figure;
% gscatter(x(1,:), x(2,:));
% hold on;
% x1 = min(x(1,:))-2:0.1:max(x(1,:))+2;
% x2 = min(x(2,:)):0.1:max(x(2,:));
% [X1,X2] = meshgrid(x1,x2);
% for c = 1:C
%     % plot the mixture models
%     F = mvnpdf([X1(:) X2(:)],mu(1:2, c)',Sigma(1:2,d*(c-1)+1:d*(c-1)+2));
%     F = reshape(F,length(x2),length(x1));
%     contour(x1,x2,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]); 
%     plot(mu(1,c),mu(2,c),'kx','LineWidth',2,'MarkerSize',10);
% end
        

%
% ite = ite-1;
% end
   
