%% Recovery step
function [a, Mu, Sigma] = GMRecovery(x, data, K, C, M, maxiter, reg, plt)

% K = N_sensor;
% M = N_recovery;
% plt = 0;
C = 1;
%
% Input:
% x: true state
% data: cell structure with 2 GMM
% { [p1] [mu1] [sigma1] [p2] [mu2] [sigma2] } 
% { prior prob. f(xn|y1:n-1)+local posterior(averaging) } 
% K: total sensors
% C: Gaussian components
% M: total sampling particles for recovery
%

% K = 1;

A1 = data{1,4}; Mu1 = data{1,5}; Sigma1 = data{1,6}; % local posterior
A2 = data{1,1}; Mu2 = data{1,2}; Sigma2 = data{1,3}; % prior prob.
P1 = GMsample(A1, Mu1, Sigma1, M/2); % draw from posterior
P2 = GMsample(A2, Mu2, Sigma2, M/2); % draw from prior
P = [P1 P2]; % particles

% ---- verify the fusion(testing, delete in final version) @1_st sensor ---- %
% if plt == 1
%     figure;
%     plot(x(1),x(2),'kx','LineWidth',2,'MarkerSize',10);
%     hold on;
%     scatter(P1(1,:),P1(2,:),'r.');
%     hold on;
%     scatter(P2(1,:),P2(2,:),'b.');
%     hold off;
%     legend('true state', 'GA-posterior', 'prior');
%     title('PDF distribution');
% end
% ---- verify the fusion(testing, delete in final version) @1_st sensor ---- %

W1 = []; % weights
W2 = [];
for p1 = 1:M/2
    w = GMprob(A1, Mu1, Sigma1, P1(:,p1))^(K-1)/...
        GMprob(A2, Mu2, Sigma2, P1(:,p1))^(K-1);
%     if w == 0
%        w = 1e-3; 
%     end
    W1 = [W1;w];
end
for p2 = 1:M/2
    w = GMprob(A1, Mu1, Sigma1, P2(:,p2))^K/...
        GMprob(A2, Mu2, Sigma2, P2(:,p2))^K;
%     if w == 0
%        w = 1e-3; 
%     end
    W2 = [W2;w];
end
% normalize
% fprintf('%e\n',max(W1))
% fprintf('%e\n',sum(W1))
% fprintf('%e\n',max(W2))
% fprintf('%e\n',sum(W2))

W1(W1==0) = 1e-50; W1(isinf(W1)) = 1e50;
W1 = W1/sum(W1);
W2(W2==0) = 1e-50; W2(isinf(W2)) = 1e50;
W2 = W2/sum(W2);
% if (sum(W1) == 0 || isnan(sum(W1)))
%     W1 = 2/M*ones(M/2,1);
% else
%     W1 = W1/sum(W1);
% end
% if (sum(W2) == 0 || isnan(sum(W2)))
%     W2 = 2/M*ones(M/2,1);
% else
%     W2 = W2/sum(W2);
% end
        
% W1= W1/sum(W1); 
% W2= W2/sum(W2);
% if plt == 1
%     figure;
%     histogram(W1,50,'BinLimits',[0.0001,0.01]);
%     title('particle weights after normalization(sampling from GA-posterior)');
%     figure;
%     histogram(W2,50,'BinLimits',[0.0001,0.01]);
%     title('particle weights after normalization(sampling from prior)');
% end
% resample
P_update = zeros(size(P,1),size(P,2));

for m = 1 : M/2
    P_update(:,m) = P1(:,find(rand <= cumsum(W1),1));   % particles with higher weights get more offsprings
end
for m = M/2+1 : M
    P_update(:,m) = P2(:,find(rand <= cumsum(W2),1));   % particles with higher weights get more offsprings
end

% ---- verify the fusion(testing, delete in final version) @1_st sensor ---- %
if plt == 1
    figure;
    plot(x(1),x(2),'kx','LineWidth',2,'MarkerSize',10);
    hold on;
    scatter(P_update(1,1:M/2),P_update(2,1:M/2),'r.');
    hold on;
    scatter(P_update(1,M/2+1:M),P_update(2,M/2+1:M),'r.');
    hold off;
    legend('true state');
    title('global posterior after recovery(particles after resampling)');
end
% ---- verify the fusion(testing, delete in final version) @1_st sensor ---- %

% GMLearn

% (to delete after) %
% GMModel1 = fitgmdist(P_update',5,'RegularizationValue',reg,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
% a1 = GMModel1.ComponentProportion';
% mu1 = GMModel1.mu';
% Sigma1 = diag_sigma(GMModel1.Sigma);
% Sigma1 = reshape(Sigma1,[size(P,1),size(P,1)*5]);
% (to delete after) %
    
% GMModel = fitgmdist(P_update',C,'RegularizationValue',1e-5,'Options',statset('MaxIter',maxiter));
GMModel = fitgmdist(P_update',C,'RegularizationValue',reg,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
% GMModel = fitgmdist(P_update',C,'RegularizationValue',0.01,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
a = GMModel.ComponentProportion';
Mu = GMModel.mu';
Sigma = diag_sigma(GMModel.Sigma);
Sigma = reshape(Sigma,[size(P,1),size(P,1)*C]);

% (to delete after) %
% plotparticles(a,Mu,Sigma,P_update,2000,x);
% (to delete after) %

debug = 1;
% GMLearn
% [a, Mu, Sigma] = GM(P, W, C);
