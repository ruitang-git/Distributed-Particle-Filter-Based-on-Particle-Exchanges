%% GM Fusion
function [a, mu, Sigma] = GMFusion(data, Eps, C, N, maxiter)
% At this file, we assume the GM statics of node k's neighboring sensors 
% are all available at node k
%
% Input: 
% data: cell data structure
% { [p(1)]  [mu(1)]  [sigma(1)] }
% { [p(nk)] [mu(nk)] [sigma(nk)] }
% { [p(Nk)] [mu(Nk)] [sigma(Nk)] }
% Eps: Metropolis weight
% N: total particles for fusion
%
% Output:
%
Eps(find(Eps==0))=[]; % delete the 0 element
Nk = size(data,1); % Nk: number of neighboring sensors at sensor k

% - constant eps - % 
% eps = 1/Nk;
% ---------------- %

d = size(data{1,2},1); % state dimension
P = []; % particle set
for nk = 1:Nk
    W = []; % weight set ; column vector
    % sampling particles
    a = data{nk,1}; Mu = data{nk,2}; Sigma = data{nk,3};
    Nnk = floor(N*Eps(nk)); % decide samping size at each node
    
    % -- constant eps ---- %
    % Nnk = floor(N*eps); % decide samping size at each node
    % -------------------- %
    
    p = GMsample(a,Mu,Sigma,Nnk); 
    % weight calculating
    for n = 1 : Nnk
        % denominator
        deno = GMprob(a, Mu, Sigma, p(:,n));
        % nominator
        no = 1;
        for nkk = 1:Nk
            a_nkk = data{nkk,1}; Mu_nkk = data{nkk,2}; Sigma_nkk = data{nkk,3};
            no = no * GMprob(a_nkk, Mu_nkk, Sigma_nkk, p(:,n))^Eps(nkk);
            
            % -- constant eps ---- %
            % no = no * GMprob(a_nkk, Mu_nkk, Sigma_nkk, p(:,n))^eps;
            % -------------------- %
        end
        w = no/deno;
%         if w == 0
%             w = 1e-3; 
%         end
        W = [W; w];
    end
    % normalize
%     if sum(W) == 0 || isnan(sum(W))
%         W = 1/length(W)*ones(1,length(W));
%     else
%         W = W/sum(W);
%     end
    W(W==0) = 1e-50; W(isinf(W)) = 1e50;
    W = W/sum(W);
    % resample
    p_update = zeros(size(p,1),size(p,2));
    for n = 1 : size(p_update,2)
        p_update(:,n) = p(:,find(rand <= cumsum(W),1));   % particles with higher weights get more offsprings
    end
    P = [P p_update];
end

% save particle_fuse P
% GMLearn
% GMModel = fitgmdist(P',C,'RegularizationValue',1e-5,'Options',statset('MaxIter',maxiter));
GMModel = fitgmdist(P',C,'RegularizationValue',1e-5,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
% GMModel = fitgmdist(P',C,'RegularizationValue',0.01,'CovarianceType','diagonal','Options',statset('MaxIter',maxiter));
a = GMModel.ComponentProportion';
mu = GMModel.mu';
Sigma = diag_sigma(GMModel.Sigma);
Sigma = reshape(Sigma,[d,d*C]);
% [a, mu, Sigma] = GM(P, W, C);