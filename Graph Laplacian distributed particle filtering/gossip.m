function ave_x = gossip(E, x, iter, plt)
%function ave_x = gossip(E, x, plt)
%
% Description:
% 'performs the Randomized Gossip Algorithm. '
%
% INPUTS:
%  E - adjacent matrix
%  x - signals at each node
%  plt - decide the plotting process (1 - yes; 0 - no)
%  iter - gossip iterations
%
% OUTPUT:
%  ave_x - average of all signals (achieved by averaging process)
%

N = size(E,1); % number of nodes
% if N == 4
%     load gossip_4
% elseif N == 20
%     load gossip_20
% elseif N == 25
%     load gossip_25
% else
% % optimization part
% cvx_begin quiet
%     variable P(N,N) nonnegative symmetric;
%     P_v = P.*E;
%     W = zeros(N);
%     for i = 1:N
%         for j = 1:N
%             e_i = zeros(N,1);
%             e_j = zeros(N,1);
%             e_i(i,1) = 1;
%             e_j(j,1) = 1;
%             w = eye(N)-0.5*(e_i-e_j)*(e_i-e_j)';
%             W = W + P_v(i,j)*w;
%         end
%     end
%     W = W/N;
%     minimize lambda_sum_largest(W,2)
%     subject to
%         sum(P_v,2) == ones(N,1);
% cvx_end
% P = full(P)./sum(full(P),2);
% end


% averaging process
x_ave = mean(x,2);
index = 1:N;
k = 1;
while iter > 0
    i = ceil(N*rand(1));
    list = find(E(i,:));
    j = list(randi(length(list)));
%     j = randsrc(1,1,[index;P(i,:)]);
    middle = (x(:,i)+x(:,j))/2;
    x(:,i) = middle;
    x(:,j) = middle;
    error = x - x_ave*ones(1,N);
    e_rg_ave(k) = sum(error.^2,'all');
    iter =  iter - 1;
    k = k+1;
end

% max consensus
x_max = max(x,[],2);
index = 1:N;
k = 1;
iter_max = 4*N;
while iter_max > 0
    i = ceil(N*rand(1));
    list = find(E(i,:));
    j = list(randi(length(list)));
%     j = randsrc(1,1,[index;P(i,:)]);
    middle = max(x(:,i),x(:,j));
    x(:,i) = middle;
    x(:,j) = middle;
    error = x - x_max*ones(1,N);
    e_rg_max(k) = sum(error.^2,'all');
    iter_max =  iter_max - 1;
    k = k+1;
end

ave_x =  x(:,1); % achieved by max consensus

% plot
if plt == 1
    figure;
    plot(log(e_rg_ave));
    xlabel("transmissions");
    ylabel("log(error)");
    title("Convergence curve of randomized Gossip(averaging)");
    figure;
    plot(log(e_rg_max));
    xlabel("transmissions");
    ylabel("log(error)");
    title("Convergence curve of randomized Gossip(maximizing)");
end
