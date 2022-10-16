function [f_hat,V] = graphLap(z, N, f, m, k, plt)
%function [f_hat,V] = graphLap(z, N, f, m, k, plt)
%
% Description:
% 'performs the graph laplacian approximation algorithm.'
%
% INPUTS:
%  z - locations of nodes (or to say, particles) 
%  N - number of nodes
%  f - signals at each node
%  m - m smallest Laplacian eigenvalues needed to be reserved
%  k - k-nearst neighbors
%  plt - decide the plotting process (1 - yes; 0 - no)
%
% OUTPUTS:
%  f_hat - approxiamtion of frequency values
%  V - graph eigen matrix
%

% plt = 0;



z = z.';
% z_va = z(:,1:4);
Idx = knnsearch(z,z,'k',k);
A = zeros(N); % adjacent matrix
for i = 1:N
    for j = 1:k
        A(i,Idx(i,j)) = 1;
        A(Idx(i,j),i) = 1;
    end
end

% plot particle cloud
particle_w = ceil(exp(sum(f,2)));
% if plt == 1
%     figure;
%     jet_color = colormap(jet(max(particle_w)));
%     gplot(A,z(:,1:2),'-'); % plot edges
% %     title('Particle Graph');
%     hold on;
%     for i = 1:N
%         selected_color = jet_color(particle_w(i),:);
%         plot(z(i,1), z(i,2), 'o','color','k','MarkerEdgeColor',selected_color,...
%          'MarkerFaceColor',selected_color,'MarkerSize',3,'LineWidth',0.01);  % plot weights by varying colors
%         hold on;
%     end
%     c.Label.String = 'weight';
%     hold off;
% end

% graph laplacian decomposition
W = dist(z'); % euclidean distance matrix
A = A.* W; % adjacent matrix in weighted graph
% A = A - eye(size(A,1));
D = diag(sum(A)); % Degree matrix
L = D-A; % laplacian matrix
[V, lambda] = eig(L);
[~,index] = sort(diag(lambda),'ascend');
V = V(:,index);


% transfer signals to frequency domain
f_hat = V'*f;

% plot
if plt == 1
    figure;
    plot(sum(abs(f_hat),2));
%     title('Frequency components of the original signals/weights');
    xlabel('j(j^{th} component)');
    ylabel('$\sum_{i=1}^N|f(j)|$');
end

% approxiamtion
f_hat = f_hat(1:m,:);
V = V(:,1:m);
% f_hat(m/2+1:end-m/2+1,:) = 0;
end
