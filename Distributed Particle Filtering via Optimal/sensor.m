function [A,cd] = sensor(n,roi)
%function [A,cd] = sensor(n, roi)
%
% Description: 
% 'generating a sensor graph '
% 
% INPUTS:
% n - number of nodes(sensors) 
% roi - communication range
%
% OUTPUTS:
% cd - coordinates of each node; size: 2*n
% A - adjacent matrix of the graph
%
%  remaining question:
%  how to generate a much more uniformly distributed network ? 
%
rng(6);
% rng(5);

r = 180; % map range r*r
n_s = ceil(sqrt(n)); % divide the space into sub-space
step = r/n_s; % range of sub-space
p = randperm(n_s^2,n);
cd = zeros(2,n); % coordinates
for i = 1:n
    if mod(p(i),n_s) == 0
        cd(1,i) = (n_s-1) * step + step * rand;
    else
        cd(1,i) = (mod(p(i),n_s) - 1) * step + step * rand;
    end
    cd(2,i) = (ceil(p(i)/n_s) - 1) * step + step * rand;
end

A = zeros(n,n); % adjacency matrix
for i = 1:n
    for j = 1:n
        if norm(cd(:,i)-cd(:,j)) < roi
            A(i,j) = 1; A(j,i) = 1;
        end
    end
end
A = A - eye(n);

rng shuffle;