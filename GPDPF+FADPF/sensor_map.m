function [A,coords] = sensor_map(N, len, width, plt)
%function [coords,A] = sensor_map(N)
%
% Description: 
% 'generating a sensor graph(grid)'
% 
% INPUTS:
% N - number of nodes(sensors)
% len - length of the map
% width - width of the map
%
% OUTPUTS:
% coords - coordinates of each node; size: 2*N
% A - adjacent matrix of the graph
%

%
plt = 0;

% calculate the coordinates
coords = [N,2];
for i = 1:sqrt(N)
    for j = 1:sqrt(N)
        coords((i-1)*sqrt(N)+j,1) = (i-1)*len/(sqrt(N)-1);
        coords((i-1)*sqrt(N)+j,2) = (j-1)*width/(sqrt(N)-1);
    end
end

% calcute the adjacent matrix
A = zeros(N);
dist_coords = dist(coords');
for i = 1:N
    for j = 1:N
        if dist_coords(i,j) <= len/(sqrt(N)-1)
            A(i,j) = 1;
        end
    end
end
A = A - eye(size(A,1));
coords = coords';

% plot
if plt == 1
    figure;
    [Xt,Yt] = gplot(A,coords','-w');
    plot(Xt,Yt,'-o','Color','k','MarkerEdgeColor','r',...
        'MarkerFaceColor','r','LineWidth',0.2);
    title(['Sensor Grid (N= ',num2str(N),')']);
end

end