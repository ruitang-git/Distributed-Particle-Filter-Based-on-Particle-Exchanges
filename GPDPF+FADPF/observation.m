%% observation function
function [h,R] = observation(x,l,flag)
%  measurement
%  range + range rate(Doppler)
%  Input
%  l: [l1,l2]^T sensor location
%  x: [x1,x2,x3,x4,x5,x6]^T target: location+velocity+acceleration
%  
%  Output:
%  h(x) = [h_range(x), h_doppler(x)]^T
%

h_range = norm(x(1:2)-l);
h_doppler = (x(3)*(x(1)-l(1))+x(4)*(x(2)-l(2)))/h_range;

q = 2.2;
h_rss = 30 - 10*log10((norm(x(1:2)-l))^q);
h_angle = atan2(x(1)-l(1),x(2)-l(2));

% h = [h_rss;h_angle];
h = [h_range;h_doppler];
% noise
% R = [0.01 0;0 0.01];
% R = [1 0;0 0.1];
R = [1 0;0 1];
v = chol(R).' * randn(2,1);
if flag == 1
    h = h + v;
end