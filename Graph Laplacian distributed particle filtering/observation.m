%% observation function
function [h,R] = observation(x,l,sigma_v,sigma_w,flag)
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
h = [h_range;h_doppler];
% noise
mu = zeros(2,1);
R = [sigma_v^2 0;0 sigma_w^2];
v = mu + chol(R).' * randn(2,1);
if flag == 1
    h = h + v;
end