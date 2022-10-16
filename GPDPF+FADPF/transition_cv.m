%% state transition model(CV)
function y = transition_cv(x,sigma_u)
% g(x) = D.x m + v
%  Input
%  t: state transition interval


v_R = sigma_u*eye(2); 
T = 1;
G = [T^2/2 0;0 T^2/2;T 0;0 T];
F = [1,0,1,0;0,1,0,1;0,0,1,0;0,0,0,1];
R = chol(v_R);
y = F * x + G * ( R.' * randn(2,1) );
