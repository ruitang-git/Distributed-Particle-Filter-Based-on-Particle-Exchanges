% script
%
% description:
% used to analyze the GLA error with varying parameter m
%

error_GLA = []; % initialize
plt = 0;
P_w_pre = P_w;
for m = 2:50:N
    if m == 2
        plt = 1;
    end
    if mod(m,100) == 2
        fprintf('%dth finished!\n',m);
    end
    [P_w_hat,V] = graphLap(x_P_update, N, P_w_pre, m, k, plt); % approximation 
    % ----- check the approximation error ------ %
    P_w = V*P_w_hat;
    error_GLA = [error_GLA norm(P_w-P_w_pre)];
    plt = 0;
end
figure;
x = 2:10:N;
plot(x,error_GLA);
xlabel('m');
ylabel('norm(W-$\hat{W}$)');
