%% fine_tuning
function [a, mu, sigma] = fine_tuning(data)
%
% Input:
% data: cell data structure
% { [p(1)]  [mu(1)]  [sigma(1)] }
% { [p(nk)] [mu(nk)] [sigma(nk)] }
% { [p(Nk)] [mu(Nk)] [sigma(Nk)] }

d = size(data{1,2},1); % dimension
C = length(data{1,1}); % number of components
Nk = size(data,1); % number of neighbors (including itself)
a = []; 
mu = [];
sigma = [];

for c  = 1:C
    a_o = data{1,1}(c);
    mu_o = data{1,2}(:,c);
    sigma_o = data{1,3}(:,(c-1)*d+1:c*d);
    a_list = a_o;
    mu_list = mu_o;
    sigma_list = sigma_o;
    for nk = 2:Nk
        kld = []; 
        C_s = C-c+1;
        for c_s = 1:C_s
            mu_s = data{nk,2}(:,c_s);
            sigma_s = data{nk,3}(:,(c_s-1)*d+1:c_s*d);
            D = KLD(mu_o, mu_s, sigma_o, sigma_s, d);
            kld = [kld D];
        end
        sn = find(kld == min(kld)); % serial number
        a_list = [a_list data{nk,1}(sn)];
        mu_list = [mu_list data{nk,2}(:, sn)];
        sigma_list = [sigma_list data{nk,3}(:,(sn-1)*d+1:sn*d)];
        data{nk,1}(sn) = [];
        data{nk,2}(:, sn) = [];
        data{nk,3}(:,(sn-1)*d+1:sn*d) = [];
    end
    
    % averaging
    a = [a; sum(a_list)/Nk];
    a_list = a_list/sum(a_list);
    mu_avg = mu_list * a_list'; % centralised mean
    sigma_avg = zeros(d); % centralised covariance
    for nk = 1:Nk
        sigma_avg = sigma_avg + a_list(nk)*(sigma_list(:,1+d*(nk-1):d*nk)+...
            (mu_avg-mu_list(:,nk))*(mu_avg-mu_list(:,nk))');
    end
    sigma_avg = sigma_avg .* eye(d);
    %
    mu = [mu mu_avg];
    sigma = [sigma sigma_avg];
end


        