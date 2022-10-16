%% convergence analysis
% cvg_analysis

% consensus无法收敛，第一第二次迭代后KLD显著下降，但随后稳步攀升

% centralised Gaussian approximation
mu_cg = mu_cen * a_cen; % centralised mean
sigma_cg = zeros(d); % centralised covariance
for c = 1:C
    sigma_cg = sigma_cg + a_cen(c)*(sigma_cen(:,1+d*(c-1):d*c)+...
        (mu_cg-mu_cen(:,c))*(mu_cg-mu_cen(:,c))');
end

% distributed Gaussian approximation(@1st sensor)
cvg_error = [];
%
%[e(@s1,t) e(@s2,t) ...]^T
%
for n_sensor = 1:N_sensor
    a_de = DATA{n_sensor,4};
    mu_de = DATA{n_sensor,5};
    sigma_de = DATA{n_sensor,6};
    mu_dg = mu_de * a_de; % centralised mean
    sigma_dg = zeros(d); % centralised covariance
    for c = 1:C
        sigma_dg = sigma_dg + a_de(c)*(sigma_de(:,1+d*(c-1):d*c)+...
            (mu_dg-mu_de(:,c))*(mu_dg-mu_de(:,c))');
    end

    % calculate KLD
    kld = KLD(mu_dg,mu_cg,sigma_dg,sigma_cg,d);
    
    cvg_error = [cvg_error; kld];
end
    