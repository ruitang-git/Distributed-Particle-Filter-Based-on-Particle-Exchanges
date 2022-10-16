% script
%
% plot_GMM
%
% description:
% plot Gaussian Mixture to check the convergence performance of 
% fusion step
%

if C==1
    figure;
    for n_sensor = 1:N_sensor
         a_k = DATA{n_sensor,4}; mu_k = DATA{n_sensor,5}; Sigma_k = DATA{n_sensor,6};
         Y=a_k*exp(-(X-mu_k(1)).^2./(2*Sigma_k(1,1)))/(2*pi*sqrt(Sigma_k(1,1)));
         plot(X,Y,'-');
         hold on;
    end
    Y=a_cen*exp(-(X-mu_cen(1)).^2./(2*sigma_cen(1,1)))/(2*pi*sqrt(sigma_cen(1,1)));
    plot(X,Y,'--');
    hold off;
elseif C==2
    figure;
    for n_sensor = 1:N_sensor
         a_k1 = DATA{n_sensor,4}(1); mu_k1 = DATA{n_sensor,5}(:,1); Sigma_k1 = DATA{n_sensor,6}(:,1:6);
         a_k2 = DATA{n_sensor,4}(2); mu_k2 = DATA{n_sensor,5}(:,2); Sigma_k2 = DATA{n_sensor,6}(:,7:12);
         Y1=a_k1*exp(-(X-mu_k1(1)).^2./(2*Sigma_k1(1,1)))/(2*pi*sqrt(Sigma_k1(1,1)));
         Y2=a_k2*exp(-(X-mu_k2(1)).^2./(2*Sigma_k2(1,1)))/(2*pi*sqrt(Sigma_k2(1,1)));
         plot(X,Y1,X,Y2,'-');
         hold on;
    end
    Y1=a_cen(1)*exp(-(X-mu_cen(1,1)).^2./(2*sigma_cen(1,1)))/(2*pi*sqrt(sigma_cen(1,1)));
    Y2=a_cen(2)*exp(-(X-mu_cen(1,2)).^2./(2*sigma_cen(1,7)))/(2*pi*sqrt(sigma_cen(1,7)));
    plot(X,Y1,'--',X,Y2,'--');
    hold off;
else
end