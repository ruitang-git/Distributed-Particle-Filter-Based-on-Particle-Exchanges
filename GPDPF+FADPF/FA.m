function x_P = FA(x_P, cd, h, R, iters_FA, iters_com, neighbors, protocol, x_true)
%%% GA
%%% Firefly Algorithm Based Fusion
% parameters
alpha = 1;
gamma = 0.03;
%alpha = 0.16;

beta0 = 0.5;
N_sensor = length(x_P);
N = size(x_P{1,1},2); % local particle size
d = size(x_P{1,1},1);
sigma_u = 0.5; t = 1;
Sigma = sigma_u^2 * ...
     [1/20*t^5 0 1/8*t^4 0 1/6*t^3 0;...
      0 1/20*t^5 0 1/8*t^4 0 1/6*t^3;...
      1/8*t^4 0 1/3*t^3 0 1/2*t^2  0;...
      0 1/8*t^4 0 1/3*t^3 0 1/2*t^2 ;...
      1/6*t^3 0 1/2*t^2 0 t 0       ;...
      0 1/6*t^3 0 1/2*t^2 0 t      ];
% Sigma = sigma_u^2 * ...
%      [1/20*t^5 0 1/8*t^4 0 ;...
%       0 1/20*t^5 0 1/8*t^4 ;...
%       1/8*t^4 0 1/3*t^3 0  ;...
%       0 1/8*t^4 0 1/3*t^3  ];
  
% figure;
% plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.'); hold on;
% plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
if strcmp(protocol, 'gossip')
    
%     for n_sensor = 1:N_sensor
%         for n = 1:N
%             x_P{n_sensor,1}(:,n) = x_P{n_sensor,1}(:,n)+alpha*chol(Sigma).'*randn(d,1);
%         end
%     end
    selected_times = zeros(1,N_sensor);
    for iter_com = 1:iters_com   
        s1 = randperm(N_sensor, 1); s2 = neighbors{s1}(randperm(length(neighbors{s1}),1));
        selected_times(s1) = selected_times(s1)+1;selected_times(s2) = selected_times(s2)+1;
        x1 = x_P{s1,1}; x2 = x_P{s2,1}; 
        I = cell(2,2); % brightness of fireflies(particles),defined using particles' weight
        for n = 1:size(x1,2)
            [hp,~] = observation(x1(:,n),cd(:,s1),0);
            I{1,1}(n) = mvnpdf(hp', h{s1,1}', R); % particle from s1, evaluted at s1;
            [hp,~] = observation(x2(:,n),cd(:,s1),0);
            I{1,2}(n) = mvnpdf(hp', h{s1,1}', R); % particle from s2, evaluted at s1;
            [hp,~] = observation(x1(:,n),cd(:,s2),0);
            I{2,1}(n) = mvnpdf(hp', h{s2,1}', R); % particle from s1, evaluted at s2;
            [hp,~] = observation(x2(:,n),cd(:,s2),0);
            I{2,2}(n) = mvnpdf(hp', h{s2,1}', R); % particle from s2, evaluted at s2;
        end
        
        %
        I{1,1} = I{1,1}.*I{2,1};
        I{2,2} = I{1,2}.*I{2,2};
        % sort(rearrange)
%         [I_copy{1,1},idx] = sort(I{1,1}); x1_copy = x1(:,idx);
%         [I_copy{2,2},idx] = sort(I{2,2}); x2_copy = x2(:,idx);
       x1_copy = x1; x2_copy = x2;
%         I_copy = I;
%         x = [x1 x2]; I_fuse = [I{1,1} I{2,2}]; I_copy = I_fuse; x_copy = x;
        %
        
%         if mod(iter_com,80) == 0
%             figure;
%             plot(x1(1,:),x1(2,:),'.'); hold on;
%             plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%             debug = 1;
%         end

        %if iter_com < 80
%             figure;
%             plot(x1(1,:),x1(2,:),'.'); hold on;
%             plot(x2_copy(1,:),x2_copy(2,:),'.'); hold on;

        % apply FA at s1
        for i = 1:N % particle from s1(operation on this particle set)
            %for j = 1:N % particle from s2 
                %if I{1,1}(i)<I{2,2}(j)
                    idx = find(I{2,2}>I{1,1}(i));
                    if ~isempty(idx)
                        j = idx(randperm(length(idx),1));
                        r = norm(x1(:,i)-x2_copy(:,j)); % distance 
    %                     x1(:,i) = x1(:,i)+beta0/(1+gamma*r^2)*(x2_copy(:,j)-x1(:,i))+alpha*chol(Sigma).'*randn(d,1);
                        x1(:,i) = x1(:,i)+beta0/(1+gamma*r^2)*(x2_copy(:,j)-x1(:,i));
    %                     noise = alpha*chol(Sigma).'*[randn(d-2,1);0.1*randn;0.1*randn]; %noise(end) = 0; noise(end-1) = 0;
                        noise = alpha*chol(Sigma).'*randn(d,1); noise(d-1:d) = 1/(1+exp(1*(selected_times(s1)-N_sensor)))*noise(d-1:d);
                        x1(:,i) = x1(:,i) + noise;
                    end
                    %x_tmp = x1(:,i); x2_tmp = x2_copy(:,j); x_tmp(end) = x2_tmp(end); x_tmp(end-1) = x2_tmp(end-1); x1(:,i) = x_tmp;
%                        x_update = x1(:,i)+beta0/(1+gamma*r^2)*(x2_copy(:,j)-x1(:,i))+alpha*chol(Sigma).'*(rand(d,1)-0.5);
%                         [hp,~] = observation(x_update,cd(:,s1),0);
%                         I_update = mvnpdf(hp', h{s1,1}', R); 
%                         if I_update > 0.1*I{1,1}(i)
%                        x1(:,i) = x_update; %I{1,1}(i) = I_update;
                     %end
                    %break;
                %end
            %end
        end

        % apply FA at s2
        for i = 1:N % particle from s2(operation on this particle set)
            %for j = 1:N % particle from s1
            idx = find(I{1,1}>I{2,2}(i)); 
            if ~isempty(idx)
                j = idx(randperm(length(idx),1));
                %if I{2,2}(i)<I{1,1}(j)
                r = norm(x1_copy(:,j)-x2(:,i)); % distance
%                     x2(:,i) = x2(:,i)+beta0/(1+gamma*r^2)*(x1_copy(:,j)-x2(:,i))+alpha*chol(Sigma).'*randn(d,1);
                x2(:,i) = x2(:,i)+beta0/(1+gamma*r^2)*(x1_copy(:,j)-x2(:,i));
%                     noise = alpha*chol(Sigma).'*[randn(d-2,1);0.1*randn;0.1*randn]; %noise(end) = 0; noise(end-1) = 0;
                noise = alpha*chol(Sigma).'*randn(d,1); noise(d-1:d) = 1/(1+exp(1*(selected_times(s2)-N_sensor)))*noise(d-1:d);
                x2(:,i) = x2(:,i) + noise;
            end
                    %x_tmp = x2(:,i); x1_tmp = x1_copy(:,j); x_tmp(end) = x1_tmp(end); x_tmp(end-1) = x1_tmp(end-1); x2(:,i) = x_tmp;
%                         x_update = x2(:,i)+beta0/(1+gamma*r^2)*(x1_copy(:,j)-x2(:,i))+alpha*chol(Sigma).'*(rand(d,1)-0.5);
%                         [hp,~] = observation(x_update,cd(:,s2),0);
%                         I_update = mvnpdf(hp', h{s2,1}', R); 
%                         if I_update > 0.1*I{2,2}(i)

%                         end
                    %break;
                %end
            %end
        end
%         if mod(iter_com,5) == 0
%             plot(x1(1,:),x1(2,:),'.'); hold on;
%             plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%             legend('before','after');
%             debug = 1;
%         end

        x_P{s1,1} = x1; x_P{s2,1} = x2;
        %if mod(iter_com,100) == 0
%         if iter_com == 100
%             figure;
%             plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.'); hold on;
%             plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%             debug = 1;
%         end
%         if iter_com == 5000
%             figure;
%             plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.'); hold on;
%             plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%             debug = 1;
%         end
    end
%     if mod(iter_com,100) == 0
%             figure;
%             plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.'); hold on;
%             plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%             de = 1;
%     end
elseif strcmp(protocol, 'broadcast')
%     for iter_com = 1:iters_com
%         for n_sensor = 1:N_sensor
%             for n = 1:N
%                 [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,n_sensor),0);
%                 I(1,n) = mvnpdf(hp', h{n_sensor,1}', R); % #neighbors*N
%                 for neighbor = 1:length(neighbors{1,n_sensor})
%                     [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,neighbors{1,n_sensor}(neighbor)),0);
%                     I(1+neighbor,n) = mvnpdf(hp', h{neighbors{1,n_sensor}(neighbor),1}', R); 
%                 end
%             end
%             for col = 2:size(I,1)
%                 I(1,:) = I(1,:).*I(col,:);
%             end
%             iters_FA_copy = iters_FA;
%             while iters_FA_copy > 0
%                 for i = 1:N % particle from n_sensor
%                     for j = 1:N % particle from n_sensor
%                         if I(1,i)<I(1,j)
%                             r = norm(x_P{n_sensor,1}(:,i)-x_P{n_sensor,1}(:,j)); % distance
%                             x_update = x_P{n_sensor,1}(:,i)+beta0/(1+gamma*r^2)*(x_P{n_sensor,1}(:,j)-x_P{n_sensor,1}(:,i))+alpha*chol(Sigma).'*randn(d,1);
%     %                         x1(:,i) = x1(:,i)+beta0*exp(-gamma*r^2)*(x2(:,j)-x1(:,i))+alpha*chol(Sigma).'*(rand(d,1)-0.5);
%                             %if I_update > I{1,1}(i)
%                                 x_P{n_sensor,1}(:,i) = x_update;
%                             %end
%                         end
%                     end
%                 end
%                 iters_FA_copy = iters_FA_copy - 1;
%             end
%         end
%         figure;
%         plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.'); hold on;
%         plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%         debug = 1;
%     end
    for iter_com = 1:iters_com
        for n_sensor = 1:N_sensor
            for neighbor = 1:length(neighbors{1,n_sensor})
                x_P{neighbors{1,n_sensor}(neighbor),1} = ...
                    [x_P{neighbors{1,n_sensor}(neighbor),1} x_P{n_sensor ,1}(:,1:N)];
            end
        end
        for n_sensor = 1:N_sensor
             I{n_sensor,1} = zeros(1,size(x_P{n_sensor,1},2));
             for n = 1:size(x_P{n_sensor,1},2)
                [hp,~] = observation(x_P{n_sensor,1}(:,n),cd(:,n_sensor),0);
                I{n_sensor,1}(n) = mvnpdf(hp', h{n_sensor,1}', R);
             end
             iters_FA_copy = iters_FA;
%                               figure;
%                 plot(x_P{n_sensor,1}(1,:),x_P{n_sensor,1}(2,:),'.'); hold on;
%                 plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
             while iters_FA_copy > 0
                 gamma = 0.03;
%                                  figure;
%                 plot(x_P{n_sensor,1}(1,1:N),x_P{n_sensor,1}(2,1:N),'.'); hold on;
                 [I{n_sensor,1}(1:N),idx] = sort(I{n_sensor,1}(1:N)); x_P{n_sensor,1}(:,1:N) = x_P{n_sensor}(:,idx);
                 for i = 1:N-1 % particle from n_sensor
                     j = i+randperm(N-i,1);
                     %for j = 1:N % particle from n_sensor
                         if I{n_sensor,1}(i)<I{n_sensor,1}(j)
                             r = norm(x_P{n_sensor,1}(:,i)-x_P{n_sensor,1}(:,j)); % distance
                             x_update = x_P{n_sensor,1}(:,i)+beta0/(1+gamma*r^2)*(x_P{n_sensor,1}(:,j)-x_P{n_sensor,1}(:,i))+alpha*chol(Sigma).'*randn(d,1);
    %                         x1(:,i) = x1(:,i)+beta0*exp(-gamma*r^2)*(x2(:,j)-x1(:,i))+alpha*chol(Sigma).'*(rand(d,1)-0.5);
                            %if I_update > I{1,1}(i)
                             [hp,~] = observation(x_update,cd(:,n_sensor),0);
                             I_update = mvnpdf(hp', h{n_sensor,1}', R); 
                             x_P{n_sensor,1}(:,i) = x_update; I{n_sensor,1}(i) = I_update;
                            %end
                         end
                     %end
                 end
                 [I{n_sensor,1}(1:N),idx] = sort(I{n_sensor,1}(1:N)); x_P{n_sensor,1}(:,1:N) = x_P{n_sensor}(:,idx);
%                 plot(x_P{n_sensor,1}(1,1:N),x_P{n_sensor,1}(2,1:N),'.'); hold on;
%                 plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%                 legend('before','after');
%                 debug = 1;
                 gamma = 0.03;
                 
%                  figure;
%                 plot(x_P{n_sensor,1}(1,N+1:end),x_P{n_sensor,1}(2,N+1:end),'.'); hold on;
                 [I{n_sensor,1}(N+1:end),idx] = sort(I{n_sensor,1}(N+1:end)); x_P{n_sensor,1}(:,N+1:end) = x_P{n_sensor}(:,N+idx);
                 for i = N+1:size(x_P{n_sensor,1},2)-1 % particle from n_sensor
                     j = i+randperm(size(x_P{n_sensor,1},2)-i,1);
                     %for j = N+1:size(x_P{n_sensor,1},2) % particle from n_sensor
                         if I{n_sensor,1}(i)<I{n_sensor,1}(j)
                             r = norm(x_P{n_sensor,1}(:,i)-x_P{n_sensor,1}(:,j)); % distance
                             x_update = x_P{n_sensor,1}(:,i)+beta0/(1+gamma*r^2)*(x_P{n_sensor,1}(:,j)-x_P{n_sensor,1}(:,i))+alpha*chol(Sigma).'*randn(d,1);
    %                         x1(:,i) = x1(:,i)+beta0*exp(-gamma*r^2)*(x2(:,j)-x1(:,i))+alpha*chol(Sigma).'*(rand(d,1)-0.5);
                            %if I_update > I{1,1}(i)
                             [hp,~] = observation(x_update,cd(:,n_sensor),0);
                             I_update = mvnpdf(hp', h{n_sensor,1}', R); 
                             x_P{n_sensor,1}(:,i) = x_update; I{n_sensor,1}(i) = I_update;
                            %end
                         end
                     %end
                 end
                 [I{n_sensor,1}(N+1:end),idx] = sort(I{n_sensor,1}(N+1:end)); x_P{n_sensor,1}(:,N+1:end) = x_P{n_sensor}(:,N+idx);
                 iters_FA_copy = iters_FA_copy - 1;

%                 plot(x_P{n_sensor,1}(1,N+1:end),x_P{n_sensor,1}(2,N+1:end),'.'); hold on;
%                 plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%                                 legend('before','after');
%                 debug = 1;
             end
             
             [~, idx1] = sort(I{n_sensor,1}(1:N),'descend');
             [~, idx2] = sort(I{n_sensor,1}(N+1:end),'descend');
             idx = [idx1(1:N/2) idx2(1:N/2)];
             x_P{n_sensor,1} = x_P{n_sensor,1}(:, idx);
        end
%         figure;
%         plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.'); hold on;
%         plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
%         debug = 1;
    end
end
% figure;
% plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.'); hold on;
% plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
% figure; plot(x_P{1,1}(1,:),x_P{1,1}(2,:),'.');hold on;plot(x_P{2,1}(1,:),x_P{2,1}(2,:),'.');hold on;plot(x_P{3,1}(1,:),x_P{3,1}(2,:),'.');hold on;plot(x_P{4,1}(1,:),x_P{4,1}(2,:),'.');hold on;plot(x_P{5,1}(1,:),x_P{5,1}(2,:),'.');hold on;
% plot(x_P{6,1}(1,:),x_P{6,1}(2,:),'.');hold on;plot(x_P{7,1}(1,:),x_P{7,1}(2,:),'.');hold on;plot(x_P{8,1}(1,:),x_P{8,1}(2,:),'.');hold on;plot(x_P{9,1}(1,:),x_P{9,1}(2,:),'.');hold on;
% plot(x_true(1),x_true(2),'kx','LineWidth',2,'MarkerSize',10); hold off;
debug = 1;

