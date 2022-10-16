% function plotparticles(a1,mu1,sigma1,a,mu,sigma,benchmark,N)
% function plotparticles(a,mu,sigma,benchmark,N,x)
function plotparticles(a1,mu1,sigma1,a,mu,sigma,x,N)
p = GMsample(a,mu,sigma,N);
p1 = GMsample(a1,mu1,sigma1,N);
figure;
plot(p(1,:),p(2,:),'b.');
hold on;
plot(p1(1,:),p1(2,:),'r.');
hold on;
plot(x(1),x(2),'kx','LineWidth',2,'MarkerSize',10);
hold off
end