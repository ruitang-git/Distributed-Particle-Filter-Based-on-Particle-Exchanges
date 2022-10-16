% calculate KL divergence
function D = KLD(mu1,mu2,sigma1,sigma2,d)
D = 0.5*(log(det(sigma2)/det(sigma1))+trace(inv(sigma2)*sigma1)-d+...
    (mu1-mu2)'*inv(sigma2)*(mu1-mu2));
