function L = ay_gmm_likelihood(mode,param,ptc_p)

% calculate likelihood
L = zeros(size(ptc_p,1),1);
for i=1:param.n_mix
    mu  = param.Model{i}.M;
    cov = 0.5*(param.Model{i}.S+param.Model{i}.S');
    %cov(2,1)=cov(1,2);
    L  = L + mvnpdf(ptc_p,mu,cov) * param.Model{i}.W; 
end
% total likelihood plus penalty term
if mode==1  % BIC
    ns = size(ptc_p,1);
    L  = -2*sum(log(L))+2*log(ns)*(param.n_mix + param.n_mix *(2*length(mu)+(length(mu)*(length(mu)-1))/2));
end
if mode==2% AIC
    
    L  = -2*sum(log(L))+2*(param.n_mix + param.n_mix *(2*length(mu)+(length(mu)*(length(mu)-1))/2));
end
if mode==3 % liklihood
end