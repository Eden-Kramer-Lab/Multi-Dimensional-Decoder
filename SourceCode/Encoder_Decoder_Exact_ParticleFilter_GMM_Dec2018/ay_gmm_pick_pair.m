function [ind_a,ind_b,Ls]=ay_gmm_pick_pair(Param,ptc_p)
% this funnction returns two possible mixtures which might be merged to one
min_eps = realmin('double');
Ps = zeros(size(ptc_p,1),Param.n_mix);
Ws = zeros(size(ptc_p,1),Param.n_mix);
for m=1:Param.n_mix
    Ps(:,m) = mvnpdf(ptc_p,Param.Model{m}.M,0.5*(Param.Model{m}.S+Param.Model{m}.S'));
    Ws(:,m) = Param.Model{m}.W;
end

Ls    = realmax;
ind_a = -1;
ind_b = -1;
for i=1:Param.n_mix
    for j=i+1:Param.n_mix
        % rescale-factor for weights
        scale = 1/(1-Param.Model{i}.W-Param.Model{j}.W);
        % rescale weights
        tWs   = scale*Ws; tWs(:,i)=0;tWs(:,j)=0;
        % calculate rescaled deviance
        temp  = -2*sum(log(min_eps+sum(tWs.*Ps,2)));
        % we remove pairs which gives the minimum change in the likelihood
        % estimate
        if temp < Ls
            Ls    = temp;
            ind_a = i;
            ind_b = j;
        end
    end
end

