function [inHPD,SErr,MSErr,AreaInd,MX] =ay_performance_estimate_2d(mode,RealXY,L,PV,dE,HPD)
dXs = dE(1);
dYs = dE(2);
% likelihood
if mode==1
    Lo  = L;
else
    Lo  = gmm_likelihood(PV,L);
end
Lo  = Lo/sum(Lo);
% put in an array format
[Ls,ind] = sort(Lo,'descend');
% find the HPD range
temp_ind  = find(cumsum(Ls)>=HPD);
if ~isempty(temp_ind)
    end_ind   = temp_ind(1);
else
    end_ind   = 1;
end
AreaInd = end_ind * dE(1)* dE(2);
% find distance
dE = (PV(ind(1:end_ind),1)-RealXY(1)).^2 + (PV(ind(1:end_ind),2)-RealXY(2)).^2;
dE = sqrt(dE);
% find minimum
inHPD = 0;
if min(dE) <= sqrt(dXs*dXs+dYs*dYs)
   inHPD = 1;
end
MSErr = dE(1);
%%---
MX    = PV'*Lo;
SErr  = sqrt((MX(1)-RealXY(1)).^2 + (MX(2)-RealXY(2)).^2);



    function ML = gmm_likelihood(ptc_p,param)
        ML= zeros(size(ptc_p,1),1);
        for z=1:param.n_mix
            mu  = param.Model{z}.M;
            cov = 0.5 * (param.Model{z}.S+param.Model{z}.S');
            ML  = ML + mvnpdf(ptc_p,mu,cov) * param.Model{z}.W; 
        end
    end
end