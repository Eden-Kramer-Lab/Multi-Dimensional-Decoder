function [inHPD,SErr,MSErr,PSErr,AreaInd,MX] = ay_performance_estimate_xx(RealXY,L,PV,dE,HPD,PXY)
%% Input
% mode = 1 means L is the real posterior
% mode = 2 means L carries mixture model
% L: either Posterior or Mixture Information
% Xs, Ys are the samples L is calculated or being calculates
% HPD, is the the HPD range
%% Output
% inHPD is either 1 or 0
% SErr: Root Squared Error between RealXY and Estimate
% MSErr: Root Squared Error between RealXY and Estimate on its maximum
dXs = dE(1);
dYs = dE(2);
dVx = dE(3);
dVy = dE(4);
% likelihood
Lo  = gmm_likelihood(PV,L);
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
% find distance
dEE = (PV(ind(1:end_ind),1)-RealXY(1)).^2 + (PV(ind(1:end_ind),2)-RealXY(2)).^2 + ...
     (PV(ind(1:end_ind),3)-RealXY(3)).^2 + (PV(ind(1:end_ind),4)-RealXY(4)).^2;
dEE = sqrt(dEE);
% find minimum
inHPD = 0;
if min(dEE) <= sqrt(dXs*dXs+dYs*dYs+dVx*dVx+dVy*dVy)
   inHPD = 1;
end
MSErr = dEE(1);
MX    = PV'*Lo;
SErr  = sqrt((MX(1)-RealXY(1)).^2 + (MX(2)-RealXY(2)).^2 +(MX(3)-RealXY(3)).^2 +(MX(4)-RealXY(4)).^2);
%%---
Lo = gmm_likelihood_xy(PXY,L);
Lo  = Lo/sum(Lo);
[Ls,ind] = sort(Lo,'descend');
temp_ind  = find(cumsum(Ls)>=HPD);
if ~isempty(temp_ind)
    end_ind   = temp_ind(1);
else
    end_ind   = 1;
end
AreaInd = end_ind * dE(1) * dE(2);
MX     = PXY'*Lo;
PSErr  = sqrt((MX(1)-RealXY(1)).^2 + (MX(2)-RealXY(2)).^2 );


    function ML = gmm_likelihood(ptc_p,param)
        ML= zeros(size(ptc_p,1),1);
        for z=1:param.n_mix
            mu  = param.Model{z}.M;
            cov = 0.5 * (param.Model{z}.S+param.Model{z}.S');
            ML  = ML + mvnpdf(ptc_p,mu,cov) * param.Model{z}.W; 
        end
    end

    function ML = gmm_likelihood_xy(ptc_p,param)
        ML= zeros(size(ptc_p,1),1);
        for z=1:param.n_mix
            mu  = param.Model{z}.M(1:2);
            cov = 0.5 * (param.Model{z}.S+param.Model{z}.S');
            cov = cov(1:2,1:2);
            ML  = ML + mvnpdf(ptc_p,mu,cov) * param.Model{z}.W; 
        end
    end
end