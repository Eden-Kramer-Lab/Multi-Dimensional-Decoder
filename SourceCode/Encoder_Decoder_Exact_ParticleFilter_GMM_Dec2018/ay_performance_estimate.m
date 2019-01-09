function [inHPD,SErr,MSErr] =ay_performance_estimate(mode,RealXY,L,Xs,Ys,HPD)
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
dXs = Xs(2)-Xs(1);
dYs = Ys(2)-Ys(1);

if mode==1    % Real Posterior
   Lo = L; 
end
if mode==2    % Mixture Model
    Lo = zeros(length(Xs),length(Ys));
    for i=1:length(Xs)
        for j=1:length(Ys)
            Lo(i,j)=gmm_likelihood([Xs(i) Ys(j)],L);
        end
    end
end
% put in an array format
Lo  = Lo/sum(Lo(:));
[n_row,n_col]= size(Lo);
Ls  = reshape(Lo,[],1);
[Ls,ind] = sort(Ls,'descend');
ind_c    = 1+floor((ind-1)./n_row); 
ind_r    = ind - (ind_c-1).*n_row;
% find the HPD range
temp_ind  = find(cumsum(Ls)>=HPD);
if ~isempty(temp_ind)
    end_ind   = temp_ind(1);
else
    end_ind   = 1;
end
% find distance
dE = (Xs(ind_r(1:end_ind))-RealXY(1)).^2 + (Ys(ind_c(1:end_ind))-RealXY(2)).^2;
dE = sqrt(dE);
% find minimum
inHPD = 0;
if min(dE) <= sqrt(dXs*dXs+dYs*dYs)
   inHPD = 1;
end
MSErr = dE(1);
SErr  = sqrt((Xs*sum(Lo,2)-RealXY(1)).^2 + (sum(Lo,1)*Ys'-RealXY(2)).^2);



    function ML = gmm_likelihood(ptc_p,param)
        ML= 0;
        for z=1:param.n_mix
            mu  = param.Model{z}.M;
            cov = param.Model{z}.S;
            tcov = 0.5 * (cov+cov');
            ML  = ML + mvnpdf(ptc_p,mu,tcov) * param.Model{z}.W; 
        end
    end
end