function [Ps,Pw] = ay_gen_particle_mem(Ns,Mo,Prior)
global MarkIntensity
global MarkInfo

MIN_DEV = 1e-6;
%% This function generates partciles on spike time
%  Type, defines how mixture weight is updated
%  Ns,   defines number of samples
%  Mo,   is the observed mark information
%  Prior,is the mixture model prior
%  Mark, is the data used for online kernel

%% It defines prior points per each mixture - s points
% number of points correspond to number of mixure
sp = zeros(Prior.n_mix,length(Prior.Model{1}.M));
for s=1:Prior.n_mix
        sp(s,:)= Prior.Model{s}.M;
end

%% Calculate likelihood of sp points
% here, we have a 2D example here
% here, mark is the cell index
lp = ones(size(sp,1),1);
for i=1:length(lp)
    % calculate likelihood of each points
    lp(i) = ay_point_likelihood_mem(sp(i,:),Mo);
end


%% Changes Priors of Models
for i=1:length(lp)
    lp(i) = max(realmin,lp(i) * Prior.Model{i}.W); 
end
Ws = lp/max(realmin,sum(lp));

%% Generate Particles
Mu = [];
Sigma = [];
for i=1:length(lp)
    Mu = [Mu;Prior.Model{i}.M];
    Sigma(:,:,i) = Prior.Model{i}.S;
end
Ps = ay_mvgmmrnd(Mu,Sigma,Ws,Ns);

%% Generate Particles
ptcParam = Prior;
for i=1:length(lp)
    ptcParam.Model{i}.W = Ws(i);
end
Pw = ones(Ns,1);
for i=1:Ns
    % calculate likelihood given mixture
    Lb = ay_gmm_likelihood(3,ptcParam,Ps(i,:));
    % calculate likelihood given prior
    La = ay_gmm_likelihood(3,Prior,Ps(i,:));
    % calculate prior
    Lt = ay_point_likelihood_mem(Ps(i,:),Mo);
    % define W
    Pw(i) = (La*Lt)/max(realmin,Lb);
end
Pw = Pw/max(realmin,sum(Pw));


%% Resample Partciles
% we add this step to avoid bad samples
cPw= cumsum(Pw);
cPs= zeros(size(Ps));
for i=1:Ns
     [~,ind] = min(abs(cPw-rand));
     cPs(i,:)= Ps(ind,:)+MIN_DEV * randn(size(Ps(1,:)));
end
Ps = cPs;
Pw = ones(Ns,1)/Ns;

    
end
    
    