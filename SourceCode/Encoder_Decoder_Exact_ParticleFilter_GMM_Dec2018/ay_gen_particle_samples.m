function cPs = ay_gen_particle(Ps,TransProc,Mo,Mark)
%% This function generates partciles on spike time
%  Type, defines how mixture weight is updated
%  Ns,   defines number of samples
%  Mo,   is the observed mark information
%  Prior,is the mixture model prior
%  Mark, is the data used for online kernel
MIN_DEV = 1e-9;

%% Draw Samples from Ps and TransProcess
%% Generate Particles
MU = TransProc.A*Ps';
Po = mvnrnd(MU',TransProc.S);
Ns = size(Ps,1);
%% Calculate likelihood of sp points
% here, we have a 2D example here
% here, mark is the cell index
lp = ones(Ns,1);
for i=1:Ns
    % calculate likelihood of each points
    lp(i) = ay_point_likelihood(Po(i,:),Mo,Mark,1);
end
Ws = lp/max(realmin,sum(lp));

%% Resample Partciles
cPw = cumsum(Ws);
cPs = zeros(size(Po));
for i=1:Ns
     [~,ind] = min(abs(cPw-rand));
     cPs(i,:)= Po(ind,:)+MIN_DEV * randn(size(Po(1,:)));
end

end
    
    