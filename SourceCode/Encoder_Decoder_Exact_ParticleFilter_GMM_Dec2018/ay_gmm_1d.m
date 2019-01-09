function L =ay_gmm_1d(MixModel,Xs,dir)
%% run filter
L   = zeros(length(Xs),1);
for n=1:MixModel.n_mix
     S = 0.5*(MixModel.Model{n}.S + MixModel.Model{n}.S');
     M = MixModel.Model{n}.M;
     Md = M(dir);
     Sd = S(dir,dir);
     L = L + MixModel.Model{n}.W * normpdf(Xs',Md,sqrt(Sd));
end


% % copy Prior to Post
% Imag = zeros(length(Xs),length(Ys));
% % Fill vectors
% for i=1:length(Xs)
%     for j=1:length(Ys)
%         for n=1:MixModel.n_mix
%             S = 0.5*(MixModel.Model{n}.S+MixModel.Model{n}.S');
%             M = MixModel.Model{n}.M;
%             L = mvnpdf([Xs(i) Ys(j)],M,S);
%             Imag(i,j)=Imag(i,j) + MixModel.Model{n}.W * L;
%         end
%     end
% end
