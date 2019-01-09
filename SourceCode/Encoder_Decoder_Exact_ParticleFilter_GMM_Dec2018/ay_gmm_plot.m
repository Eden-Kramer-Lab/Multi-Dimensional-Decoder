function Imag =ay_gmm_plot(MixModel,Xs,Ys)
%% run filter
[tx,ty]=meshgrid(Xs,Ys);
 sx  = reshape(tx,numel(tx),1);
 sy  = reshape(ty,numel(ty),1);
 sxy = [sx sy];
 L   = zeros(length(sx),1);
 
 for n=1:MixModel.n_mix
     S = 0.5*(MixModel.Model{n}.S+MixModel.Model{n}.S');
     M = MixModel.Model{n}.M;
     L = L + MixModel.Model{n}.W * mvnpdf(sxy,M,S);
 end
 Imag =reshape(L,size(tx))';

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
