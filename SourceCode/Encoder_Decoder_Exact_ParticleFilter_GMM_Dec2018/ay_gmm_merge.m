function eParam = ay_gmm_merge(mode,Param,ptc_p,iter,extra)
COV_DEGENERACY = 1e2;
MIN_EIG        = 1e-3;
ADD_EIG        = 1e-3;


dim   = size(ptc_p,2);
% mode==1, adding a new component
if mode==1
   n_mix        = Param.n_mix;
   eParam.n_mix = n_mix+1;
   eParam.iter  = iter;
   for i=1:n_mix
       eParam.Model{i}.M  = Param.Model{i}.M;
       eParam.Model{i}.S  = Param.Model{i}.S;
       eParam.Model{i}.W  = Param.Model{i}.W; 
   end
   % update new parameter
   Ps = zeros(length(ptc_p),1);
   for i=1:eParam.n_mix-1
        tS   = 0.5*(eParam.Model{i}.S+eParam.Model{i}.S');
        Ps   = Ps + mvnpdf(ptc_p,eParam.Model{i}.M,tS)* eParam.Model{i}.W;
   end
   % candidate initial
   tPs    = median(Ps);
   ind    = find(Ps>tPs);
   M      = sum(ptc_p(ind,:))/max(realmin,length(ind));
   ptc_mp = ptc_p(ind,:)-repmat(M,length(ind),1);
   S      = (sum(Ps(ind))/sum(Ps))*(ptc_mp'*ptc_mp)/max(realmin,length(ind)); 
   alpha  = 0.5; 
   m_eig = min(eig(S));
   if m_eig < MIN_EIG ||  cond(S)>COV_DEGENERACY
         eParam.Model{eParam.n_mix}.S = S + (MIN_EIG-m_eig)*eye(size(S));
   else
         eParam.Model{eParam.n_mix}.S = S;
   end
   eParam.Model{eParam.n_mix}.M = M ;
   % main loop
   for i=1:eParam.iter
       % likelihood given new mixture
       tS   = 0.5*(eParam.Model{end}.S+eParam.Model{end}.S');
       Pt = mvnpdf(ptc_p,eParam.Model{end}.M,tS);
       % find ht
       ht    = (alpha * Pt) ./ (alpha * Pt + (1-alpha) * Ps); 
       alpha = sum(ht)/length(ht);
       % calculate new mean
       Mx  = sum(repmat(ht,1,dim).*ptc_p)/sum(ht);
       % calculate new variance
       ptc_mp =  ptc_p-repmat(Mx,size(ptc_p,1),1);
       Sx     = ((repmat(ht,1,dim).*ptc_mp)'*ptc_mp)/max(realmin,sum(ht)); 
       m_eig = min(eig(Sx));
       if m_eig < MIN_EIG ||  cond(Sx)>COV_DEGENERACY
       % update the last one
            eParam.Model{end}.S = Sx + (MIN_EIG-m_eig)*eye(size(Sx));
       else
           eParam.Model{end}.S = Sx;
       end
       eParam.Model{end}.M = Mx;
   end
   % updater new term
   for i=1:eParam.n_mix-1
       eParam.Model{i}.W  = (1-alpha) * Param.Model{i}.W; 
   end
   eParam.Model{end}.W = alpha;
end
% mode==2, dropping component
if mode==2
   % This is the original number of mixtures 
   n_mix = Param.n_mix;
   % This is the new 
   eParam.n_mix = n_mix-1;
   eParam.iter  = iter;
   scale = 1/(1-Param.Model{extra(1)}.W-Param.Model{extra(2)}.W);
   ind   = 1;
   for i=1:n_mix
       if i~=extra(1) && i~=extra(2)
            eParam.Model{ind}.M  = Param.Model{i}.M;
            eParam.Model{ind}.S  = Param.Model{i}.S;
            eParam.Model{ind}.W  = Param.Model{i}.W * scale; 
            ind = ind+1;
       end
   end
   alpha = Param.Model{extra(1)}.W + Param.Model{extra(2)}.W;
   eParam.Model{eParam.n_mix}.M = (Param.Model{extra(1)}.M*Param.Model{extra(1)}.W + Param.Model{extra(2)}.M*Param.Model{extra(2)}.W)/alpha;
   eParam.Model{eParam.n_mix}.S = (Param.Model{extra(1)}.S*Param.Model{extra(1)}.W^2 + Param.Model{extra(2)}.S*Param.Model{extra(2)}.W^2)/alpha^2;
   m_eig = min(eig(eParam.Model{eParam.n_mix}.S));
   if m_eig< MIN_EIG ||  cond(eParam.Model{eParam.n_mix}.S)>COV_DEGENERACY
         eParam.Model{eParam.n_mix}.S = eParam.Model{eParam.n_mix}.S + (MIN_EIG-m_eig)*eye(size(eParam.Model{eParam.n_mix}.S));
   end
   % update new parameters
   Ps = zeros(length(ptc_p),1);
   for i=1:eParam.n_mix-1
       tS   = 0.5*(eParam.Model{i}.S+eParam.Model{i}.S');
       temp = mvnpdf(ptc_p,eParam.Model{i}.M,tS);
       Ps   = Ps + temp * eParam.Model{i}.W;
   end
   for i=1:eParam.iter
       % likelihood given new mixture
       tS   = 0.5*(eParam.Model{end}.S+eParam.Model{end}.S');
       Pt = mvnpdf(ptc_p,eParam.Model{end}.M,tS);
       % find ht
       ht    = (alpha * Pt) ./ max(realmin,alpha * Pt + (1-alpha) * Ps); 
       alpha = sum(ht)/length(ht);
       % calculate new mean
       Mx  = sum(repmat(ht,1,dim).*ptc_p)/sum(ht);
       % calculate new variance
       ptc_mp =  ptc_p-repmat(Mx,size(ptc_p,1),1);
       Sx     = ((repmat(ht,1,dim).*ptc_mp)'*ptc_mp)/max(realmin,sum(ht)); 
       % update the last one
       m_eig = min(eig(Sx));
       if m_eig< MIN_EIG ||  cond(Sx)>COV_DEGENERACY
            eParam.Model{end}.S = Sx+ (MIN_EIG-m_eig)*eye(size(Sx));
       else
           eParam.Model{end}.S  = Sx; 
       end
       eParam.Model{end}.M = Mx;
   end
   % updater new term
   for i=1:eParam.n_mix-1
       eParam.Model{i}.W  = (1-alpha) * eParam.Model{i}.W; 
   end
   eParam.Model{end}.W = alpha;
end
