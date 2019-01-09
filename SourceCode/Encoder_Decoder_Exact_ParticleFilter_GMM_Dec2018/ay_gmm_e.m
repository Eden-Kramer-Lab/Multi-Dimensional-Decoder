function Param=ay_gmm_e(ptc_p,n_mix,n_iter)
COV_DEGENERACY = 1e2;
MIN_EIG        = 1e-3;
% data dimension
dim = size(ptc_p,2);
% create param
Param.n_mix  = n_mix;
Param.iter   = n_iter;

part_len = size(ptc_p,1)/Param.n_mix;
for m=1:Param.n_mix
    % minimum index
    ind_a = floor(max(1,(m-1)* part_len + 1));
    % maximum index
    ind_b = floor(min(m* part_len,size(ptc_p,1)));
    % calculate
    Param.Model{m}.M = mean(ptc_p(ind_a:ind_b,:));
    Param.Model{m}.S = cov(ptc_p(ind_a:ind_b,:));
    Param.Model{m}.S = Param.Model{m}.S + 1e-2 *eye(size(Param.Model{m}.S)); 
    Param.Model{m}.W = 1/Param.n_mix; 
end
% for gmm
%% run GMM
% first, normalize L
% GMM loop
for i=1:Param.iter
    % Calculate each point probability
    Ps = zeros(size(ptc_p,1),Param.n_mix);
    for m=1:Param.n_mix
        tS   = 0.5*(Param.Model{m}.S+Param.Model{m}.S');
        temp    = mvnpdf(ptc_p,Param.Model{m}.M,tS);
        Ps(:,m) = temp * Param.Model{m}.W;
    end
    Ps = Ps./repmat(max(realmin,sum(Ps,2)),1,n_mix);
    % Find Ws
    Ws = zeros(Param.n_mix,1);
    for m=1:Param.n_mix
        Ws(m) = sum(Ps(:,m));
    end
    Ws = Ws/max(realmin,sum(Ws));
    
    % Find Ms
    Mx  = zeros(Param.n_mix,dim);
    Sxx = zeros(Param.n_mix,dim,dim);
    for m=1:Param.n_mix
        % scale
        scale     = max(realmin,sum(Ps(:,m)));
        % mean
        Mx(m,:)   = sum(repmat(Ps(:,m),1,dim).*ptc_p)/scale;
        % covariance
        ptc_mp    = ptc_p-repmat(Mx(m,:),size(ptc_p,1),1);
        Sxx(m,:,:)= ((repmat(Ps(:,m),1,dim).*ptc_mp)'*ptc_mp)/scale;
    end
    % check valid as well
    valid_cov = ones(Param.n_mix,1); 
    w_scale   = 1;
    for m=1:Param.n_mix
        Param.Model{m}.M = Mx(m,:);
        Param.Model{m}.S = squeeze(Sxx(m,:,:));
        if cond(Param.Model{m}.S)> COV_DEGENERACY || min(eig(Param.Model{m}.S))< MIN_EIG
            valid_cov(m)=0;
            w_scale = w_scale - Ws(m);
        end
        Param.Model{m}.W = Ws(m); 
    end
    if Param.n_mix > 1
        if sum(1-valid_cov)== Param.n_mix
            tParam   = [];
            % eigen value
            m_eig = min(eig(Param.Model{1}.S));
            tParam.Model{1}.S = Param.Model{1}.S  + (MIN_EIG-m_eig)*eye(size(Param.Model{1}.S));
            % mean 
            tParam.Model{1}.M = Param.Model{1}.M*Param.Model{1}.W;
            for m=2:Param.n_mix 
                tParam.Model{1}.M = tParam.Model{1}.M + (Param.Model{m}.M*Param.Model{m}.W);
            end
            % weight
            tParam.Model{1}.W = 1;
            tParam.n_mix = 1;
            n_mix        = 1;
        else
            tParam   = [];
            % update
            ind_fill = 0;
            for m=1:Param.n_mix 
                if valid_cov(m)
                    ind_fill = ind_fill+1;
                    tParam.Model{ind_fill}.M = Param.Model{m}.M;
                    tParam.Model{ind_fill}.S = Param.Model{m}.S;
                    tParam.Model{ind_fill}.W = Param.Model{m}.W/max(realmin,w_scale);
                end
            end
            n_mix  = ind_fill;
            tParam.n_mix = ind_fill;
        end
        Param = tParam;
    else
        m_eig = min(eig(Param.Model{1}.S));
        if m_eig< MIN_EIG
            Param.Model{1}.S = Param.Model{1}.S  + (MIN_EIG-m_eig)*eye(size(Param.Model{1}.S));
        end
    end
    
    
end

end


