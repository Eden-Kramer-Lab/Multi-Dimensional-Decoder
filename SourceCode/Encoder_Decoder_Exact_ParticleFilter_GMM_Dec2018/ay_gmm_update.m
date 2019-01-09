function [Post,Update] = ay_gmm_update(Prior,Mark,steps)
max_spread  = 1000;
r_drop      = 1e-1;
n_xy        = 3;
%% This function updates GMMs over non-spiking period
%  Prior,is the mixture model prior
%  Mark, is the data used for online kernel
% Copy Prior to Posterior
Post = Prior;
% Run update Rule
Update = 0;
    % update
    Ws = zeros(Post.n_mix,1);
    for i =1:Post.n_mix
        temp_steps = steps;
        [G,H]= ay_grd_hessian(Prior.Model{i}.M,Mark,temp_steps);
        for xy=1:(n_xy-1)
            temp_steps = steps / 2^xy;
            [Gt,Ht]= ay_grd_hessian(Prior.Model{i}.M,Mark,temp_steps);
            G = G + Gt;
            H = H + Ht;
        end
        G = G /n_xy;
        H = H /n_xy;
        % find optimal r to build a psd covariance
        if all(eig(H)>0)
             Post.Model{i}.S = pinv(pinv(Prior.Model{i}.S)+H);
             Post.Model{i}.S = 0.5*(Post.Model{i}.S+Post.Model{i}.S'); 
             Post.Model{i}.M = (Prior.Model{i}.M' - Post.Model{i}.S *G)';
        else
             Post.Model{i}.S = Prior.Model{i}.S;
             Post.Model{i}.M = Prior.Model{i}.M;
        end
        Ws(i) = Post.Model{i}.W * weight_update(Prior.Model{i}.M,Prior.Model{i}.S,Post.Model{i}.M,Post.Model{i}.S);


        
        
%         r    = 1.0;
%         stop = true;
%         while stop
%            temp = pinv(Prior.Model{i}.S)+r*H;
%            if all(eig(temp)>eps)
%                 stop = false;
%             else
%                 r = max(0,r - r_drop);
%                 if r==0
%                     stop = false;
%                 end
%             end
%         end
%         if r>0
%              Post.Model{i}.S = pinv(pinv(Prior.Model{i}.S)+r*H);
%              Post.Model{i}.S = 0.5*(Post.Model{i}.S+Post.Model{i}.S'); 
%              % update mean
%              A = Prior.Model{i}.M' - Post.Model{i}.S *(G-0.5*(1-r)*H*Prior.Model{i}.M');
%              B = eye(size(Post.Model{i}.S))+0.5*(1-r)*Post.Model{i}.S*H;
%              Post.Model{i}.M = (pinv(B)*A)';
%          else
%              Post.Model{i}.S = Prior.Model{i}.S;
%              Post.Model{i}.M = Prior.Model{i}.M;
%          end
%          Ws(i) = Post.Model{i}.W * weight_update(Prior.Model{i}.M,Prior.Model{i}.S,Post.Model{i}.M,Post.Model{i}.S);
        
    end
    % update weights
    Ws = Ws/max(realmin,sum(Ws));
    for i =1:Post.n_mix
        Post.Model{i}.W = Ws(i);
        if max(eig(Post.Model{i}.S))> max_spread
            Update = 1;
        end
    end
    
    %% weight update
    function L = weight_update(pre_mean,pre_cov,post_mean,post_cov)
        % get likelihood first - nominator
        Mo = zeros(length(Mark.Cell),1);
        La = max(realmin,ay_point_likelihood(pre_mean,Mo,Mark));
        % get denominator
        Lb = mvnpdf(pre_mean,pre_mean,pre_cov);
        % calculate the rate
        Lc = mvnpdf(pre_mean,post_mean,post_cov);
        % calculate L
        L = (La*Lb)/Lc;
    end

end 
%     %% rate at a point
%     function [La,G,H] = pp_rate(pp,mo,hx,hy,hvx,hvy)
%         % generate points around pp and calculate its likelihood
%         La = zeros(3,3,3,3);
%         for ix=-1:1,
%            for iy=-1:1,
%               for ivx=-1:1,
%                  for ivy=-1:1,
%                      pp_e= [pp(1)+ix*hx  pp(2)+iy*hy  pp(3)+ivx*hvx  pp(4)+ivy*hvy];
%                      La(ix+2,iy+2,ivx+2,ivy+2) = log(eps+ay_point_likelihood(pp_e,mo,Mark,num_time_step));
%                  end
%               end
%            end
%         end
%         % Calculate gradient around pp
%         G(1,1) = (La(3,2,2,2)-La(1,2,2,2))/(2*hx);
%         G(2,1) = (La(2,3,2,2)-La(2,1,2,2))/(2*hy);
%         G(3,1) = (La(2,2,3,2)-La(2,2,1,2))/(2*hvx);
%         G(4,1) = (La(2,2,2,3)-La(2,2,2,1))/(2*hvy);
%         
%         % Calculate Hessian around pp
%         H(1,1) = (La(3,2,2,2)+La(1,2,2,2)-2*La(2,2,2,2))/hx^2;
%         H(2,2) = (La(2,3,2,2)+La(2,1,2,2)-2*La(2,2,2,2))/hy^2;
%         H(3,3) = (La(2,2,3,2)+La(2,2,1,2)-2*La(2,2,2,2))/hvx^2;
%         H(4,4) = (La(2,2,2,3)+La(2,2,2,1)-2*La(2,2,2,2))/hvy^2;
%         
%         H(1,2) = (La(3,3,2,2)-La(3,1,2,2)-La(1,3,2,2)+La(1,1,2,2))/(4*hx*hy);  H(2,1)=H(1,2);
%         H(1,3) = (La(3,2,3,2)-La(3,2,1,2)-La(1,2,3,2)+La(1,2,1,2))/(4*hx*hvx); H(3,1)=H(1,3);
%         H(1,4) = (La(3,2,2,3)-La(3,2,2,1)-La(1,2,2,3)+La(1,2,2,1))/(4*hx*hvy); H(4,1)=H(1,4);
%         
%         H(2,3) = (La(2,3,3,2)-La(2,1,3,2)-La(2,3,1,2)+La(2,1,1,2))/(4*hy*hvx);  H(3,2)=H(2,3);
%         H(2,4) = (La(2,3,2,3)-La(2,1,2,3)-La(2,3,2,1)+La(2,1,2,1))/(4*hy*hvy);  H(4,2)=H(2,4);
%         
%         H(3,4) = (La(2,2,3,3)-La(2,2,1,3)-La(2,2,3,1)+La(2,2,1,1))/(4*hvx*hvy); H(4,3)=H(3,4);
%         
%     end
    