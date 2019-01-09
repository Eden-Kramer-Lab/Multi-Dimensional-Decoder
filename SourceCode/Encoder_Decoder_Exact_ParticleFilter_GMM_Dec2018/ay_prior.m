function cur_prior = ay_prior(prev_post,trans,num_time_step)
if nargin == 2
    num_time_step = 1;
end
%% one-step update given a mixture model
% note, just mean and covariance changes per mixture
cur_prior = prev_post;
for i=1:cur_prior.n_mix
    Mt = prev_post.Model{i}.M;
    St = prev_post.Model{i}.S;
    for n=1:num_time_step
        % update mean
        Mt = (trans.A * Mt')';
        % update covariance
        St = trans.A * St * trans.A' + trans.S;
    end
    cur_prior.Model{i}.M = Mt;
    cur_prior.Model{i}.S = St;
end



