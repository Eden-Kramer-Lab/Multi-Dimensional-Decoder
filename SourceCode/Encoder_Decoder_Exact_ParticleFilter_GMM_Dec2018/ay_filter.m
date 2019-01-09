function Post =ay_filter(L,Prior,Xs,Ys,Trans,num_time_step)
%% run filter
% copy Prior to Post
Post     = Prior;
% Fill vectors
TransMap      = zeros(length(Xs)*length(Ys),2);
ReshapedPrior = zeros(length(Xs)*length(Ys),1);
ind = 1;
for i=1:length(Xs)
    for j=1:length(Ys)
        TransMap(ind,:)= (Trans.A^num_time_step) * [Xs(i);Ys(j)];
        ReshapedPrior(ind) = Prior(i,j);
        ind = ind + 1;
    end
end
T = eye(size(Trans.A));
for t=1:num_time_step-1
    T = T+Trans.A^t;
end
% run for each i and j
for i=1:length(Xs)
    for j=1:length(Ys)
        % x,y
        m_trans = mvnpdf([Xs(i) Ys(j)],TransMap,T*Trans.S);
        Post(i,j) = L(i,j)*sum(m_trans.*ReshapedPrior);
    end
end
Post = Post./sum(sum(Post));
