A=load('Result.txt');
A=A(1:800,:);
ind_a = find(A(:,3));

ind_b = find(A(:,3)==0);

OVERALL=  [mean(A(:,6))   mean(A(:,7))  mean(A(:,8))  sqrt(mean(A(:,9).^2))   sqrt(mean(A(:,10).^2))  sqrt(mean(A(:,11).^2))    sqrt(mean(A(:,12).^2))   sqrt(mean(A(:,13).^2))  sqrt(mean(A(:,14).^2))]
NON_SPIKE =  [mean(A(ind_b,6))   mean(A(ind_b,7))  mean(A(ind_b,8))  sqrt(mean(A(ind_b,9).^2))   sqrt(mean(A(ind_b,10).^2))  sqrt(mean(A(ind_b,11).^2))    sqrt(mean(A(ind_b,12).^2))   sqrt(mean(A(ind_b,13).^2))  sqrt(mean(A(ind_b,14).^2))]
SPIKE =  [mean(A(ind_a,6))   mean(A(ind_a,7))  mean(A(ind_a,8))  sqrt(mean(A(ind_a,9).^2))   sqrt(mean(A(ind_a,10).^2))  sqrt(mean(A(ind_a,11).^2))    sqrt(mean(A(ind_a,12).^2))   sqrt(mean(A(ind_a,13).^2))  sqrt(mean(A(ind_a,14).^2))]