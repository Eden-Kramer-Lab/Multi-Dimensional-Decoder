function [L,ML] = ay_likelihood(Mo,Mark,num_time_step,mx,my)
if nargin==2
    num_time_step = 1;
end

dT  = Mark.dT * num_time_step;
dTN = (Mark.Path.T(end)-Mark.Path.T(1))/length(Mark.Path.T);

temp    = zeros(length(Mark.Cell),length(mx),length(my));
% norm_st = sqrt(det(2*pi*Mark.Kernel.St));
norm_st  = 1;

for i=1:length(mx)
     for j=1:length(my)
          % den
          temp_1  = mvnpdf([Mark.Path.X Mark.Path.Y],[mx(i) my(j)],Mark.Kernel.St);
          temp_1(isnan(temp_1))=0;
          den     = max(realmin,sum(temp_1)*norm_st);
          % num
          for m=1:length(Mark.Cell)
                if ~isempty(Mark.Cell{m}.X)
                    temp_1 = mvnpdf([Mark.Cell{m}.X' Mark.Cell{m}.Y'],[mx(i) my(j)],Mark.Kernel.Sm);
                    temp_1(isnan(temp_1))=0;
                    num    = sum(temp_1);
                    temp(m,i,j)= num/(den*dTN);
                else
                    temp(m,i,j)= NaN;
                end
          end
     end
end
  
%% likelihood
L   = ones(length(mx),length(my));
for m=1:length(Mark.Cell)
        temp_a = squeeze(temp(m,:,:));
        if ~isnan(sum(sum(temp_a)))
            L = (L.*(exp(-dT*temp_a).*((dT*temp_a).^Mo(m))))/factorial(Mo(m));
        end
end
ML  = temp;

