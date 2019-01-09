function L = ay_point_likelihood(Point,Mo,Mark,num_time_step)

%% Likelihood of a point
if nargin==3
    num_time_step = 1;
end

dT    = Mark.dT * num_time_step;
dTN   = (Mark.Path.T(end)-Mark.Path.T(1))/length(Mark.Path.T);
%%%----
%norm_st = sqrt(det(2*pi*Mark.Kernel.St));
norm_st = 1;

temp_1  = mvnpdf([Mark.Path.X Mark.Path.Y],[Point(1) Point(2)],Mark.Kernel.St);
temp_1(isnan(temp_1))=0;
den     = max(realmin , sum(temp_1)*norm_st);
               
% num
temp = zeros(length(Mark.Cell),1);
for m=1:length(Mark.Cell)
    if ~isempty(Mark.Cell{m}.X)
        % [m i j]
        temp_1 = mvnpdf([Mark.Cell{m}.X' Mark.Cell{m}.Y'],[Point(1) Point(2)],Mark.Kernel.Sm);
        temp_1(isnan(temp_1))= 0;
        num    = sum(temp_1);
        temp(m)= num/(den*dTN);
    else
        temp(m)= NaN;
    end
end
%% likelihood
L   = 1;
for m=1:length(Mark.Cell)
   if ~isnan(temp(m))
        L = (L.*(exp(-dT*temp(m)).*((dT*temp(m)).^Mo(m))))/factorial(Mo(m));
   end
end


