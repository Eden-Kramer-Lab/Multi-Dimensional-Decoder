function [Grd,Hs] = ay_grd_hessian(Point,Mark,steps)
if length(Point)==2
    %% we assume denominator is firly flat
    Grd = realmin*ones(2,1);
    Hs  = realmin*eye(2,2);
    hx  = steps(1);
    hy  = steps(2);

    La = zeros(3,3);
    mo = zeros(length(Mark.Cell),1);
    for ix=-1:1
       for iy=-1:1
          pp_e = [Point(1)+ix*hx  Point(2)+iy*hy];
          La(ix+2,iy+2) = -log(max(realmin,ay_point_likelihood(pp_e,mo,Mark,1)));
       end
    end

    Grd(1,1) = (La(3,2)-La(1,2))/(2*hx);
    Grd(2,1) = (La(2,3)-La(2,1))/(2*hy);

    % Calculate Hsessian around pp
    Hs(1,1) = (La(3,2)+La(1,2)-2*La(2,2))/hx^2;
    Hs(2,2) = (La(2,3)+La(2,1)-2*La(2,2))/hy^2;

    Hs(1,2) = (La(3,3)-La(3,1)-La(1,3)+La(1,1))/(4*hx*hy);
    Hs(2,1) = Hs(1,2);
else
        %% we assume denominator is firly flat
    Grd = realmin*ones(4,1);
    Hs  = realmin*ones(4,4);
    hx  = steps(1);
    hy  = steps(2);
    hvx = steps(3);
    hvy = steps(4);

    La = zeros(3,3,3,3);
    mo = zeros(length(Mark.Cell),1);
    for ix=-1:1
       for iy=-1:1
          for ivx=-1:1
             for ivy=-1:1
                 pp_e= [Point(1)+ix*hx  Point(2)+iy*hy  Point(3)+ivx*hvx  Point(4)+ivy*hvy];
                 La(ix+2,iy+2,ivx+2,ivy+2) = -log(max(realmin,ay_point_likelihood(pp_e,mo,Mark,1)));
             end
          end
       end
    end


    Grd(1,1) = (La(3,2,2,2)-La(1,2,2,2))/(2*hx);
    Grd(2,1) = (La(2,3,2,2)-La(2,1,2,2))/(2*hy);
    Grd(3,1) = (La(2,2,3,2)-La(2,2,1,2))/(2*hvx);
    Grd(4,1) = (La(2,2,2,3)-La(2,2,2,1))/(2*hvy);

    % Calculate Hsessian around pp
    Hs(1,1) = (La(3,2,2,2)+La(1,2,2,2)-2*La(2,2,2,2))/hx^2;
    Hs(2,2) = (La(2,3,2,2)+La(2,1,2,2)-2*La(2,2,2,2))/hy^2;
    Hs(3,3) = (La(2,2,3,2)+La(2,2,1,2)-2*La(2,2,2,2))/hvx^2;
    Hs(4,4) = (La(2,2,2,3)+La(2,2,2,1)-2*La(2,2,2,2))/hvy^2;

    Hs(1,2) = (La(3,3,2,2)-La(3,1,2,2)-La(1,3,2,2)+La(1,1,2,2))/(4*hx*hy);
    Hs(2,1) = Hs(1,2);

    Hs(1,3) = (La(3,2,3,2)-La(3,2,1,2)-La(1,2,3,2)+La(1,2,1,2))/(4*hx*hvx);
    Hs(3,1)=Hs(1,3);

    Hs(1,4) = (La(3,2,2,3)-La(3,2,2,1)-La(1,2,2,3)+La(1,2,2,1))/(4*hx*hvy);
    Hs(4,1)=Hs(1,4);

    Hs(2,3) = (La(2,3,3,2)-La(2,1,3,2)-La(2,3,1,2)+La(2,1,1,2))/(4*hy*hvx);
    Hs(3,2)=Hs(2,3);

    Hs(2,4) = (La(2,3,2,3)-La(2,1,2,3)-La(2,3,2,1)+La(2,1,2,1))/(4*hy*hvy);
    Hs(4,2)=Hs(2,4);

    Hs(3,4) = (La(2,2,3,3)-La(2,2,1,3)-La(2,2,3,1)+La(2,2,1,1))/(4*hvx*hvy);
    Hs(4,3)=Hs(3,4);

   

end
   
