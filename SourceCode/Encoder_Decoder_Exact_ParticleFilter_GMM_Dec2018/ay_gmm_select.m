function Param = ay_gmm_select(mode,MaxMix,IterA,IterB,Ps,Method)
    if mode==1
        % Run with maximum number
        Param = ay_gmm_e(Ps,MaxMix,IterA);
        if Param.n_mix>1
            L = ay_gmm_likelihood(1,Param,Ps);
            loop = 1;
            while loop
               if  Param.n_mix > 1 
                  % pick pairs
                  [ind_a,ind_b]= ay_gmm_pick_pair(Param,Ps);
                  % merge
                  tParam = ay_gmm_merge(2,Param,Ps,IterB,[ind_a ind_b]);
                  % check likelihood
                  tL    = ay_gmm_likelihood(Method,tParam,Ps);
                  if tL < L 
                      L = tL;
                      Param = tParam;
                  else
                      loop = 0;
                  end
               else
                  loop = 0;
               end
            end
        end
    else
        % Run with maximum number
        Param = ay_gmm_e(Ps,Pw,1,IterA);
        L     = ay_gmm_likelihood(1,Param,Ps,Pw);
        loop  = 1;
        while loop
           if Param.n_mix < MaxMix 
              % merge
              tParam = ay_gmm_merge(1,Param,Ps,IterB);
              % check likelihood
              tL     = ay_gmm_likelihood(Method,tParam,Ps);
              if tL < L 
                  L = tL;
                  Param = tParam;
              else
                  loop = 0;
              end
           else
              loop = 0;
           end
        end
    end
end    