function [Cd_DS_3D] = UnstCD(Dvis,Cd_Ind,Values_360,Values_360r,aE,r,Cd0)
        
        R_InRange=(r>0.8*r(end));        % boolean (TRUE for outer section) 
        % DRAG 3D
        
        % complete unsteady drag coefficient
        
        % No rotational correction for outer sections to avoid over
        % prediction (drag reduces from 0.8R - R during rotation)
        Cd_2St=interp1(Values_360.Alpha,Values_360.Cd,aE(R_InRange,:),'spline');   
        rr=r(~R_InRange)'.*ones(size(aE(~R_InRange,:)));
        
        %%
        
        Cd_3St=interp2(r,Values_360r.Alpha,Values_360r.Cd,rr,aE(~R_InRange,:),'spline'); 
        Cd_St = vertcat(Cd_3St, Cd_2St);
        
        % complete unsteady drag coefficient
        Cd_DS_3D=Cd_Ind+Cd_St+(Cd_St-Cd0).*Dvis;
        
end