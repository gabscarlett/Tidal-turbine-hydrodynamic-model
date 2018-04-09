function [Values_360, Values_360r] = PreProcessor1(aoa,Cl_2d,Cd_2d,Cn_2d,Clin,LinRange,B,r,c,az)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Gabriel Scarlett August 2017 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % NOTE:
        % Non-rotational:   function of angle of attack only.
        % Rotational:       function of both angle of attack and blade location.

        % OUTPUTS
        % non-rotational static coefficients with deep stall applied.
        % non-rotational separation point.
        % rotational coefficients with deep stall applied.
        % Rotational separation point.
        % Angle of attack in {-pi,pi}

        % Deep stall  non-rotational values (-180 <-> 180)
        [Values_360.Alpha,Values_360.Cl, Values_360.Cd, Values_360.Cn] = VitExt(aoa',Cl_2d',Cd_2d');
        
        Values_360.F=sep_point(Values_360.Alpha,az,Values_360.Cn,Clin,LinRange); % Non-rotational separation point (-180 <-> 180)

        
        
        % Rotational flow augmentation  
        
        F=sep_point(aoa,az,Cn_2d,Clin,LinRange); % for rotational solution
        [Cl_3d,Cd_3d,~,~] = stall_delay(B,r,c,aoa,F,Cl_2d,Cd_2d,Clin,az);
        

        % Deep stall  rotational values (-180 <-> 180)
        for K =1:length(r)
        [Values_360r.Alpha, Values_360r.Cl(K,:), Values_360r.Cd(K,:),Values_360r.Cn(K,:)] = VitExt(aoa',Cl_3d(:,K)',Cd_3d(:,K)');  
        end
        
        
        Values_360r.F=sep_point(Values_360r.Alpha,az,Values_360r.Cn,Clin,LinRange); % Rotational separation point (-180 <-> 180)
        
        Values_360r.Alpha=Values_360r.Alpha';
        Values_360r.F=Values_360r.F';
        Values_360r.Cd=Values_360r.Cd';
        Values_360r.Cl=Values_360r.Cl';
        Values_360r.Cn=Values_360r.Cn';
  
end