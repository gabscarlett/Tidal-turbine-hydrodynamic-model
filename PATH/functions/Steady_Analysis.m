function [P_s,T_s, Cl_S3d, Cd_S3d, CN_S3d, CT_S3d, AoA_s, MZs, MXs] = Steady_Analysis(U0,TSR,pitch,Blades,t,rho,omega,NBsec)


        % STEADY INDUCTION FACTORS
        load TGL_blade_profile
        r=linspace(rad(1),rad(end),NBsec);
        B=interp1(rad,B,r,'PCHIP')+pitch; c=interp1(rad,c,r,'PCHIP');
        
        load S814_static_data
        
        % PREPROCESSOR
        [~, Values_360r] = PreProcessor1(aoa,Cl_2d,Cd_2d,Cn_2d,Clin,LinRange,B,r,c,az);

        % STEADY INDUCTION FACTORS
        
        [a_s3d,ap_s3d] = BEM_steady3D(U0,TSR,Blades,r,c,B,Values_360r);
        % Remove end points
        a_s3d=a_s3d(2:end-1);
        ap_s3d=ap_s3d(2:end-1);
        r=r(2:end-1);c=c(2:end-1);B=B(2:end-1);
        Values_360r.Cl=Values_360r.Cl(:,2:end-1);Values_360r.Cd=Values_360r.Cd(:,2:end-1);
        

        % EQUIVALENT STEADY ANGLE OF ATTACK
        AoA_s=90-atan2d(r.*omega.*(1+ap_s3d'),U0.*(1-a_s3d'))-(rad2deg(B)).*ones(size(t'));
        AoA_s=AoA_s';

        % Steady (rotational) coefficients
        rr=r'.*ones(size(AoA_s));
        Cl_S3d=interp2(r,Values_360r.Alpha,Values_360r.Cl,rr,deg2rad(AoA_s),'spline'); 
        Cd_S3d=interp2(r,Values_360r.Alpha,Values_360r.Cd,rr,deg2rad(AoA_s),'spline');
        
        % Flow angle
        PHI=deg2rad(AoA_s)+B';
       
        % relative velocity
        Wrel_s = sqrt((U0.*(1-a_s3d')).^2+(r.*omega.*(1+ap_s3d')).^2).*ones(size(t'));
        
        % Dynamic pressure
        FF_s=0.5.*rho.*c'.*Wrel_s'.^2;
    
        % Moments, thrust, power
        [MZs,MXs,T_s,P_s]=loads(FF_s,Cl_S3d,Cd_S3d,PHI,r,omega); % STEADY
        
        % tangential and thrust force coefficients
        CN_S3d=Cl_S3d.*cos(PHI)+Cd_S3d.*sin(PHI);
        CT_S3d=Cl_S3d.*sin(PHI)-Cd_S3d.*cos(PHI);
end