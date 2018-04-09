function [U_axial,U_theta,a,a_p,t,steps,Tr_steps]=InitialConditions(data,TSR,U0,ZTb,Tr,r,B,c,dt,Values_360r)
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Gabriel Scarlett August 2017 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THIS FUNCTION DETERMINES THE VELOCITIES AT EACH BLADE AS IT ROTATES USING
% THE REDAPT FLOW DATA.

% INITIAL CONDITIONS FOR THE INDUCTION FACTORS ARE DETERMINED USING THE
% QUASI-STEADY BEM OUTPUT WITH THE BLADE ONE FLOW FIELD

        % Flip values to pass to be BEM 
        Values_360r.Alpha=Values_360r.Alpha';
        Values_360r.F=Values_360r.F';
        Values_360r.Cd=Values_360r.Cd';
        Values_360r.Cl=Values_360r.Cl';
        Values_360r.Cn=Values_360r.Cn';

%%%%%%%%%%%%%%%%%%%%%%% LOAD REDAPT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         U1 = double(data.adcp1.U);                     % streamwise ReDAPT velocity
         W1 = double(data.adcp1.W);                     % depthwise ReDAPT velocity
         depth=size(U1,2)-1;                            % Depth
         z0=-depth-1+ZTb;                               % hub position
         TIME=2*length(U1)-2;                           % time of sample
         t_ReD=linspace(0,TIME,length(U1(:,1)));        % time index
         d_ReD=linspace(-depth,0,length(U1(1,:)));      % space (in z) index
         Blades=3;

         Tr_steps=floor(TIME/Tr);                       % number of rotations in the time sample
         T_end=Tr_steps*Tr;                             % total time in terms for rotations

         steps=round(Tr/dt);                            % steps per iteration
         total_steps=Tr_steps*steps;                    % total steps in sample

         t=linspace(0,T_end,total_steps);               % time array
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%   Time step blade rotation, determine all velocity components and induction factors   %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
         
%          % Preallocation
             U_axial_1 = zeros(length(r), length(t));
             U_theta_1 = zeros(length(r), length(t));
             U_axial_2 = zeros(length(r), length(t));
             U_theta_2 = zeros(length(r), length(t));
             U_axial_3 = zeros(length(r), length(t));
             U_theta_3 = zeros(length(r), length(t));
             
             a = zeros(length(r), length(t));
             a_p = zeros(length(r), length(t));

             omega = abs(U0*TSR/r(end));      % rotational speed of blades (rad/s)         
             psi=omega.*t;                    % azimuthal position over one rotation
         
             % PHASE LAG
             p1=0;
             p2=deg2rad(120);
             p3=deg2rad(240);
 
 
    parfor i=1:length(psi)
        
        %%%%%%%%%%%%%%%%%
        %%%% BLADE_1 %%%% 
        %%%%%%%%%%%%%%%%%
        
        z_blade_1 = r'.*sin(psi(i)-p1);
        z_1 = z0 + z_blade_1;                       % depth of blade
        
        % look up redapt data   
        ut1=interp2(d_ReD,t_ReD,U1,z_1',t(i),'cubic'); % ADCP 1 streamwise velocity
        wt1=interp2(d_ReD,t_ReD,W1,z_1',t(i),'cubic'); % ADCP 1 depthwise velocity
        
        U_axial_1(:,i)=ut1;                           % Axial velocity component
        U_theta_1(:,i) = omega.*r +wt1.*cos(psi(i)-p1);  % Tangential velocity
        
        TSR_unsteady_1=U_theta_1(:,i)./U_axial_1(:,i);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % induction factors calculation based on blade 1 only %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [a(:,i),a_p(:,i)] = BEM_TIME_UNSTEADY(U0,TSR,TSR_unsteady_1,Blades,r,c,B,Values_360r);
        % SET END CONDITION
        
        
        %%%%%%%%%%%%%%%%%
        %%%% BLADE_2 %%%% 
        %%%%%%%%%%%%%%%%%
        
        z_blade_2 = r'.*sin(psi(i)-p2);
        z_2 = z0 + z_blade_2;                       % depth of blade
        
        % look up value from redapt data   
        ut2=interp2(d_ReD,t_ReD,U1,z_2',t(i),'cubic'); % ADCP 1 streamwise velocity
        wt2=interp2(d_ReD,t_ReD,W1,z_2',t(i),'cubic'); % ADCP 1 depthwise velocity
        
        U_axial_2(:,i)=ut2;                              % Axial velocity component
        U_theta_2(:,i) = omega.*r +wt2.*cos(psi(i)-p2);  % Tangential velocity

        %%%%%%%%%%%%%%%%%
        %%%% BLADE_3 %%%%
        %%%%%%%%%%%%%%%%%
        
        z_blade_3 = r'.*sin(psi(i)-p3);
        z_3 = z0 + z_blade_3;                       % depth of blade
        
        % look up redapt data   
        ut3=interp2(d_ReD,t_ReD,U1,z_3',t(i),'cubic'); % ADCP 1 streamwise velocity
        wt3=interp2(d_ReD,t_ReD,W1,z_3',t(i),'cubic'); % ADCP 1 depthwise velocity
        
        U_axial_3(:,i)=ut3;                              % Axial velocity component
        U_theta_3(:,i) = omega.*r +wt3.*cos(psi(i)-p3);  % Tangential velocity
        
        
        i
    end
        U_axial(:,:,1)=U_axial_1;
        U_axial(:,:,2)=U_axial_2;
        U_axial(:,:,3)=U_axial_3;

        U_theta(:,:,1)=U_theta_1;
        U_theta(:,:,2)=U_theta_2;
        U_theta(:,:,3)=U_theta_3;

        % Set end conditions
        a(1,:)=0.990;
        a(100,:)=0.990;
        a_p(1,:)=0.001;
        a_p(100,:)=0.001;

    end