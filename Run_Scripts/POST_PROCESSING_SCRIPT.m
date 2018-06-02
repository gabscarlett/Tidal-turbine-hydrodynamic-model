
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                 %%%
%%%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%
%%%       %%%%%% Written by Gabriel Scarlett August 2017 %%%%%%     %%%
%%%       %%%%%%     The University of Edinburgh, UK     %%%%%%     %%%
%%%       %%%%%%         School of Engineering           %%%%%%     %%%
%%%       %%%%%%       Institute for Energy Systems      %%%%%%     %%%
%%%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%
%%%                                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ASSIGN PATHS TO FUNCTIONS AND DATA: WORK PC
% path(path,genpath('\Users\s1040865\Dropbox\PhD\Modelling\Programs\Matlab\Working_Folder_April_2018\PATH\functions')); % Functions
% path(path,genpath('\Users\s1040865\Dropbox\PhD\Modelling\Programs\Matlab\Working_Folder_April_2018\PATH\data')); % Input data
% path(path,genpath('\Users\s1040865\Dropbox\PhD\Modelling\Programs\Matlab\Working_Folder_April_2018\PATH\Saved_Simulation_data')); % Saved simulation data

% ASSIGN PATHS TO FUNCTIONS AND DATA: HOME W520
path(path,genpath('\Users\gabsc\Dropbox\PhD\Modelling\Programs\Matlab\Working_Folder_April_2018\PATH\functions')); % Functions
path(path,genpath('\Users\gabsc\Dropbox\PhD\Modelling\Programs\Matlab\Working_Folder_April_2018\PATH\data')); % Input data
path(path,genpath('\Users\gabsc\Dropbox\PhD\Modelling\Programs\Matlab\Working_Folder_April_2018\PATH\Saved_Simulation_data')); % Saved simulation data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST PROCESSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ITERATES OVER 3 BLADES
% % USES MEASURED ONSET FLOW DATA FROM REDAPT AND INDUCTION FACTORS FROM COUPLED MODEL


% DETERMINED BENDING MOMENTS, POWER AND THRUST COEFFICIENTS FOR THE ROTOR
% DETERMINED LOCAL MEAN SECTIONAL LOAD, MOMENT AND TORQUE COEFFICIENTS
% COMPARES UNSTEADY, QUASI-STEADY AND STEADY PREDICTION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INITIALISATION CALLS
% FOLLOWING CALLS TO PRE-PROCESS DATA
%
% 1: CALL TO PreProcessor1       - Extrapolate static aerofoil data
%                                  to deep stall region.
%                                - Apply rotational augmentation
%                                  correction
%                                - Compute the point of trailing edge separation
%                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLADE LOADS
% FOLLOWING CALLS MADE 
%
% 2:  CALL TO wag                 - Determine unsteady linear Cl for each blade(2D array)
% 3:  CALL TO DS_2D               - Determine 2D unsteady non-linear Cl and Cd  for each blade (2D array)
% 4:  CALL TO DS_3D               - Determine 3D unsteady non-linear Cl and for each blade (2D array)
% 5:  CALL TO UnstCD              - Determine 3D unsteady non-linear Cd and for each blade (2D array)
% 6:  CALL TO loads               - Determine the power, thrust, and moments for each blade (UNSTEADY)
% 7:  CALL TO loads               - Determine the power, thrust, and moments for each blade (QUASI-STEADY)
% 8:  CALL TO Steady              - Determine steady state components 
% 9:  CALL TO GRAPHS              - Produced graphs for Renewable Energy paper 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TURBINE SPECIFICATIONS AND OPERATING CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

graph_settings
clear, clc, close all

%%%%%%%%%%%%%%%%%%%%%% Operating conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TSR=3.5;                   % tip speed ratio
% file_sim = 'ReDAPT_Unsteady_TSR_3p5';
% Pitch = 1.2;               % pitch angle (deg) 

% TSR=4.0;                   % tip speed ratio
% file_sim = 'ReDAPT_Unsteady_TSR_4';
% Pitch = 0.2;              % pitch angle (deg) 

TSR=4.5;                   % tip speed ratio
file_sim = 'ReDAPT_Unsteady_TSR_4p5';
Pitch = 0.1;               % pitch angle (deg) 

load(file_sim)

U0=2.77;                % streamwise current (m/s)
% Hh:                   % distance from bed to hub centre

g=9.81;                 % acceleration due to gravity (m/s^2)
rho = 1025;             % density of sea water (kg/m^3)


%%%%%%%%%%%%%%%%%%%% Turbine specifications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_turb ='TGL_TURBINE';
load(file_turb)

% THIS IS THE FULL SCALE TGL BLADE GIVEN BY GRETTON
% Blades:               % number of blades
pitch=deg2rad(Pitch);   % operational pitch applied (degree)

% Discretise blade
NBsec = 100; % number of blade sections
r=linspace(rad(1),rad(end),NBsec);
r=r(2:end-1); % remove end values as there is no flow there due to tip losses
R_InRange=(r>0.8*r(end));        % boolean (TRUE for outer section) 
% Interpolate twist and chord
B=interp1(rad,B,r,'PCHIP')+pitch; c=interp1(rad,c,r,'PCHIP');

R=(rad(end));               % radius of blade
mu=r./R;                    % normalised radial position
A=pi*R^2;                   % swept area (m^2)
Aeff=pi*(R^2-r(1)^2);       % effective area (m^2)
omega = abs(U0*TSR/R);      % rotational speed of blades (rad/s)
Tr = (2*pi)/omega;          % period of rotation (s)

%%%%%%%%%%%%%%%%%%%%%%% LOAD REDAPT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         file_flow = 'ReDAPT_FlowSample';
         load(file_flow)


         U1 = double(data.adcp1.U);  % streamwise ReDAPT velocity
         W1 = double(data.adcp1.W);  % depthwise ReDAPT velocity
         depth=size(U1,2)-1;         % water depth 
         z0=-depth-1+Hh;            % hub depth (m)
         TIME=2*length(U1)-2;
         t_ReD=linspace(0,TIME,length(U1(:,1))); % time index
         d_ReD=linspace(-depth,0,length(U1(1,:))); % space (in z) index
         
         % Temporal data
         
         Tr_steps=round(TIME/Tr);           % number of rotations in the time sample
         T_end=Tr_steps*Tr;                 % total time in terms for rotations
         dt=Tr/72;
         steps=round(Tr/dt);                % steps per iteration
         total_steps=Tr_steps*steps;        % total steps in sample

         t=linspace(0,T_end,total_steps);   % time array
         omega = abs(U0*TSR/R);             % rotational speed of blades (rad/s)         
         psi=omega.*t;                      % azimuthal position over one rotation
         
         % PHASE LAG
         phase=[0,deg2rad(120),deg2rad(240)];
        
        
        
        %% PREPROCESSOR
        % ROTATONAL AUGMENTATION / STALL DELAY / SEPARATION POINT
        
        file_foil = 'S814_static_data';
        load(file_foil)

        [Values_360, Values_360r] = PreProcessor1(aoa,Cl_2d,Cd_2d,Cn_2d,Clin,LinRange,B,r,c,az);

        file_ds ='S814_DS_parameters';
         %% LOAD UNSTEADY INDUCTION FACTORS DETERMINED BY COUPLED MODEL
       
        Axial_induction(1,:)=[];
        Axial_induction(end,:)=[];
        
        Tangential_induction(1,:)=[];
        Tangential_induction(end,:)=[];

        % Add steps for first rotation to align array size
        Lr=length(r);
       
        a=[ones(Lr,steps).*Axial_induction(:,1),Axial_induction]; % concatenate first rotation
        ap=[ones(Lr,steps).*Tangential_induction(:,1),Tangential_induction]; % concatenate first rotation
        
            if length(psi) > length(a) % then need to concatenate the last rotation, too
                
                    a=[a,ones(Lr,steps).*Axial_induction(:,end)]; % concatenate last rotation
                    ap=[ap,ones(Lr,steps).*Tangential_induction(:,end)]; % concatenate last rotation
            end


for i=1:Blades
        
        %% LOOPS OVER EACH EACH BLADE
        
        i
        z_blade = r'.*sin(psi-phase(i));
        z = z0 + z_blade;                       % depth of blade
        
        tt=t.*ones(size(z));
        % look up redapt data   
        ut1=interp2(d_ReD,t_ReD,U1,z,tt,'cubic');      % ADCP 1 streamwise velocity
        wt1=interp2(d_ReD,t_ReD,W1,z,tt,'cubic');      % ADCP 1 depthwise velocity
        
        U_axial=ut1;                                   % Axial velocity component
        U_theta = omega*r' +wt1.*cos(psi-phase(i));    % Tangential velocity
        
        % Relative velocity
        Wrel(:,:,i) = sqrt((U_axial.*(1-a)).^2+(U_theta.*(1+ap)).^2);
        
        % Angle of attack
        AoA(:,:,i) = 90-atan2d(U_theta.*(1+ap),U_axial.*(1-a))-(rad2deg(B)');
        
        % DEAL WITH NANS
        AoA(:,:,i) = fillmissing(AoA(:,:,i),'linear',2,'EndValues','nearest');
        Wrel(:,:,i) = fillmissing(Wrel(:,:,i),'linear',2,'EndValues','nearest');
        
        
        % CALL INDICIAL SOLUTION
        % pass AoA history to indicial load model
        [Cl_us(:,:,i),Cl_c,Cl_nc,Ds,aE(:,:,i)] = wag(c,dt,U0,deg2rad(AoA(:,:,i)),az,Clin);
        
        % CALL DYNAMIC STALL SOLUTION
        % NON-ROTATIONAL SOLUTION
        [~,~,Cl_DS_2D(:,:,i), Dvis_2d, Cd_Ind_2d, ff_2d(:,:,i), fff_2d(:,:,i), VortexTracker_2d(:,:,1)] =.....
            DS_2D(Values_360,Cl_us(:,:,i),Cl_c,Cl_nc,Ds,aE(:,:,i),deg2rad(AoA(:,:,i)),file_ds);
        
        % DRAG NON-ROTATIONAL
        
        %Cd_2dSt=interp1(aoa,Cd_2d,aE(:,:,i),'spline');
        Cd_2dSt=interp1(Values_360.Alpha,Values_360.Cd,aE(:,:,i),'spline');
        Cd_DS_2D(:,:,i)=Cd_Ind_2d+Cd_2dSt+(Cd_2dSt-Cd0).*Dvis_2d; 
                    
        % CALL DYNAMIC STALL SOLUTION
        % ROTATIONAL AUGMENTATION SOLUTION
        [~,~,Cl_DS_3D(:,:,i), Dvis, Cd_Ind, ff_3d(:,:,i), fff_3d(:,:,i), VortexTracker_3d(:,:,1)] =.....
            DS_3D(B,c,Values_360r,r,Cl_us(:,:,i),Cl_c,Cl_nc,Ds,aE(:,:,i),deg2rad(AoA(:,:,i)),file_ds);        
                
        % DRAG ROTATIONAL 
        [Cd_DS_3D(:,:,i)] = UnstCD(Dvis,Cd_Ind,Values_360,Values_360r,aE(:,:,i),r,Cd0);
    
        
        % Quasi-steady (2D) coefficients
        % Interpotation with the 360 degree values
        Cl_QS2d(:,:,i)=interp1(Values_360r.Alpha,Values_360.Cl,deg2rad(AoA(:,:,i)),'spline');
        Cd_QS2d(:,:,i)=interp1(Values_360r.Alpha,Values_360.Cd,deg2rad(AoA(:,:,i)),'spline');
        
       
        
        % Quasi-steady (3D) coefficients
        % Interpotation with the rotatonal 360 degree values
        rr=r'.*ones(size(AoA(:,:,i)));
        Cl_QS3d(:,:,i)=interp2(r,Values_360r.Alpha,Values_360r.Cl,rr,deg2rad(AoA(:,:,i)),'spline');
        Cd_QS3d(:,:,i)=interp2(r,Values_360r.Alpha,Values_360r.Cd,rr,deg2rad(AoA(:,:,i)),'spline');
        Cd_QS3d(R_InRange,:,i)=Cd_QS2d(R_InRange,:,i); % VALUES NEAR THE TIP = 2D VALUE
        
        
        % Flow angle
        phi=deg2rad(AoA(:,:,i))+B';
        
        % Dynamic pressure
        FF=0.5.*rho.*c'.*Wrel(:,:,i).^2;
        
        % Moments, thrust, power
        [MZ(:,i),MX(:,i),T(:,i),P(:,i)]=loads(FF,Cl_DS_3D(:,:,i),Cd_QS3d(:,:,i),phi,r,omega); % UNSTEADY
        
        [MZqs(:,i),MXqs(:,i),Tqs(:,i),Pqs(:,i)]=loads(FF,Cl_QS3d(:,:,i),Cd_QS3d(:,:,i),phi,r,omega); % Quasi-steady


        % NORMAL AND TANGENTIAL FORCES TO THE ROTOR PLANE
        
        CN_Ds3d(:,:,i)=Cl_DS_3D(:,:,i).*cos(phi)+Cd_DS_3D(:,:,i).*sin(phi);
        CT_Ds3d(:,:,i)=Cl_DS_3D(:,:,i).*sin(phi)-Cd_DS_3D(:,:,i).*cos(phi);
        
        CN_Ds2d(:,:,i)=Cl_DS_2D(:,:,i).*cos(phi)+Cd_DS_2D(:,:,i).*sin(phi);
        CT_Ds2d(:,:,i)=Cl_DS_2D(:,:,i).*sin(phi)-Cd_DS_2D(:,:,i).*cos(phi);
        
        CN_Qs3d(:,:,i)=Cl_QS3d(:,:,i).*cos(phi)+Cd_QS3d(:,:,i).*sin(phi);
        CT_Qs3d(:,:,i)=Cl_QS3d(:,:,i).*sin(phi)-Cd_QS3d(:,:,i).*cos(phi);
        
        CN_Qs2d(:,:,i)=Cl_QS2d(:,:,i).*cos(phi)+Cd_QS2d(:,:,i).*sin(phi);
        CT_Qs2d(:,:,i)=Cl_QS2d(:,:,i).*sin(phi)-Cd_QS2d(:,:,i).*cos(phi);

end
        

        %% STEADY ANALYSIS
        
        [P_s, T_s, Cl_S3d, Cd_S3d, CN_S3d, CT_S3d, AoA_s, MZs, MXs, f] = Steady(U0,TSR,pitch,t,rho,omega,NBsec,file_turb,file_foil);


        %% SEPERATION POINT

        % Quasi-steady (non-rotational) separation point
        f_2d=interp1(Values_360.Alpha,Values_360.F,deg2rad(AoA),'spline');
        
        InRange1=f_2d>1;
        f_2d(InRange1)=0.99;
        
        % Quasi-steady (rotational) separation point
        rrr=r'.*ones(size(AoA));
        f_3d=interp2(r,Values_360r.Alpha,Values_360r.F,rrr,deg2rad(AoA),'spline');
        
        InRange2=f_3d>0.99;
        f_3d(InRange2)=1;
        
       %% CALL GRAPH FUNCTION TO MAKE PLOTS 
       % FOR (RENEWABLE ENERGY JOURNAL 2018)
        
        
         GRAPHS(r,R,t,Tr,A,U0,P,P_s,T,T_s,Cl_DS_3D,Cl_QS3d,Cd_DS_3D,Cd_QS3d,CN_Ds3d,....
             CN_Qs3d,CT_Ds3d,CT_Qs3d,Cl_S3d,Cd_S3d,CN_S3d,CT_S3d,ff_3d,f_3d,f,aoa,AoA,AoA_s,Values_360r,MZ,MZs,MZqs,MX,MXs,MXqs)

