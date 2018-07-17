% Written by G.T. Scarlett May 2018

% Calculate CP and CT for a range of tip-speed ratios
% Plot CP and CT against tip-speed ratio
% Compare with Aerodyn Predication


clear
clc
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% OPERATION PARAMATERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=9.81;                 % acceleration due to gravity (m/s^2)
rho = 1025;             % density of sea water (kg/m^3)
U0=2.77;                % streamwise current (m/s)
Blades=3;               % number of blades

Pitch=0.1;              % operational pitch angle (deg)
pitch=deg2rad(Pitch);   % operational pitch angle (rad)

%%%%%%%%%%%%%%%%%%%% Turbine specifications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load TGL_BLADE_PROFILE % THIS IS THE FULL SCALE TGL BLADE GIVEN BY GRETTON

load TGL_TURBINE % THIS IS THE FULL SCALE TGL BLADE GIVEN BY GRETTON

NBsec = 100; % number of blade sections

r=linspace(rad(1),rad(end),NBsec);
% cubic spline for blade twist and chord
B=interp1(rad,B,r,'PCHIP')+pitch;
c=interp1(rad,c,r,'PCHIP');

R=(r(end));                     % radius of blade
A=pi*(R^2-r(1)^2);              % swept area
dr=r(2)-r(1);                   % increment size

% loop TSR
TSR=0.5:0.25:8;               % tip speed ratio


load S814_360_BEM_data % AEROFOIL DATA CORRECTED EXTRAPOLATED FOR DEEP STALL


parfor i=1:length(TSR)
   
omega = abs(U0*TSR(i)/R);      % rotational speed of blades (rad/s)



%% CALL BEM
%%%%%%%%%%%%%%%%%%%%%%%% determine induction factors %%%%%%%%%%%%%%%%%%%%%%

% a:    axial induction factor zero yaw
% ap:   tangential induction factor
% Ct:   tangential force coefficient
% Cn:   normal force coefficient


% CALL BEM FUNCTION
[a,ap,Ct,Cn] = BEM_2D(U0,TSR(i),Blades,r,c,B,Alpha,CL,CD);

W=sqrt((U0.*(1-a(1:end-1))).^2 + (r(1:end-1).*omega.*(1+ap(1:end-1))).^2);
q=0.5*rho*W.^2;
Y=Ct.*q.*c(1:end-1);
X=Cn.*q.*c(1:end-1);
T=Y.*r(1:end-1);

%%%%%  SUM TORQUE AND THRUST AS AVERAGE BETWEEN 2 ELEMENTS

Q=0;
Thr=0;

for ii=1:length(T)-1
    Q=Q+0.5*(T(ii)+T(ii+1)).*dr;
    Thr=Thr+0.5*(X(ii)+X(ii+1)).*dr;
end

% CP
Pout(i)=3*Q*omega;
Pin=0.5*rho*A*U0.^3;
CP(i)=Pout(i)/Pin;
% CT
Tout(i)=3*Thr;
Tin=0.5*rho*A*U0.^2;
CT(i)=Tout(i)/Tin;

i
end


%% PLOT RESULTS WITH AERODYN PRESICTION

% CALL GRAPH SETTINGS
graph_settings

% DRAW VERTICAL LINE SHOWING CP MAX
y1=linspace(0,max(CP),10);
x1=ones(size(y1))*4.5;
y2=linspace(0,max(CT),10);
x2=ones(size(y1))*4.5;

% LOAD AERODYN SIMULATION RESULTS
load('ADyn_CP_CT')

% PLOT
figure(1)
ax1=subplot(2,1,1); % CP
plot(ax1,TSR,CP,'ko',TSR_Adyn,CP_Adyn,'k:',x1,y1,'r--')
ylabel('$$C_P$$')
xlabel('$$\lambda$$')
grid off
axis([0 8 0 0.6])


ax2=subplot(2,1,2); % CT
plot(ax2,TSR,CT,'ko',TSR_Adyn,CT_Adyn,'k:',x2,y2,'r--')
ylabel('$$C_T$$')
xlabel('$$\lambda$$')
legend('Implementation','AeroDyn','Location','best')
legend boxoff
grid off
axis([0 8 0 1])

