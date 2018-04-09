function [Alpha, Cl_360, Cd_360, Cn_360] = VitExt(aoa,Cl_nd,Cd_nd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Gabriel Scarlett August 2017 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VITERNA EXTRAPOLATION METHOD

% Function to extend lift and drag coefficients through full 360 degree

% Inputs
% aoa - angle of attack
% Cl_nd 2d or 2d Cl values in given - angle of attack range
% Cd_nd 2d or 2d Cd values in given - angle of attack range

% Outputs
% Alpha - in range -pi to pi
% Cl_360 2d or 2d Cl values in Alpha range
% Cd_360 2d or 2d Cd values in Alpha range

%% %%%%%%%%%%%%%%%%% VITERNA EXTRAPOLATION INTERNAL FUNCTION %%%%%%%%%%%%%%

AR=10;
CdMax=1.11+0.018*AR;
A1=0.5*CdMax;
B1=CdMax;

alpha_high=aoa(end);
Cl_high=Cl_nd(end);
Cd_high=Cd_nd(end);

A2=(Cl_high -CdMax*sin(alpha_high)*cos(alpha_high))*(sin(alpha_high)/cos(alpha_high)^2);
B2=Cd_high - CdMax*sin(alpha_high)^2/cos(alpha_high);


funCl= @(ALPHA)A1*sin(2*ALPHA)+A2*cos(ALPHA).^2./sin(ALPHA);
funCd= @(ALPHA)B1.*sin(ALPHA).^2+B2.*cos(ALPHA);
Cl_adj=0.7;

%%                                                      TOWARDS +180 DEGREE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% aoa high <-> 90 degree %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=20;

Alpha_h1=linspace(alpha_high,pi/2,n);
Alpha_h1(1)=[];
Cl_360_h1 = funCl(Alpha_h1);
Cd_360_h1=funCd(Alpha_h1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 90 degree <-> 180 degree - aoa high %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Alpha_h2=linspace(pi/2,pi-alpha_high,n);
Alpha_h2(1)=[];
Cl_360_h2=-Cl_adj*funCl(pi-Alpha_h2);
Cd_360_h2=funCd(pi-Alpha_h2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 180 degree - aoa high <-> 180 degree %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Alpha_h3=linspace(pi-alpha_high,pi,n);
Alpha_h3(1)=[];
Cl_360_h3=(Alpha_h3-pi)/alpha_high*Cl_high*Cl_adj; % linear variation
Cd_360_h3=funCd(pi-Alpha_h3);

%%                                                      TOWARDS -180 DEGREE

alpha_low=aoa(1);
Cl_low=Cl_nd(1);
Cd_low=Cd_nd(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -alpha_high <-> alpha_low %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if alpha_low <= -alpha_high
            Alpha_L1 = [];
            Cl_360_L1 = [];
            Cd_360_L1 = [];
           
            alpha_low_max = alpha_low;
else
            
            
            Alpha_L1 = linspace(-alpha_high,alpha_low,n);
            Alpha_L1(1)=[];
            Alpha_L1(end)=[];
            % Note: this is done slightly differently than AirfoilPrep for better continuity
            Cl_360_L1= -Cl_high*Cl_adj + (Alpha_L1+alpha_high)/(alpha_low+alpha_high)*(Cl_low+Cl_high*Cl_adj);
            Cd_360_L1 = Cd_low + (Alpha_L1-alpha_low)/(-alpha_high-alpha_low)*(Cd_high-Cd_low);
            alpha_low_max = -alpha_high;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -90 degree <-> -aoa high %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Alpha_L2 = linspace(-pi/2,alpha_low_max,n);
        Alpha_L2(1) =[];
        Cl_360_L2 = -Cl_adj*funCl(-Alpha_L2);
        Cd_360_L2 = funCd(-Alpha_L2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -180 + alpha_high  <-> -90 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        Alpha_L3 = linspace(-pi+alpha_high,-pi/2,n);
        Alpha_L3(1) =[];
        Cl_360_L3 = Cl_adj*funCl(Alpha_L3+pi);
        Cd_360_L3 = funCd(Alpha_L3+pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -180 degree <-> -180 + alpha_high %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        
             
        Alpha_L4 = linspace(-pi,-pi+alpha_high,n); 
        Cl_360_L4 = (Alpha_L4+pi)/alpha_high*Cl_high*Cl_adj;  % linear variation
        Cd_360_L4 = funCd(Alpha_L4+pi);
        
%%                                                      CONCATENATE

Alpha = horzcat(Alpha_L4,Alpha_L3,Alpha_L2,Alpha_L1,aoa,Alpha_h1,Alpha_h2,Alpha_h3);
Cl_360= horzcat(Cl_360_L4,Cl_360_L3,Cl_360_L2,Cl_360_L1,Cl_nd,Cl_360_h1,Cl_360_h2,Cl_360_h3);
Cd_360= horzcat(Cd_360_L4,Cd_360_L3,Cd_360_L2,Cd_360_L1,Cd_nd,Cd_360_h1,Cd_360_h2,Cd_360_h3);

% some logic to avoid negative drag
InRange=Cd_360>0;
Cd_360=Cd_360.*InRange + (1-InRange).*Cd_low;

% Final check to avoid duplicates as this GREATLY screws things up!
[~,ia] = unique(Alpha);   % find unique values
isun=false(size(Alpha));  % create logical array
isun(ia)=true;            % assign true to the unique values
% now use logical indexing to remove non unique elements
Alpha=Alpha(isun);
Cl_360=Cl_360(isun);
Cd_360=Cd_360(isun);

% normal force coefficient

Cn_360=Cl_360.*cos(Alpha)+Cd_360.*sin(Alpha); % normal coefficient

end