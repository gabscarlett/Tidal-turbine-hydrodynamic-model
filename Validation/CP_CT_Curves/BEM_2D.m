function [a_m,a_pm,CT_m,CN_m] = BEM_2D(U0,TSR,Nb,r,c,B,Alpha,Cl,Cd)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Written by Gabriel Scarlett May 2017   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BEM IMPLEMENTATION OF NING 2014, WIND ENERGY
% CONVERGES ON A SINGLE PARAMETER (PHASE ANGLE - PHI)
% TO SOLVE THE RESIDUAL EQUATION, CALLS FUNCTION ZERO (BRENT ALGORITHM)

% THIS VERSION LOOPS THROUGH BLADE SECTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEADY BEM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OUTPUTS INDUCTION FACTORS USING STATIC DATA

% INPUTS

% input: U0           % streamwise current (m/s)
% input: TSR          % tip speed ratio
% input: Nb           % number of blades
% input: r            % radius (m)
% input: c            % chord (m)
% input: B            % twist (rad)
% input: alpha        % LOOK-UP: angle of attack (rad)
% input: Cl           % LOOK-UP: lift coefficient
% input: Cd           % LOOK-UP: drag coefficient

% input: Values_360       % Cl, Cd, alpha look-up table for each blade section

Rt = r(end);                 % Turbine radius at tip
Rh=0.75*r(1);                % Radius of the hub
omega = abs(U0*TSR/Rt);      % rotational speed of blades (rad/s)

% Preallocation
a_m=ones(1,length(r));
a_pm=ones(1,length(r));

syms PHI

for i=1:length(r)-1
    
%               local blade values

                TSr=omega*r(i)/U0;       % local tip speed ratio
                SLr=Nb*c(i)/(2*pi*r(i));  % local solidity
            
%               calculate initial inflow angle
                
                a1=2+pi*TSr*SLr;
                a2=4-4*pi*TSr*SLr;
                a3=pi*TSr^2 *SLr*(8*B(i)+pi*SLr);
                
                a=0.25*(a1-sqrt(abs(a2+a3))); % initial guess Moriarty
                ap=0;
                phi=atan((1-a)/(TSr*(1+ap)));  % inflow angle
                
                
            % initial values
            Err=1;
            Ep=1E-9;

            j=0;
            
while (max(Err) > Ep)
  j=1+j;
                aoa=phi-B(i);

%               Look up Cl and Cd with calculated aoa

                CL=interp1(Alpha,Cl,aoa,'spline'); 
                CD=interp1(Alpha,Cd,aoa,'spline');

%               calculate body forces
 
                CN=CL*cos(phi)+CD*sin(phi);
                CT=CL*sin(phi)-CD*cos(phi);
                
%               tip and hub losses
 
                ftip=0.5*Nb*(Rt-r(i))/(r(i)*abs(sin(phi)));
                Ftip=(2/pi)*acos(exp(-ftip));
                
                fhub=0.5*Nb*(r(i)-Rh)/(Rh*abs(sin(phi)));
                Fhub=(2/pi)*acos(exp(-fhub));

                F=Ftip*Fhub; % total losses

%               convenience parameters
                
                K=SLr*CN/(4*F*sin(phi)^2);
                
                K_p=SLr*CT/(4*F*sin(phi)*cos(phi));

                
                
                
%%              Determine where the solution lies

    
    if K <= 2/3

%               METHOD 1 - momentum theory

                a=K/(1+K);
    else

%               METHOD 2 - Glauret's empirical method

                y1=2*F*K-((10/9) - F);
                y2=2*F*K-F*((4/3) - F);
                y3=2*F*K-((25/9) - 2*F);
                
                if abs(y3) <1E-6
                    a=1-1/(2*sqrt(y2));
                else
                    a=(y1-sqrt(y2))/y3;
                end  
    end
    
               a_p=K_p/(1-K_p);
               
               
               % residual equation %% SOLVED USING BRENTS METHOD
                
              
               R= @(PHI) sin(PHI)/(1-a) - cos(PHI)*(1-K_p)/TSr; % evaluate residual in terms of PHI only
               
               Rpb=@(PHI) sin(PHI)*(1-K)-cos(PHI)*(1-K_p)/TSr;
               
               if R(Ep)*R(0.5*pi) < 0 
                   
                  phi_0= zero ( Ep, 0.5*pi, eps, eps, R);
                   
               elseif Rpb(-0.25*pi)*Rpb(-Ep) < 0
                 
                   phi_0= zero (-0.25*pi,-Ep, eps, eps, Rpb);
               else
               

                   phi_0= zero (0.5*pi,pi-Ep, eps, eps, Rpb);
               end
               
               
               Err=R(phi); % error (from current iteration)

               phi=phi_0; % update new value of flow angle for next iteration
               
               
               % lets break out of the loop if no convergence after 200
               % iterations
               if j>200
                   break
               end
               
end
% SAVE CONVERGED VALUES

a_m(i)=a;
a_pm(i)=a_p;
CT_m(i)=CT;
CN_m(i)=CN;
end

% SET END BOUNDARY

% root
a_m(1)=0.990;
a_pm(1)=0.001;

% tip
a_m(end)=0.990;
a_pm(end)=0.001;

end




               

