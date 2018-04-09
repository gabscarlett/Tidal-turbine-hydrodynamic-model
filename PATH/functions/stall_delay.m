function [Cl_3d,Cd_3d,Cn_3d,Ct_3d] = stall_delay(B,r,c,alpha,F,Cl_2d,Cd_2d,Clin,az)
%% applies stall-delay to 2D coefficients as a function of angle of atack and radial coordinate

% uses the method of Lindenburg C, 2004, for inboard blade sections
% applies correction of Lindenburg C, 2003, to Cl for outboard blade sections

% empirical values (for S809)
b=1.6; % b1
d=0.25; % b2

%%          CALCULATE COEFFICIENTS
% stall delayed to greater angle of attack due to rotation
delta_aRot= d/(2*pi)*b.*(c./r).*cos(alpha+B).^2;
alpha_rot=alpha+delta_aRot;

% 3d lift coefficient for r < 0.8R
Cl_3d_in=Cl_2d+b.*(c./r).*cos(alpha+B).^2.*((1-F).^2.*cos(alpha_rot)+d.*cos(alpha_rot -az))-Clin.*delta_aRot;

% 3d drag coefficient for all sections
Cd_3d=Cd_2d+b.*(c./r).*cos(alpha+B).^2.*(1-F).^2.*sin(alpha_rot);


% planform area outward from each radial section

%Preallocation
S=ones(1,length(r));

for i=1:length(r)-1
    
S(i)=trapz(r(i:end),c(i:end));   
end
% S(end+1)=1 (non-trivial)

span=r(end)-r; % outward span at each radial section

AR=span.^2./S; % aspect ratio of outer blade from each radial section.

Clp=Clin.*(alpha-az); % lift coefficient in potential flow

% 3d lift coefficient for  r >= 0.8 R
Cl_3d_out=Cl_2d-cos(alpha+B).^2.*exp(-1.5.*AR).*(Clp-Cl_2d).*Cl_2d./Clp; 

InRange=(r>=0.8*r(end)); % boolean array inner; inner == False, outer == True

Cl_3d=Cl_3d_in.*(1-InRange) + Cl_3d_out.*InRange; % use boolean array to form final 3d lift coefficient

Cn_3d=Cl_3d.*cos(alpha)+Cd_3d.*sin(alpha); % normal coefficient
Ct_3d=Cl_3d.*sin(alpha)-Cd_3d.*cos(alpha); % tangential coefficient
end

