function [Cn_DS, Ct_DS, Cl_DS, Dvis, Cd_Ind, ff, fff, VortexTracker] = DS_3D(B, c, Values_360r, r, Cp, Cn_c, Cn_nc, ds, aE, at)

%       DYNAMIC STALL MODULE FOR NON-UNIFORM FORCING

%       Written G.Scarlett, The University of Edinburgh, April 2017

%       Dynamic stall model following Sheng Glabraith and Cotton

%       Mach effects are ignored.


%%              LOAD DYNAMIC STALL DATA

% determine static separation point for simulation
rr=r'.*ones(size(at));
f=interp2(r,Values_360r.Alpha,Values_360r.F,rr,at,'spline');  % spline outside loop
% !!!! Avoid negative values !!!
InRange=(f>0);
f=f.*InRange;

faE=interp2(r,Values_360r.Alpha,Values_360r.F,rr,(aE),'spline');
InRange=(faE>0);
faE=faE.*InRange;

load S814_DS_parameters

% angle shift to account for rotational stall delay

delta_aRot_ass= 0.25/(2*pi)*1.6.*(c./r).*cos(ass+B).^2;
ass=ass+delta_aRot_ass;

delta_aRot_ads0= 0.25/(2*pi)*1.6.*(c./r).*cos(ads0+B).^2;
ads0=ads0+delta_aRot_ads0;


da1=ads0-ass;
                
a_cr = ads0';
                
%       Initial values
        
        tau_Prev=0;
        Da_Prev=0;
        Df_Prev=0;
        
%%      *********************** Part 1 ************************************ 

%       Time step through angle of attach history 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ela=length(at(1,:)); % original length of AoA
% lets run many times for numerical convergence
for i=1:3
at=horzcat(at,at); % make vector with of 16*AoA length
f=horzcat(f,f); % make vector with of 16*AoA length
faE=horzcat(faE,faE); % make vector with of 16*AoA length
aE=horzcat(aE,aE); % make vector with of 16*AoA length
Cp=horzcat(Cp,Cp); % same to attached flow solution
Cn_c=horzcat(Cn_c,Cn_c); % circular component of attached solution
Cn_nc=horzcat(Cn_nc,Cn_nc); % non-circular component of attached solution
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ela2= length(at(1,:));

% Preallocation
% Preallocation
Cn_DS=ones(size(at));
Ct_DS=ones(size(at));
ff=ones(size(at));
fff=ones(size(at));
Dvis=ones(size(at));
aInd=ones(size(at));

for i=2:ela2
    
%%       Calculate lagged angle of attack
    
%       Deficiency function
    
        Da=Da_Prev.*exp(-ds'./Ta)+(at(:,i)-at(:,i-1)).*exp(-ds'./(2*Ta)); 
    
        a_p = at(:,i) - Da; 

            
%       Determine ff from look up table using lagged angle
          
         ff(:,i)=interp2(r',Values_360r.Alpha,Values_360r.F,r',(a_p-da1'),'linear'); % linear inside loop for speed!!

%       Determine fff same as original LB model, but with Tv
% 
        Df=Df_Prev.*exp(-ds'./Tv)+(ff(:,i)-ff(:,i-1)).*exp(-ds'./(2*Tv));

%       Apply lag to separation point

        fff(:,i)=ff(:,i)-Df; 
        
%       Calculate lift due to trailing edge separation

        Cn_Ts=Cp(:,i).*((1+sqrt(fff(:,i)))/2).^2; 
        

   
%%      *********************** Part 2 ************************************ 

%                             VORTEX LIFT
 

%       Evaluate vortex tracking time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InRange2= (abs(a_p) >= a_cr); % stalled if true
tau_0 = InRange2.*(tau_Prev + ds');  % Update vortex passage time        
        
%       Calculate vortex shape function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InRange3=(tau_0 > 0) & (tau_0 <= Tv);
Vx3=(sin(pi.*tau_0./(2*Tv))).^(1.5); % first vortex

InRange4= (tau_0 > Tv);
Vx4 = (cos(pi.*(tau_0 - Tv)./TvL)).^2; % subsequent shedding


Vx=InRange3.*Vx3+(InRange4.*(1-InRange3)).*Vx4;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InRange5 =(abs(at(:,i)) < abs(at(:,i-1))); % if pitching down

Vx = ~InRange5.*Vx;
VortexTracker(:,i)=Vx;  
    
%       Calculate vortex lift contribution
        Cnv =B1.*(fff(:,i)-f(:,i)).*Vx;
    
    
%%                    TOTAL NORMAL AND TANGENTIAL CONTRIBUTION
                
%      ************************** NORMAL COMPONENT ************************

                            Cn_DS(:,i) = Cn_Ts + Cnv; 
                
%      ************************ TANGENTIAL COMPONENT **********************

                            Ct_DS(:,i)=eta.*Clin.*(aE(:,i)-az).^2.*(sqrt(ff(:,i))-E0);
                            
%      ************************** DRAG COMPONENT **************************
                            
                            Dvis(:,i)= (0.5*(1-sqrt(ff(:,i)))).^2 - (0.5*(1-sqrt(faE(:,i)))).^2;
                            aInd(:,i)=(at(:,i)-aE(:,i));
    
    
%%                         Update values

                            Da_Prev=Da;
                            tau_Prev=tau_0;
                            Df_Prev=Df;

end

% lets slice the second last period to ensure
% a converged solution with no end points

ff =ff(:,end-2*ela+1:end-ela); % f and ff can be output for post processing
fff =fff(:,end-2*ela+1:end-ela);

Cn_DS =Cn_DS(:,end-2*ela+1:end-ela);
Ct_DS =Ct_DS(:,end-2*ela+1:end-ela);
Dvis =Dvis(:,end-2*ela+1:end-ela);
aInd =aInd(:,end-2*ela+1:end-ela);
at =at(:,end-2*ela+1:end-ela);
VortexTracker=VortexTracker(:,end-2*ela+1:end-ela);

% Lift coefficient
Cl_DS=Cn_DS.*cos(at)+Ct_DS.*sin(at);

% Induced drag term
Cd_Ind=aInd.*Cl_DS;

end