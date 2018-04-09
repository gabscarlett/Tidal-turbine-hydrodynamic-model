function [Cl,Cl_c,Cl_nc,Ds,aE] = wag(c,Dt,U,a,az,Clin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Gabriel Scarlett August 2017 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % ATTACHED FLOW SOLUTION FOR GENERAL FORCING
 
% Indicial lift using numerical recursive method from Leishman, 2006

%%%%%%%%%%%%%%%%%%%%%%%%%

% Wagner (Jones, W.P)
A1=0.165;
A2=0.335;
b1=0.041;
b2=0.32;

% Reduced time (s)
Ds=2*U*Dt./c;

% Previous values
X_Prev=0;
Y_Prev=0;

ela = length(a(1,:));

% lets run many times and take second last solution
for i=1:4
a=horzcat(a,a); % make vector 16 time the length of AoA
end

at = diff(a,1,2)/Dt;   % first time derivative of AoA

ela2= length(a(1,:));

% Preallocation
aE=ones(size(a));
Cl_c=ones(size(a));
Cl_nc=ones(size(a));

for i=2:ela2-1
   
    
% Exact integration method
X=X_Prev.*exp(-b1*Ds') +A1*(a(:,i)-a(:,i-1))./Ds'.*(1-exp(-b1*Ds'))/b1;
Y=Y_Prev.*exp(-b2*Ds') +A2*(a(:,i)-a(:,i-1))./Ds' .*(1-exp(-b2.*Ds'))/b2;

aE(:,i)=a(:,i) -X -Y;                 % Effective angle of attack

Cl_c(:,i)=Clin.*(aE(:,i)-az);          % Circulatory lift component

Cl_nc(:,i)=0.5*Clin.*at(:,i).*c'/(2*U); % Added mass of AoA oscillations


% Update previous values
X_Prev=X;
Y_Prev=Y;

end

% lets slice the second last period
% ela+1 because we stop at ela2 - 1
aE=aE(:,end-2*ela+1:end-ela);
Cl_c =Cl_c(:,end-2*ela+1:end-ela);
Cl_nc =Cl_nc(:,end-2*ela+1:end-ela);

% interpolate the endpoint
aE(:,end)=aE(:,end-1)+(aE(:,end-2)-aE(:,end-3));
Cl_c(:,end)=Cl_c(:,end-1)+(Cl_c(:,end-2)-Cl_c(:,end-3));
Cl_nc(:,end)=Cl_nc(:,end-1)+(Cl_nc(:,end-2)-Cl_nc(:,end-3));

Cl=Cl_c + Cl_nc;                            % Total linear unsteady lift coefficient
end


