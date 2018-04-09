function [MY,MX,T,P]=loads(F,Cl,Cd,phi,r,omega)

    % Written by Gabriel.T.Scarlett July 2017
    
    % THIS FUNCTION CALLS simps
    % simps: numerical integration of the distributed blade loads
    
    % Caluculates moments, thrust and power for a single blade
    % Determines the root and edgewise bending moments
    % Determines the total thrust forcing acting on a
    % single blade
    % Determines the power contribution of the blade

    CN=Cl.*cos(phi)+Cd.*sin(phi);      % Normal to rotor plane
    CT=Cl.*sin(phi)-Cd.*cos(phi);      % Tangential to rotor plane
    
    FN=F.*CN;                          % distributed normal load N/m
    FT=F.*CT;                          % distributed tangential load N/m
    
    MY=simps(r,FN.*r');                % root bending (N m)
    MX=simps(r,FT.*r');                % edgewise bending (N m)
    
    T=simps(r,FN);                     % total thrust force (N) 
    P=MX.*omega;                       % output power  (W)
end
