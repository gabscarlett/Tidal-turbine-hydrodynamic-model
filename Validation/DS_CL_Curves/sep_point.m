function [f]=sep_point(alpha_f,az,Cn,Clin,LinRange)
%% Create separation point look up table using normal coefficient

%                       KIRCHHOFF THEORY INVERTED TO SOLVE FOR F
%                       MAX F SET TO 1
%                       LINEAR LOW ANGLE VALUES SET TO 1
%                       OUT OF RANGE VALUES REMOVED


%       Use this script to load in data and read static normal coefficient
%       (Cn) vs angle of attack (alpha) and linear curve (Clin) the range
%       of the linear lift curve (LinRange) and the angle of zero lift (az).
%       OK to use Cl if Cn is not easily available


% determine f through inversion of the Kirchhoff equation
% see Thwaits

F=(2*sqrt(abs(Cn./(Clin.*(alpha_f-az))))-1).^2;

% Set values in linear 2D range to 1
InRange1=((alpha_f>=LinRange(1))& alpha_f<=LinRange(2));
F=F.*(~InRange1)+InRange1;

% locate and set 3D range to 1
InRange2= F>1;
f=F.*~InRange2 +InRange2;