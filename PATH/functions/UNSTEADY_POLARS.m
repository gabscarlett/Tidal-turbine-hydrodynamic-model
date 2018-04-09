function [ValuesOut]=UNSTEADY_POLARS(ValuesIn,Na)

% THIS FUNCTION TAKES LIFT AND DRAG COEFFIICENTS FROM THE UNSTEADY DYNAMIC
% STALL MODEL. THE VALUES ARE BINNED AND MEAN VALUES CALCULATED SO AS TO
% FORM DISCREATE LIFT AND DRAG LOOK-UP TABLES WITH ANGLE OF ATTACK.


% INPUTS FROM UNSTEADY SIMULATION
% ValuesIn:
% Cl_DS_3d
% Cd_DS_3d
% AoA

%Na = number of bins


%% BIN DATA AND TAKE AVERAGE

Nr=size(ValuesIn.AoA_In,1); % blade sections

% Prealocation
Cl_M=zeros(Nr,Na);
Cd_M=zeros(Nr,Na);
Alpha_M=zeros(Nr,Na);

% ClSec=zeros(1,Na);
% CdSec=zeros(1,Na);

% loop through each blade section
for j=1:Nr

ClSec=ValuesIn.Cl_In(j,:);
CdSec=ValuesIn.Cd_In(j,:);

[I,EDGES] = discretize(ValuesIn.AoA_In(j,:), Na); % bin alpha values 

% loop through each bin
for i=1:Na

INDEX=find(I==i);                % find the indexes in each bin
Cl_M(j,i)=mean(ClSec(INDEX));    % take the mean of the lift coefficient for each bin
Cd_M(j,i)=mean(CdSec(INDEX));    % take the mean of the drag coefficient for each bin
Alpha_M(j,i)=mean(EDGES(i:i+1)); % take the mean alpha value between the edges of each bin (Na+1 bins gives Na values)
end
%% interpolate for the zero entries (bins where mean(0) ocurred creating NaN values)

nanx = isnan(Cl_M(j,:));
t    = 1:numel(Cl_M(j,:));

Cl_M(j,nanx) = interp1(t(~nanx), Cl_M(1,~nanx), t(nanx));
Cd_M(j,nanx) = interp1(t(~nanx), Cd_M(1,~nanx), t(nanx));
end

% Set negative values of drag to zero
InRange =(Cd_M>0);
Cd_M=Cd_M.*InRange;

%% NOW BEST FIT THE DATA

% Prealocation
CL_M=zeros(Nr,Na);
CD_M=zeros(Nr,Na);

w = ones(size(1,Na)); w([1 end]) = 100;

for k=1:Nr

    % lift
[~,CL_M(k,:),~] = spaps(Alpha_M(k,:),Cl_M(k,:),1, w, 3); % SMOOTHING SPLINE

    % drag
[~,CD_M(k,:),~] = spaps(Alpha_M(k,:),Cd_M(k,:),1.e-2, w, 3); % SMOOTHING SPLINE
end

%%  SAVE VALUES AS STRUCTURED ARRRAY

ValuesOut.aoa=deg2rad(Alpha_M);
ValuesOut.Cl=CL_M;
ValuesOut.Cd=CD_M;

end
