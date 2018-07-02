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

% This script plots the static and dynamic lift coefficients at a location
% on the blade for both the rotational and non-rotational cases. The
% measured static and dynamic data is taken from the OSU wind tunnel tests.

% FOR MEASURED DATA SEE: 

%      Reuss Ramsay, Hoffmann, Gregorek - 1995
%      Effects of grit roughness and pitch oscillations on the S809 airfoil
     
clear
close all
clc
graph_settings 

% LOAD DATA

file_turb ='TGL_BLADE_PROFILE';
load(file_turb)                         % load turbine blade profile
r=rad; % radial position on blade

file_foil = 'S814_static_data';
load(file_foil)                 % Load airfoil data

file_ds ='S814_DS_parameters';


% UNSTEADY PITCHING MOTION
 
k=0.091;     % reduced frequency
U0 = 2;      % freestream velocity
W=k*2*U0./c; % frequency
T=2*pi./W;   % Period
 
% Time array
n=2;     
steps=5000;
T_steps=steps/n;
t=linspace(0,round(n.*max(T)),steps);
dt=t(2);
 
% AoA time history 
a0 = deg2rad(10.75);    % amplitude of AoA oscillations
am=deg2rad(13.8);       % mean lift value
at= a0.*sin(W'.*t)+am;  % harmonic pitching motion

%% ============================ PRE-PROCESSOR =========================== %
 
        % ROTATONAL AUGMENTATION / STALL DELAY / SEPARATION POINT
[Values_360, Values_360r] = PreProcessor1(aoa,Cl_2d,Cd_2d,Cn_2d,Clin,LinRange,B,r,c,az);

%% Attached flow solution

[Cl_US,Cl_c,Cl_nc,Ds,aE] = wag(c,t(2),U0,at,az,Clin);
     
%% Dynamic stall without rotation
     
[~,~,Cl_DS_2d,~,~,~,~,~] = DS_2D(Values_360, Cl_US, Cl_c, Cl_nc, Ds, aE, at, file_ds);
     
%% Dynamic stall with rotation

[~, ~, Cl_DS_3d, ~, ~, ~, ~, ~] = DS_3D(B, c, Values_360r, r, Cl_US, Cl_c, Cl_nc, Ds, aE, at, file_ds);
 
%% plots

Alpha=rad2deg(aoa);
j=7; % location where r = 0.47R;

tt=steps-T_steps; % plot last period

load S814_OSU_k_091_a0_10 % LOAD WIND TUNNEL DATA FROM OSU TESTS

close all

% PLOT CL vs ALPHA
figure(1)
plot(Alpha,Cl_2d,'ko-',a_14,Cl_14,'kd',rad2deg(at(j,tt:end)),Cl_DS_2d(j,tt:end),'k',rad2deg(Values_360r.Alpha),Values_360r.Cl(:,j),'r+-',rad2deg(at(j,tt:end)),Cl_DS_3d(j,tt:end),'r')
xlabel('$$\alpha$$ [deg]')
ylabel('$$C_L$$')
legend('Measured static', 'Measured dynamic','Predicted dynamic','Predicted static (rotational)','Predicted dynamic (rotational)','location','best')
legend boxoff
axis([0 30 0.5 2.5])

