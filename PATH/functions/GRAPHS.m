function [] = GRAPHS(r,R,t,Tr,A,U0,P,P_s,T,T_s,Cl_DS_3D,Cl_QS3d,Cd_DS_3D,Cd_QS3d,CN_Ds3d,CN_Qs3d,CT_Ds3d,CT_Qs3d,Cl_S3d,Cd_S3d,CN_S3d,CT_S3d,ff_3d,f_3d,f,aoa,AoA,AoA_s,Values_360r,MZ,MZs,MZqs,MX,MXs,MXqs)

% THIS FUNCTION MAKES GRAPHS FOR RENEWABLE ENERGY 2017 JOURNAL PAPER

% calls function:
% graph_settings - sets the font and latex compiler for graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PLOT CP, CT, CD, CL, AOA, f    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho=1025;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POWER AND THRUST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TOTAL POWER
Pin=0.5*rho*A*U0.^3; % total power available to rotor
Pout=sum(P,2); 
CP=Pout./Pin;
mean_CP=ones(size(CP))*mean(CP);

Pout_s=3*P_s;
CP_s=Pout_s/Pin;

% percentage difference from mean - unsteady
POWER_DIFF=(mean_CP(1)-CP_s(1))/CP_s(1);

% TOTAL THRUST
Tin=0.5*rho*A*U0.^2;
Tout=sum(T,2); 
CT=Tout./Tin;
mean_CT=ones(size(CT))*mean(CT);
Tout_s=3*T_s;
CT_s=Tout_s/Tin;

% percentage difference from steady and unsteady
THRUST_DIFF=(mean_CT(1)-CT_s(1))/CT_s(1);


%% FIGURE: 
figure(1)
% CP
ax1 = subplot(2,1,1);
plot(ax1,t/Tr,CP,'k',t/Tr,mean_CP,'k--',t/Tr,CP_s,'r:')
ylabel('$$C_P$$')
xlabel('$$t/T_r$$')
axis([0 10 0 1])

% CT
ax2 = subplot(2,1,2);
plot(ax2,t/Tr,CT,'k',t/Tr,mean_CT,'k--',t/Tr,CT_s,'r:')
ylabel('$$C_T$$')
xlabel('$$t/T_r$$')
legend('Unsteady','Unsteady mean','Steady','Location','best')
legend boxoff
axis([0 10 0 1]) 

%print('CP_CT','-depsc2','-r1000');

    
%% mean CL, CD, CT and CMx

        % MEAN VALUES
% % CL
 
Mean_Clds3D=mean(Cl_DS_3D(:,:,1),2);
Mean_Clqs3D=mean(Cl_QS3d(:,:,1),2);

 
% % CD

Mean_Cdds3D=mean(Cd_DS_3D(:,:,1),2);
Mean_Cdqs3D=mean(Cd_QS3d(:,:,1),2);

% % C_thrust
 
Mean_CNds3D=mean(CN_Ds3d(:,:,1),2);
Mean_CNqs3D=mean(CN_Qs3d(:,:,1),2);
 
% % C_tan
 
Mean_CTds3D=mean(CT_Ds3d(:,:,1),2);
Mean_CTqs3D=mean(CT_Qs3d(:,:,1),2);

%% FIGURE: 
figure(2)
ax1=subplot(2,2,1);
plot(ax1,r/R,Mean_Clds3D,'k',r/R,Mean_Clqs3D,'k--',r/R,Cl_S3d(:,1),'r:')
ylabel('$$\bar{C_L}$$')
ax2=subplot(2,2,2);
plot(ax2,r/R,Mean_Cdds3D,'k',r/R,Mean_Cdqs3D,'k--',r/R,Cd_S3d(:,1),'r:')
ylabel('$$\bar{C_D}$$')
legend('Unsteady','Quasi-steady','Steady','Location','best');
legend boxoff

ax3=subplot(2,2,3);
plot(ax3,r/R,Mean_CNds3D,'k',r/R,Mean_CNqs3D,'k--',r/R,CN_S3d(:,1),'r:')
ylabel('$$\bar{C_T}$$')
xlabel('$$r/R$$')
ax4=subplot(2,2,4);
plot(ax4,r/R,Mean_CTds3D.*r'/R,'k',r/R,Mean_CTqs3D.*r'/R,'k--',r/R,CT_S3d(:,1).*r'/R,'r:')
ylabel('$$\bar{C_Q}$$')
xlabel('$$r/R$$')
%print('mean_cl_cd','-depsc2','-r1000');

% plot mean sectional angle of attack
aSteady=mean(AoA_s,2);
aUnsteady=mean(AoA(:,:,1),2);

graph_settings
figure(500)
plot(r/R,aSteady,'r',r/R,aUnsteady,'k')
ylabel('$$\alpha$$ [deg]')
xlabel('$$r/R$$')
legend('Steady','Unsteady mean') 
%% % MULTI PLOT SHOWING THE EFFECTS OF UNSTEADINESS AT 3 BLADE LOCATIONS
X=95;
XX=30;
XXX=1;
len1 =0;
len2 =5;
BLADE=3;

figure(3)
% f at 3 sections
ax1=subplot(3,3,1);
plot(ax1,t/Tr,ff_3d(X,:,BLADE),'k',t/Tr,f_3d(X,:,BLADE),'k--',t/Tr,f(X,:),'r:')
ylabel('$$f$$')
legend('Unsteady','Quasi-steady','Steady','Location','best');
legend boxoff
axis([len1 len2 0 1.5])

ax2=subplot(3,3,2);
plot(ax2,t/Tr,ff_3d(XX,:,BLADE),'k',t/Tr,f_3d(XX,:,BLADE),'k--',t/Tr,f(XX,:),'r:')
axis([len1 len2 0 1.5])

ax3=subplot(3,3,3);
plot(ax3,t/Tr,ff_3d(XXX,:,BLADE),'k',t/Tr,f_3d(XXX,:,BLADE),'k--',t/Tr,f(XXX,:),'r:')
axis([len1 len2 0 1.5])

% AoA at 3 sections
ax4=subplot(3,3,4);
plot(ax4,t/Tr,AoA(X,:,BLADE),'k',t/Tr,AoA_s(X,:),'r:')
ylabel('$$\alpha$$ \rm [deg]')
axis([len1 len2 0 30])

ax5=subplot(3,3,5);
plot(ax5,t/Tr,AoA(XX,:,BLADE),'k',t/Tr,AoA_s(XX,:),'r:')
axis([len1 len2 0 30])

ax6=subplot(3,3,6);
plot(ax6,t/Tr,AoA(XXX,:,BLADE),'k',t/Tr,AoA_s(XXX,:),'r:')
axis([len1 len2 0 30])

% Cl time history at 3 sections
ax7=subplot(3,3,7);
plot(ax7,t/Tr,Cl_DS_3D(X,:,BLADE),'k',t/Tr,Cl_QS3d(X,:,BLADE),'k--',t/Tr,Cl_S3d(X,:),'r:')
ylabel('$$C_L$$')
xlabel('\it t/Tr') 
axis([len1 len2 0 3])

ax8=subplot(3,3,8);
plot(ax8,t/Tr,Cl_DS_3D(XX,:,BLADE),'k',t/Tr,Cl_QS3d(XX,:,BLADE),'k--',t/Tr,Cl_S3d(XX,:),'r:')
xlabel('\it t/Tr')
axis([len1 len2 0 3])

ax9=subplot(3,3,9);
plot(ax9,t/Tr,Cl_DS_3D(XXX,:,BLADE),'k',t/Tr,Cl_QS3d(XXX,:,BLADE),'k--',t/Tr,Cl_S3d(XXX,:),'r:')
xlabel('\it t/Tr') 
axis([len1 len2 0 3])

%print('Unsteady_effects','-depsc2','-r1000');

%% Hysterisis plot

% lift and drag
Alpha=rad2deg(Values_360r.Alpha);
figure(25)

ax1=subplot(2,3,1);
plot(ax1,AoA(X,1:200,1),Cl_DS_3D(X,1:200,1),'k',Alpha,Values_360r.Cl(:,X),'r:')
ylabel('$$C_L$$ \rm [-]')
legend('Unsteady','Static curve','Location','best');
legend boxoff
axis([0 25 0.5 3])

ax2=subplot(2,3,2);
plot(ax2,AoA(XX,1:200,1),Cl_DS_3D(XX,1:200,1),'k',Alpha,Values_360r.Cl(:,XX),'r:')
axis([0 25 0.5 3])

ax3=subplot(2,3,3);
plot(ax3,AoA(XXX,1:200,1),Cl_DS_3D(XXX,1:200,1),'k',Alpha,Values_360r.Cl(:,XXX),'r:')
axis([0 25 0.5 3])

ax4=subplot(2,3,4);
plot(ax4,AoA(X,1:200,1),Cd_DS_3D(X,1:200,1),'k',Alpha,Values_360r.Cd(:,X),'r:')
xlabel('$$\alpha$$ \rm [deg]')
ylabel('$$C_D$$ \rm [-]')
axis([0 25 -0.04 0.6])

ax5=subplot(2,3,5);
plot(ax5,AoA(XX,1:200,1),Cd_DS_3D(XX,1:200,1),'k',Alpha,Values_360r.Cd(:,XX),'r:')
xlabel('$$\alpha$$ \rm [deg]')
axis([0 25 -0.04 0.6])

ax6=subplot(2,3,6);
plot(ax6,AoA(XXX,1:200,1),Cd_DS_3D(XXX,1:200,1),'k',Alpha,Values_360r.Cd(:,XXX),'r:')
xlabel('$$\alpha$$ \rm [deg]')
axis([0 25 -0.04 0.6])
%print('Hst','-depsc2','-r1000');


%% CMz and CMy plot

M_IN=0.5*A*rho*U0.^2*R;

CMz=MZ(:,1)/M_IN;
CMx=MX(:,1)/M_IN;
Tdum=0:0.5:20;
mean_cmz=mean(CMz).*ones(size(Tdum));
mean_cmx=mean(CMx).*ones(size(Tdum));

CMz_qs=MZqs(:,1)/M_IN;
CMx_qs=MXqs(:,1)/M_IN;
mean_cmz_qs=mean(CMz_qs).*ones(size(t));
mean_cmx_qs=mean(CMx_qs).*ones(size(t));

CMz_s=MZs/M_IN;
CMx_s=MXs/M_IN;

p_diff_cmz=(mean_cmz(1)-CMz_s(1))/CMz_s(1);
p_diff_cmx=(mean_cmx(1)-CMx_s(1))/CMx_s(1);


p_diff_cmz2=(mean_cmz(1)-mean_cmz_qs(1))/mean_cmz_qs(1);
p_diff_cmx2=(mean_cmx(1)-mean_cmx_qs(1))/mean_cmx_qs(1);

figure(112)
ax4=subplot(2,1,1);
plot(ax4,t/Tr,CMz,'k',Tdum,mean_cmz,'k-o',t/Tr,CMz_qs,'b-.',t/Tr,mean_cmz_qs,'b--',t/Tr,CMz_s,'r:','MarkerSize',4)
xlabel('\it t/Tr')
ylabel('$$C_{My}$$')
axis([0 5 0.1 0.3])
ax2=subplot(2,1,2);
plot(ax2,t/Tr,CMx,'k',Tdum,mean_cmx,'k-o',t/Tr,CMx_qs,'b-.',t/Tr,mean_cmx_qs,'b--',t/Tr,CMx_s,'r:','MarkerSize',4 )
xlabel('\it t/Tr')
ylabel('$$C_{Mx}$$')
axis([0 5 0 0.2])
legend('Unsteady','Unsteady mean','Quasi-steady','Quasi-steady mean','Steady','location','Best')
legend boxoff
%print('CMB','-depsc2','-r1000');


end