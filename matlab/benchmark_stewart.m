clear all
close all

addpath('.\functions')

%time settings
animate = 0; %animation on/off
dt      = 0.001; %step time
t_end   = 60;   %total time
t       = 0:dt:t_end;
N       = length(t);

trans_on  = 1; %translation on/off
rot_on    = 1; %rotation on/off

q   = zeros(6,N); %pose
qd  = zeros(6,N); %velocity
qdd = zeros(6,N); %acceleration

q0 = [0 0 1.5 0 0 0]';  %inital position
qd0 = zeros(6,1);

%set random translation/rotation
for i = 1:N
    
    r = [trans_on*0.5*sin(0.25*t(i)),trans_on*0.5*sin(0.3*t(i)),1.5+trans_on*0.4*sin(0.4*t(i))]';       %translation
    phi = [rot_on*1*10/180*pi*sin(0.1*t(i)),0.5*1*rot_on*(-1)*20/180*pi*sin(0.2*t(i)),0.5*rot_on*1*30/180*pi*sin(0.3*t(i))]';   %rotation
    q(:,i) = [r;phi];
   
    %velocity
    qd(:,i)=(q(:,i)-q0)/dt;         
    q0 = q(:,i);
   
    %acceleration
    qdd(:,i)=(qd(:,i)-qd0)/dt;
    qd0 = qd(:,i);

end
% 
 tic
 [L_bsk,Ld_bsk,Ldd_bsk,b]             = StewartIK_bsk(q,qd,qdd,t,animate);  %IK using basic spatial kinematic
 IK_bsk_time = toc
 tic
 [q_bsk,qd_bsk,qdd_bsk]               = StewartFK_bsk(q0,L_bsk,Ld_bsk,Ldd_bsk,t,animate);  %FK using basic spatial kinematic
 FK_bsk_time = toc
 tic
 [D,Dd,Ddd,Dvel,Dacc,q_cga,Vb,Vbd]    = StewartFK_CGA(q0,L_bsk,Ld_bsk,Ldd_bsk,dt,t,animate,b);    %FK using conformal geometric algebra
 FK_cga_time = toc
  tic
 [L_cga,Ld_cga,Ldd_cga]            = StewartIK_CGA(D,Vb,Vbd,t,animate); %IK using conformal geometric algebra
 IK_cga_time = toc



