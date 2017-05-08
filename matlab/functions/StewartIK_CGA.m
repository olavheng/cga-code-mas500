function [L,Ld,Ldd] = StewartIK_CGA(DB,qd,qdd,t,sim)
%function [L,Ld,Ldd] = StewartIK_CGA(DB,qd,qdd,t,sim)
%
% This function is a inverse kinematic algortithm for a 6-SPS stewart
% platform. It applies Conformal Geometric Algebra and the property of the 
% inner product for distrance measure. From given pose, veocity and
% acceleration of the platform, the leg-length, leg-speed and
% leg-acceleration are returned.
% 
% Created by Olav Heng
% INPUT:
% DB        : The displacement versor [x0 x1 x2 x3 y0 y1 y2 y3]' for the stewart platfrom, size:(8,n)
% qd        : The velocity screw [x y z phi theta psi]' for the stewart platfrom body frame, size:(6,n)
% qdd       : The acceleration screw [x y z phi theta psi]' for the stewart platfrom body frame, size:(6,n)
% t         : Time for each step n, required for simulation, size:(1,n) , set 0 if not used.
% sim       : Simulation on/off: 1/0
%
% OUTPUT
% L         : Leg-length from given pose (6,n)
% Ld        : Leg extraction speed from pose and velocity. (6,n)
% Ldd       : Leg extraction acceleration from pose, velocity and
% acceleration of the platform (6,n)

N = length(DB);

%Stewart parameter
jointLoc = [15,105,135,225,255,345]/180*pi; %location of the base points "a"
a0       = [2,0,0]';

Az = @(psi)[cos(psi),-sin(psi), 0;
            sin(psi), cos(psi), 0;
            0       , 0       , 1];

%distribute spherical joints "a" in the base, and "bp" on the platform
a   = zeros(3,6);
bp  = zeros(3,6);
for i = 1:length(jointLoc)
    
    a(:,i) = Az(jointLoc(i)) * a0;
    
    if i == 6 %shifting joint position in bp so that [a1,a2...] and [bp1,bp2...] corrospond to each other
    bp(:,1)   = Az(60/180*pi) * a(:,i);
    else
    bp(:,i+1) = Az(60/180*pi) * a(:,i);
    end

end
%%%Define sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = quat2rotm(D) %Dual quaternions to homogenious transformation matrix
    X0 = D(1,1);
    X1 = D(2,1);
    X2 = D(3,1);
    X3 = D(4,1);
    Y0 = D(5,1);
    Y1 = D(6,1);
    Y2 = D(7,1);
    Y3 = D(8,1);

    H      = zeros(4,4);
    H(1,:) = [X1^2-X2^2-X3^2+X0^2, 2*(X1*X2-X3*X0), 2*(X1*X3+X2*X0), 2*(X1*Y0-Y1*X0-X2*Y3+Y2*X3)];
    
    H(2,:) = [2*(X1*X2+X3*X0), -X1^2+X2^2-X3^2+X0^2, 2*(-X1*X0+X2*X3), 2*(X1*Y3-Y1*X3+X2*Y0-Y2*X0)];
    
    H(3,:) = [2*(X1*X3-X2*X0), 2*(X1*X0+X2*X3), -X1^2-X2^2+X3^2+X0^2, 2*(-X1*Y2+Y1*X2+X3*Y0-Y3*X0)];
    
    H(4,:) = [0, 0, 0, (X1^2+X2^2+X3^2+X0^2)];
    
    H = H / H(4,4);
end

function PlotAxes(H,axisLength,lineWidth) %a function for ploting axes in figures
% ** Input **
% H         : 4 by 4 homogeneous transformation matrix
% axisLegth : length of the axises
% lineWidth : width of the axis lines


% Coordinate center
u = H(1:3,4);

v = u + axisLength*H(1:3,1);
plot3([u(1) v(1)],[u(2) v(2)],[u(3) v(3)],'r','LineWidth',lineWidth);
hold on

v = u + axisLength*H(1:3,2);
plot3([u(1) v(1)],[u(2) v(2)],[u(3) v(3)],'g','LineWidth',lineWidth);

v = u + axisLength*H(1:3,3);
plot3([u(1) v(1)],[u(2) v(2)],[u(3) v(3)],'b','LineWidth',lineWidth);
end

function y =  sandwich(D,b) %geometric product sandwich function
 
     X0 = D(1,1);
     X1 = D(2,1);
     X2 = D(3,1);
     X3 = D(4,1);
     Y0 = D(5,1);
     Y1 = D(6,1);
     Y2 = D(7,1);
     Y3 = D(8,1);

     b1 = b(1,1);
     b2 = b(2,1);
     b3 = b(3,1);
     b5 = 1;

     p1 = (X0^2*b1+2.*X0*X2*b3-2.*X0*X3*b2-2.*X0*Y1*b5+X1^2*b1+2.*X1*X2*b2+2.*X1*X3*b3+2.*X1*Y0*b5-X2^2*b1-2.*X2*Y3*b5-X3^2*b1+2.*X3*Y2*b5);%/(X0^2+X1^2+X2^2+X3^2);
     p2 = (X0^2*b2-2.*X0*X1*b3+2.*X0*X3*b1-2.*X0*Y2*b5-X1^2*b2+2.*X1*X2*b1+2.*X1*Y3*b5+X2^2*b2+2.*X2*X3*b3+2.*X2*Y0*b5-X3^2*b2-2.*X3*Y1*b5);%/(X0^2+X1^2+X2^2+X3^2);
     p3 = (X0^2*b3+2.*X0*X1*b2-2.*X0*X2*b1-2.*X0*Y3*b5-X1^2*b3+2.*X1*X3*b1-2.*X1*Y2*b5-X2^2*b3+2.*X2*X3*b2+2.*X2*Y1*b5+X3^2*b3+2.*X3*Y0*b5);%/(X0^2+X1^2+X2^2+X3^2);
     
     y = [p1 p2 p3]';

end

function Ld = LD_cga(D,V,L,a,bp) %constraint equation for actuator speed


    a1 = a(1);
    a2 = a(2);
    a3 = a(3);
    a5 = 1;

    b1 = bp(1);
    b2 = bp(2);
    b3 = bp(3);
    b5 = 1;

    v1 = V(1,1);
    v2 = V(2,1);
    v3 = V(3,1);
    w1 = V(4,1);
    w2 = V(5,1);
    w3 = V(6,1);

    X0 = D(1,1);
    X1 = D(2,1);
    X2 = D(3,1);
    X3 = D(4,1);
    Y0 = D(5,1);
    Y1 = D(6,1);
    Y2 = D(7,1);
    Y3 = D(8,1);

     Ld = -1/L*(-1.0*a2*X0^2*b3*w1+1.0*a2*b5*X2^2*v2-1.00*a5*X0^2*b1*v1+1.00*a5*X2^2*b3*v3+1.00*a5*X3^2*b2*v2+1.00*a5*X1^2*b2*v2...
         +1.00*a5*X3^2*b1*v1+1.00*a5*X1^2*b3*v3-1.00*a5*X2^2*b2*v2-1.00*a5*X3^2*b3*v3-1.00*a5*X0^2*b3*v3-1.00*a5*X1^2*b1*v1...
         -1.00*a5*X0^2*b2*v2+1.00*a5*X2^2*b1*v1+1.0*a2*X1^2*b1*w3-1.0*a2*X3^2*b3*w1+1.0*a2*b5*X1^2*v2+1.0*a2*X0^2*b5*v2...
         +1.0*a2*X1^2*b3*w1-1.0*a2*X3^2*b1*w3-1.0*a2*X2^2*b1*w3+1.0*a2*X0^2*b1*w3+1.0*a2*X2^2*b3*w1+1.0*a2*b5*X3^2*v2...
         -1.0*a1*X2^2*b2*w3+1.0*a1*X2^2*b5*v1+1.0*a1*X3^2*b2*w3-1.0*a1*X1^2*b3*w2+1.0*a1*X3^2*b3*w2+1.0*a1*X1^2*b5*v1...
         -1.0*a1*X2^2*b3*w2-1.0*a1*X0^2*b2*w3+1.0*a1*X0^2*b3*w2+1.0*a1*X1^2*b2*w3+1.0*a1*X0^2*b5*v1+1.0*a1*b5*X3^2*v1...
         -1.0*a3*X1^2*b2*w1+1.0*a3*X3^2*b5*v3+1.0*a3*X0^2*b5*v3+1.0*a3*X2^2*b5*v3-1.0*a3*X3^2*b2*w1-1.0*a3*X1^2*b1*w2...
         +1.0*a3*X2^2*b2*w1+1.0*a3*X0^2*b2*w1+1.0*a3*X3^2*b1*w2+1.0*a3*b5*X1^2*v3+1.0*a3*X2^2*b1*w2-1.0*a3*X0^2*b1*w2...
         +2*a3*X2*X1*b1*w1-2*a3*X2*X1*b2*w2-2*a3*Y2*X0*b5*w1+2*a3*X0*Y1*b5*w2+2*a3*X2*Y0*b5*w1+2*a3*X2*Y3*b5*w2-2*a3*X3*X1*b3*w2...
         -2*a3*X0*X1*b3*w1+2*a3*X0*X3*b1*w1+2*a3*X0*X3*b2*w2-2*a3*Y2*X3*b5*w2-2*a3*X3*Y1*b5*w1-2*a3*Y0*X1*b5*w2+2*a3*X1*Y3*b5*w1...
         +2.0*a5*Y1*X0*b5*v1-2.0*a5*Y1*X2*b5*v3+2.0*a5*Y1*X3*b5*v2+2.0*a5*b5*Y3*v3*X0-2.0*a5*b5*Y3*v2*X1+2.0*a5*b5*Y3*v1*X2...
         -2.0*a5*Y0*X1*b5*v1-2.0*a5*X0*b1*v2*X3+2.0*a5*X0*b2*v1*X3+2.0*a5*b5*Y2*v2*X0+2.0*a5*b5*Y2*v3*X1-2.0*a5*b5*Y2*v1*X3...
         -2.0*a5*X3*X1*b1*v3-2.0*a5*X3*X1*b3*v1-2.0*a5*X3*X2*b2*v3-2.0*a5*X3*X2*b3*v2-2.0*a5*X0*b2*v3*X1+2.0*a5*X0*b3*v2*X1...
         +2.0*a5*X0*b1*v3*X2-2.0*a5*X0*b3*v1*X2-2.0*a5*X1*X2*b1*v2-2.0*a5*X1*X2*b2*v1-2.0*a5*Y0*X2*b5*v2-2.0*a5*Y0*X3*b5*v3...
         -2*a1*X2*X1*b1*w3+2*a2*Y0*X1*b5*w3-2*a3*X2*X0*b3*w2+2*a3*X2*X3*b3*w1-2*a1*X2*X0*b1*w2-2*a1*X0*X3*b1*w3+2*a1*X2*X3*b2*w2...
         -2*a1*X2*X3*b3*w3+2*a1*Y1*w3*X3*b5+2*a1*X3*Y0*b5*w2-2*a1*Y2*X1*b5*w2-2*a1*Y3*X1*b5*w3-2*a1*X2*Y0*b5*w3+2*a1*X2*Y1*b5*w2...
         +2*a1*Y2*X0*b5*w3-2*a1*Y3*X0*b5*w2+2*a1*X1*X3*b1*w2+2*a1*X1*X0*b2*w2+2*a1*X1*X0*b3*w3+2*a2*Y2*w1*X1*b5-2*a2*X0*X3*b2*w3...
         -2*a2*X3*X1*b1*w1+2*a2*X3*X1*b3*w3-2*a2*X3*X2*b2*w1-2*a2*X0*X1*b2*w1+2*a2*X0*X2*b1*w1+2*a2*X0*X2*b3*w3+2*a2*X1*X2*b2*w3...
         -2*a2*Y3*X2*b5*w3-2*a2*X3*Y0*b5*w1+2*a2*Y2*w3*X3*b5-2*a2*Y1*X2*b5*w1-2*a2*X0*Y1*b5*w3+2*a2*Y3*X0*b5*w1);
end
%  
function Ldd = LDD_cga(D,V,Vd,L,LD,a,bp) %constaint equation for actuator acceleration


    a1 = a(1);
    a2 = a(2);
    a3 = a(3);
    a5 = 1;

    b1 = bp(1);
    b2 = bp(2);
    b3 = bp(3);
    b5 = 1;

    v1 = V(1,1);
    v2 = V(2,1);
    v3 = V(3,1);
    w1 = V(4,1);
    w2 = V(5,1);
    w3 = V(6,1);

    v1d = Vd(1,1);
    v2d = Vd(2,1);
    v3d = Vd(3,1);
    w1d = Vd(4,1);
    w2d = Vd(5,1);
    w3d = Vd(6,1);

     X0 = D(1,1);
     X1 = D(2,1);
     X2 = D(3,1);
     X3 = D(4,1);
     Y0 = D(5,1);
     Y1 = D(6,1);
     Y2 = D(7,1);
     Y3 = D(8,1);
 
     Ldd = 1/L*(-LD^2-(-2*a1*X0*X1*b3*w1*w2+2*a1*X2*Y1*b5*w2d-2.0*a1*X1*b2*X2*w2^2-1.0*a1*X2^2*b5*w3*v2-2.0*a1*X3*b3*X1*w2^2 ...
         -2*a1*X1*b1*w3d*X2-2*a1*X3*b3*w3d*X2-2*a1*X2*Y0*b5*w3d+2*a1*X1*b3*w3d*X0+1.00*a5*b3*w1*X0^2*v2+1.00*a5*b2*w1*X1^2*v3...
         -1.00*a5*b2*w1*X0^2*v3+1.00*a5*b3*w1*X3^2*v2-1.00*a5*b2*w3*X3^2*v1-1.00*a5*b1*w2*X2^2*v3-1.00*a5*b2*w1*X2^2*v3...
         -1.00*a5*b1*w3*X0^2*v2+1.00*a5*b1*w2*X1^2*v3+1.00*a5*b3*w2*X1^2*v1-1.00*a5*b3*w1*X2^2*v2+1.00*a5*b2*w3*X0^2*v1...
         +2*a2*X3*b5*w2*w1*Y2-2*a2*w2*X2*b5*w1*Y3+2*a2*w2*X2*b5*w3*Y1-2*a3*X0*X3*b2*w1*w3+2*a3*X1*X2*b1*w2*w3+2*a3*X1*X2*b2*w1*w3...
         +2*a3*X2*X3*b3*w2*w3+2*a3*w2*X1*b5*w3*Y3-2*a3*X3*b5*w2*w3*Y1+2*a3*w3*X3*b5*w1*Y2+2*a3*w1*X1*b5*w3*Y0+2*a3*X2*b5*w3*w2*Y0...
         -2*a3*X0*X1*b3*w2*w3-2.0*a5*b3*w1*X2*v3*X3+2.0*a5*b2*w1*X3*v2*X2+2*a1*X1*X3*b1*w1*w3+2*a1*X2*X3*b2*w1*w3-2*a3*w3*X0*b5*w2*Y2...
         +2*a3*X1*X3*b3*w1*w3+2*a3*X0*X2*b3*w1*w3+2*a3*X0*X3*b1*w2*w3+2*a2*w3*X3*b5*w2*Y0+2*a2*X1*X3*b1*w2*w3+2*a2*X1*X3*b3*w1*w2...
         +2*a2*X0*X1*b2*w2*w3-2*a2*X0*X2*b1*w2*w3-2*a2*X0*b5*w2*w3*Y3-2*a2*X1*b5*w3*w2*Y2+2*a2*X0*X2*b3*w1*w2+2*a2*X1*X2*b2*w1*w2...
         -2*a2*X0*X3*b2*w1*w2+2*a2*X2*X3*b2*w2*w3-2*a1*X1*Y3*b5*w3d+1.0*a1*X0^2*b5*w2*v3+1.0*a1*X0^2*b2*w1*w2-1.0*a1*X3^2*b2*w1*w2...
         +2*a1*X3*Y0*b5*w2d+2*a1*X3*Y1*b5*w3d-2.0*a1*X3*Y2*b5*w2^2-1.0*a3*X3^2*b2*w2*w3-1.0*a3*X0^2*b5*v1*w2+1.0*a3*X3^2*b5*v2*w1...
         +1.0*a3*X0^2*b5*w1*v2-1.0*a3*X3^2*b5*w2*v1+1.0*a3*X0^2*b2*w2*w3+1.0*a3*X1^2*b1*w1*w3-1.0*a3*X1^2*b2*w2*w3-1.0*a3*X3^2*b1*w1*w3...
         -2*a3*X1*b2*w2d*X2+2*a3*X1*b1*w1d*X2-2.0*a3*X1*b1*X3*w2^2-2*a3*Y2*X3*b5*w2d-2*a3*Y2*X0*b5*w1d-2.0*a3*X1*b1*X3*w1^2 ...
        -2*a3*X3*Y1*b5*w1d+2*a3*Y1*X0*b5*w2d-2.0*a3*X2*Y1*b5*w2^2+2*a3*X3*b1*w1d*X0+2.0*a3*X2*b1*X0*w1^2+2.0*a3*X2*b1*X0*w2^2 ...
         +2*a3*X3*b3*w1d*X2-2*a3*X2*b3*w2d*X0+2*a3*X3*b2*w2d*X0-2.0*a3*X3*b2*X2*w1^2-2.0*a3*X3*b2*X2*w2^2-2*a3*X1*Y0*b5*w2d+2*a3*X1*Y3*b5*w1d...
         +2.0*a3*X1*Y2*b5*w1^2+2.0*a3*X1*Y2*b5*w2^2+2*a2*X3*Y2*b5*w3d-2.0*a2*X1*Y3*b5*w1^2-2.0*a2*X3*b3*X2*w3^2-2*a2*X3*b2*w1d*X2...
         -2*a2*X3*b2*w3d*X0-2*a2*X2*Y1*b5*w1d-2*a2*X0*Y1*b5*w3d-2.0*a2*X3*b3*X2*w1^2-2.0*a2*X1*Y3*b5*w3^2+2*a2*X1*Y0*b5*w3d...
        +2.0*a2*X3*Y1*b5*w3^2+2.0*a2*X3*Y1*b5*w1^2-2*a2*X3*Y0*b5*w1d-2.0*a2*Y0*X2*w3^2*b5-2.0*a2*Y0*X2*b5*w1^2-2*a2*Y3*X2*b5*w3d...
        +2*a2*X1*b2*w3d*X2-2.0*a2*X1*b1*X2*w1^2-2.0*a2*X1*b1*X2*w3^2+2*a2*Y3*X0*b5*w1d-1.0*a2*w1*X1^2*b5*v3-1.0*a2*X0^2*b5*w1*v3...
        +1.0*a2*X2^2*b5*w3*v1+1.0*a2*X0^2*b5*v1*w3-2.0*a5*b5*Y0*v3d*X3-2.0*a5*b5*Y0*v1d*X1-2.0*a5*b5*Y0*v2d*X2-2.0*a5*X2*X1*b1*v2d...
        -2.0*a5*X2*b2*v1d*X1+2.0*a5*Y3*X2*b5*v1d+2.0*a5*Y3*X0*b5*v3d-2.0*a5*Y3*X1*b5*v2d+1.00*a5*b2*w3*X2^2*v1+1.00*a5*b3*w2*X2^2*v1...
        +1.00*a5*b1*w3*X3^2*v2-1.00*a5*b3*w1*X1^2*v2-1.00*a5*b3*w2*X3^2*v1-2.0*a3*X2*Y1*b5*w1^2-2.0*a3*X1*b2*X0*w2^2-2*a3*X1*b3*w2d*X3...
       -2*a3*X1*b3*w1d*X0-2.0*a3*X3*Y0*b5*w1^2-2.0*a3*X1*b2*X0*w1^2+2*a3*X2*Y3*b5*w2d+2*a3*X2*Y0*b5*w1d-1.0*a3*w2*X1^2*b5*v1...
        +1.0*a3*w1*X2^2*b5*v2-1.0*a3*w2*X2^2*b5*v1+1.0*a3*w1*X1^2*b5*v2+1.0*a3*X0^2*b1*w1*w3-1.0*a3*X2^2*b1*w1*w3+1.0*a3*X2^2*b2*w2*w3...
        -2*a1*X1*Y2*b5*w2d-1.0*a1*X1^2*b5*w3*v2+2.0*a1*X2*Y3*b5*w3^2-2.0*a1*X1*Y0*b5*w3^2+1.0*a1*w2*X2^2*b5*v3+1.0*a1*X3^2*b5*w2*v3...
        +2.0*a1*X3*b2*X0*w3^2-2*a1*X3*b1*w3d*X0+2.0*a1*X2*Y3*b5*w2^2-2.0*a1*X2*b3*X0*w2^2+2*a1*w1*X2*b5*w3*Y1+2*a2*w1*X1*b5*w2*Y0...
        +2.0*a5*b3*w3*X2*v1*X3-2.0*a5*b2*w2*X0*v1*X1-2.0*a5*b3*w3*X0*v1*X1-2.0*a5*X1*Y3*b5*v3*w1-2.0*a5*X3*Y0*b5*v1*w2-2.0*a5*X3*Y1*b5*v1*w3...
        -2.0*a5*X0*Y2*b5*v1*w3+2.0*a5*X0*Y3*b5*v1*w2-2*a3*X0*b5*w1*w3*Y1-2*a3*X2*b5*w3*w1*Y3+2.0*a5*X1*Y3*b5*v1*w3+2.0*a5*X2*Y1*b5*v2*w1...
       +2.0*a5*X2*Y3*b5*v2*w3+2.0*a5*b3*w1*X0*v3*X1-2.0*a5*b3*w3*X3*v2*X1+2.0*a5*b1*w3*X0*v1*X3+2.0*a5*b3*w2*X3*v3*X1-2.0*a5*b1*w2*X1*v1*X3...
       +2.0*a5*X3*Y0*b5*v2*w1-2.0*a5*X3*Y2*b5*v2*w3-2.0*a5*X1*Y0*b5*v2*w3-2.0*a5*X1*Y2*b5*v2*w1-2.0*a5*b2*w3*X2*v2*X1-2.0*a5*b1*w1*X3*v3*X0...
        +2.0*a5*b2*w3*X3*v2*X0-1.0*a2*X2^2*b5*v3*w1+1.0*a2*X1^2*b1*w1*w2-1.0*a2*X1^2*b3*w2*w3-1.0*a2*X2^2*b1*w1*w2-1.0*a2*X3^2*b1*w1*w2...
        -1.0*a2*X2^2*b3*w2*w3+1.0*a2*X3^2*b3*w2*w3+1.0*a2*X0^2*b3*w2*w3+1.0*a2*X0^2*b1*w1*w2-1.0*a2*w1*X3^2*b5*v3+1.0*a2*X1^2*b5*w3*v1...
        +1.0*a2*w3*X3^2*b5*v1+2*a2*X2*b1*w1d*X0+2.0*a1*X3*b2*X0*w2^2-1.0*a1*w3*X3^2*b5*v2-1.0*a1*X1^2*b2*w1*w2-1.0*a1*X1^2*b3*w1*w3...
        +1.0*a1*w2*X1^2*b5*v3-2.0*a1*X1*b2*X2*w3^2-2.0*a1*X1*Y0*b5*w2^2+1.0*a1*X2^2*b2*w1*w2-2.0*a1*X3*b3*X1*w3^2-2.0*a1*X3*Y2*b5*w3^2 ...
        +2*a1*Y2*X0*b5*w3d-2*a1*Y3*X0*b5*w2d+2.0*a1*Y1*X0*w3^2*b5+1.0*a1*X3^2*b3*w1*w3+1.0*a1*X0^2*b3*w1*w3-1.0*a1*X2^2*b3*w1*w3...
        +2*a1*X3*b1*w2d*X1+2*a1*X3*b2*w2d*X2-1.0*a1*w3*X0^2*b5*v2+2.0*a1*Y1*X0*w2^2*b5-2.0*a1*X2*b3*X0*w3^2+2*a1*X1*b2*w2d*X0...
        -2*a1*X2*b1*w2d*X0+2.0*a3*X0*Y3*b5*w1^2+2.0*a3*X0*Y3*b5*w2^2-2.0*a3*X3*Y0*b5*w2^2-2.0*a2*X3*b1*X0*w1^2-2.0*a2*X3*b1*X0*w3^2 ...
        +2*a2*X2*b3*w3d*X0+2*a2*X1*Y2*b5*w1d+2.0*a2*Y2*X0*w1^2*b5+2.0*a2*Y2*X0*w3^2*b5-2*a2*X1*b2*w1d*X0+2*a2*X3*b3*w3d*X1...
        +2.0*a2*X1*b3*X0*w1^2+2.0*a2*X1*b3*X0*w3^2-2*a2*X3*b1*w1d*X1+1.00*a5*b1*w2*X0^2*v3-1.00*a5*b1*w3*X1^2*v2-1.00*a5*b2*w3*X1^2*v1...
        +1.00*a5*b2*w1*X3^2*v3-1.00*a5*b3*w2*X0^2*v1+1.00*a5*b1*w3*X2^2*v2-1.00*a5*b1*w2*X3^2*v3+2.0*a5*X2*X0*b1*v3d-2.0*a5*X2*X0*b3*v1d...
        -2.0*a5*X2*X3*b3*v2d-2.0*a5*X2*b2*v3d*X3-2.0*a5*Y1*X2*b5*v3d+2.0*a5*Y1*X3*b5*v2d+2.0*a5*Y1*X0*b5*v1d+2.0*a5*X3*X0*b2*v1d...
        -2.0*a5*X3*b1*v2d*X0-2.0*a5*b5*Y2*v1d*X3+2.0*a5*b5*Y2*v3d*X1+2.0*a5*b5*Y2*v2d*X0+2.0*a5*X0*b3*v2d*X1-2.0*a5*X0*b2*v3d*X1...
        -2.0*a5*X3*b1*v3d*X1-2.0*a5*X3*X1*b3*v1d-2.0*a5*X2*Y0*b5*v3*w1-2.0*a5*X2*Y3*b5*v3*w2+2.0*a5*X1*Y2*b5*v1*w2-2.0*a5*b2*w2*X3*v3*X0...
        -2.0*a5*X0*Y1*b5*v3*w2+2.0*a5*X0*Y2*b5*v3*w1+2.0*a5*X0*Y1*b5*v2*w3-2.0*a5*X0*Y3*b5*v2*w1+2.0*a5*X1*Y0*b5*v3*w2+2.0*a5*X2*Y0*b5*v1*w3...
        -2.0*a5*X2*Y1*b5*v1*w2+2.0*a5*b2*w2*X2*v3*X1-2.0*a5*b3*w3*X2*v2*X0+2.0*a5*b1*w3*X1*v1*X2-2.0*a5*b1*w1*X2*v2*X0...
        -2.0*a5*b1*w1*X2*v3*X1+2.0*a5*b2*w1*X1*v2*X0-2.0*a5*b2*w2*X2*v1*X3+2.0*a5*X3*Y1*b5*v3*w1+2.0*a5*X3*Y2*b5*v3*w2+2.0*a5*b1*w2*X0*v1*X2...
        +2.0*a5*b3*w2*X2*v3*X0+2.0*a5*b1*w1*X3*v2*X1+2*a1*X1*X2*b1*w1*w2+2*a1*X0*X1*b2*w1*w3+2*a1*X2*X3*b3*w1*w2-2*a1*X0*X2*b1*w1*w3...
        -2*a1*w1*X3*b5*w2*Y1-2*a1*X1*b5*w3*w1*Y2+2*a1*w2*X1*b5*w1*Y3-2*a1*X0*b5*w2*w1*Y2+2*a1*X0*X3*b1*w1*w2+2*a1*w1*X3*b5*w3*Y0...
        +2*a1*w1*X2*b5*w2*Y0-2*a1*w3*X0*b5*w1*Y3-2*a2*X0*b5*w1*w2*Y1-1.00*a5*X0^2*b1*v1d-1.00*a5*X2^2*b5*v1^2-1.00*a5*X0^2*b5*v3^2 ...
        -1.00*a5*X1^2*b5*v2^2+1.00*a5*X3^2*b2*v2d-1.00*a5*X1^2*b5*v3^2-1.00*a5*X3^2*b5*v1^2+1.00*a5*X3^2*b1*v1d-1.00*a5*X0^2*b2*v2d...
        +1.00*a5*X1^2*b3*v3d-1.00*a5*X1^2*b5*v1^2-1.00*a5*X2^2*b5*v2^2-1.00*a5*X3^2*b5*v3^2-1.00*a5*X0^2*b5*v1^2-1.00*a5*X2^2*b5*v3^2 ...
        -1.00*a5*X3^2*b5*v2^2+1.00*a5*X1^2*b2*v2d+1.00*a5*X2^2*b1*v1d-1.00*a5*X0^2*b3*v3d-1.00*a5*X3^2*b3*v3d-1.00*a5*X2^2*b2*v2d...
        -1.00*a5*X1^2*b1*v1d-1.00*a5*X0^2*b5*v2^2+1.00*a5*X2^2*b3*v3d+1.00*a3*X1^2*b3*w2^2+1.0*a3*X3^2*b1*w2d+1.0*a3*X3^2*b5*v3d...
        -1.0*a3*X0^2*b1*w2d+1.00*a3*X2^2*b3*w1^2+1.00*a3*X2^2*b3*w2^2-1.0*a3*X3^2*b2*w1d-1.00*a3*X0^2*b3*w2^2-1.00*a3*X0^2*b3*w1^2 ...
        +1.0*a1*X3^2*b2*w3d-1.0*a1*X2^2*b2*w3d+1.0*a1*X0^2*b3*w2d+1.00*a1*X2^2*b1*w3^2+1.00*a1*X3^2*b1*w3^2-1.0*a1*X1^2*b3*w2d...
        +1.00*a1*X3^2*b1*w2^2-1.0*a1*X0^2*b2*w3d+1.0*a3*X0^2*b5*v3d+1.0*a3*X0^2*b2*w1d+1.00*a3*X1^2*b3*w1^2+1.0*a3*X2^2*b1*w2d...
        -1.0*a3*X1^2*b2*w1d-1.0*a3*X1^2*b1*w2d-1.00*a3*X3^2*b3*w1^2-1.00*a3*X3^2*b3*w2^2+1.0*a3*X2^2*b2*w1d+1.0*a3*X2^2*b5*v3d...
        +1.0*a3*b5*X1^2*v3d+1.0*a2*b5*X1^2*v2d+1.0*a2*X0^2*b1*w3d+1.0*a2*X2^2*b3*w1d-1.0*a2*X0^2*b3*w1d+1.00*a2*X1^2*b2*w1^2 ...
        +1.00*a2*X1^2*b2*w3^2-1.0*a2*X2^2*b1*w3d+1.0*a2*X2^2*b5*v2d+1.0*a2*b5*X3^2*v2d-1.00*a2*X0^2*b2*w1^2-1.00*a2*X0^2*b2*w3^2 ...
        +1.0*a2*X1^2*b3*w1d-1.0*a2*X3^2*b1*w3d+1.0*a2*X1^2*b1*w3d-1.00*a2*X2^2*b2*w1^2-1.0*a2*X3^2*b3*w1d+1.00*a2*X3^2*b2*w3^2 ...
        +1.00*a2*X3^2*b2*w1^2+1.0*a2*X0^2*b5*v2d-1.00*a2*X2^2*b2*w3^2+1.0*a1*X3^2*b3*w2d-1.00*a1*X1^2*b1*w2^2-1.00*a1*X1^2*b1*w3^2 ...
        +1.0*a1*b5*X2^2*v1d+1.0*a1*b5*X3^2*v1d-1.00*a1*X0^2*b1*w2^2-1.00*a1*X0^2*b1*w3^2+1.0*a1*b5*X1^2*v1d+1.0*a1*X1^2*b2*w3d...
        +1.00*a1*X2^2*b1*w2^2+1.0*a1*X0^2*b5*v1d-1.0*a1*X2^2*b3*w2d));
end


%calculate actuator length,speed,acceleration
%DB   = zeros(8,N);
b    = zeros(3,6,N);
L    = zeros(6,N);
Ld   = zeros(6,N);
Ldd  = zeros(6,N);

for i = 1:N

    %transform rotation and translation into dual quaternion  
   % DB(:,i)    = qToDualQuat(q(:,i));
    % current position "r" of platform frame
    %r          = [ q(1,i) , q(2,i) , q(3,i)]';
     X0 = DB(1,i);
     X1 = DB(2,i);
     X2 = DB(3,i);
     X3 = DB(4,i);
     Y0 = DB(5,i);
     Y1 = DB(6,i);
     Y2 = DB(7,i);
     Y3 = DB(8,i);
     
    r = [2*(X1*Y0-Y1*X0-X2*Y3+Y2*X3),2*(X1*Y3-Y1*X3+X2*Y0-Y2*X0),2*(-X1*Y2+Y1*X2+X3*Y0-Y3*X0)]'/(X1^2+X2^2+X3^2+X0^2); %platform frame position
    %body frame velocity to inertia frame velocity 
    VB(1:3,1)  = [ qd(1,i) , qd(2,i) , qd(3,i) ]' - cross( [qd(4,i) qd(5,i) qd(6,i)]' , r );
    VB(4:6,1)  = qd(4:6,i);
    %body frame acceleration to inertia frame acceleration 
    VBd(1:3,1) = [qdd(1,i),qdd(2,i),qdd(3,i)]'-(cross([qdd(4,i) qdd(5,i) qdd(6,i)]',r)+cross([qd(4,i) qd(5,i) qd(6,i)]',[qd(1,i) qd(2,i) qd(3,i)]'));
    VBd(4:6,1) = qdd(4:6,i);
    
    %length,speed and acceleration of linear actuators
    for j = 1:6
    b(:,j,i)    = sandwich(DB(:,i),bp(:,j));           %platform joint location "b" from the dual quaternion transformation
    
    b4          = 0.5 * (b(1,j,i)^2 + b(2,j,i)^2 + b(3,j,i)^2); %CGA representation "b4" from b
    
    L(j,i)      = sqrt(2*abs(a(1,j)*b(1,j,i) + a(2,j)*b(2,j,i) + a(3,j)*b(3,j,i) - ((a(1,j)^2 + a(2,j)^2 + a(3,j)^2)/2) - b4));     %length
    Ld(j,i)     = LD_cga(DB(:,i),VB,L(j,i),a(:,j),bp(:,j));                                                                         %speed
    Ldd(j,i)    = LDD_cga(DB(:,i),VB,VBd,L(j,i),Ld(j,i),a(:,j),bp(:,j));                                                            %acceleration
    
    end


end
%% Animation figure settings
 if sim == true,
%clf
h = figure;
%scrsz = get(0, 'ScreenSize'); 
%h=figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
%h.Position = [380 200 1000 800];
axis equal, grid on
view(25,20)
%view(90,0)
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2,2]*1.5);
ylim([-2,2]*1.5);
zlim([0,3]);



tStart = cputime();
tFactor = 1;
frame = 0;
while ishandle(h)

    % Calculate elapsed time
    tElapsed = (cputime() - tStart)*tFactor;
    
    [~,i] = min(abs(t-tElapsed));
    title(['t = ',num2str(t(i),'%.2f'),' s']);
    
    
    % Platfrom coordiante
%     H01(1:3,1:3) = Axyz(phi(:,i));
%     H01(1:3,4) = r(:,i);
    

    H = quat2rotm(DB(:,i));
    
    % Transform objets
    hold on
        cla
        % Base joints
        plot3(a(1,:),a(2,:),a(3,:),'*r')
        
        % platform joints
        plot3(b(1,:,i),b(2,:,i),b(3,:,i),'*b')
      
        %platform frame
        for k=1:6
        line([a(1,k),b(1,k,i)],[a(2,k),b(2,k,i)],[a(3,k),b(3,k,i)],'LineWidth',2) 

        end
        
        for k=1:6
      
            if k ==6
        line([b(1,k,i),b(1,1,i)],[b(2,k,i),b(2,1,i)],[b(3,k,i),b(3,1,i)],'Color','red','LineWidth',2)  
            else
        line([b(1,k,i),b(1,k+1,i)],[b(2,k,i),b(2,k+1,i)],[b(3,k,i),b(3,k+1,i),],'Color','red','LineWidth',2) 
            end
        end
        
        % Platform coordinate system
        PlotAxes(H,0.3,2)
       
        % Base coordinate system
        %PlotAxes(eye(4,4),0.3,2)
     
        
    hold off

    if i == length(t)
        break;
    end
    
    drawnow;
    frame = frame + 1;
end

fprintf('FPS = %f \n',frame/t(i));
 end
end