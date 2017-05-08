function [D,Dd,Ddd,Dvel,Dacc,q_cga,Vb,Vbd] = StewartFK_CGA(q0,L,Ld,Ldd,dt,t,sim,b_temp)
%function [D,Dd,Ddd,Dvel,Dacc,q_cga,Vb,Vbd] = StewartFK_CGA(q0,L,Ld,Ldd,dt,t,sim,b_temp)
% 
% This function is a forward kinematic algortithm for a 6-SPS stewart
% platform. It applies Conformal Geometric Algebra for defining constaint
% equation. From given leg-length, leg-speed and leg-acceleration, the pose
% velocity and acceleration of the platform frame is returned.
% 
% Created by Olav Heng.
% INPUT:
% q0        : Initial pose components [x y z phi theta psi]' for the stewart platfrom, size:(6,1)
% L         : Leg-length of stewart platform, size:(6,n)
% Ld        : Leg extraction speed, size:(6,n)
% Ldd       : Leg extraction acceleration, size:(6,n)
% dt        : Required for integration of Dd and Ddd, to obtain Dvel and Dacc, set 0 if not required.
% t         : Time for each step n, required for simulation, size:(1,n), set 0 if not used.
% sim       : Simulation on/off: 1/0
% b_temp    : Template used for verification in simulation plot. size:(3,6,n), set 0 if not used. 
%
% OUTPUT
% D         : Pose of the platform returned as displacement versor, size:(8,n)
% Dd        : The first derivative displacement versor, size:(8,n)
% Ddd       : The second derivative displacement versor, size:(8,n)
% Dvel      : Displacement versor derived from integration of Dd, size:(8,n)
% Dacc      : Displacement versor derived from double integration of Ddd, size:(8,n)
% q_cga     : The pose components [x y z phi theta psi]', derived from D, size:(6,n)
% Vb        : The velocity components of [x y z phi theta psi]' for the stewart platfrom body frame, size:(6,n)
% Vbd       : The acceleration components of [x y z phi theta psi]' for the stewart platfrom body frame, size:(6,n)

global a bp
N = length(L); 

%Stewart parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jointLoc = [15,105,135,225,255,345] / 180 * pi;         %location of the base points "a"
a0       = [2,0,0]';

%rotation around z-axis for distribution of joints.
Az = @(psi)[cos(psi),-sin(psi), 0;
            sin(psi), cos(psi), 0;
            0       , 0       , 1];
        
a   = zeros(3,6); %base joints
bp  = zeros(3,6);  
%distribute spherical joints "a" in the base, and "bp" on the platform
for i = 1:length(jointLoc)
    
    a(:,i) = Az(jointLoc(i)) * a0;
    
    if i == 6                  %shifting joint position in bp so that [a1,a2...] and [bp1,bp2...] corrospond to each other
    bp(:,1) = Az(60/180*pi) * a(:,i);
    else
    bp(:,i+1) = Az(60/180*pi) * a(:,i);
    end
    
end

%"a4" CGA representation of sphere with radius "L" and center "a"
s4 = zeros(6,N);
for i = 1 : N
        for j = 1:6
            s4(j,i) = (a(1,j)^2 + a(2,j)^2 + a(3,j)^2) / 2 - 1/2 * L(j,i)^2;
        end
end

%define sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function D = qToDualQuat(q) %rotational and translational vector to dual quaternions (displacement versor)

    x = q(1,1);
    y = q(2,1);
    z = q(3,1);
    roll  = q(4,1);
    pitch = q(5,1);
    yaw   = q(6,1);


    c1 = cos(roll/2);
    c2 = cos(pitch/2);
    c3 = cos(yaw/2);
    s1 = sin(roll/2);
    s2 = sin(pitch/2);
    s3 = sin(yaw/2);

    X0 = c1*c2*c3+s1*s2*s3;
    X1 = s1*c2*c3+c1*s2*s3;
    X2 = c1*s2*c3-s1*c2*s3; 
    X3 = c1*c2*s3+s1*s2*c3;

     Y0 = 0.5*(X3*z+X2*y+X1*x);
     Y1 = -0.5*(X0*x-X2*z+X3*y);
     Y2 = -0.5*(X0*y+X1*z-X3*x);
     Y3 = -0.5*(X0*z-X1*y+X2*x);


    D = ([X0 X1 X2 X3 Y0 Y1 Y2 Y3]');%/sqrt(X0^2+X1^2+X2^2+X3^2)+2*(X0*Y0+X1*Y1+X2*Y2+X3*Y3));
end

function E = rotm2eulerXYZ(H)

H = H./H(4,4);

if H(1,3) < 1
    
    if H(1,3) > -1

        thetaY = asin(H(1,3));
        thetaX = atan2(-H(2,3),H(3,3));
        thetaZ = atan2(-H(1,2),H(1,1));
        
    else
        
        thetaY = -pi/2;
        thetaX = -atan2(-H(2,1),H(2,2));
        thetaZ = 0;
        
    end
else
        thetaY = pi/2;
        thetaX = atan2(H(2,1),H(2,2));
        thetaZ = 0;
end

E =[H(1,4) H(2,4) H(3,4) thetaX thetaY thetaZ ]';
end

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

function Dn = normalizeD(D) %normalizing the dual quaternion
     Dabs = sqrt(D(1,1)^2 + D(2,1)^2 + D(3,1)^2 + D(4,1)^2) ;
     Dn   = [D(1,1) D(2,1) D(3,1) D(4,1) D(5,1) D(6,1) D(7,1) D(8,1)]' / Dabs;
end

function J = Jacobi(D) %jacobian matrix

    % Rename states
    X0 = D(1);
    X1 = D(2);
    X2 = D(3);
    X3 = D(4);
    Y0 = D(5);
    Y1 = D(6);
    Y2 = D(7);
    Y3 = D(8);

    J = zeros(6,6);
        for joint = 1:6

            b1 = bp(1,joint);
            b2 = bp(2,joint);
            b3 = bp(3,joint);
            b5 = 1;

            a1 = a(1,joint);
            a2 = a(2,joint);
            a3 = a(3,joint);
            a5 = 1;

         J(joint,1:6) = [2.0*b5*Y3*X2*a5-2.0*X1*b2*X2*a5+2.0*b5*Y1*X0*a5+X3^2*b5*a1+X2^2*b5*a1+b5*X0^2*a1+X2^2*b1*a5+X3^2*b1*a5+X1^2*b5*a1-X1^2*b1*a5+2.0*X0*b2*X3*a5-2.0*Y2*X3*b5*a5-2.0*Y0*X1*b5*a5-X0^2*b1*a5-2.0*X0*b3*X2*a5-2.0*X1*b3*X3*a5,...
                        -2.0*X1*X2*b1*a5-2.0*b5*Y3*X1*a5+2.0*X1*b3*X0*a5-2.0*Y0*X2*b5*a5-2.0*X3*b3*X2*a5-X2^2*b2*a5+X1^2*b2*a5+2.0*Y2*X0*b5*a5-2.0*X0*b1*X3*a5+X3^2*b2*a5+X2^2*b5*a2+b5*X3^2*a2+X1^2*b5*a2+b5*X0^2*a2-X0^2*b2*a5+2.0*b5*Y1*X3*a5,...
                        2.0*X0*b1*X2*a5-2.0*X1*b2*X0*a5+2.0*Y2*X1*b5*a5-2.0*X3*b2*X2*a5-2.0*b5*Y1*X2*a5-2.0*Y0*X3*b5*a5+X0^2*b5*a3+X1^2*b5*a3-X3^2*b3*a5+X2^2*b3*a5+X2^2*b5*a3+X1^2*b3*a5+2.0*b5*Y3*X0*a5-2.0*X1*X3*b1*a5+X3^2*b5*a3-X0^2*b3*a5,...
                        -2.*X1*b2*X0*a2-2.*Y1*X2*b5*a2-2.*X0*b3*X1*a3-2.*Y1*X3*b5*a3+2.*Y2*X1*b5*a2-2.*X0*Y2*b5*a3+2.*X2*b1*X1*a3+2.*Y3*X0*b5*a2-2.*X3*b1*X1*a2-X1^2*b2*a3+X0^2*b2*a3+X2^2*b2*a3-X3^2*b2*a3-X0^2*b3*a2+2.*X3*X0*b1*a3+2.*X2*b3*X3*a3+2.*X2*Y0*b5*a3-2.*X3*b2*X2*a2+X2^2*b3*a2+X1^2*b3*a2-X3^2*b3*a2+2.*X2*X0*b1*a2+2.*Y3*X1*b5*a3-2.*X3*Y0*b5*a2,...
                        2.*Y0*X3*b5*a1-2.*X3*Y2*b5*a3+2.*Y1*X0*b5*a3+2.*X0*b2*X1*a1-X2^2*b3*a1+X3^2*b3*a1+X3^2*b1*a3-X0^2*b1*a3-X1^2*b3*a1-X1^2*b1*a3+X2^2*b1*a3+X0^2*b3*a1-2.*X2*b3*X0*a3+2.*X3*b2*X2*a1-2.*Y3*X0*b5*a1+2.*Y3*X2*b5*a3-2.*Y0*X1*b5*a3-2.*X3*b3*X1*a3-2.*X2*b1*X0*a1-2.*X2*b2*X1*a3-2.*X1*Y2*b5*a1+2.*X2*Y1*b5*a1+2.*X3*X0*b2*a3+2.*X3*b1*X1*a1,...
                        -2.*X2*Y3*b5*a2-2.*X2*b1*X1*a1-2.*Y1*X0*b5*a2-X2^2*b2*a1+X3^2*b2*a1-X0^2*b2*a1-X3^2*b1*a2-2.*X3*b2*X0*a2+2.*Y0*X1*b5*a2+2.*X3*b3*X1*a2-2.*X3*b1*X0*a1+2.*X2*X1*b2*a2+2.*X0*b3*X1*a1+2.*X3*Y2*b5*a2+X0^2*b1*a2-X2^2*b1*a2+X1^2*b2*a1+X1^2*b1*a2-2.*Y3*X1*b5*a1+2.*X3*Y1*b5*a1+2.*X2*b3*X0*a2+2.*X0*Y2*b5*a1-2.*Y0*X2*b5*a1-2.*X3*b3*X2*a1];
        end
end

function Dd = Ddot_b(V,D) %First derivative from velocity Vb.

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

    X0d =-.5*(X1*w1+X2*w2+X3*w3); %Id
    X1d = .5*(X0*w1+X2*w3+X3*w2); %e23 
    X2d = .5*(X0*w2-X1*w3-X3*w1); %e13  
    X3d = .5*(X0*w3+X1*w2+X2*w1); %e12 


    Y0d = .5*(X1*v1+X2*v2+X3*v3-Y1*w1-Y2*w2-Y3*w3); %e1^e2^e3^einf
    Y1d =-.5*(X0*v1-X2*v3+X3*v2+Y0*w1-Y2*w3+Y3*w2); %e1^einf
    Y2d =-.5*(X0*v2+X1*v3-X3*v1+Y0*w2+Y1*w3-Y3*w1); %e2^einf
    Y3d =-.5*(X0*v3-X1*v2+X2*v1-Y0*w3-Y1*w2+Y2*w1); %e3^einf

Dd = [X0d X1d X2d X3d Y0d Y1d Y2d Y3d]';
end

function Ddd = Dddot_b(V,Vd,D,Dd) %Second derivative of D from velocity Vb, and acceleration Vbd

    vd1 = Vd(1,1);
    vd2 = Vd(2,1);
    vd3 = Vd(3,1);
    wd1 = Vd(4,1);
    wd2 = Vd(5,1);
    wd3 = Vd(6,1);

     Xd0 = Dd(1,1);
     Xd1 = Dd(2,1);
     Xd2 = Dd(3,1);
     Xd3 = Dd(4,1);
     Yd0 = Dd(5,1);
     Yd1 = Dd(6,1);
     Yd2 = Dd(7,1);
     Yd3 = Dd(8,1);

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

    X0dd =-.5*(Xd1*w1+Xd2*w2+Xd3*w3)-.5*(X1*wd1+X2*wd2+X3*wd3); %Id
    X1dd = .5*(Xd0*w1+Xd2*w3+Xd3*w2)+.5*(X0*wd1+X2*wd3+X3*wd2); %e23 
    X2dd = .5*(Xd0*w2-Xd1*w3-Xd3*w1)+.5*(X0*wd2-X1*wd3-X3*wd1); %e13  
    X3dd = .5*(Xd0*w3+Xd1*w2+Xd2*w1)+.5*(X0*wd3+X1*wd2+X2*wd1); %e12 

    Y0dd = .5*(Xd1*v1+Xd2*v2+Xd3*v3-Yd1*w1-Yd2*w2-Yd3*w3)+.5*(X1*vd1+X2*vd2+X3*vd3-Y1*wd1-Y2*wd2-Y3*wd3); %e1^e2^e3^einf
    Y1dd =-.5*(Xd0*v1-Xd2*v3+Xd3*v2+Yd0*w1-Yd2*w3+Yd3*w2)-.5*(X0*vd1-X2*vd3+X3*vd2+Y0*wd1-Y2*wd3+Y3*wd2); %e1^einf
    Y2dd =-.5*(Xd0*v2+Xd1*v3-Xd3*v1+Yd0*w2+Yd1*w3-Yd3*w1)-.5*(X0*vd2+X1*vd3-X3*vd1+Y0*wd2+Y1*wd3-Y3*wd1); %e2^einf
    Y3dd =-.5*(Xd0*v3-Xd1*v2+Xd2*v1-Yd0*w3-Yd1*w2+Yd2*w1)-.5*(X0*vd3-X1*vd2+X2*vd1-Y0*wd3-Y1*wd2+Y2*wd1); %e3^einf

    Ddd = [X0dd X1dd X2dd X3dd Y0dd Y1dd Y2dd Y3dd]';
end

function Fi = Phi(D,S4) %constraint equations Phi for the interative fsolver.

% Rename states
    X0 = D(1,1);
    X1 = D(2,1);
    X2 = D(3,1);
    X3 = D(4,1);
    Y0 = D(5,1);
    Y1 = D(6,1);
    Y2 = D(7,1);
    Y3 = D(8,1);

    Fi = zeros(8,1);

    p5      = X0^2 + X1^2 + X2^2 + X3^2-1;
    p6      = 2.*(X0*Y0 + X1*Y1 + Y2*X2 + Y3*X3);
    Fi(7,1) = p6;
    Fi(8,1) = p5;

    b5      = 1;
        for joint = 1:6

                b1 = bp(1,joint);
                b2 = bp(2,joint);
                b3 = bp(3,joint);
                b4 = (bp(1,joint)^2+bp(2,joint)^2+bp(3,joint)^2)/2;
                a1 = a(1,joint);
                a2 = a(2,joint);
                a3 = a(3,joint);
                a4 = S4(joint,1);
                s5 = 1;

            Fi(joint,1) = -1.0*X3^2*b4*s5-2.0*Y3^2*b5*s5-1.0*X3^2*b5*a4-X1^2*b3*a3-1.0*X0^2*b4*s5-X1^2*b2*a2+X1^2*b1*a1-1.0*X1^2*b5*a4-1.0*X1^2*b4*s5-1.0*X2^2*b5*a4...
                          -1.0*X2^2*b4*s5-X3^2*b2*a2-X3^2*b1*a1+X3^2*b3*a3-2.0*Y1^2*b5*s5-2.0*Y0^2*b5*s5-2.0*Y2^2*b5*s5-X2^2*b3*a3-X2^2*b1*a1+X2^2*b2*a2+X0^2*b1*a1...
                          +X0^2*b3*a3+X0^2*b2*a2-1.0*X0^2*b5*a4+2.*X3*X1*b1*a3+2.*X3*X2*b2*a3-2.0*X3*Y0*b3*s5-2.*X1*Y2*b5*a3-2.*X0*Y1*b5*a1+2.0*X0*Y1*b1*s5+2.0*Y1*X2*b3*s5...
                          -2.0*Y1*X3*b2*s5+2.*X1*Y0*b5*a1+2.*X3*Y0*b5*a3+2.*X2*Y0*b5*a2-2.*X3*Y1*b5*a2+2.*X2*Y1*b5*a3+2.0*X0*Y3*b3*s5+2.0*Y3*X1*b2*s5-2.0*Y3*X2*b1*s5...
                          -2.0*X2*Y0*b2*s5+2.*X1*X2*b2*a1-2.0*X1*Y0*b1*s5+2.*X3*X2*b3*a2+2.0*X0*Y2*b2*s5-2.0*Y2*X1*b3*s5+2.0*Y2*X3*b1*s5+2.*X3*Y2*b5*a1+2.*X1*Y3*b5*a2...
                          -2.*X0*Y3*b5*a3-2.*X2*Y3*b5*a1+2.*X3*X1*b3*a1+2.*X0*X3*b1*a2+2.*X1*X2*b1*a2+2.*X0*X1*b2*a3-2.*X0*X1*b3*a2-2.*X0*X2*b1*a3+2.*X0*X2*b3*a1-2.*X0*Y2*b5*a2...
                          -2.*X0*X3*b2*a1;

        end
end

function Vd = Phi_acc(D,V,L,Ld,Ldd) %rigth hand side acceleration constaints

    v1 = V(1,1);
    v2 = V(2,1);
    v3 = V(3,1);
    w1 = V(4,1);
    w2 = V(5,1);
    w3 = V(6,1);

    X0 = D(1);
    X1 = D(2);
    X2 = D(3);
    X3 = D(4);
    Y0 = D(5);
    Y1 = D(6);
    Y2 = D(7);
    Y3 = D(8);

    Vd = zeros(6,1);
    Id = (-X2^2-X3^2-X1^2-X0^2);

    for joint = 1:6
        b1 = bp(1,joint);
        b2 = bp(2,joint);
        b3 = bp(3,joint);
        b5 = 1;

        a1 = a(1,joint);
        a2 = a(2,joint);
        a3 = a(3,joint);
        a5 = 1;

        s4dd = Ld(joint)^2+L(joint)*Ldd(joint);

        PddS = -2*w3*X0*b5*w1*Y1*a3-2*w2*X0*b5*w3*Y2*a3+2*w3*X1*b5*w1*Y0*a3+2*w3*X1*b5*w2*Y3*a3+2*w1*X3*b5*w3*Y2*a3...
            -2*X0*X1*b3*w2*w3*a3+2*X1*X3*b3*w1*w3*a3+2*w2*Y0*b5*w3*X2*a3-2*w3*Y1*b5*w2*X3*a3-2*w1*X2*b5*w3*Y3*a3...
            +2.0*X2*Y1*b5*v2*w1*a5+2.0*X2*Y3*b5*v2*w3*a5+2.0*b1*w3*X2*v1*X1*a5-2.0*b2*w2*X1*v1*X0*a5-2.0*b3*w3*X1*v1*X0*a5...
            +2.0*b2*w1*X2*v2*X3*a5+2.0*b2*w3*X0*v2*X3*a5-2.0*X1*Y2*b5*v2*w1*a5+2.0*X2*Y0*b5*v1*w3*a5-2.0*X2*Y1*b5*v1*w2*a5...
            -2.0*b2*w2*X0*v3*X3*a5-2.0*b2*w3*X1*v2*X2*a5+2*X1*X3*b3*w1*w2*a2-2*w1*X0*b5*w2*Y1*a2+2*w3*X2*b5*w2*Y1*a2...
            +2*w2*X3*b5*w3*Y0*a2+2*w2*X3*b5*w1*Y2*a2-2*w3*X0*b5*w2*Y3*a2-2.0*X0*Y3*b5*v2*w1*a5+2.0*X1*Y0*b5*v3*w2*a5....
            -2.0*X1*Y3*b5*v3*w1*a5-2.0*X3*Y0*b5*v1*w2*a5-2.0*X3*Y1*b5*v1*w3*a5+2*w2*X1*b5*w1*Y0*a2-2.0*X3*b3*X1*w2^2*a1...
            -X2^2*b3*w2*w3*a2-X2^2*b3*w1*w3*a1-2.0*X2*X1*b1*w3^2*a2-X2^2*b1*w1*w2*a2-w1*X1^2*b5*v3*a2+w1*X2^2*b5*v2*a3...
            -2.0*Y0*X1*b5*w2^2*a1+2.0*Y3*X0*w1^2*b5*a3-2.0*Y0*X2*b5*w1^2*a2+w2*X0^2*b5*v3*a1-b1*w3*X1^2*v2*a5+w1*X1^2*b5*v2*a3...
            -2.0*X2*X1*b1*w1^2*a2-w1*X3^2*b5*v3*a2+2.0*Y1*X0*w2^2*b5*a1+X2^2*b2*w2*w3*a3-b2*w3*X1^2*v1*a5-2.0*X3*b1*X0*w1^2*a2...
            +w3*X1^2*b5*v1*a2-X3^2*b2*w1*w2*a1+2*w2*X2*b5*w1*Y0*a1-2*w3*X1*b5*w1*Y2*a1-2*w1*X3*b5*w2*Y1*a1-2*w2*X0*b5*w1*Y2*a1...
            -2*w1*X0*b5*w3*Y3*a1+2.0*X3*Y1*b5*v3*w1*a5+2.0*X3*Y2*b5*v3*w2*a5-2.0*X0*Y1*b5*v3*w2*a5+2.0*X0*Y2*b5*v3*w1*a5...
            -2.0*X1*Y0*b5*v2*w3*a5+b1*w3*X2^2*v2*a5-2.0*X3*Y2*b5*w2^2*a1+2.0*Y3*X2*b5*w2^2*a1-2.0*X2*b2*X3*w2^2*a3...
            +w2*X2^2*b5*v3*a1+X1^2*b1*w1*w2*a2-w1*X2^2*b5*v3*a2+b1*w2*X0^2*v3*a5-2.0*X3*X1*b1*w1^2*a3-w2*X2^2*b5*v1*a3...
            -b2*w1*X2^2*v3*a5-b2*v1*X3^2*w3*a5+2.0*Y2*X0*w3^2*b5*a2-2.0*X3*X1*b1*w2^2*a3+X2^2*b2*w1*w2*a1-2.0*X0*b2*X1*w1^2*a3...
            -X3^2*b1*w1*w2*a2-2.0*Y0*X1*b5*w3^2*a1-b3*w2*X0^2*v1*a5+X0^2*b3*w2*w3*a2+w3*X3^2*b5*v1*a2-b2*w1*X0^2*v3*a5...
            -w2*X1^2*b5*v1*a3-2.0*X2*b2*X1*w3^2*a1+2.0*b1*w3*X3*v1*X0*a5-2.0*X0*Y2*b5*v1*w3*a5+2*X0*X1*b2*w1*w3*a1...
            +2.0*b3*w2*X1*v3*X3*a5-2.0*b1*w2*X3*v1*X1*a5+2.0*X1*Y3*b5*v1*w3*a5+2.0*b2*w1*X0*v2*X1*a5+2.0*b3*w1*X1*v3*X0*a5...
            -2.0*b1*w1*X0*v3*X3*a5+2.0*X0*Y1*b5*v2*w3*a5+2*X0*X2*b3*w1*w3*a3+2*X0*X3*b1*w2*w3*a3-2*X0*X3*b2*w1*w3*a3...
            +2*X1*X2*b1*w2*w3*a3+2*X1*X2*b2*w1*w3*a3+2*X2*X3*b3*w2*w3*a3-2.0*b1*w1*X1*v3*X2*a5+2.0*b2*w2*X1*v3*X2*a5...
            -X0^2*b5*v3^2*a5+X2^2*b3*w1^2*a3-X1^2*b5*v3^2*a5-X3^2*b5*v1^2*a5-X2^2*b2*w3^2*a2-X0^2*b2*w1^2*a2-X0^2*b1*w3^2*a1...
            -X1^2*b5*v1^2*a5-X2^2*b5*v3^2*a5+X3^2*b1*w3^2*a1-X0^2*b1*w2^2*a1+X2^2*b3*w2^2*a3-X2^2*b5*v1^2*a5-X0^2*b5*v1^2*a5...
            -X3^2*b5*v3^2*a5+X3^2*b1*w2^2*a1-X0^2*b5*v2^2*a5-X3^2*b5*v2^2*a5-X1^2*b5*v2^2*a5-X0^2*b3*w2^2*a3+X1^2*b2*w3^2*a2...
            +X1^2*b3*w1^2*a3+X3^2*b2*w1^2*a2+X2^2*b1*w2^2*a1-X0^2*b3*w1^2*a3+X3^2*b2*w3^2*a2+X1^2*b2*w1^2*a2+X2^2*b1*w3^2*a1...
            -X2^2*b5*v2^2*a5+X1^2*b3*w2^2*a3-2.0*X3*b3*X2*w1^2*a2-X3^2*b2*w2*w3*a3-2.0*X2*b2*X3*w1^2*a3-2.0*Y0*X2*b5*w3^2*a2...
            -b1*w2*X3^2*v3*a5+2.0*Y2*X0*w1^2*b5*a2-X3^2*b1*w1*w3*a3+2.0*X1*Y2*b5*w2^2*a3+2.0*Y3*X2*b5*w3^2*a1+b1*v3*X1^2*w2*a5...
            +X0^2*b3*w1*w3*a1-2.0*Y3*X1*b5*w1^2*a2+w1*X3^2*b5*v2*a3+w2*X3^2*b5*v3*a1+w3*X2^2*b5*v1*a2-2.0*X3*Y2*b5*w3^2*a1...
            -2.0*X2*Y1*b5*w2^2*a3-w3*X1^2*b5*v2*a1+2.0*X2*b1*X0*w2^2*a3-2.0*Y0*X3*b5*w2^2*a3-X1^2*b2*w2*w3*a3+b2*w1*X3^2*v3*a5...
            -2.0*X2*b3*X0*w2^2*a1+2.0*X2*b1*X0*w1^2*a3-b3*v1*X3^2*w2*a5+2.0*Y1*X0*w3^2*b5*a1-w2*X0^2*b5*v1*a3+X0^2*b1*w1*w2*a2...
            +w2*X1^2*b5*v3*a1+w3*X0^2*b5*v1*a2+X3^2*b3*w2*w3*a2+2.0*X3*b2*X0*w2^2*a1+X3^2*b3*w1*w3*a1+2.0*X3*Y1*b5*w3^2*a2...
            +X0^2*b2*w1*w2*a1-2.0*Y3*X1*b5*w3^2*a2+b3*v2*X0^2*w1*a5-X1^2*b2*w1*w2*a1-2.0*X0*b2*X1*w2^2*a3+2.0*X1*b3*X0*w3^2*a2...
            -b1*v2*X0^2*w3*a5+2.0*X1*Y2*b5*w1^2*a3-2.0*X3*b1*X0*w3^2*a2+b3*w2*X1^2*v1*a5-2.0*X2*b2*X1*w2^2*a1...
            +2.0*X1*b3*X0*w1^2*a2+b2*w3*X0^2*v1*a5-2.0*Y0*X3*b5*w1^2*a3-2.0*X3*b3*X1*w3^2*a1-w3*X0^2*b5*v2*a1+b2*w3*X2^2*v1*a5...
            -X3^2*b3*w1^2*a3-X3^2*b3*w2^2*a3-X1^2*b1*w3^2*a1-X0^2*b2*w3^2*a2-X2^2*b2*w1^2*a2-X1^2*b1*w2^2*a1...
            -2.0*X2*b3*X0*w3^2*a1+2.0*Y3*X0*w2^2*b5*a3-b3*w1*X1^2*v2*a5-w1*X0^2*b5*v3*a2+X1^2*b1*w1*w3*a3-2.0*X3*b3*X2*w3^2*a2...
            -2.0*X2*Y1*b5*w1^2*a3+v2*X0^2*b5*w1*a3+b3*w1*X3^2*v2*a5-w3*X3^2*b5*v2*a1+b2*v3*X1^2*w1*a5-X1^2*b3*w2*w3*a2...
            +2.0*X3*Y1*b5*w1^2*a2+X0^2*b1*w1*w3*a3+X0^2*b2*w2*w3*a3-w3*X2^2*b5*v2*a1+2.0*X3*b2*X0*w3^2*a1+b3*w2*X2^2*v1*a5...
            -b3*w1*X2^2*v2*a5-X1^2*b3*w1*w3*a1-b1*w2*X2^2*v3*a5-X2^2*b1*w1*w3*a3-v1*X3^2*b5*w2*a3+b1*w3*X3^2*v2*a5...
            -2*w2*X1*b5*w3*Y2*a2-2*w1*X2*b5*w2*Y3*a2+2.0*X3*Y0*b5*v2*w1*a5-2.0*X3*Y2*b5*v2*w3*a5-2.0*b1*w1*X0*v2*X2*a5...
            -2.0*b3*w3*X0*v2*X2*a5-2.0*b2*w2*X3*v1*X2*a5-2*X0*X1*b3*w1*w2*a1-2*X0*X2*b1*w1*w3*a1+2*X0*X3*b1*w1*w2*a1...
            +2*X2*X3*b2*w1*w3*a1+2*X2*X3*b3*w1*w2*a1+2*X1*X2*b1*w1*w2*a1+2*X1*X3*b1*w1*w3*a1+2*X2*X3*b2*w2*w3*a2...
            +2*X1*X2*b2*w1*w2*a2+2*X0*X1*b2*w2*w3*a2-2*X0*X2*b1*w2*w3*a2+2*X0*X2*b3*w1*w2*a2-2*X0*X3*b2*w1*w2*a2...
            +2*X1*X3*b1*w2*w3*a2+2*w3*X2*b5*w1*Y1*a1+2.0*b3*w3*X3*v1*X2*a5+2.0*b3*w2*X0*v3*X2*a5+2.0*b1*w1*X1*v2*X3*a5...
            -2.0*b3*w3*X1*v2*X3*a5+2.0*b1*w2*X2*v1*X0*a5-2.0*b3*w1*X3*v3*X2*a5+2.0*X1*Y2*b5*v1*w2*a5+2.0*X0*Y3*b5*v1*w2*a5...
            -2.0*X2*Y0*b5*v3*w1*a5-2.0*X2*Y3*b5*v3*w2*a5+2*w1*X3*b5*w3*Y0*a1+2*w1*X1*b5*w2*Y3*a1;


        Vd(joint,1) = Id*(s4dd)-PddS;
    end
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



%solving displacement versor D, body velocity Vb, body acceleration Vbd using CGA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.Algorithm = 'levenberg-marquardt';             %solver for fsolve-function
opts.Display = 'off';

D       = zeros(8,N);
D0      = qToDualQuat(q0);      %initial position
J       = zeros(6,6,N);
Phi_vel = zeros(6,1);
V       = zeros(6,N);
Vb      = zeros(6,N);
Vd      = zeros(6,N);
Vbd     = zeros(6,N);

for i = 1:N
    
    funPhi = @(M) Phi(M,s4(:,i));   %function for solver
    %Displacement versor "D"
    if i == 1;
        D(:,i) = fsolve(funPhi,D0,opts); %interative solver
    else
       D(:,i) = fsolve(funPhi,D(:,i-1),opts);   %interative solver
    end
    
    %Dual quaternion
    X0 = D(1,i);
    X1 = D(2,i);
    X2 = D(3,i);
    X3 = D(4,i);
    Y0 = D(5,i);
    Y1 = D(6,i);
    Y2 = D(7,i);
    Y3 = D(8,i);
    
    %Define Jacobian
    J(:,:,i) = Jacobi(D(:,i));                                          
   
    %velocity screw "V" in interia frame 
    for j = 1:6
    a4d = L(j,i)*Ld(j,i);
    Phi_vel(j,1) = -(X0^2+X1^2+X2^2+X3^2)*a4d;
    end
    V(:,i) = J(:,:,i)\Phi_vel;                           
    
    %current position "r" of platform frame
    r = [2*(X1*Y0-Y1*X0-X2*Y3+Y2*X3),2*(X1*Y3-Y1*X3+X2*Y0-Y2*X0),2*(-X1*Y2+Y1*X2+X3*Y0-Y3*X0)]'/(X1^2+X2^2+X3^2+X0^2);
    %velocity screw "Vb" in body frame
    Vb(1:3,i) = [V(1,i),V(2,i),V(3,i)]'+cross([V(4,i) V(5,i) V(6,i)]',r);% ; V(4:6,i)];      
    Vb(4:6,i) = V(4:6,i);
                                
    %Acceleration screw "Vd" in inertia frame
    Vd(:,i) = J(:,:,i)\Phi_acc(D(:,i),V(:,i),L(:,i),Ld(:,i),Ldd(:,i));   
    %Acceleration screw "Vbd" in body frame
    Vbd(1:3,i) = [Vd(1,i),Vd(2,i),Vd(3,i)]'+cross([Vd(4,i) Vd(5,i) Vd(6,i)]',r)+cross([Vb(4,i) Vb(5,i) Vb(6,i)]',[Vb(1,i) Vb(2,i) Vb(3,i)]');% ; Vd(4:6,1)];
    Vbd(4:6,i) = Vd(4:6,i);
end

H     = zeros(4,4,N); %homogenius transformation 
q_cga = zeros(6,N);   %pose of platform frame in 3D space

%converting displacement versor "D" to homogenius transformation "H" and pose "q_cga" of platform frame in 3D space
for i = 1:N    
H(1:4,1:4,i)  = quat2rotm(D(:,i));  
q_cga(:,i)    = rotm2eulerXYZ(H(:,:,i));
end   
%%%%Validation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if sim == true,

%integration of Dd for validation of solution Vb
Dvel      = D(:,1);       %displacment versor from integration of Dd
Dd        = zeros(8,1);   %derivative of D
Dvel(:,1) = D(:,1);       %inital

    for i = 1:N
        Dd(:,i)     = Ddot_b(Vb(:,i),Dvel(:,i));
        Dvel(:,i+1) = normalizeD(Dvel(:,i) + Dd(:,i)*dt); %normalized for each step to maintain D*Dr=1 
    end

%integration of Ddd for validation of solution Vbd
Dacc      = zeros(8,N); %displacment versor from integration of Ddd
Ddd       = zeros(8,N); %second derivative of D
Dda(:,1)  = zeros(8,1); %derivative of D from integration of Ddd
Dacc(:,1) = D(:,1);     %inital

    for i = 1:N  
        Ddd(:,i+1)  = Dddot_b(Vb(:,i),Vbd(:,i),Dacc(:,i),Dda(:,i));
        Dda(:,i+1)  = (Dda(:,i)+(Ddd(:,i)*dt));
        Dacc(:,i+1) = normalizeD(Dacc(:,i)+(Dda(:,i))*dt+0.5*(Ddd(:,i)*dt^2)); %normalized for each step to maintain D*Dr=1
    end
 
   %% Animation figure settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clf
h = figure;
%scrsz = get(0, 'ScreenSize'); 
%h=figure('Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
%h.Position = [380 200 1000 800];
axis equal, grid on
view(25,20)
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2,2]*1.5);
ylim([-2,2]*1.5);
zlim([0,3]);

tStart = cputime();
tFactor = 1;
frame = 0;

bsol_CGA  = zeros(3,6,N);
bsol_CGAv = zeros(3,6,N);
bsol_CGAa = zeros(3,6,N);
while ishandle(h)

    % Calculate elapsed time
    tElapsed = (cputime() - tStart)*tFactor;
    
    [~,i] = min(abs(t-tElapsed));
    title(['t = ',num2str(t(i),'%.2f'),' s']);
    
    
    % Platfrom coordiante
    
     % H  = quat2rotm(D(:,i));
      Hv = quat2rotm(Dvel(:,i));
      Ha = quat2rotm(Dacc(:,i));
      
    for j = 1:6   
%         bsol_CGA(:,j,i) = H*[bp(:,j);1];
%         bsol_CGAv(:,j,i) = Hv*[bp(:,j);1];
%         bsol_CGAa(:,j,i) = Ha*[bp(:,j);1];
        bsol_CGA(:,j,i)  = sandwich(D(:,i),bp(:,j));
        bsol_CGAv(:,j,i) = sandwich(Dvel(:,i),bp(:,j));
        bsol_CGAa(:,j,i) = sandwich(Dacc(:,i),bp(:,j));
    end
  
        % Transform objets
        hold on
        cla
        % Base joints
        plot3(a(1,:),a(2,:),a(3,:),'*r')
        
        % PLatfrom joints
       plot3(b_temp(1,:,i),b_temp(2,:,i),b_temp(3,:,i),'*b','LineWidth',2)
%        plot3(bsol_CGA(1,:,i),bsol_CGA(2,:,i),bsol_CGA(3,:,i),'*r')
%        plot3(bsol_CGAv(1,:,i),bsol_CGAv(2,:,i),bsol_CGAv(3,:,i),'*g')
%        plot3(bsol_CGAa(1,:,i),bsol_CGAa(2,:,i),bsol_CGAa(3,:,i),'*m')
%       
        for k=1:6
        line([a(1,k),bsol_CGA(1,k,i)],[a(2,k),bsol_CGA(2,k,i)],[a(3,k),bsol_CGA(3,k,i)],'LineWidth',2) 
        line([a(1,k),bsol_CGAv(1,k,i)],[a(2,k),bsol_CGAv(2,k,i)],[a(3,k),bsol_CGAv(3,k,i)],'LineWidth',2)
        line([a(1,k),bsol_CGAa(1,k,i)],[a(2,k),bsol_CGAa(2,k,i)],[a(3,k),bsol_CGAa(3,k,i)],'LineWidth',2)
        end
        
        for k=1:6
      
            if k ==6
      l  =  line([bsol_CGA(1,k,i),bsol_CGA(1,1,i)],[bsol_CGA(2,k,i),bsol_CGA(2,1,i)],[bsol_CGA(3,k,i),bsol_CGA(3,1,i)],'Color','red','LineWidth',2,'DisplayName','D')  ;
      lv =  line([bsol_CGAv(1,k,i),bsol_CGAv(1,1,i)],[bsol_CGAv(2,k,i),bsol_CGAv(2,1,i)],[bsol_CGAv(3,k,i),bsol_CGAv(3,1,i)],'Color','green','LineWidth',2,'DisplayName','D');
      la =   line([bsol_CGAa(1,k,i),bsol_CGAa(1,1,i)],[bsol_CGAa(2,k,i),bsol_CGAa(2,1,i)],[bsol_CGAa(3,k,i),bsol_CGAa(3,1,i)],'Color','magenta','LineWidth',2,'DisplayName','D');
            else
        line([bsol_CGA(1,k,i),bsol_CGA(1,k+1,i)],[bsol_CGA(2,k,i),bsol_CGA(2,k+1,i)],[bsol_CGA(3,k,i),bsol_CGA(3,k+1,i),],'Color','red','LineWidth',2) 
        line([bsol_CGAv(1,k,i),bsol_CGAv(1,k+1,i)],[bsol_CGAv(2,k,i),bsol_CGAv(2,k+1,i)],[bsol_CGAv(3,k,i),bsol_CGAv(3,k+1,i),],'Color','green','LineWidth',2)
        line([bsol_CGAa(1,k,i),bsol_CGAa(1,k+1,i)],[bsol_CGAa(2,k,i),bsol_CGAa(2,k+1,i)],[bsol_CGAa(3,k,i),bsol_CGAa(3,k+1,i),],'Color','magenta','LineWidth',2)
            end
        end
        
        % Platfrom coordinate
         PlotAxes(H(:,:,i),0.3,2)
         PlotAxes(Hv,0.3,2)
         PlotAxes(Ha,0.3,2)
         
        % Base coordinate system
        %PlotAxes(eye(4,4),0.3,2)
     
    
       
    hold off

    if i == length(t)
        break;
    end
    
    drawnow;
    frame = frame + 1;
end
 h = legend('D','D from $\dot{D}$','D from $\ddot{D}$');
 set(h,'Interpreter','latex')
 
fprintf('FPS = %f \n',frame/t(i));
else
    %q_cga = sym('q_cga');
    Dd = sym('Dd');
    Ddd = sym('Ddd');
    Dacc = sym('Dacc');
    Dvel = sym('Dvel');
end   
end