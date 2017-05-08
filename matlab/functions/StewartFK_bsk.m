function [q,qd,qdd] = StewartFK_bsk(q0,L,LD,LDD,t,sim)
%function [q,qd,qdd] = StewartFK_bsk(q0,L,LD,LDD,t,sim)
% 
% This function is a forward kinematic algortithm for a 6-SPS stewart
% platform. It applies 3D Spatial vector kinematic for defining constaint
% equation. From given leg-length, leg-speed and leg-acceleration, the pose
% velocity and acceleration of the platform frame is returned.
% 
% Created by Olav Heng.
% INPUT:
% q0        : Initial pose components [x y z phi theta psi]' for the stewart platfrom, size:(6,1)
% L         : Leg-length of stewart platform, size:(6,n)
% Ld        : Leg extraction speed, size:(6,n)
% Ldd       : Leg extraction acceleration, size:(6,n)
% t         : Time for each step n, required for simulation, size:(1,n)
% sim       : Simulation on/off: 1/0
%
% OUTPUT
% q         : The pose components [x y z phi theta psi]', size:(6,n)
% qd        : The velocity components of [x y z phi theta psi]' for the stewart platfrom body frame, size:(6,n)
% qdd       : The acceleration components of [x y z phi theta psi]' for the stewart platfrom body frame, size:(6,n)

global a bp

N = length(L);

%Stewart parameter
jointLoc = [15,105,135,225,255,345]/180*pi; %location of the base points "a"
a0 = [2,0,0]';

Az = @(psi)[cos(psi),-sin(psi), 0;
            sin(psi), cos(psi), 0;
            0       , 0       , 1];

%distribute spherical joints "a" in the base, and "bp" on the platform
a  = zeros(3,6);
bp = zeros(3,6);
for i = 1:length(jointLoc)
    
    a(:,i) = Az(jointLoc(i))*a0;
    
        if i == 6 %shifting joint position in bp so that [a1,a2...] and [bp1,bp2...] corrospond to each other
        bp(:,1) = Az(60/180*pi)*a(:,i);
        else
        bp(:,i+1) = Az(60/180*pi)*a(:,i);
        end

end

%%%define sub-function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ws = scew(w)

    ws = zeros(3,3);

    ws(1,:) = [0 -w(3,1) w(2,1)];
    ws(2,:) = [w(3,1) 0 -w(1,1)];
    ws(3,:) = [-w(2,1) w(1,1) 0];
    
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

function A = Axyz(phi)

    rx = phi(1,1);
    ry = phi(2,1);
    rz = phi(3,1);

    A = eye(3,3);

    A(1,1:3) = [cos(ry)*cos(rz), -cos(ry)*sin(rz), sin(ry)];
    A(2,1:3) = [cos(rx)*sin(rz) + cos(rz)*sin(rx)*sin(ry), cos(rx)*cos(rz) - sin(rx)*sin(ry)*sin(rz), -cos(ry)*sin(rx)];
    A(3,1:3) = [sin(rx)*sin(rz) - cos(rx)*cos(rz)*sin(ry), cos(rz)*sin(rx) + cos(rx)*sin(ry)*sin(rz),  cos(rx)*cos(ry)];
end

function J = Jacobi(q)

    % Rename states
    r = q(1:3,1);
    phi = q(4:6,1);

    J = zeros(6,6);
    
        for joint = 1:6
            b = Axyz(phi)*bp(:,joint);
            d = r + b - a(:,joint);

        J(joint,1:6) = [2*d', 2*d'*(scew([1 0 0]')*b),2*d'*(scew([0 1 0]')*b),2*d'*(scew([0 0 1]')*b)];


        end

end

function Fi = Phi(q,L)

    % Rename states
    r   = q(1:3,1);
    phi = q(4:6,1);
    Axyz(phi);

    Fi = zeros(6,1);
    
            for joint = 1:6
                d = r + Axyz(phi)*bp(:,joint) - a(:,joint);

                Fi(joint,1) = d'*d - L(joint,1)^2;
            end

end

function acc  = acc(q,qd,L,Ld,Ldd)
    
    % Rename states
    r     = q(1:3,1);
    phi   = q(4:6,1);
    rd    = qd(1:3,1);
    phid  = qd(4:6,1);

    acc = zeros(6,1);

        for joint = 1:6

        b   = Axyz(phi)  *  bp(:,joint);
        d   = r + b - a(:,joint);

        dd  = rd + scew(phid)*b;
        ddd = (scew(phid)  *  scew(phid)) * b;


        acc(joint,1) =   -2 * (d'*ddd) -   2 * (dd'*dd)  +  2 * Ld(joint,1)^2   +   2 * L(joint,1) * Ldd(joint,1);

        end
end





% Solve Forward Kinematics
q = zeros(6,N);
qd = zeros(6,N);
qdd = zeros(6,N);
J = zeros(6,6,N);
opts.Algorithm = 'levenberg-marquardt';
opts.Display = 'off';
tic
for i = 1:N
    funPhi = @(qc) Phi(qc,L(:,i));
   
    if i == 1;
        q(:,i) = fsolve(funPhi,q0,opts);
    else
        q(:,i) = fsolve(funPhi,q(:,i-1),opts);
    end
    
     J(:,:,i) = Jacobi([q(1:3,i)' q(4:6,i)']');
     
     qd(:,i) = J(:,:,i) \ (2* (L(:,i) .* LD(:,i)));
    
     qdd(:,i) = J(:,:,i) \ acc(q(:,i),qd(:,i),L(:,i),LD(:,i),LDD(:,i));

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

bsol = zeros(4,6,N);
H = eye(4,4);

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
    

    H(1:3,1:3) = Axyz(q(4:6,i));
    H(1:3,4)   = q(1:3,i);
    for j = 1:6
    bsol(:,j,i) = H*[bp(:,j);1];
    end
    % Transform objets
    hold on
        cla
        % Base joints
        plot3(a(1,:),a(2,:),a(3,:),'*r')
        
        % platform joints
        plot3(bsol(1,:,i),bsol(2,:,i),bsol(3,:,i),'*b')
      
        %platform frame
        for k=1:6
        line([a(1,k),bsol(1,k,i)],[a(2,k),bsol(2,k,i)],[a(3,k),bsol(3,k,i)],'LineWidth',2) 

        end
        
        for k=1:6
      
            if k ==6
        line([bsol(1,k,i),bsol(1,1,i)],[bsol(2,k,i),bsol(2,1,i)],[bsol(3,k,i),bsol(3,1,i)],'Color','red','LineWidth',2)  
            else
        line([bsol(1,k,i),bsol(1,k+1,i)],[bsol(2,k,i),bsol(2,k+1,i)],[bsol(3,k,i),bsol(3,k+1,i),],'Color','red','LineWidth',2) 
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