% +=====================================================================+
% | Course: ME 400 - Mechanical Vibrations
% | Title: Project
% | Developed by: Owen Glascoe
% +=====================================================================+
clear
close all
clc
%% [SET UP]
% +---------------------------------------------------------------------+
% | Inputs
% +---------------------------------------------------------------------+
m_b = 1;             % Mass in kg
base = 0.2;
height = 0.01;
J = 1;%1/12*(base^4+height^4);  % Polar Moment in kg*m^2
k1 = 1;              % Spring 1 Constant in N/m
k2 = 1;              % Spring 2 Constant in N/m
l1 = 2;              % Distance from CG in m
l2 = 1;              % Distance from CG in m
c = zeros(2,2);      % Damping Constant in Ns/m
F_0 = 1;             % Excitation Force in N
omega_0 = 0;         % Excitation Frequency in rad/s
z_0 = 0;             % Initial Displacement in m
theta_0 = 0;         % Initial Rotation in rad
v_i = [0;0];         % Initial Velocity in m/s
x_i = [z_0;theta_0]; % Initial Displacement
%% [MODAL]
% +---------------------------------------------------------------------+
% | Mass in kg
% +---------------------------------------------------------------------+
m11 = m_b;                
m12 = 0;                  
m21 = 0;                  
m22 = J;                  
m = [m11 m12; m21 m22];   
% +---------------------------------------------------------------------+
% | Stiffness in N/m
% +---------------------------------------------------------------------+
k11 = k1+k2;              
k12 = k2*l2 - k1*l1;      
k21 = k2*l2 - k1*l1;      
k22 = k1*l1^2 + k2*l2^2;  
k = [k11 k12; k21 k22];   
% +---------------------------------------------------------------------+
% | Mode Shape & Natural Frequencies in rad/s
% +---------------------------------------------------------------------+
[psi,omega_n_2] = eig(k,m);      % Eigenvalue/Eigenvector problem
omega_n = sqrt(diag(omega_n_2)); % Natural Frequency in rad/s
omega_n1 = omega_n(1);           % Natural Frequency in rad/s
omega_n2 = omega_n(2);           % Natural Frequency in rad/s
% +---------------------------------------------------------------------+
% | Modal Matrices
% +---------------------------------------------------------------------+
M = transpose(psi)*m*psi; % Modal Mass Matrix
K = transpose(psi)*k*psi; % Modal Stiffness Matrix
%% [PLOTTING]
% +---------------------------------------------------------------------+
% | Plot the 2DOF Time History of X1 and X2
% +---------------------------------------------------------------------+  
tspan =  [0 20];
[t,X] = ode45(@(t,X) two_DOF_forced_vibes(t,X,m,c,k,F_0,omega_0), tspan, [x_i; v_i]);
figure (1)
yyaxis left
plot(t,X(:,1),'LineWidth',1.25)
ylabel('Displacement (m)');
yyaxis right
plot(t,X(:,2),'LineWidth',1.25)
ylabel('Displacement (rad)');
xlabel('Time (s)');
title('Time History of Z & \theta');
legend('Z','\theta')
grid on;


%% [Animation]
% +---------------------------------------------------------------------+
% | 2DOF Animation
% +---------------------------------------------------------------------+
animate = 1;
width = 0.025;
figure (2)
if (animate == 1)
    ax1 = -l1+cos(X(:,2));
    ay1 = X(:,1)-l1*X(:,2);
    ax2 = l2+cos(X(:,2));
    ay2 = X(:,1)+l2*X(:,2);
    cx1 = cos(X(:,2));
    cz = X(:,1);
    h = plot([ax1(1),ax2(1)],[ay1(1),ay2(1)],'-k',[ax1(1), ax1(1)],[ay1(1),2],'-r',[ax2(1), ax2(1)],[ay2(1),2],'-r',[cx1(1)-width,cx1(1)+width],[cz(1)-width,cz(1)+width],'-b','LineWidth',3);
    xlabel('x, m')
    ylabel('z, m')
    axis equal
    ylim([-0.5,2])
    xlim([-2,3])
    for i = 2:length(ax1)
        set(h(1),'XData',[ax1(i),ax2(i)],'YData',[ay1(i),ay2(i)]);
        set(h(2),'XData',[ax1(i),ax1(i)],'YData',[ay1(i),2]);
        set(h(3),'XData',[ax2(i),ax2(i)],'YData',[ay2(i),2]);
        set(h(4),'XData',[cx1(i)-width,cx1(i)+width],'YData',[cz(i)-width,cz(i)+width]);
        drawnow
        pause(1/30);
    end
end

%% [FUNCTIONS]
% +---------------------------------------------------------------------+
% | 2DOF Forced Harmonic Vibration
% +---------------------------------------------------------------------+
function dXdt = two_DOF_forced_vibes(t, x, m, c, k, F_0, omega)
    % x(1): Displacement of Mass 1
    % x(2): Displacement of Mass 2    
    % x(3): Velocity of Mass 1
    % x(4): Velocity of Mass 2
    dXdt = zeros(4,1);                      % Preallocate
    temp = zeros(4,4);                      % Preallocate
    m_inv = inv(m);                         % Inverse of Mass Matrix
    temp(1:2, 3:4) = eye(2);                
    temp(3:4, 1:2) = -m_inv*k;               
    temp(3:4, 3:4) = -m_inv*c;              
    Fa = zeros(4,1);                        % Preallocate
    Fa(3:4,1) = m_inv*[F_0*cos(omega*t); 0]; % Applied Force
    Fa_0(3:4,1) = m_inv*[0*cos(omega*t); 0]; % Applied Force
    
%{
 if t==0
        dXdt = temp*x + Fa;
    else
        dXdt = temp*x + Fa_0;
    end 
%}

    
    dXdt = temp*x + Fa;                     % Solution
end