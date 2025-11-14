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
J = 1;               % Polar Moment in kg*m^2
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
tspan = [0 20];   
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

figure (2)
if (animate == 1)
    ax1 = 0-cos(X(:,2));
    ay1 = X(:,1)-l1*X(:,2);
    ax2 = 0+cos(X(:,2));
    ay2 = X(:,1)+l2*X(:,2);
    cz = X(:,1);
    h = plot([ax1(1),ax2(1)],[ay1(1),ay2(1)],'-k',[ax1(1), ax1(1)],[0,ay1(1)],'-r',[ax2(1), ax2(1)],[0,ay2(1)],'-r','LineWidth',3);
    p = plot(ax1(1)+l1,cz(1),'-b');
    xlabel('x, m')
    ylabel('y, m')
    axis equal
    ylim([-0.5,2])
    xlim([-2,2])
    for i = 2:length(ax1)
        set(h(1),'XData',[ax1(i),ax2(i)],'YData',[ay1(i),ay2(i)]);
        set(h(2),'XData',[ax1(i),ax1(i)],'YData',[0,ay1(i)]);
        set(h(3),'XData',[ax2(i),ax2(i)],'YData',[0,ay2(i)]);
        %p = plot(ax1(i)+l1,cz(i),'-b');
        drawnow
        pause(1/30);
    end
end

%% [FUNCTIONS]
% +---------------------------------------------------------------------+
% | 2DOF Forced Harmonic Vibration
% +---------------------------------------------------------------------+
function dXdt = two_DOF_forced_vibes(t, x, m, c, k, F0, omega)
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
    Fa(3:4,1) = m_inv*[F0*cos(omega*t); 0]; % Applied Force
    dXdt = temp*x + Fa;                     % Solution
end








%{
syms J theta theta_ddot l_1 l_2 k_1 k_2 z m z_ddot

d=8.04;%mm
J=pi()*d^4/32; %solid cylinder shaft

%EOMs
eqn(1) = J*theta_ddot+(l_2*k_2-l_1*k_1)*z+(k_1*l_1^2+k_2*l_2^2)*theta==0;
eqn(2) = m*z_ddot+k_1*(z-l_1*theta)+k_2*(z-l_2*theta)==0;

z=solve(eqn,z); 
%}
