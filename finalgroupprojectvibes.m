%% Vibrations Final Project Group 4

clear 
close all 
clc

syms J theta theta_ddot l_1 l_2 k_1 k_2 z m z_ddot

d=8.04;%mm
J=pi()*d^4/32; %solid cylinder shaft

%EOMs
eqn(1) = J*theta_ddot+(l_2*k_2-l_1*k_1)*z+(k_1*l_1^2+k_2*l_2^2)*theta==0;
eqn(2) = m*z_ddot+k_1*(z-l_1*theta)+k_2*(z-l_2*theta)==0;

z=solve(eqn,z);