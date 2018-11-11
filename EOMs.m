clc
clear
close all


syms phi_1 phi_2 phi_3; %angles between the directions of the p-joints and x0 axis
syms m_slider m_MP J_MP;
syms initial_Px_1 initial_Px_2 initial_Px_3; %initial x-coordinates of P joint (A_j) in the BP frame
syms initial_Py_1 initial_Py_2 initial_Py_3; %initial y-coordinates of P joint (A_j) in the BP frame
syms x_c y_c theta_c; %our cartesian variables
syms x_c_dot y_c_dot theta_c_dot;
syms a; %side of triangle
syms u_1 u_2 u_3; %displacements of P joints from initial locations

phi = [phi_1; phi_2; phi_3]; 
initial_Px = [initial_Px_1; initial_Px_2; initial_Px_3]; 
initial_Py = [initial_Py_1; initial_Py_2; initial_Py_3]; 
p = [x_c; y_c; theta_c]; 
p_dot = [x_c_dot; y_c_dot; theta_c_dot]; 
M_s = [m_slider 0 0; 0 m_slider 0; 0 0 m_slider];
u = [u_1; u_2; u_3]; 

% computing D
%d = zeros(3,1);
for i=1:3
    d(i) = (p(1)-initial_Px(i))*cos(phi(i)) + (p(2)-initial_Py(i))*sin(phi(i)) + sqrt(3)/3*a*cos(p(3)-pi/6);
end
D = [d(1) 0 0; 0 d(2) 0; 0 0 d(3)];

% computing K
%K = zeros(3,3);
for j=1:3
    K(j, 1) = p(1) - initial_Px(j) + sqrt(3)/3*a*cos(p(3)+phi(j)-pi/6);
    K(j, 2) = p(2) - initial_Py(j) + sqrt(3)/3*a*sin(p(3)+phi(j)+pi/6);
    K(j, 3) = -sqrt(3)/3*a*(-(p(1) - initial_Px(j))*sin(p(3)+phi(j)-pi/6)+(p(2) - initial_Py(j))*cos(p(3)+phi(j)-pi/6)+u(j)*sin(p(3)-pi/6));
end

KE = transpose((inv(D)*K*p_dot))*M_s*(inv(D)*K*p_dot);
% diff(transpose((inv(D)*K*p_dot))*M_s*(inv(D)*K*p_dot), x_c);
% diff(transpose((inv(D)*K*p_dot))*M_s*(inv(D)*K*p_dot), y_c);
% diff(transpose((inv(D)*K*p_dot))*M_s*(inv(D)*K*p_dot), theta_c);


%%%%%%%%%% Design Parameters %%%%%%%%%%
phi_dsgn = [0; 2/3*pi; -pi/3];
m_slider_dsgn = 0.5;
R = 0.5;
r = 0.2;

%assuming in figure on pg8 that B3A1=B1A2=B2A3=sqrt(3)*r
initial_Px_dsgn = [-R*sqrt(3)/2+sqrt(3)*r; R*sqrt(3)/2-sqrt(3)/2*r; -r*sqrt(3)/2];
initial_Py_dsgn = [0; -R/2+3*r/2; R-3*r/2];

KE = transpose((inv(D)*K*p_dot))*M_s*(inv(D)*K*p_dot);

for i=1:3
    KE = subs(KE, phi(i), phi_dsgn(i));
    KE = subs(KE, initial_Px(i), initial_Px_dsgn(i));
    KE = subs(KE, initial_Py(i), initial_Py_dsgn(i));
end

KE = subs(KE, m_slider, m_slider_dsgn);
vpa(KE)
