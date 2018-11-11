clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=0.3;   %m,radius of the BP
r=0.1; %m,radius of the MP
l=R-r;  %m,length of the limb rod
d0=sqrt(3)*r;  %m,initial inertia position of the P joint
d=l/10; %m,diameter of the limb rod
h=d;    %m,height of the MP

rho_r=2.83*1000;    %kg/m^3,density of the limb rod,aluminium alloy
rho_MP=rho_r;        %kg/m^3,density of the MP

m_s=0.5;    %kg,mass of the slider
I_s=0;      %kgm^2,second moment of inertia of the slider w.r.t. its c.o.m.
m_r=pi/4*rho_r*d^2*l;   %kg,mass of the limb rod
I_r=m_r/12*(l^2+3/4*d^2);         %kgm^2,second moment of inertia of the limb rod w.r.t. its c.o.m.
m_MP=3*sqrt(3)/4*rho_MP*h*r^2;  %kg,mass of the MP
I_MP=3*sqrt(3)/16*rho_MP*h*r^4; %kgm^2,second moment of inertia of the MP triangle w.r.t. its c.o.m.
m_l=5;  %kg,mass of the load
I_l=m_l*r^2/2;  %kgm^2,second moment of inertia of the load w.r.t. its c.o.m.

%directions of the P joint
Phi=[0,2/3*pi,-2/3*pi];
for j=1:3
   e(:,j)=[cos(Phi(j));sin(Phi(j))];
end

%coordinates of the BP vertex B_j
for j=1:3
   xy_B(:,j)=R*[cos(Phi(j)-5*pi/6);sin(Phi(j)-5*pi/6)];
end

%%%%%%%%%%%%%%%%%%%%%%%% Desired Trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=10;   %Hz,shaking frequecy
Circle=2;
N=1000; %sampling number in one cycle
T=1/f/N*Circle;    %s,sampling interval
A_x=5/1000; %m,amplitude of the translation in x axis
A_y=0/1000; %m,amplitude of the translation in y a;   
A_r=0*pi/180;   %rad,amplitude of the rotation aroud z axis

%Cartesian variables
for i=1:N
   t(i)=T*i;
   X(i)=A_x*sin(2*pi*f*t(i)+pi/3);
   Y(i)=A_y*sin(2*pi*f*t(i));
   Theta(i)=A_r*sin(2*pi*f*t(i)-pi/3);
end
t_v=t(1:N-1);   t_a=t(1:N-2);

%Cartesian velocities
X_dot=diff(X,1,2)/T;
Y_dot=diff(Y,1,2)/T;
Theta_dot=diff(Theta,1,2)/T;
for i=1:N-1
	p_dot(:,i)=[Theta_dot(i);X_dot(i);Y_dot(i)]; 
end

%%%%%%%%%%%%%%%%%%%%%%%% Joint variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial joint variables
for j=1:3
    u(j,1)=0;
    Vphi_A(j,1)=pi/6;
    Vphi_M(j,1)=5*pi/6;
    gamma(:,j,1)=[u(j,1);Vphi_A(j,1);Vphi_M(j,1)]; 
end

E=[0,-1;1,0];
%joint-displacement & joint-rate array
for i=1:N-1 %jacobian
    for j=1:3
       p_AM(:,j,i)=[l*cos(Vphi_A(j,i)+Phi(j));l*sin(Vphi_A(j,i)+Phi(j))];   %vector A_jM_j
%       p_M(:,j,i)=[r*cos(Vphi_A(j,i)+Vphi_M(j,i)+Phi(j)-pi/6);r*sin(Vphi_A(j,i)+Vphi_M(j,i)+Phi(j)-pi/6)]; %vector M_jC%
       p_M(:,j,i)=[r*cos(pi+Theta(i)+Phi(j)-pi/6);r*sin(pi+Theta(i)+Phi(j)-pi/6)]; %vector M_jC
       p_A(:,j,i)=p_AM(:,j,i)+p_M(:,j,i);   %vector A_jC
       %Jacobian of the jth limb at ith moment
       J(:,:,j,i)=[0,1,1;e(:,j),E*p_A(:,j),E*p_M(:,j)];   J_in(:,:,j,i)=inv(J(:,:,j,i));
       %joint-rate array of the jth limb at ith moment
       gamma_dot(:,j,i)=J_in(:,:,j,i)*p_dot(:,i);
 end
    %joint-displacement array at the (i+1)th moment
    gamma(:,:,i+1)=gamma(:,:,i)+gamma_dot(:,:,i)*T;
    
    %joint displacements at the (i+1)th moment
	u(:,i+1)=gamma(1,:,i+1);
	Vphi_A(:,i+1)=gamma(2,:,i+1);
	Vphi_M(:,i+1)=gamma(3,:,i+1);
    
end

%joint velocities at the ith moment
u_dot(:,:)=gamma_dot(1,:,:);
Vphi_A_dot(:,:)=gamma_dot(2,:,:);
Vphi_M_dot(:,:)=gamma_dot(3,:,:);

%joint-acceleration array
gamma_ddot=diff(gamma,2,3)/T^2;

%joint accelerations 
u_ddot(:,:)=gamma_ddot(1,:,:);
Vphi_A_ddot(:,:)=gamma_ddot(2,:,:);
Vphi_M_ddot(:,:)=gamma_ddot(3,:,:);

%%%%%%%%%%%%%%%%%%%% Twist-shaping matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N-1
   for j=1:3
       Us(:,:,j)=[0,0,0;e(:,j),zeros(2,1),zeros(2,1)];
       Ur(:,:,j)=[0,1,0;e(:,j),p_AM(:,j,i)/2,zeros(2,1)];
   end
   U1(:,:,i)=[Us(:,:,1),zeros(3),zeros(3)];
   U2(:,:,i)=[zeros(3),Us(:,:,2),zeros(3)];
   U3(:,:,i)=[zeros(3),zeros(3),Us(:,:,3)];
   U4(:,:,i)=[Ur(:,:,1),zeros(3),zeros(3)];
   U5(:,:,i)=[zeros(3),Ur(:,:,2),zeros(3)];
   U6(:,:,i)=[zeros(3),zeros(3),Ur(:,:,3)];
   U7(:,:,i)=[J(:,:,1,i),zeros(3),zeros(3)];
   U(:,:,i)=[U1(:,:,i); U2(:,:,i); U3(:,:,i); U4(:,:,i); U5(:,:,i); U6(:,:,i); U7(:,:,i)];
   Jin(:,:,i)=[J_in(:,:,1,i);J_in(:,:,2,i);J_in(:,:,3,i)];
   D(:,:,i)=zeros(3);
   for j=1:3
       K(j,1,i)=r*sin(Vphi_M(j,i)-pi/6);
       K(j,2,i)=cos(Vphi_A(j,i)+Phi(j));
       K(j,3,i)=sin(Vphi_A(j,i)+Phi(j));
       D(j,j,i)=cos(Vphi_A(j,i));
   end
   Th(:,:,i)=Jin(:,:,i)*inv(K(:,:,i))*D(:,:,i);
   TS(:,:,i)=U(:,:,i)*Th(:,:,i);
end
%Twist-shaping derivatives
TS_dot=diff(TS,1,3)/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inertia dyad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N-1
    M1(:,:,i)=[0,zeros(1,2);zeros(2,1),m_s*eye(2)];
    M2(:,:,i)=[0,zeros(1,2);zeros(2,1),m_s*eye(2)];
    M3(:,:,i)=[0,zeros(1,2);zeros(2,1),m_s*eye(2)];
    for j=1:3
        p_Cr(:,j,i)=xy_B(:,j)+(d0+u(j,i))*e(:,j)+p_AM(:,j,i)/2;
        Iz(j+3,i)=I_r+m_r*p_Cr(:,j,i)'*p_Cr(:,j,i);
    end
    M4(:,:,i)=[Iz(4,i),zeros(1,2);zeros(2,1),m_r*eye(2)];
    M5(:,:,i)=[Iz(5,i),zeros(1,2);zeros(2,1),m_r*eye(2)];
    M6(:,:,i)=[Iz(6,i),zeros(1,2);zeros(2,1),m_r*eye(2)];
    m7=m_MP+m_l;
    I7=I_MP+I_l;
    Iz(7,i)=I7;%+m7*(X(i)^2+Y(i)^2);
    M7(:,:,i)=[Iz(7,i),zeros(1,2);zeros(2,1),m7*eye(2)];
    M(:,:,i)=[M1(:,:,i),zeros(3,18);
              zeros(3),M2(:,:,i),zeros(3,15);
              zeros(3,6),M3(:,:,i),zeros(3,12);
              zeros(3,9),M4(:,:,i),zeros(3,9);
              zeros(3,12),M5(:,:,i),zeros(3,6);
              zeros(3,15),M6(:,:,i),zeros(3);
              zeros(3,18),M7(:,:,i)];
end

%%%%%%%%%%%%%%%%%%%%%%%% Inverse Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N-2
    I(:,:,i)=TS(:,:,i)'*M(:,:,i)*TS(:,:,i); %inertia matrix
    C(:,:,i)=TS(:,:,i)'*M(:,:,i)*TS_dot(:,:,i); %Coriolis and centrifugal matrix
    f_u(:,i)=I(:,:,i)*u_ddot(:,i)+C(:,:,i)*u_dot(:,i); %Motor forces
end

%%%%%%%%%%%%%%%%%%%%%%%%% Motor Power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N-2
   for j=1:3
      P(j,i)=u_dot(j,i)*f_u(j,i); 
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)   %%%%% Desired Trajectory %%%%%
subplot(3,1,1);
plot(t,X*1000,'LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('$x_C$ (mm)','interpreter','latex','FontSize',14,'FontWeight','bold');
title('Cartesian Displacement','FontSize',16,'FontWeight','bold');
xlim([0,1/f*Circle]);
grid on;
subplot(3,1,2);
plot(t,Y*1000,'LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('$y_C$ (mm)','interpreter','latex','FontSize',14,'FontWeight','bold');
xlim([0,1/f*Circle]);
grid on;
subplot(3,1,3);
plot(t,Theta/pi*180,'LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('$\theta$ (degree)','interpreter','latex','FontSize',14,'FontWeight','bold');
grid on;

figure(2)   %%%%% Motor displacements %%%%%
plot(t,u(1,:)*1000,t,u(2,:)*1000,':',t,u(3,:)*1000,'-.','LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('Displacement (mm)','FontSize',14,'FontWeight','bold');
title('Motor Displacement','FontSize',16,'FontWeight','bold');
legend('$u_1$','$u_2$','$u_3$','interpreter','latex','FontWeight','bold');
xlim([0,1/f*Circle]);
grid on;

figure(3)   %%%%% Motor velocities %%%%%
plot(t_v,u_dot(1,:),t_v,u_dot(2,:),':',t_v,u_dot(3,:),'-.','LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('Velocity (m/s)','FontSize',14,'FontWeight','bold');
title('Motor Velocity','FontSize',16,'FontWeight','bold');
legend('$\dot{u}_1$','$\dot{u}_2$','$\dot{u}_3$','interpreter','latex','FontWeight','bold');
xlim([0,1/f*Circle]);
grid on;

figure(4)   %%%%% Motor accelerations %%%%%
plot(t_a,u_ddot(1,:),t_a,u_ddot(2,:),':',t_a,u_ddot(3,:),'-.','LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('Acceleration ($m/s^2$)','interpreter','latex','FontSize',14,'FontWeight','bold');
title('Motor Acceleration','FontSize',16,'FontWeight','bold');
legend('$\ddot{u}_1$','$\ddot{u}_2$','$\ddot{u}_3$','interpreter','latex','FontWeight','bold');
xlim([0,1/f*Circle]);
grid on;
hold on

figure(5)   %%%%% Motor Force %%%%%
plot(t_a,f_u(1,:),t_a,f_u(2,:),':',t_a,f_u(3,:),'-.','LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('Force ($N$)','interpreter','latex','FontSize',14,'FontWeight','bold');
title('Motor Forces','FontSize',16,'FontWeight','bold');
legend('$f_1$','$f_2$','$f_3$','interpreter','latex','FontWeight','bold');
xlim([0,1/f*Circle]);
grid on;

figure(6)   %%%%% Motor Power %%%%%
plot(t_a,P(1,:),t_a,P(2,:),':',t_a,P(3,:),'-.','LineWidth',2);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('Power ($W$)','interpreter','latex','FontSize',14,'FontWeight','bold');
title('Motor Power','FontSize',16,'FontWeight','bold');
legend('$P_1$','$P_2$','$P_3$','interpreter','latex','FontWeight','bold');
xlim([0,1/f*Circle]);
grid on;



