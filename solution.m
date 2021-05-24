%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% This script solves the response x(t) for ODE given in %%%%%%%% 
%%%%% Q4.47, and compares the results with the response for %%%%%%%%
%%%%% ODE given in Q14.27.                                  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Step function is chosen as F(t)=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Ramp Function is chosen as F(t)=rt; where r = 1;%%%%%%%%%%%%%%
%%%%% Responses are plotted with time 0<t<20;%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initial conditions for forced vibration are zero;%%%%%%%%%%%%%
%%%%% Initial conditions for free vibration is x(0)=1;v(0)=1.%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solution 
clear 
global m c k r mu0 mu1 mu2 mu3 step x_initial
m = 1;
c = 0.5;
k = 1;
mu0 = 0;
mu1 = 0.01;
mu2 = 0.1;
mu3 = 1;
%------Define the shape for step function and ramp function--------
step = 1;
r = 1;
%------Define initial condition of homogeneous solution------------
% initial displacement:
x_initial(1) = 1;
% initial velocity:
x_initial(2) = 1;

%------Solve x(t) using the defined functions--------
[t1,x1] = ode45(@step_mu0,[0 20],[0 0]);
[t2,x2] = ode45(@step_mu1,[0 20],[0 0]);
[t3,x3] = ode45(@step_mu2,[0 20],[0 0]);
[t4,x4] = ode45(@step_mu3,[0 20],[0 0]);
[t5,x5] = ode45(@ramp_mu0,[0 20],[0 0]);
[t6,x6] = ode45(@ramp_mu1,[0 20],[0 0]);
[t7,x7] = ode45(@ramp_mu2,[0 20],[0 0]);
[t8,x8] = ode45(@ramp_mu3,[0 20],[0 0]);

%-------plot the results---------------------
figure(1)
plot(t1,x1(:,1))
hold on
plot(t2,x2(:,1))
hold on
plot(t3,x3(:,1))
hold on
plot(t4,x4(:,1))
hold off
legend('Step Function Response mu = 0','Step Function Response mu = 0.01',...
'Step Function Response mu = 0.1','Step Function Response mu = 1')
xlabel('Time t')
ylabel('x(t)')
title('Step Function Response')

figure(2)
plot(t5,x5(:,1))
hold on
plot(t6,x6(:,1))
hold on
plot(t7,x7(:,1))
hold on
plot(t8,x8(:,1))
hold off

legend('Ramp Function Response mu = 0','Ramp Function Response mu = 0.01',...
'Ramp Function Response mu = 0.1','Ramp Function Response mu = 1')
xlabel('Time t')
ylabel('x(t)')
title('Ramp Function Response')

%------solve responses for free vibration-------------- 
[t21,x21] = ode45(@H_mu0,[0 20],x_initial);
[t22,x22] = ode45(@H_mu1,[0 20],x_initial);
[t23,x23] = ode45(@H_mu2,[0 20],x_initial);
[t24,x24] = ode45(@H_mu3,[0 20],x_initial);

%------plot free vibration responses-------------------
figure(3)
plot(t21,x21(:,1))
hold on
plot(t22,x22(:,1))
hold on
plot(t23,x23(:,1))
hold on
plot(t24,x24(:,1))
hold off
legend('Homogeneous Response mu = 0','Homogeneous Response mu = 0.01',...
'Homogeneous Response mu = 0.1','Homogeneous Response mu = 1')
xlabel('Time t')
ylabel('x(t)')
title('Response (Homogeneous equation)')

end



%-----Define Step Function ODE with different mu values---------
function xdot1 = step_mu0(t,x)
global m c k mu0 step
xdot1 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu0*x(1)^3 + step)];
end

function xdot2 = step_mu1(t,x)
global m c k mu1 step
xdot2 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu1*x(1)^3 + step)];
end

function xdot3 = step_mu2(t,x)
global m c k mu2 step
xdot3 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu2*x(1)^3 + step)];
end

function xdot4 = step_mu3(t,x)
global m c k mu3 step
xdot4 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu3*x(1)^3 + step)];
end

%-----Define Ramp Function ODE with different mu values---------
function xdot5 = ramp_mu0(t,x)
global m c k r mu0
xdot5 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu0*x(1)^3 + r*t)];
end

function xdot6 = ramp_mu1(t,x)
global m c k r mu1 
xdot6 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu1*x(1)^3 + r*t)];
end

function xdot7 = ramp_mu2(t,x)
global m c k r mu2
xdot7 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu2*x(1)^3 + r*t)];
end

function xdot8 = ramp_mu3(t,x)
global m c k r mu3
xdot8 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu3*x(1)^3 + r*t)];
end


%-----Define Homogeneous ODE with different mu values---------
function xdot9 = H_mu0(t,x)
global m c k mu0
xdot9 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu0*x(1)^3)];
end

function xdot10 = H_mu1(t,x)
global m c k mu1 
xdot10 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu1*x(1)^3)];
end

function xdot11 = H_mu2(t,x)
global m c k mu2
xdot11 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu2*x(1)^3)];
end

function xdot12 = H_mu3(t,x)
global m c k mu3
xdot12 = [x(2); (1/m)*(-c*x(2)-k*x(1)-mu3*x(1)^3)];
end