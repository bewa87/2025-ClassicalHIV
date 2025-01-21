% Numerical Simulation Of Primary HIV-Infection

% Step 1: Initial Conditions And Problem Parameters
% Compare Stafford et al. 2000/ Alizon & Magnus 2012

% Initial Conditions

T_zero   = 10;
T_i_zero = 0;
V_zero   = 0.000001;

% Problem Parameters

r     = 0.17;
beta  = 0.00065;
d     = 0.01;
delta = 0.39;
pig   = 850;
c     = 3;

% Step 2: Equilibrium Point

T_star   = (c*delta)/(beta*pig);
T_i_star = (beta*pig*r-c*d*delta)/(beta*delta*pig);
V_star   = (beta*pig*r-c*d*delta)/(beta*c*delta);

% Step 3: Problem Definition

T            = 10;
h            = 0.0001;
t            = 0:h:T;
t            = t';
T_rk4        = zeros(length(t),1);
T_i_rk4      = zeros(length(t),1);
V_rk4        = zeros(length(t),1);
T_rk4(1)     = T_zero;
T_i_rk4(1)   = T_i_zero;
V_rk4(1)     = V_zero;

% Step 3: Comparison Solution Of Runge-Kutta-4

f1 = @(x,y,z) r - beta*x*z - d*x;
f2 = @(x,y,z) beta*x*z - delta*y;
f3 = @(x,y,z) pig*y - c*z;

for j = 1:1:(length(t)-1)
   k1_x         = f1(T_rk4(j),T_i_rk4(j),V_rk4(j));
   k1_y         = f2(T_rk4(j),T_i_rk4(j),V_rk4(j));
   k1_z         = f3(T_rk4(j),T_i_rk4(j),V_rk4(j));
   k2_x         = f1(T_rk4(j) + (h/2)*k1_x, T_i_rk4(j) + (h/2)*k1_y, V_rk4(j) + (h/2)*k1_z);
   k2_y         = f2(T_rk4(j) + (h/2)*k1_x, T_i_rk4(j) + (h/2)*k1_y, V_rk4(j) + (h/2)*k1_z);
   k2_z         = f3(T_rk4(j) + (h/2)*k1_x, T_i_rk4(j) + (h/2)*k1_y, V_rk4(j) + (h/2)*k1_z);
   k3_x         = f1(T_rk4(j) + (h/2)*k2_x, T_i_rk4(j) + (h/2)*k2_y, V_rk4(j) + (h/2)*k2_z);
   k3_y         = f2(T_rk4(j) + (h/2)*k2_x, T_i_rk4(j) + (h/2)*k2_y, V_rk4(j) + (h/2)*k2_z);
   k3_z         = f3(T_rk4(j) + (h/2)*k2_x, T_i_rk4(j) + (h/2)*k2_y, V_rk4(j) + (h/2)*k2_z);
   k4_x         = f1(T_rk4(j) + h*k3_x, T_i_rk4(j) + h*k3_y, V_rk4(j) + h*k3_z);
   k4_y         = f2(T_rk4(j) + h*k3_x, T_i_rk4(j) + h*k3_y, V_rk4(j) + h*k3_z);
   k4_z         = f3(T_rk4(j) + h*k3_x, T_i_rk4(j) + h*k3_y, V_rk4(j) + h*k3_z);
   T_rk4(j+1)   = T_rk4(j)   + (h/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x);
   T_i_rk4(j+1) = T_i_rk4(j) + (h/6)*(k1_y + 2*k2_y + 2*k3_y + k4_y);
   V_rk4(j+1)   = V_rk4(j)   + (h/6)*(k1_z + 2*k2_z + 2*k3_z + k4_z);
endfor

##figure(4)
##
##subplot(3,1,1);
##plot(t(1:1000:end),T_rk4(1:1000:end),'linestyle','-','marker','x');
##title('Load of Target CD4-cells');
##xlabel('t [days]');
##ylabel('T(t)','interpreter','tex');
##
##subplot(3,1,2);
##plot(t(1:1000:end),T_i_rk4(1:1000:end),'linestyle','-','marker','x','color','red');
##title('Load of Infected Target CD4-Cells');
##xlabel('t [days]');
##ylabel('T_{i}(t)','interpreter','tex');
##
##subplot(3,1,3);
##plot(t(1:1000:end),V_rk4(1:1000:end),'linestyle','-','marker','x','color','black');
##title('Load of Viral Particles');
##xlabel('t [days]');
##ylabel('V(t)','interpreter','tex');
##ylim([0, 2000]);

% Step 4: Grid Refinement For NSFDM

h = [0.1 0.05 0.025 0.0125];

% Step 4.1: Solve NSFDM For h(1)

t01            = 0:h(1):T;
T_nsfdm01      = zeros(length(t01),1);
T_i_nsfdm01    = zeros(length(t01),1);
V_nsfdm01      = zeros(length(t01),1);
T_nsfdm01(1)   = T_zero;
T_i_nsfdm01(1) = T_i_zero;
V_nsfdm01(1)   = V_zero;

% Loop

for j = 1:1:(length(t01)-1)
   T_nsfdm01(j+1)   = (T_nsfdm01(j) + h(1)*r)/(1 + h(1)*beta*V_nsfdm01(j) + h(1)*d);
   T_i_nsfdm01(j+1) = (T_i_nsfdm01(j) + h(1)*beta*T_nsfdm01(j+1)*V_nsfdm01(j))/(1 + h(1)*delta);
   V_nsfdm01(j+1)   = (V_nsfdm01(j) + h(1)*pig*T_i_nsfdm01(j+1))/(1 + h(1)*c);
endfor

err_x_01  = max(abs(T_rk4(1:1000:end) - T_nsfdm01));
err_y_01  = max(abs(T_i_rk4(1:1000:end) - T_i_nsfdm01));
err_z_01  = max(abs(V_rk4(1:1000:end) - V_nsfdm01));
err_01    = max([err_x_01 err_y_01 err_z_01]);

% Step 4.2: Solve NSFDM For h(2)

t02            = 0:h(2):T;
T_nsfdm02      = zeros(length(t02),1);
T_i_nsfdm02    = zeros(length(t02),1);
V_nsfdm02      = zeros(length(t02),1);
T_nsfdm02(1)   = T_zero;
T_i_nsfdm02(1) = T_i_zero;
V_nsfdm02(1)   = V_zero;

% Loop

for j = 1:1:(length(t02)-1)
   T_nsfdm02(j+1)   = (T_nsfdm02(j) + h(2)*r)/(1 + h(2)*beta*V_nsfdm02(j) + h(2)*d);
   T_i_nsfdm02(j+1) = (T_i_nsfdm02(j) + h(2)*beta*T_nsfdm02(j+1)*V_nsfdm02(j))/(1 + h(2)*delta);
   V_nsfdm02(j+1)   = (V_nsfdm02(j) + h(2)*pig*T_i_nsfdm02(j+1))/(1 + h(2)*c);
endfor

err_x_02  = max(abs(T_rk4(1:500:end) - T_nsfdm02));
err_y_02  = max(abs(T_i_rk4(1:500:end) - T_i_nsfdm02));
err_z_02  = max(abs(V_rk4(1:500:end) - V_nsfdm02));
err_02    = max([err_x_02 err_y_02 err_z_02]);

% Step 4.3: Solve NSFDM For h(3)

t03            = 0:h(3):T;
T_nsfdm03      = zeros(length(t03),1);
T_i_nsfdm03    = zeros(length(t03),1);
V_nsfdm03      = zeros(length(t03),1);
T_nsfdm03(1)   = T_zero;
T_i_nsfdm03(1) = T_i_zero;
V_nsfdm03(1)   = V_zero;

% Loop

for j = 1:1:(length(t03)-1)
   T_nsfdm03(j+1)   = (T_nsfdm03(j) + h(3)*r)/(1 + h(3)*beta*V_nsfdm03(j) + h(3)*d);
   T_i_nsfdm03(j+1) = (T_i_nsfdm03(j) + h(3)*beta*T_nsfdm03(j+1)*V_nsfdm03(j))/(1 + h(3)*delta);
   V_nsfdm03(j+1)   = (V_nsfdm03(j) + h(3)*pig*T_i_nsfdm03(j+1))/(1 + h(3)*c);
endfor

err_x_03  = max(abs(T_rk4(1:250:end) - T_nsfdm03));
err_y_03  = max(abs(T_i_rk4(1:250:end) - T_i_nsfdm03));
err_z_03  = max(abs(V_rk4(1:250:end) - V_nsfdm03));
err_03    = max([err_x_03 err_y_03 err_z_03]);

% Step 4.4: Solve NSFDM For h(4)

t04            = 0:h(4):T;
T_nsfdm04      = zeros(length(t04),1);
T_i_nsfdm04    = zeros(length(t04),1);
V_nsfdm04      = zeros(length(t04),1);
T_nsfdm04(1)   = T_zero;
T_i_nsfdm04(1) = T_i_zero;
V_nsfdm04(1)   = V_zero;

% Loop

for j = 1:1:(length(t04)-1)
   T_nsfdm04(j+1)   = (T_nsfdm04(j) + h(4)*r)/(1 + h(4)*beta*V_nsfdm04(j) + h(4)*d);
   T_i_nsfdm04(j+1) = (T_i_nsfdm04(j) + h(4)*beta*T_nsfdm04(j+1)*V_nsfdm04(j))/(1 + h(4)*delta);
   V_nsfdm04(j+1)   = (V_nsfdm04(j) + h(4)*pig*T_i_nsfdm04(j+1))/(1 + h(4)*c);
endfor

err_x_04  = max(abs(T_rk4(1:125:end) - T_nsfdm04));
err_y_04  = max(abs(T_i_rk4(1:125:end) - T_i_nsfdm04));
err_z_04  = max(abs(V_rk4(1:125:end) - V_nsfdm04));
err_04    = max([err_x_04 err_y_04 err_z_04]);

##% Step 4.5: Solve NSFDM For h(5)
##
##t05            = 0:h(5):T;
##T_nsfdm05      = zeros(length(t05),1);
##T_i_nsfdm05    = zeros(length(t05),1);
##V_nsfdm05      = zeros(length(t05),1);
##T_nsfdm05(1)   = T_zero;
##T_i_nsfdm05(1) = T_i_zero;
##V_nsfdm05(1)   = V_zero;
##
##% Loop
##
##for j = 1:1:(length(t05)-1)
##   T_nsfdm05(j+1)   = (T_nsfdm05(j) + h(5)*r)/(1 + h(5)*beta*V_nsfdm05(j) + h(5)*d);
##   T_i_nsfdm05(j+1) = (T_i_nsfdm05(j) + h(5)*beta*T_nsfdm05(j+1)*V_nsfdm05(j))/(1 + h(5)*delta);
##   V_nsfdm05(j+1)   = (V_nsfdm05(j) + h(5)*pig*T_i_nsfdm05(j+1))/(1 + h(5)*c);
##endfor
##
##err_x_05  = max(abs(T_rk4(1:125:end) - T_nsfdm05));
##err_y_05  = max(abs(T_i_rk4(1:125:end) - T_i_nsfdm05));
##err_z_05  = max(abs(V_rk4(1:125:end) - V_nsfdm05));
##err_05    = max([err_x_05 err_y_05 err_z_05]);

##% Step 4.6: Solve NSFDM For h(6)
##
##t06            = 0:h(6):T;
##T_nsfdm06      = zeros(length(t06),1);
##T_i_nsfdm06    = zeros(length(t06),1);
##V_nsfdm06      = zeros(length(t06),1);
##T_nsfdm06(1)   = T_zero;
##T_i_nsfdm06(1) = T_i_zero;
##V_nsfdm06(1)   = V_zero;
##
##% Loop
##
##for j = 1:1:(length(t06)-1)
##   T_nsfdm06(j+1)   = (T_nsfdm06(j) + h(6)*r)/(1 + h(6)*beta*V_nsfdm06(j) + h(6)*d);
##   T_i_nsfdm06(j+1) = (T_i_nsfdm06(j) + h(6)*beta*T_nsfdm06(j+1)*V_nsfdm06(j))/(1 + h(6)*delta);
##   V_nsfdm06(j+1)   = (V_nsfdm06(j) + h(6)*pig*T_i_nsfdm06(j+1))/(1 + h(6)*c);
##endfor
##
##err_x_06  = max(abs(T_rk4(1:125:end) - T_nsfdm06));
##err_y_06  = max(abs(T_i_rk4(1:125:end) - T_i_nsfdm06));
##err_z_06  = max(abs(V_rk4(1:125:end) - V_nsfdm06));
##err_06    = max([err_x_06 err_y_06 err_z_06]);

% Step 5: Error Plot

err_h = h;
err   = [err_01 err_02 err_03 err_04];

figure(4)
loglog(err_h,err,'linestyle','-','marker','x');
title('Convergence Plot For Non-Standard Finite-Difference-Method');
xlabel('h');
ylabel('error');