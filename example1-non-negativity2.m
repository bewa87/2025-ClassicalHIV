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

T            = 200;
h            = 0.1;
t            = 0:h:T;
t            = t';
T_ee         = zeros(length(t),1);
T_rk2        = zeros(length(t),1);
T_nsfdm      = zeros(length(t),1);
T_i_ee       = zeros(length(t),1);
T_i_rk2      = zeros(length(t),1);
T_i_nsfdm    = zeros(length(t),1);
V_ee         = zeros(length(t),1);
V_rk2        = zeros(length(t),1);
V_nsfdm      = zeros(length(t),1);
T_ee(1)      = T_zero;
T_rk2(1)     = T_zero;
T_nsfdm(1)   = T_zero;
T_i_ee(1)    = T_i_zero;
T_i_rk2(1)   = T_i_zero;
T_i_nsfdm(1) = T_i_zero;
V_ee(1)      = V_zero;
V_rk2(1)     = V_zero;
V_nsfdm(1)   = V_zero;

% Step 3.1: Solve Explicit Eulerian Scheme

for j = 1:1:(length(t)-1)
  T_ee(j+1)   = T_ee(j)   + h*(r - beta*V_ee(j)*T_ee(j) - d*T_ee(j));
  T_i_ee(j+1) = T_i_ee(j) + h*(beta*V_ee(j)*T_ee(j) - delta*T_i_ee(j));
  V_ee(j+1)   = V_ee(j)   + h*(pig*T_i_ee(j) - c*V_ee(j));
endfor

% Step 3.2: Solve RK2 Scheme

f1 = @(x,y,z) r - beta*x*z - d*x;
f2 = @(x,y,z) beta*x*z - delta*y;
f3 = @(x,y,z) pig*y - c*z;

for j = 1:1:(length(t)-1)
   k1_x         = f1(T_rk2(j),T_i_rk2(j),V_rk2(j));
   k1_y         = f2(T_rk2(j),T_i_rk2(j),V_rk2(j));
   k1_z         = f3(T_rk2(j),T_i_rk2(j),V_rk2(j));
   k2_x         = f1(T_rk2(j) + h*k1_x, T_i_rk2(j) + h*k1_y, V_rk2(j) + h*k1_z);
   k2_y         = f2(T_rk2(j) + h*k1_x, T_i_rk2(j) + h*k1_y, V_rk2(j) + h*k1_z);
   k2_z         = f3(T_rk2(j) + h*k1_x, T_i_rk2(j) + h*k1_y, V_rk2(j) + h*k1_z);
   T_rk2(j+1)   = T_rk2(j)   + (h/2)*(k1_x+k2_x);
   T_i_rk2(j+1) = T_i_rk2(j) + (h/2)*(k1_y+k2_y);
   V_rk2(j+1)   = V_rk2(j)   + (h/2)*(k1_z+k2_z);
endfor

% Step 3.3: Solve NSFDM Scheme

for j = 1:1:(length(t)-1)
   T_nsfdm(j+1)   = (T_nsfdm(j) + h*r)/(1 + h*beta*V_nsfdm(j) + h*d);
   T_i_nsfdm(j+1) = (T_i_nsfdm(j) + h*beta*T_nsfdm(j+1)*V_nsfdm(j))/(1 + h*delta);
   V_nsfdm(j+1)   = (V_nsfdm(j) + h*pig*T_i_nsfdm(j+1))/(1 + h*c);
endfor

% Step 4: Plotting Of Solutions

figure(2)

subplot(3,1,1);
plot(t,T_ee);
title('Load of Target CD4-cells (Explicit Eulerian Scheme, h = 0.1)');
xlabel('t [days]');
ylabel('T(t)','interpreter','tex');
yticks([0 4 8 12]);

subplot(3,1,2);
plot(t,T_rk2);
title('Load of Target CD4-cells (Second-Order Runge-Kutta Scheme, h = 0.1)');
xlabel('t [days]');
ylabel('T(t)','interpreter','tex');
yticks([0 4 8 12]);

subplot(3,1,3);
plot(t,T_nsfdm);
title('Load of Target CD4-cells (Non-Standard Finite-Difference-Method, h = 0.1)');
xlabel('t [days]');
ylabel('T(t)','interpreter','tex');
yticks([0 4 8 12]);