beep off;
addpath('Config')
addpath('Config\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

%% Physical and geometrical parameters

dynFlag = 'H-T'; % Heat transfer desired
flag_control = 'Free'; % Specifies the OCP bounds
FullyBounded = 1; % Bound pressure, mass-flow rate and restrain mass trajectory

[V2, L2, D2, V1, L1, A_int2, D1, A_tube, Ta, P1, P2, T1, T2, ...
    u2, u1, rho1, rho2, m1, m2, h_conv, h1, ...
    h_in, h_out, S_in, S_out, c_CFRP, c_metal, m_CFRP, m_metal, T_wall, ...
    k_CFRP, k_metal, t_CFRP, t_metal] = Parameters(dynFlag,1);

u0 = 0; % Starts with valve closed
x_init = [m2;u2;m1;u1;Ta;Ta;u0]; % Real initial condition

% All tables needed for the interpolation
[D, hx, U, hy, T, P, H, K_ratio] = Tables(V2);

%% Simulation parameters

% Since we optimize for time, total time is not known, we just have to
% choose a certain number of time steps
N = 100; % Nb. of multiple-shooting intervals
N_sub = 20; % Nb. of sub-integration intervals

%% Optimisation parameters
import casadi.*

% Defines optimal control problem bounds
[lb,ub] = OCP_Bounds(N,dynFlag,V2,flag_control);

% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opti = casadi.Opti(); % Optimisation problem definition

% The 6 system states + the free control input + time-step
% M = [x1(1) ... x1(Nr) x2(1) ... x2(Nr) ... ... x6(1) ... x6(Nr) u(1) ... u(Nr) h]
M = opti.variable(7*N + 1,1);
epsilon = opti.variable(1,1); % Slack for mass reference
slack_P = opti.variable(1,1); % Slack for pressure limit

% Optimisation parameters (used as parameters so that they can be changed
% later in the code without having to run the jacobian and hessian
% constructions when we change them
r = opti.parameter(); % Mass reference
Q_eps = opti.parameter(); % Mass slack cost coeff.
Q_t = opti.parameter(); % Time cost coeff.
UA_IN = opti.parameter(); % Internal heat transfer coefficient
UA_OUT = opti.parameter(); % Outer heat transfer coefficient
UA_int = opti.parameter(); % Inside heat transfer coefficient
T_max = opti.parameter(); % Max. temp. allowed
P_max = opti.parameter(); % Max. pressure allowed, with slack
C_d = opti.parameter(); % Discharging coeff.
A_pipe = opti.parameter(); % Piping area
Q_m = opti.parameter(); % Mass ref. cost coeff.
Q_p = opti.parameter(); % Pressure slack cost coeff.
max_slack_P = opti.parameter(); % Bound for pressure slack

% Reference mass : r. Thus, set initial mass flow rate as the one th;at
% allow to reach our target within a time of N*dt_init
h_init = ub(7*N + 1); % Starts as high as possible
ref_init = 0.267; % Reference initialisation
ref_mass_flow_rate = ref_init/(N*h_init);
U_init = ref_mass_flow_rate*ones(N,1); % Free input initialisation

% Initial optimiser values
X0 = [m2*ones(N,1);u2*ones(N,1);m1*ones(N,1);u1*ones(N,1);Ta*ones(2*N,1); ...
    U_init;h_init];

% Unscaled variable (for computational speed and accuracy)
Mus = M.*X0;

h = Mus(7*N + 1); % Time-step
h_sub = h/N_sub; % Sub-integration time step
input = Mus(6*N + 1:7*N); % Free control input vector

% Path variables (to store properties)
T2 = MX.zeros(N,1);
h1 = MX.zeros(N,1); 
P1 = MX.zeros(N,1);
P2 = MX.zeros(N,1);
Q_IN = MX.zeros(N,1);
Q_OUT = MX.zeros(N,1);
Q_int = MX.zeros(N,1);
q_dot_max = MX.zeros(N,1);

% Dynamics definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dynamics parametrisation

param = [c_CFRP c_metal k_CFRP k_metal h_in h_out m_CFRP m_metal ...
                t_CFRP t_metal D2 L2 Ta UA_IN UA_OUT UA_int];

x_dot = Dynamics(param,0); % Dynamics Casadi function

% Construct the appropriate interpolant objects : Temperature, specific
% enthalpy, pressure and specific heat ratio
inter_opts = struct('inline',true); % Allows jacobian computation w.r.t. coeff.

% Grid vectors
xgrid = D;
ygrid = U;

% Interpolant objects definition
Temp = casadi.interpolant('Temp','bspline',{xgrid, ygrid},T(:),inter_opts);
Enth = casadi.interpolant('Enth','bspline',{xgrid, ygrid},H(:),inter_opts);
Pressure = casadi.interpolant('Pressure','bspline',{xgrid, ygrid},P(:),inter_opts);
Gamma = casadi.interpolant('Gamma','bspline',{xgrid, ygrid},K_ratio(:),inter_opts);

%% Optimisation problem

disp('====================================================================')
disp(['Just began jacobian and gradient computation at time : ',num2str(month(datetime)), ...
    '/',num2str(day(datetime)),' : ',num2str(hour(datetime)),'h',num2str(minute(datetime)), ...
    ' ',num2str(second(datetime)),' s'])
disp('====================================================================')

tic;

% Constraints on M instead of Mus in the aim of ensuring that they are
% respected (otherwise they might be too small and thus ignored)
opti.subject_to( {(lb./X0 <= M) , (M <= ub./X0)} );

% Enforce initial conditions
opti.subject_to( x_init == Mus(1:N:7*N) ); 

for i = 1:N - 1
        
    x_final = Mus(i:N:6*N); % Current states for sub-integration start
        
    T2(i) = Temp([x_final(1)/V2 x_final(2)]);
    opti.subject_to( T2(i) <= T_max ); % Temp. constraint
    P2(i) = Pressure([x_final(1)/V2 x_final(2)]);
    
    % Control input constraints (non-negative + limited by hydrogen
    % expansion rate)
    gamma = Gamma([x_final(3)/V1 x_final(4)]);
    P1(i) = Pressure([x_final(3)/V1 x_final(4)]);
    rho1 = x_final(3)/V1;
    
    opti.subject_to( input(i) >= 0); % Input non-negativity
    q_dot_max(i) = C_d*(rho1*gamma*P1(i))^(1/2)*A_pipe*(2/(gamma + 1))^((gamma + 1)/(2*(gamma - 1)));
    
    if FullyBounded
        
        opti.subject_to( P2(i) <= P_max + slack_P); % pressure constraint
        opti.subject_to( input(i) <= q_dot_max(i) ); % Hydrogen max. expansion rate
        
    end
    
    % Start sub-integration
    for k = 1:N_sub
        
        x = x_final;
        
        % Parameters and lookup-tables evaluations
        
        rho2 = x(1)/V2;
        U2 = x(2);
        rho1 = x(3)/V1;
        U1 = x(4);
        T2_val = Temp([rho2 U2]);
        v1 = input(i)/rho1/A_tube;
        h1(i) = Enth([rho1 U1]);
        
        % RK4 sub-integration
        
        k1 = x_dot(x, input(i), T2_val, h1(i), v1);
        k1 = k1(1:6);
        k2 = x_dot(x + h_sub/2*k1, input(i), T2_val, h1(i), v1);
        k2 = k2(1:6);
        k3 = x_dot(x + h_sub/2*k2, input(i), T2_val, h1(i), v1);
        k3 = k3(1:6);
        k4 = x_dot(x + h_sub*k3, input(i), T2_val, h1(i), v1);
        k4 = k4(1:6);
        x_final = x + h_sub/6*(k1 + 2*k2 + 2*k3 + k4);
        
    end
    
    % Dynamics continuity constraints
    opti.subject_to( Mus(i + 1:N:6*N) == x_final );
    
    % Evaluate heat transfer
    HT = x_dot(Mus(i:N:6*N), input(i), T2(i), h1(i), v1);
    Q_IN(i) = HT(7);
    Q_OUT(i) = HT(8);
    Q_int(i) = HT(9);
    
end

% Final properties and constraints evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho2 = Mus(N,1)/V2;
U2 = Mus(2*N,1); 
T2(N) = Temp([rho2 U2]);

opti.subject_to( T2(N) <= T_max ); % Final temp. constraint
P2(N) = Pressure([rho2 U2]);
opti.subject_to( input(N) == 0 ); % Last input is zero
opti.subject_to( P2(N) <= P_max + slack_P); % Final pressure constraint

% Pressure slack constraints
opti.subject_to( { (0 <= slack_P) , (slack_P <= max_slack_P) } );

rho1 = Mus(3*N,1)/V1;
U1 = Mus(4*N,1); 
P1(N) = Pressure([rho1 U1]);

gamma = Gamma([rho1 U1]);
q_dot_max(N) = C_d*(rho1*gamma*P1(N))^(1/2)*A_pipe*(2/(gamma + 1))^((gamma + 1)/(2*(gamma - 1)));

% Enforces stabilization at reference (and no overshooting). Mass must be
% within 1 % at most of desired reference

if FullyBounded
    
    opti.subject_to( Mus(N-9:N,1) <= r*ones(10,1));
    opti.subject_to( Mus(N-9:N,1) >= r*ones(10,1) - epsilon);
    opti.subject_to( { (epsilon >= 0) , (epsilon <= 0.005*r) } ); 
    
end

opti.subject_to( Mus(N) >= r );

% Evaluate heat transfer
v1 = input(N)/rho1/A_tube;
HT = x_dot(Mus(N:N:6*N), input(N), T2(N), h1(N), v1);
Q_IN(N) = HT(7);
Q_OUT(N) = HT(8);
Q_int(N) = HT(9);

% Cost variables
time = MX.sym('t',1);
Qeps = MX.sym('Q_eps',1);
Qt = MX.sym('Qt',1);
u_var = MX.sym('u_var',N);
eps = MX.sym('eps',1);
mass = MX.sym('mass',N);
ref = MX.sym('ref',1);
Qm = MX.sym('Qm',1);
Qp = MX.sym('Qp',1);
slackP = MX.sym('slackP',1);

if ~FullyBounded
    
    obj = time;
    objective = Function('obj',{time},{obj});
    opti.minimize(objective(h));
    
else
    
    % Minimize time, error w.r.t. reference and control effort to attain ref.
    obj = Qt*time^2 + eps*Qeps*eps + ...
        (mass - ref)'*Qm*(mass - ref) + slackP^2*Qp;
    objective = Function('obj',{time,eps,Qeps,Qt,mass,ref,Qm,slackP,Qp},{obj});
    opti.minimize(objective(h,epsilon,Q_eps,Q_t,Mus(1:N),r,Q_m,slack_P,Q_p));

end

t_der1 = toc;
disp(['Just finished jacobian and hessian computation. Duration : ',num2str(t_der1),' s'])
disp('====================================================================')

%% Solver setup

p_opts = struct('expand',false);
s_opts = struct('max_iter',1000);
% scaling = struct(, 'nlp_scaling_method','none' 'bound_relax_factor',,'hessian_approximation', 'exact' );

opti.solver('ipopt',p_opts,s_opts); % set numerical backend

%% Solve OCP

tic;
disp(['Just began the solving procedure at time : ',num2str(month(datetime)), ...
    '/',num2str(day(datetime)),' : ',num2str(hour(datetime)),'h',num2str(minute(datetime)), ...
    ' ',num2str(second(datetime)),' s'])
disp('====================================================================')

opti.set_initial(Mus,X0); % Set optimiser initial values

% Define all parameters
opti.set_value(r, 0.267);
opti.set_value(Q_eps,1e4);
opti.set_value(Q_t,1/lb(7*N + 1)^2);
opti.set_value(UA_IN,300);
opti.set_value(UA_OUT,4.5);
opti.set_value(UA_int,250);
opti.set_value(T_max,273+65);
opti.set_value(P_max,30000);
opti.set_value(A_pipe,pi*0.00635^2/4);
opti.set_value(C_d,1);
opti.set_value(Q_m,1/ref_init^2);
opti.set_value(max_slack_P,3000);
opti.set_value(Q_p,1/3000);

solution = opti.solve();   % Actually solve the optimisation problem

t_sim2 = toc;
disp('====================================================================')
disp(['Finished solving in ',num2str(t_sim2),' [s], plotting'])
disp('====================================================================')

%% Plotting
%% Find P and T with tables

% Retrieve unscaled solution
x = solution.value(Mus);

% Retrieve all quantities for plotting
rho1 = x(2*N + 1:3*N)/V1; % [kg/m3], density vector
T1 = zeros(N,1); 
P1 = solution.value(P1);
rho2 = x(1:N)/V2; % [kg/m3], density vector
T2 = solution.value(T2); 
P2 = solution.value(P2);
h1 = solution.value(h1);
h = x(7*N + 1,1);
ref = solution.value(r);
h_in_val = solution.value(h_in);
h_out_val = solution.value(h_out);
T_max_val = solution.value(T_max);
P_max_val = solution.value(P_max);
Q_IN_val = solution.value(Q_IN);
Q_OUT_val = solution.value(Q_OUT);
Q_int_val = solution.value(Q_int);
q_dot_max = solution.value(q_dot_max);

u = solution.value(input);

for i = 1:N
    
    T1(i) = refpropm('T','D',rho1(i),'U',x(3*N + i),'hydrogen');
    
end

% Time vector for plotting
time = linspace(0,N*h,N);

%% Results

% Close figures if already opened
close(findobj('Type','Figure','number',1))
close(findobj('Type','Figure','number',2))
close(findobj('Type','Figure','number',3))
close(findobj('Type','Figure','number',5))
close(findobj('Type','Figure','number',6))

% Set plotting general preferences
set(0, 'DefaultLineLineWidth', 1.4);
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

fig1 = figure(1);
title('Tank 2 states trajectories')

delete(findall(gcf,'type','annotation'))

subplot(3,2,1)
p1 = plot(time,x(1:N));
hold on, grid on
p2 = yline(ref,'r--');
xlabel('Time [s]')
ylabel('Mass in tank 2 [kg]');
ylim([0.8*x_init(1) 1.2*ref])
legend([p1 p2],'$m_2$','$m_{target}$');

subplot(3,2,2)
p2 = plot(time,x(N + 1:2*N));
grid on
xlabel('Time [s]')
ylabel('Specific internal energy in tank 2 [J/kg]')
legend([p2],'$u_2$')

subplot(3,2,3)
p1 = plot(time,T2);
hold on, grid on
p2 = yline(T_max_val,'r--');
xlabel('Time [s]')
ylabel('Temperature in tank 2 [K]')
ylim([290 T_max_val*1.1])
legend('$T_2$','$T_{max}$')

subplot(3,2,4)
p2 = plot(time,P2);
grid on
xlabel('Time [s]')
ylabel('Pressure in tank 2 [kPa]');
legend([p2],'$P_2$');

subplot(3,2,5)
p2 = stairs(time,u,'linewidth',1.2);
hold on, grid on
plot(time, q_dot_max)
xlabel('Time [s]')
ylabel('Mass flow rate [kg/s]');
legend('$u_{control}(k+1) = u^{\ast}(k+1)$', '$\dot{m}_{max}$',...
    'location','best');
annotation('textbox', [0.55, 0.22, 0.35, 0.12], 'string', ...
    sprintf('$\\Delta t = %.7f \\ [s]$ ', ...
    h),'Interpreter','latex');
annotation('textbox', [0.55, 0.1, 0.35, 0.12], 'string', ...
    sprintf('$m_2(N) = %.5f \\ [kg]$ \n$\\epsilon_{m_2} = %.5f$ %s (error)',x(N), ...
    100*(ref - x(N))/(ref - x_init(1)),'\%'),'Interpreter','latex');

error = floor(100*(ref - x(N))/(ref - x_init(1)));

fig2 = figure(2);
title('Tank 1 states trajectories')

subplot(2,2,1)
p1 = plot(time,x(2*N + 1:3*N));
grid on
xlabel('Time [s]')
ylabel('Mass in tank 1 [kg]');
legend([p1],'$m_1$');

subplot(2,2,2)
p2 = plot(time,x(3*N + 1:4*N));
grid on
xlabel('Time [s]')
ylabel('Specific internal energy in tank 1 [J/kg]')
legend([p2],'$u_1$')

subplot(2,2,3)
p1 = plot(time,T1);
grid on
xlabel('Time [s]')
ylabel('Temperature in tank 1[K]')
legend([p1],'$T_1$')

subplot(2,2,4)
p2 = plot(time,P1);
grid on
xlabel('Time [s]')
ylabel('Pressure in tank 1 [kPa]');
legend([p2],'$P_1$');

fig3 = figure(3);
title('Temperatures and heat transfers evolution within tank 2 lining')

yyaxis left
plot(time, x(4*N + 1:5*N))
hold on, grid on
plot(time, x(5*N + 1:6*N))
xlabel('Time [s]')
ylabel('Temperature \ [K]')

yyaxis right

plot(time, Q_IN_val)
hold on
plot(time, Q_OUT_val)
plot(time, Q_int_val)
ylabel('Heat flux $\dot{Q} \ [W]$')
legend('$T_{alu}$', '$T_{CFRP}$','$\dot{Q}_{IN}$','$\dot{Q}_{OUT}$', ...
    '$\dot{Q}_{int}$','location','best')

fig5 = figure(5);
title('Control input and tank 2 mass evolution')

yyaxis left
stairs(time,u,'linewidth',1.4);
hold on, grid on
plot(time, q_dot_max)
xlabel('Time [s]')
ylabel('Mass flow rate $\dot{m} = u \ [kg/s]$');
ylim([0 1.1*max(q_dot_max)])

yyaxis right
p1 = plot(time,x(1:N));
hold on
p2 = yline(ref,'r--');
ylabel('Mass [kg]');
ylim([0.8*x_init(1) 1.4*ref])

legend('$u_{control}(k+1) = u^{\ast}(k+1)$', ...
    'Mass in tank 2 $m_2$','$m_{target}$', 'location', 'best', ...
    'FontSize',10);

fig6 = figure(6);
title('Tank 2 gas temperature and pressure evolution')

yyaxis left
plot(time,T2);
hold on, grid on
yline(T_max_val,'r--','linewidth',1.1)
xlabel('Time [s]')
ylabel('$Temperature \ [K]$');
ylim([280 1.1*max(T2)])

yyaxis right
plot(time,P2);
ylabel('$Pressure \ [kPa]$');
ylim([0 1.2*max(P2)])

legend('$T_2$', '$T_{max}$', '$P_2$','location','northwest');