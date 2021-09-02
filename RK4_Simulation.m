% close all;
addpath('Config')
addpath('Config\casadi-windows-matlabR2016a-v3.5.5')
beep off;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% RK4 simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code aims at performing a simple Runge-Kutta 4 simulation of the   %
% system dynamics. Several preferences can be set :                       %
% - 'dynFlag' : simply specifies whether heat transfers should accounted  %
% for or not. Keep to 'H-T' to have them ON                               %
% - 'VarHT' : specifies whether heat transfers according to the Nusselt   %
% correlation should used instead of constant heat transfer coefficients  %
% - 'switching' : specifies whether switching the internal heat transfer  %
% coefficient is desired                                                  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General plotting preferences

set(0, 'DefaultLineLineWidth', 1.2);
set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%% Physical and geometrical parameters

dynFlag = 'H-T'; % Include heat transfer
VarHT = 1; % Time-varying heat transfer ?
switching = 1; % Internal heat transfer coefficient switching

% Simulation parameters retrieved from configuration file
[V2, L2, D2, V1, L1, A_int2, D1, A_tube, Ta, P1, P2, T1, T2, ...
    u2, u1, rho1, rho2, m1, m2, h_conv, h1, ...
    h_in, h_out, S_in, S_out, c_CFRP, c_metal, m_CFRP, m_metal, T_wall, ...
    k_CFRP, k_metal, t_CFRP, t_metal] = Parameters(dynFlag, 1);

x_init = [m2;u2;m1;u1;Ta;Ta]; % Real initial condition

% Use interpolated tables
[D, hx, U, hy, T, P, H] = Tables(V2);

%% Simulation parameters

% Choose the totale simulation time. Time-step will be adapted to 0.05 [s].
time = 100; % [s]
hr = 0.05; % [s]
Nr = time/hr; % Nb. of simulation steps

%% Simulation parameters 

% Choose the simulation input to use :
% input = [linspace(0,0.0005,Nr/2) 0.0005*ones(1,Nr/2)]; % [kg/s]
% input = linspace(0.005,0.0001,Nr);
input = 0.00165*ones(Nr,1);

% Dynamics definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Heat transfer parameters for calculating Re and Ra numbers
g = 9.81; % gravity acceleration constant , [m/s^2]
% For saving the heat transfer coefficients history
UA_IN_hist = zeros(Nr,1);
UA_OUT_hist = zeros(Nr,1);
UA_int_hist = zeros(Nr,1);

% Dynamics parameters
UA_IN = 300;
UA_int = 250;
UA_OUT = 4.5;
% Parameters for the dynamics function evaluation
param = [c_CFRP c_metal k_CFRP k_metal h_in h_out m_CFRP m_metal ...
                t_CFRP t_metal D2 L2 Ta UA_IN UA_OUT UA_int];

%% Integration

disp('====================================================================')
disp(['Integration starts at : ',num2str(month(datetime)), ...
    '/',num2str(day(datetime)),' : ',num2str(hour(datetime)),'h',num2str(minute(datetime)), ...
    ' ',num2str(second(datetime)),' s'])
disp('====================================================================')
tic;

% State vector : [x1(1) ... x1(Nr) x2(1) ... x2(Nr) ... ... x6(1) ... x6(Nr)]
X = zeros(6*Nr,1);
% Initialise simulation
X(1:Nr:6*Nr) = x_init;

for i = 1:Nr - 1
        
    % If switching internal heat transfer simulation is desired
    if switching
            
        if i*hr > 60 input(i) = 0; UA_IN = 250; end
        
    end
    
    % Actual states
    x = X(i:Nr:6*Nr);
    
    % Temperature evaulation
    T2(i) = BilinearInterpolation(T,x(1)/V2,x(2),D,U);
    
    % Pressure evaluation
    P2(i) = BilinearInterpolation(P,x(1)/V2,x(2),D,U);

    % Parameters for heat transfers in tank 2
    
    if VarHT % If time-varying heat transfer desired
        
        % Volumetric expansion coefficient [1/K]
        beta = refpropm('B','D',x(1)/V2,'U',x(2),'hydrogen'); 
        % Constant pressure specific heat capacity [J/(kg K)]
        Cp = refpropm('C','D',x(1)/V2,'U',x(2),'hydrogen'); 
        % Fluid thermal conductivity [W/(m K)]
        lambda = refpropm('L','D',x(1)/V2,'U',x(2),'hydrogen'); 
        % Absolute (or dynamic) viscosity
        mu = refpropm('V','D',x(1)/V2,'U',x(2),'hydrogen'); 
        
        % Rayleigh number
        RaL = g*beta*(T2(i) - x(5))*Cp*(x(1)/V2)^2*L2^3/(mu*lambda);
        % Reynolds number
        ReD = 4*input(i)/(mu*pi*D2);
        % Nusselt number 
        Nu = 0.56*ReD^(0.67) + 0.104*RaL^(0.352);
        % Internal mixed heat transfer coefficient
        h_in = Nu*lambda/L2;
        
        param = [c_CFRP c_metal k_CFRP k_metal h_in h_out m_CFRP m_metal ...
                t_CFRP t_metal D2 L2 Ta UA_IN UA_OUT UA_int];
        
    end
    
    % Fluid velocity at tank 2 inlet
    v1 = input(i)/(x(1)/V2)/A_tube;
    % Enthalpy conservation due to isenthalpic expansion
    h1(i) = BilinearInterpolation(H,x(3)/V1,x(4),D,U);

    % RK4 integration

    k1 = Dynamics_Sim(x, input(i), T2(i), h1(i), v1, param, VarHT);
    k1 = k1(1:6); % 3 last output being the heat transfers
    k2 = Dynamics_Sim(x + hr/2*k1, input(i), T2(i), h1(i), v1, param, VarHT);
    k2 = k2(1:6);
    k3 = Dynamics_Sim(x + hr/2*k2, input(i), T2(i), h1(i), v1, param, VarHT);
    k3 = k3(1:6);
    k4 = Dynamics_Sim(x + hr*k3, input(i), T2(i), h1(i), v1, param, VarHT);
    k4 = k4(1:6);
    X(i + 1:Nr:6*Nr) = x + hr/6*(k1 + 2*k2 + 2*k3 + k4);
    
    % Just retrieving the heat transfer coefficients history
    HTcoeff = Dynamics_Sim(x, input(i), T2(i), h1(i), v1, param, VarHT);
    UA_IN_hist(i) = HTcoeff(7);
    UA_OUT_hist(i) = HTcoeff(8);
    UA_int_hist(i) = HTcoeff(9);
    
    disp(['Iteration =================== ',num2str(i),' / ',num2str(Nr - 1)])
end

% Retrieve last temp. and pressure values
T2(Nr) = BilinearInterpolation(T,X(Nr)/V2,X(2*Nr),D,U);
P2(Nr) = BilinearInterpolation(P,X(Nr)/V2,X(2*Nr),D,U);

t_der1 = toc;
disp(['Integration completed. Duration : ',num2str(t_der1),' s'])
disp('====================================================================')

%% Results

% Time vector for plotting
time = linspace(0,(Nr - 1)*hr,Nr);

fig1 = figure(1);

yyaxis left
p3 = plot(time,X(1:Nr),'-');
xlabel('Time [s]')
ylabel('Mass [kg]');
% ylim([0.0028 0.0045])

yyaxis right
p4 = plot(time,input,'-');
ylabel('Mass flow rate [kg/s]');
% ylim([2.69e6 2.73e6])

title('Tank 2 temperature and pressure evolutions')
legend([p3 p4],'$m_2$', '$\dot{m} = u$','location','best');

fig2 = figure(3);

yyaxis left
plot(time,T2,'linewidth',1.2);
hold on, grid on
plot(time,X(4*Nr + 1:5*Nr),'linewidth',1.2)
plot(time,X(5*Nr + 1:6*Nr),'linewidth',1.2)
xlabel('Time [s]')
ylabel('$Temperature \ [K]$');

yyaxis right
plot(time,P2,'linewidth',1.2);
xlabel('Time [s]')
ylabel('$Pressure \ [kPa]$');
ylim([1500 1.03*max(P2)])

title('Tank 2 gas and lining temperatures and gas pressure')
legend('$T_2$','$T_{alu}$','$T_{CFRP}$','$P_2$','location','best');