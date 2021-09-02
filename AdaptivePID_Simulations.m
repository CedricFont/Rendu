close all;
beep off;
addpath('Config')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%% Train optimal PID controller %%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program simulates the adaptive I-PID. First, the adaptive PID is   %
% loaded from an external file, and then a different heat transfer        %
% coefficient is chosen in order to stimulate it.                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Physical and geometrical parameters

dynFlag = 'H-T';
flag_control = 'I-PID';
[V2, L2, D2, V1, L1, A_int2, D1, A_tube, Ta, P1, P2, T1, T2, ...
    u2, u1, rho1, rho2, m1, m2, h_conv, h1, ...
    h_in, h_out, S_in, S_out, c_CFRP, c_metal, m_CFRP, m_metal, T_wall, ...
    k_CFRP, k_metal, t_CFRP, t_metal] = Parameters(dynFlag,0);
x_init = [m2;u2;m1;u1;Ta;Ta]; % Real initial condition
[D, hx, U, hy, T, P, H, K_ratio] = Tables(V2);

switching_time = 20;
t_Ki = 4;

%% Simulation parameters and initialisation

% Controller definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First input is zero since nothing has been measured yet

K_param = [switching_time, t_Ki];

% File from which adaptive controller will be retrieved
filename = 'Config\Adaptive_IPID.xls';
N = 100; % Nb. of multiple shooting time-step
N_sub = 20; % Nb. of sub-integration time-steps
h = 1; % Here, the time-step is fixed because this is the maximum time-step 
% found for the lowest heat transfer coefficient
h_sub = h/N_sub; % Sub-int. time-step
time = N*h;

% Non-parametric models associated with different PID controllers
T_model = zeros(3,N);
K_PID = zeros(3,4);

for i = 1:11
    
    K_PID(i,:) = readmatrix(filename,'Sheet',strcat('Sheet',num2str(i)),'Range','F2:I2');
    T_model(i,:) = readmatrix(filename,'Sheet',strcat('Sheet',num2str(i)),'Range',strcat('S2:S',num2str(N + 1)));

end

% Path variables (to store properties)
T2 = zeros(N,1);
h1 = zeros(N,1); 
input = zeros(N,1);
Q_IN = zeros(N,1);
Q_OUT = zeros(N,1);
Q_int = zeros(N,1);
X = zeros(6*N,1);

r = 0.267; % Reference

% Dynamics definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncertain parameters
UA_OUT = 4.5;
UA_int = 250;
UA_IN = 250;
% Can be used for testing several heat transfer coeff.
% UA_IN = 250:10:350; 
            
T_max = 338; % Temp. constraint

%% Optimisation problem

PlotAll = 0; % To have lots of plots
adaptive = 1; % Simulate with adaptive feature ON ?

% Choose nominal controller
% Here, UA_IN = [300,250,350]
nb_models = size(T_model,1);
middleK = 6; % Middle model as symmetry point
K = K_PID(6,:); % Starting model
lambda = 0.8; % RLS forgetting factor
RLS_error = zeros(nb_models,1);
RLS_save = zeros(nb_models,N);
hist = 9; % Hysteresis time
model = zeros(N,1); % Save models history
changes = []; % Save models changes
current = middleK; % Current models
HasChanged = 0;
Nsim = length(UA_IN);

% If several heat transfers are desired to be tested sequentially, all
% plots will stand on the same figure
for s = 1:Nsim
    
    param = [c_CFRP c_metal k_CFRP k_metal h_in h_out m_CFRP m_metal ...
                t_CFRP t_metal D2 L2 Ta UA_IN(s) UA_OUT UA_int];

    x_final = x_init;

    disp(['============= Integration started =============']);

    for i = 1:N - 1
        
        X(i:N:6*N) = x_final;

        T2(i) = BilinearInterpolation(T,x_final(1,1)/V2,x_final(2,1),D,U);

        input(i + 1) = ControlLaw(X(1:N), i, K, h, r, K_param, flag_control);
        
        param = [c_CFRP c_metal k_CFRP k_metal h_in h_out m_CFRP m_metal ...
                t_CFRP t_metal D2 L2 Ta UA_IN UA_OUT UA_int];

        for k = 1:N_sub

            x = x_final;

            % Parameters and lookup-tables evaluations

            rho2 = x(1)/V2;
            U2 = x(2);
            rho1 = x(3)/V1;
            U1 = x(4);
            T2_val = BilinearInterpolation(T,rho2,U2,D,U);
            v1 = input(i)/rho1/A_tube;
            h1(i) = BilinearInterpolation(H,rho1,U1,D,U);

            % RK4 integration
            k1 = Dynamics_Sim(x, input(i), T2_val, h1(i), v1, param, 0);
            k1 = k1(1:6);
            k2 = Dynamics_Sim(x + h_sub/2*k1, input(i), T2_val, h1(i), v1, param, 0);
            k2 = k2(1:6);
            k3 = Dynamics_Sim(x + h_sub/2*k2, input(i), T2_val, h1(i), v1, param, 0);
            k3 = k3(1:6);
            k4 = Dynamics_Sim(x + h_sub*k3, input(i), T2_val, h1(i), v1, param, 0);
            k4 = k4(1:6);
            x_final = x + h_sub/6*(k1 + 2*k2 + 2*k3 + k4);

        end
        
        if adaptive
        
            % Controller adaptation
            RLS_error = lambda*RLS_error + (T_model(:,i) - repmat(T2(i),nb_models,1)).^2;
            RLS_save(:,i) = RLS_error;
            [~,model(i)] = min(RLS_error);

            if i > hist
                
                
                for j = 1:nb_models

                    if model(i - hist:i)' == j*ones(hist,1)

                        K = K_PID(middleK - (j - middleK),:);
                        if current ~= j
                            changes = [changes ; i];
                        end
                        current = j;

                    end

                end

            end
            
        end

        % Evaluate heat transfer
        HT = Dynamics_Sim(X(i:N:6*N), input(i), T2_val, h1(i), v1, param, 0);
        Q_IN(i) = HT(7);
        Q_OUT(i) = HT(8);
        Q_int(i) = HT(9);

        disp(['Iteration ================== ',num2str(i),'/',num2str(N)]);

    end

    X(N:N:6*N,1) = x_final;

    rho2 = X(N,1)/V2;
    U2 = X(2*N,1); 
    T2(N) = BilinearInterpolation(T,rho2,U2,D,U);

    %% Plotting routine

    rho1 = X(2*N + 1:3*N)/V1; % [kg/m3], density vector
    T1 = zeros(N,1); 
    P1 = zeros(N,1);
    rho2 = X(1:N)/V2; % [kg/m3], density vector
    P2 = zeros(N,1);
    u = zeros(N,1);


    for i = 1:N

        T1(i) = BilinearInterpolation(T,rho1(i),X(3*N + i),D,U);
        P1(i) = BilinearInterpolation(P,rho1(i),X(3*N + i),D,U);
        P2(i) = BilinearInterpolation(P,rho2(i),X(N + i),D,U);

    end

    time = linspace(0,N*h,N);
    x = X;
    %% Results plotting

    linewidth = 1.2;
    set(0, 'DefaultLineLineWidth', linewidth);
    set(0,'defaultTextInterpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    delete(findall(gcf,'type','annotation'))
    
    error = floor(100*(r - x(N))/(r - x_init(1)));
    
    if PlotAll
        
        fig1 = figure(1);
        
        subplot(3,2,1)
        p1 = plot(time,x(1:N));
        p2 = yline(r,'r--');
        xlabel('Time [s]')
        ylabel('Mass in tank 2 [kg]');
        ylim([0.8*x_init(1) 1.2*r])
        legend([p1 p2],'$m_2$','$m_{target}$');

        subplot(3,2,2)
        p2 = plot(time,x(N + 1:2*N));
        xlabel('Time [s]')
        ylabel('Specific internal energy in tank 2 [J/kg]')
        legend([p2],'$u_2$')

        subplot(3,2,3)
        p1 = plot(time,T2);
        p2 = yline(T_max,'r--');
        xlabel('Time [s]')
        ylabel('Temperature in tank 2 [K]')
        ylim([290 T_max*1.1])
        legend('$T_2$','$T_{max}$')

        subplot(3,2,4)
        p2 = plot(time,P2);
        xlabel('Time [s]')
        ylabel('Pressure in tank 2 [kPa]');
        legend([p2],'$P_2$');

        subplot(3,2,5)
        p2 = stairs(time,input,'linewidth',1.2);
        xlabel('Time [s]')
        ylabel('Mass flow rate [kg/s]');
        legend('$u_{control}(k+1) = K_{P} \cdot e(k) + K_I \cdot \sum_{i=k-n_I}^{k}e(i) \delta t + K_D \cdot \frac{y(k)-y(k - 1)}{\delta t}$', ...
            'location','best');
        annotation('textbox', [0.55, 0.22, 0.35, 0.12], 'string', ...
            sprintf('$K_p = %.7f$ \n$K_{I1} = %.7f$ \n$K_d = %.7f$ \n$K_{I2} = %.7f$ \n$\\Delta t = %.7f \\ [s]$ ', ...
            Kp,Ki1,Kd,Ki2,h),'Interpreter','latex');
        annotation('textbox', [0.55, 0.1, 0.35, 0.12], 'string', ...
            sprintf('$m_2(N) = %.5f \\ [kg]$ \n$\\epsilon_{m_2} = %.5f$ %s (error)',x(N), ...
            100*(r - x(N))/(r - x_init(1)),'\%'),'Interpreter','latex');

        fig2 = figure(2);

        subplot(2,2,1)
        p1 = plot(time,x(2*N + 1:3*N));
        xlabel('Time [s]')
        ylabel('Mass in tank 1 [kg]');
        legend([p1],'$m_1$');

        subplot(2,2,2)
        p2 = plot(time,x(3*N + 1:4*N));
        xlabel('Time [s]')
        ylabel('Specific internal energy in tank 1 [J/kg]')
        legend([p2],'$u_1$')

        subplot(2,2,3)
        p1 = plot(time,T1);
        xlabel('Time [s]')
        ylabel('Temperature in tank 1[K]')
        legend([p1],'$T_1$')

        subplot(2,2,4)
        p2 = plot(time,P1);
        xlabel('Time [s]')
        ylabel('Pressure in tank 1 [kPa]');
        legend([p2],'$P_1$');
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    blue = [0,0,1];
    green = [0,0.5,0];
    red = [0.8,0,0];
    orange = [1,0.549,0];
    
    color1 = blue;
    color2 = green;
    
    line1 = '-';
    line2 = '--';
    line3 = '.';

    fig3 = figure(3);
    title('Temperatures and heat-transfers within the tank lining')

    hold on
    yyaxis left
    p1 = plot(time, x(4*N + 1:5*N), line1, 'color', color2);
    hold on, grid on
    p2 = plot(time, x(5*N + 1:6*N), line2, 'color', color2);
    xlabel('Time [s]')
    ylabel('Temperature \ [K]')

    yyaxis right

    hold on
    p3 = plot(time, Q_IN, line1, 'color', color1);
    hold on, grid on
    p4 = plot(time, Q_OUT, line2, 'color', color1);
    p5 = plot(time, Q_int, line3, 'color', color1);
    ylabel('Heat flux $\dot{Q} \ [W]$')
    ax = gca;
    ax.YAxis(1).Color = color1;
    ax.YAxis(2).Color = color2;

    fig5 = figure(5);
    title('Control input and system output (mass)')

    hold on
    yyaxis left
    p6 = stairs(time,input,line1,'linewidth',linewidth, 'color', color2);
    grid on, hold on
    if ~isempty(changes)
        for j = 1:length(changes)

            p18 = xline(changes(j)*h, 'k-','linewidth',1.2);

        end
    end
    xlabel('Time [s]')
    ylabel('Mass flow rate $\dot{m} = u \ [kg/s]$');
    ylim([0 1.4*max(input)])

    hold on
    yyaxis right
    p7 = plot(time,x(1:N), line1, 'color', color1);
    hold on, grid on
    p8 = yline(r,'r--','linewidth',linewidth);
    ylabel('Mass [kg]');
    ylim([0.8*x_init(1) 1.4*r])
    ax = gca;
    ax.YAxis(1).Color = color1;
    ax.YAxis(2).Color = color2;

    fig6 = figure(6);
    title('Tank internal pressures and temperatures')

    hold on
    yyaxis left
    p9 = plot(time,T2,line1, 'color', color1);
    hold on, grid on
    p10 = yline(T_max,'r--','linewidth',1.2);
    xlabel('Time [s]')
    ylabel('$Temperature \ [K]$');
    ylim([290 1.03*max(T2)])

    hold on
    yyaxis right
    p11 = plot(time,P2, line1, 'color', color2);
    grid on
    ylabel('$Pressure \ [kPa]$');
    ylim([0 1.2*max(P2)])
    ax = gca;
    ax.YAxis(1).Color = color1;
    ax.YAxis(2).Color = color2;
    

end

legend([p1 p2 p3 p4 p5],'$T_{alu}$', '$T_{CFRP}$','$\dot{Q}_{IN}$','$\dot{Q}_{OUT}$', ...
        '$\dot{Q}_{int}$','location','best')
legend([p6 p7 p8],'$u_{control}(k+1) = \Big ( K_{I1} \cdot \sum_{i=1}^{k}e(i) \delta t \Big )_{k = 1}^{t_1} + \Big ( K_{P} \cdot e(k) + K_{I2} \cdot \sum_{i=k-n_I}^{k}e(i) \delta t + K_D \cdot \frac{y(k)-y(k - 1)}{\delta t} \Big )_{k=t_1}^{N}$', ...
        'Mass in tank 2 $m_2$','$m_{target}$', 'location', 'best', ...
        'FontSize',6,'location','best');
legend([p9 p10 p11],'$T_2$', '$T_{max}$', '$P_2$','location','northwest','location','best');

legend('$T_{alu}$', '$T_{CFRP}$','$\dot{Q}_{IN}$','$\dot{Q}_{OUT}$', ...
        '$\dot{Q}_{int}$','location','best')
legend('$u_{control}(k+1) = \Big ( K_{I1} \cdot \sum_{i=1}^{k}e(i) \delta t \Big )_{k = 1}^{t_1} + \Big ( K_{P} \cdot e(k) + K_{I2} \cdot \sum_{i=k-n_I}^{k}e(i) \delta t + K_D \cdot \frac{y(k)-y(k - 1)}{\delta t} \Big )_{k=t_1}^{N}$', ...
        'Controller switch','Mass in tank 2 $m_2$','$m_{target}$', 'location', 'best', ...
        'FontSize',6,'location','best');
legend('$T_2$', '$T_{max}$', '$P_2$','location','northwest','location','best');

