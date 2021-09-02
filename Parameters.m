function [V2, L2, D2, V1, L1, A_int2, D1, A_tube, Ta, P1, P2, T1, T2, ...
    u2, u1, rho1, rho2, m1, m2, h_conv, h1, ...
    h_in, h_out, S_in, S_out, c_CFRP, c_metal, m_CFRP, m_metal, T_wall, ...
    k_CFRP, k_metal, t_CFRP, t_metal] = Parameters(flag, base)

    addpath('Config\REFPROP')
    
    % Tanks specifications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Tank 2 : 13L, 250gr @ 300 bar 
    V2 = 13*1e-3; % [m3]
    L2 = 0.8;
    D2 = sqrt(V2/L2/pi*4);
    A_int2 = D2*pi*L2;
    % Tank 1 : 10 times first tank
    if base
        
        T_base = 293;
        V1 = 10*V2;
        P1 = 30000; % [kPa], tank 1 initial internal pressure
        
    else 
        
        T_base = 298;
        V1 = 10*V2;
        P1 = 30000; % [kPa], tank 1 initial internal pressure
        
    end
    L1 = 5;
    D1 = sqrt(V1/L1/pi*4);
    % Pipe geometry
    A_tube = 0.00635^2/4*pi; % Tube of 2 [cm] of diameter
    Ta = T_base; % [K], external temp.
    
    % Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    P2 = 200; % [kPa], tank 2 initial internal pressure
    T1 = T_base; % [K], tank 1 initial temp.
    T2 = T_base + 1; % [K], tank 2 initial temp.
    u2 = refpropm('U','T',T2,'P',P2,'hydrogen'); % Initial int. energy in tank 2
    u1 = refpropm('U','T',T1,'P',P1,'hydrogen'); % Initial int. energy in tank 2
    rho1 = refpropm('D','T',T1,'P',P1,'hydrogen'); % [kg/m3], tank 1
    rho2 = refpropm('D','T',T2,'P',P2,'hydrogen'); % [kg/m3], tank 2
    h1 = refpropm('H','T',T1,'P',P1,'hydrogen'); % [J/kg]
    m1 = rho1*V1; % [kg], tank 1 initial mass
    m2 = rho2*V2; % [kg], tank 2 initial mass
    
    % Heat transfer parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(flag,'H-T')
        
        t_CFRP = 10e-3; % Tank CFRP thickness
        t_metal = 5e-3; % Aluminium liner thickness
        k_CFRP = 0.55; % [W/mK], thermal conductivity 
        k_metal = 180;
        alpha_CFRP = 0.45e-6; % [m^2/s], thermal diffusivity 
        alpha_metal = 74.4e-6; % 
        rho_CFRP = 1530; % [kg/m3], density
        rho_metal = 2700;
        c_CFRP = k_CFRP/(rho_CFRP*alpha_CFRP); % [J/kgK], spec. heat capacity
        c_metal = k_metal/(rho_metal*alpha_metal);
        
        T_wall = T_base; % Initial tank wall temperature
        h_in = 10; % [W/m^2*K], inner convection coefficient
        h_out = 5; % [W/m^2*K], outer convection coefficient
        V_metal = pi/4*((D2 + 2*t_metal)^2 - D2^2);
        V_CFRP = pi/4*((D2 + 2*t_metal + 2*t_CFRP)^2 - (D2 + 2*t_metal)^2);
        m_metal = rho_metal*V_metal;
        m_CFRP = rho_CFRP*V_CFRP;
        S_in = pi*D2*L2; 
        S_out = pi*(D2 + 2*t_metal + 2*t_CFRP)*L2;
        
        h_conv = 0;
        
    else
        
        h_conv = 0; % [W/m2]
        [h_in, h_out, S_in, S_out, c_CFRP, m_CFRP, T_wall] = deal(0);
        
    end

end