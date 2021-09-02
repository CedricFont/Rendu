function [lb,ub] = OCP_Bounds(N,flag,V2,flag_control) % Fully numerical, no need for a gradient

    T_min = 200;
    T_max = 380;
    Pmin = 100;
    Pmax = 60000;
    % Min and max allowable internal energy in tank 2
    u_min = refpropm('U','T',T_min,'P',Pmax,'hydrogen'); % [J/kg]
    u_max = refpropm('U','T',T_max,'P',Pmin,'hydrogen');  
    m_min1 = 0; % [kg]
    m_min2 = 0;
    m_max1 = 20;
    m_max2 = 1;
    Kp_min = 0;
    Kp_max = 1;
    Ki_min = 0;
    Ki_max = 1;
    Kd_min = 0;
    Kd_max = 1;
    K_min = [Kp_min;Ki_min;Kd_min];
    K_max = [Kp_max;Ki_max;Kd_max];
    
    if strcmp(flag,'H-T') && strcmp(flag_control,'PID')
        T_min = 200;
        T_max = 370;
        lb = [m_min2*ones(N,1);u_min*ones(N,1);m_min1*ones(N,1);u_min*ones(N,1); ...
            T_min*ones(2*N,1);1e-2;K_min]; % Optimizer lower bound
        ub = [m_max2*ones(N,1);u_max*ones(N,1);m_max1*ones(N,1);u_max*ones(N,1); ...
            T_max*ones(2*N,1);1;K_max]; % Upper bound
    end
    
    if strcmp(flag,'H-T') && strcmp(flag_control,'PI')
        T_min = 200;
        T_max = 370;
        K_min = 0;
        K_max = 100;
        lb = [m_min2*ones(N,1);u_min*ones(N,1);m_min1*ones(N,1);u_min*ones(N,1); ...
            T_min*ones(2*N,1);1e-2;[K_min;K_min]]; % Optimizer lower bound
        ub = [m_max2*ones(N,1);u_max*ones(N,1);m_max1*ones(N,1);u_max*ones(N,1); ...
            T_max*ones(2*N,1);0.5;[K_max;K_max]]; % Upper bound
    end
    
    if strcmp(flag,'2tanks') && strcmp(flag_control,'PID-2')
        
        lb = [m_min*ones(N,1);u_min*ones(N,1);m_min*ones(N,1);u_min*ones(N,1); ...
            m_min*ones(N,1);u_min*ones(N,1);1;[K_min;Ki_min]]; % Optimizer lower bound
        ub = [m_max*ones(N,1);u_max*ones(N,1);m_max*ones(N,1);u_max*ones(N,1); ...
            m_max*ones(N,1);u_max*ones(N,1);20;[K_max;Ki_max]]; % Upper bound
        
    end
    
    if strcmp(flag,'H-T') && strcmp(flag_control,'I-PID')
        T_min = 220;
        T_max = 370;
        lb = [m_min2*ones(N,1);u_min*ones(N,1);m_min1*ones(N,1);u_min*ones(N,1); ...
            T_min*ones(2*N,1);1e-2;[K_min;1e-4]]; % Optimizer lower bound
        ub = [m_max2*ones(N,1);u_max*ones(N,1);m_max1*ones(N,1);u_max*ones(N,1); ...
            T_max*ones(2*N,1);1;[K_max;1e-1]]; % Upper bound
    end
    
    if strcmp(flag_control,'Free')
        
        % Min and max mass flow rate (control policy bounds)
        q_min = 0;
        q_max = 10; % Set to random for now
        
        lb = [m_min2*ones(N,1);u_min*ones(N,1);m_min1*ones(N,1);u_min*ones(N,1); ...
            T_min*ones(2*N,1);q_min*ones(N,1);1e-2]; % Optimizer lower bound
        ub = [m_max2*ones(N,1);u_max*ones(N,1);m_max1*ones(N,1);u_max*ones(N,1); ...
            T_max*ones(2*N,1);q_max*ones(N,1);1]; % Optimizer upper bound
        
    end
    
end