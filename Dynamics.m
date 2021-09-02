function f = Dynamics(param,flag)

    addpath('Config\casadi-windows-matlabR2016a-v3.5.5')
    import casadi.*

    % Parameters

    c_CFRP = param(1);
    c_metal = param(2);
    k_CFRP = param(3);
    k_metal = param(4);
    h_in = param(5);
    h_out = param(6);
    m_CFRP = param(7);
    m_metal = param(8);
    t_CFRP = param(9);
    t_metal = param(10);
    D = param(11);
    L = param(12);
    T_inf = param(13);
    UA_IN = param(14);
    UA_OUT = param(15);
    UA_int = param(16);
    R1 = D/2;
    R2 = R1 + t_metal;
    R3 = R2 + t_CFRP;
    
    if flag

        UA_IN = ( 1/(h_in*2*pi*R1*L) + log(R2/R1)/(2*pi*L*k_metal) )^(-1);
        UA_int = ( log(R2/R1)/(2*pi*L*k_metal) + log(R3/R2)/(2*pi*L*k_CFRP) )^(1);
        UA_OUT = ( log(R3/R2)/(2*pi*L*k_CFRP) + 1/(h_out*2*pi*R3*L) )^(1);
        
    end

    % Function variables

    x = MX.sym('x',6,1);
    u = MX.sym('u',1);
    T2 = MX.sym('T2',1);
    h1 = MX.sym('h1',1);
    v1 = MX.sym('v1',1);

    x1 = x(1); % Tank 2 mass
    x2 = x(2); % Tank 2 specific internal energy
    x3 = x(3); % Tank 1 mass
    x4 = x(4); % Tank 1 specific internal energy
    x5 = x(5); % T_alu
    x6 = x(6); % T_CFRP

    % Function outputs

    dx1dt = u;
    dx2dt = (-UA_IN*(T2 - x5) + u*((v1^2)/2 + h1 - x2))/x1;
    dx3dt = -u;
    dx4dt = (-u*(v1^2/2 + h1 - x4))/x3;

    Q_IN = UA_IN*(T2 - x5);
    Q_OUT = UA_OUT*(x6 - T_inf);
    Q_int = UA_int*(x5 - x6);
    dx5dt = ( Q_IN - Q_int )/(m_metal*c_metal);
    dx6dt = ( Q_int - Q_OUT )/(m_CFRP*c_CFRP);
    
    % Dynamics function definition

    f = Function('f',{x,u,T2,h1,v1},...
        {[dx1dt;dx2dt;dx3dt;dx4dt;dx5dt;dx6dt;Q_IN;Q_OUT;Q_int]});
    
end
