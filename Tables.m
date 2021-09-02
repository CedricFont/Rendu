function [D, hx, U, hy, T, P, H, K] = Tables(V2)
%% Compute table grid for T(D,U)
% Define lower and upper bounds both for D and U

% Retrieve maximum density
D_max1 = refpropm('D','T',298,'P',30000,'hydrogen');
m_max2 = 0.25; % [kg]
D_max2 = m_max2/V2;
D_max = max(D_max1,D_max2);
D_max = 2*D_max;
Nx = 100; % [kg/m3]
D = linspace(0, D_max, Nx);
hx = (D(end) - D(1))/Nx;

Tmin = 220;
Tmax = 380;
Pmin = 100; % 1 [bar]
Pmax = 35000; % 300 [bar]
Ny = 100;
% Min and max allowable internal energy in tank 2
u_min = refpropm('U','T',Tmin,'P',Pmax,'hydrogen'); % [J/kg]
u_max = refpropm('U','T',Tmax,'P',Pmin,'hydrogen'); 
U = linspace(u_min, u_max, Ny);
hy = (U(end) - U(1))/Ny;

%% Load files

T = load('T2x2');
P = load('P2x2');
H = load('H2x2');
K = load('K2x2');

T = T.T;
P = P.P;
H = H.H;
K = K.K;

end