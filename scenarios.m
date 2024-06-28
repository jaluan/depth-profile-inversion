function [d18Op,Tdgla,Z1,dT2,dZ2,dT3,dZ3,E0] = scenarios(scn)

switch scn
    case 1
%   Scn1 / synv4_1
    d18Op = 4.713;  % d18O threshold (3.5-5)
    Tdgla = 11e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.0185;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.7;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.23;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.28;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.17;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 2
% Scn2 / synv4_2
    d18Op = 4.499;  % d18O threshold (3.5-5)
    Tdgla = 11e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.0185;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.7;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.23;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.28;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.17;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 3
% Scn3 / synv4_3
    d18Op = 4.22;  % d18O threshold (3.5-5)
    Tdgla = 11e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.0185;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.7;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.23;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.28;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.17;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 4
% Scn4/ synv4_4
    d18Op = 3.814;  % d18O threshold (3.5-5)
    Tdgla = 11e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.0185;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.7;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.23;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.28;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.17;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 5
% Scn5/ synv4_5
    d18Op = 4.713;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.15;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.1;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 6
% Scn6/ synv4_6
    d18Op = 4.499;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.15;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.1;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 7
% Scn7/ synv4_7, 
    d18Op = 4.22;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.15;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.1;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 8
% Scn8/ synv4_8, 
    d18Op = 3.814;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.15;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.1;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 9
% Scn9/ synv4_9, Lang tids eksponering i ~1 m dybde inden plucking, lidt is.
    d18Op = 4.713;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.12;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.2;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = -0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)
    
    case 10
% Scn10/ synv4_10, Lang tids eksponering i ~1 m dybde inden plucking, medium is.
    d18Op = 4.499;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.12;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.2;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = -0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 11
% Scn11/ synv4_11, Lang tids eksponering i ~1 m dybde inden plucking, meget is.
    d18Op = 4.22;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.12;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.2;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = -0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 12
% Scn12/ synv4_12, Lang tids eksponering i ~1 m dybde inden plucking, meget is.
    d18Op = 3.814;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.12;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = 0.2;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 1.4;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = -0.57;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = -0.8;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 13
% Scn13/ synv4_13, Hurtig erosion ca. 415-800 ka, derefter langsom erosion, meget is.
    d18Op = 4.713;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.26;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.8;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.14;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.5;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.99;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 14
% Scn14/ synv4_14, Hurtig erosion ca. 415-800 ka, derefter langsom erosion, medium is.
    d18Op = 4.499;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.26;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.8;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.14;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.5;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.99;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 15
% Scn15/ synv4_15, Hurtig erosion ca. 415-800 ka, derefter langsom erosion, lidt is.
    d18Op = 4.22;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.26;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.8;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.14;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.5;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.99;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)

    case 16
% Scn16/ synv4_16, Hurtig erosion ca. 415-800 ka, derefter langsom erosion, lidt is.
    d18Op = 3.814;  % d18O threshold (3.5-5)
    Tdgla = 15e-3; % Time of last deglaciation [Myr] (2-15e-3)
    Z1    = 0.02;    % Depth at time of deglaciation [m] (0-0.05)
    dT2   = 0.26;  % Time change from Tdgla [Myr] (0-2 Myr)
    dZ2   = -0.8;   % Depth change during T2 [log10(m)] (-1 to 1)
    dT3   = 0.14;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
    dZ3   = 0.5;  % Depth change during T3 [log10(m)] (-1 to 1)
    E0    = 0.99;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)
end