%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%          INPUT FILE FOR SIMULATION OF STEAM REFORMER             %
%                FOR PRODUCTION OF SYNTHESIS GAS                   %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global GASCONST Ncomp RP RHOcat EPS Dp TEMPout LAMBDAst RADIUSi ...
    RADIUSo MMASS AJ AX ACTEN ENT298 ENT948 ADENT CP LAMBDA B S sumny rp...
    pin MM mpart uin Tin FRACin

%Initial data
%-------------------------------------------------------------------
Tin = 794;                % Initial temperature      [K]
pin = 29e5;               % Initial pressure         [Pa]
uin = 1.89;               % Velocity                 [m/s]


% Constants
%-------------------------------------------------------------------
GASCONST  = 8.3145E3;     % Gas constant             [J/kmoleK]
Ncomp     = 6;            % Number of components     [-]
ZP        = 30;           % Number of axial discretization points
RP        = mpart;            % Number of radial discretization points
mpart     = 50;            % Number of radial discretization points in the pellet

% Catalyst data
%-------------------------------------------------------------------
RHOcat = 2355.2;          % Density catalyst         [kgcat/m^3]
EPS    = 0.528;           % Pore fraction            [-]
Dp     = 0.0173;          % Particle diameter        [m]
rp     = Dp/2;            % Particle radius          [m]
hp     = 30000;           % Heat transfer coefficient[W/m^2K]
k_cat  = 50;              % Particle conductivity    [W/mK]
av     = 3/rp.*(1-EPS);   % Particle surface area per volume [1/m]

% Tube data
%-------------------------------------------------------------------
TEMPout  = 1100;          % Temp. outside the tube   [K]
LAMBDAst = 52;            % Heat coef. for tube metal[W/mK]
RADIUSi  = 0.051;         % Inner radius of the tube [m]
RADIUSo  = 0.066;         % Outer radius of the tube [m]
LENGTH   = 7.0;           % Tube length              [m]


% Initial massfraction of the components             [-]
%-------------------------------------------------------------------
FRACin(1) = 0.1911;       % CH4
FRACin(2) = 0.0001;       % CO
FRACin(3) = 0.0200;       % CO2
FRACin(4) = 0.0029;       % H2
FRACin(5) = 0.7218;       % H2O
FRACin(6) = 0.0641;       % N2


% Molemass of the components                         [kg/kmole]
%-------------------------------------------------------------------

MMASS(1) = 16.04;         % Molemass CH4
MMASS(2) = 28.01;         % Molemass CO
MMASS(3) = 44.01;         % Molemass CO2
MMASS(4) =  2.02;         % Molemass H2
MMASS(5) = 18.02;         % Molemass H2O
MMASS(6) = 28.01;         % Molemass N2

MM = 1/(sum(FRACin./MMASS));

% Preexponential factors for the rate constants      [kmole/kgcat h]
%-------------------------------------------------------------------
AJ(1) = 4.255E15;         % Factor for rx. 1
AJ(2) = 1.955E6 ;         % Factor for rx. 2
AJ(3) = 1.020E15;         % Factor for rx. 3


% Preexponential factors for the adsorbtion constants
%-------------------------------------------------------------------
AX(1) = 6.65E-4;          % Factor for CH4           [bar^-1]
AX(2) = 8.23E-5;          % Factor for CO            [bar^-1]
AX(3) = 6.12E-9;          % Factor for H2            [bar^-1]
AX(4) = 1.77E5;           % Factor for H2O           [-]


% Activation energies for the reactions              [J/kmole]
%-------------------------------------------------------------------
ACTEN(1) = 240.1E6;       % Activation energy for rx. 1
ACTEN(2) = 67.13E6;       % Activation energy for rx. 2
ACTEN(3) = 243.9E6;       % Activation energy for rx. 3


% Reaction enthalpies at 298K                        [J/kmole]
%-------------------------------------------------------------------
ENT298(1) = 206.1E6;      % Enthalpy for rx. 1 
ENT298(2) = -41.15E6;     % Enthalpy for rx. 2
ENT298(3) = 164.9E6;      % Enthalpy for rx. 3


% Reaction enthalpies at 948K                        [J/kmole]
%-------------------------------------------------------------------
ENT948(1) = 224.E6;       % Enthalpy for rx. 1
ENT948(2) = -37.30E6;     % Enthalpy for rx. 2
ENT948(3) = 187.5E6;      % Enthalpy for rx. 3


% Adsorption enthalpies                              [J/kmole]
%-------------------------------------------------------------------
ADENT(1) =-38.28E6;       % Enthalpy for adsorption of CH4
ADENT(2) =-70.65E6;       % Enthalpy for adsorption of CO
ADENT(3) =-82.90E6;       % Enthalpy for adsorption of H2
ADENT(4) = 88.68E6;       % Enthalpy for adsorption of H2O


% Heat capasity coefficients for the components        
%-------------------------------------------------------------------
CP(1,1) = 1.925E4;        % 1. coefficient for CH4   [J/kmoleK]
CP(1,2) = 5.213E1;        % 2. coefficient for CH4   [J/kmoleK^2]
CP(1,3) = 1.197E-2;       % 3. coefficient for CH4   [J/kmoleK^3]
CP(1,4) =-1.132E-5;       % 4. coefficient for CH4   [J/kmoleK^4]

CP(2,1) = 3.087E4;        % 1. coefficient for CO    [J/kmoleK]
CP(2,2) =-1.285E1;        % 2. coefficient for CO    [J/kmoleK^2]
CP(2,3) = 2.789E-2;       % 3. coefficient for CO    [J/kmoleK^3]
CP(2,4) =-1.272E-5;       % 4. coefficient for CO    [J/kmoleK^4]

CP(3,1) = 1.980E4;        % 1. coefficient for CO2   [J/kmoleK]
CP(3,2) = 7.344E1;        % 2. coefficient for CO2   [J/kmoleK^2]
CP(3,3) =-5.602E-2;       % 3. coefficient for CO2   [J/kmoleK^3]
CP(3,4) = 1.715E-5;       % 4. coefficient for CO2   [J/kmoleK^4]

CP(4,1) = 2.714E4;        % 1. coefficient for H2    [J/kmoleK]
CP(4,2) = 0.9274E1;       % 2. coefficient for H2    [J/kmoleK^2]
CP(4,3) =-1.381E-2;       % 3. coefficient for H2    [J/kmoleK^3]
CP(4,4) = 0.7645E-5;      % 4. coefficient for H2    [J/kmoleK^4]

CP(5,1) = 3.224E4;        % 1. coefficient for H2O   [J/kmoleK]
CP(5,2) = 0.1924E1;       % 2. coefficient for H2O   [J/kmoleK^2]
CP(5,3) = 1.055E-2;       % 3. coefficient for H2O   [J/kmoleK^3]
CP(5,4) = 0.3596E-5;      % 4. coefficient for H2O   [J/kmoleK^4]

CP(6,1) = 3.115E4;        % 1. coefficient for N2    [J/kmoleK]
CP(6,2) =-1.357E1;        % 2. coefficient for N2    [J/kmoleK^2]
CP(6,3) = 2.680E-2;       % 3. coefficient for N2    [J/kmoleK^3]
CP(6,4) =-1.168E-5;       % 4. coefficient for N2    [J/kmoleK^4]


% Conductivity coefficients for the components
%-------------------------------------------------------------------
LAMBDA(1,1) =-1.869E-3;   % 1. coefficient for CH4   [W/mK]
LAMBDA(1,2) = 8.727E-5;   % 2. coefficient for CH4   [W/mK^2]
LAMBDA(1,3) = 1.179E-7;   % 3. coefficient for CH4   [W/mK^3]
LAMBDA(1,4) =-3.614E-11;  % 4. coefficient for CH4   [W/mK^4]

LAMBDA(2,1) = 5.067E-4;   % 1. coefficient for CO    [W/mK]
LAMBDA(2,2) = 9.1025E-5;  % 2. coefficient for CO    [W/mK^2]
LAMBDA(2,3) =-3.524E-8;   % 3. coefficient for CO    [W/mK^3]
LAMBDA(2,4) = 8.199E-12;  % 4. coefficient for CO    [W/mK^4]

LAMBDA(3,1) =-7.215E-3;   % 1. coefficient for CO2   [W/mK]
LAMBDA(3,2) = 8.015E-5;   % 2. coefficient for CO2   [W/mK^2]
LAMBDA(3,3) = 5.477E-9;   % 3. coefficient for CO2   [W/mK^3]
LAMBDA(3,4) =-1.053E-11;  % 4. coefficient for CO2   [W/mK^4]

LAMBDA(4,1) = 8.099E-3;   % 1. coefficient for H2    [W/mK]
LAMBDA(4,2) = 6.689E-4;   % 2. coefficient for H2    [W/mK^2]
LAMBDA(4,3) =-4.158E-7;   % 3. coefficient for H2    [W/mK^3]
LAMBDA(4,4) = 1.562E-10;  % 4. coefficient for H2    [W/mK^4]

LAMBDA(5,1) = 7.341E-3;   % 1. coefficient for H2O   [W/mK]
LAMBDA(5,2) =-1.013E-5;   % 2. coefficient for H2O   [W/mK^2]
LAMBDA(5,3) = 1.801E-7;   % 3. coefficient for H2O   [W/mK^3]
LAMBDA(5,4) =-9.100E-11;  % 4. coefficient for H2O   [W/mK^4]

LAMBDA(6,1) = 3.919E-4;   % 1. coefficient for N2    [W/mK]
LAMBDA(6,2) = 9.966E-5;   % 2. coefficient for N2    [W/mK^2]
LAMBDA(6,3) =-5.067E-8;   % 3. coefficient for N2    [W/mK^3]
LAMBDA(6,4) = 1.504E-11;  % 4. coefficient for N2    [W/mK^4]


% Viscosity coefficients
%-------------------------------------------------------------------
B(1) = 1.00E-6;           % Coefficient for CH4      [kg/msK^0.5]
B(2) = 1.50E-6;           % Coefficient for CO       [kg/msK^0.5]
B(3) = 1.50E-6;           % Coefficient for CO2      [kg/msK^0.5]
B(4) = 0.65E-6;           % Coefficient for H2       [kg/msK^0.5]
B(5) = 1.74E-6;           % Coefficient for H2O      [kg/msK^0.5]
B(6) = 1.40E-6;           % Coefficient for N2       [kg/msK^0.5]

S(1) = 165;               % Coefficient for CH4      [K]
S(2) = 220;               % Coefficient for CO       [K]
S(3) = 220;               % Coefficient for CO2      [K]
S(4) =  67;               % Coefficient for H2       [K]
S(5) = 626;               % Coefficient for H2O      [K]
S(6) = 108;               % Coefficient for N2       [K]


% Diffusion volums
%-------------------------------------------------------------------
sumny(1) = 25.14;         % Coefficient for CH4      [-]
sumny(2) = 18.01;         % Coefficient for CO       [-]
sumny(3) = 26.90;         % Coefficient for CO2      [-]
sumny(4) =  6.12;         % Coefficient for H2       [-]
sumny(5) = 13.10;         % Coefficient for H2O      [-]
sumny(6) = 18.50;         % Coefficient for N2       [-]


