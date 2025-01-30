winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15      % total simulation time
%f = 230e12;         % frequency of the wave 230 THz\
f = 230e12;
lambda = c_c/f;     % wavelength

% Grid setup
xMax{1} = 20e-6;
nx{1} = 200;
ny{1} = 0.75*nx{1};


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0;   % magnetic permeability

% defines the inclusion region
epi{1} = ones(nx{1},ny{1})*c_eps_0;     % free space everywhere
%epi{1}(125:150,55:95)= c_eps_0*11.3;    % region of higher permittivity

% additional inclusion
%epi{1}(50:75,20:45)= c_eps_0*35;    % region of higher permittivity

% Adding a weird lightning bolt inclusion (for fun!)

% Create a lightning bolt shape as a mask
lightning_mask = zeros(nx{1}, ny{1});  % Start with a blank grid

% Define lightning bolt branches (you can adjust these as you like)
lightning_mask(70:110, 40:60) = 1;  % First line (horizontal)
lightning_mask(110:130, 50:70) = 1;  % Second line (diagonal)
lightning_mask(130:150, 60:80) = 1;  % Third line (diagonal)

% Add jaggedness for a more "lightning bolt" effect
lightning_mask(90:100, 60:65) = 1;  % Branching out a bit
lightning_mask(120:130, 80:85) = 1;  % Another branch

% Now apply the lightning bolt permittivity
epi{1}(lightning_mask == 1) = c_eps_0 * 50;  % Set permittivity in the lightning bolt region


% initialize for conductivity of materials
sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1};             % spacial step size
dt = 0.25*dx/c_c;               % time step size
nSteps = round(tSim/dt*2);      % number of steps for the simulation
yMax = ny{1}*dx;                % 
nsteps_lamda = lambda/dx        % Gid points fit in one wavelength

movie = 1;      % animation of simulation
% define plot values
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

% Boundary conditions 
bc{1}.NumS = 2;                     % Number of sources at the boundary
bc{1}.s(1).xpos = nx{1}/(4) + 1;    % Position of the source in the x-direction
bc{1}.s(1).type = 'ss';             % Type of source (sinusoidal)
bc{1}.s(1).fct = @PlaneWaveBC;      % Source behaviour at the boundary

% second source
bc{1}.s(2).ypos = ny{1}/4 + 1;      % Position of the source in the y-direction
bc{1}.s(2).type = 'ss';              % Type of source (sinusoidal)
bc{1}.s(2).fct = @PlaneWaveBC;      % Source behaviour at the boundary

% Variables to control the plane wave parameters
% mag = -1/c_eta_0;
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
%st = 15e-15;
st = -0.05;     % pulse width of source wave
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
% call to create plane wave
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

% Parameters for the second source (propagating in the y-direction)
mag_y = 1;                           % Amplitude of the wave (y-direction)
phi_y = 0;                           % Phase of the wave (y-direction)
omega_y = f * 2 * pi;                % Angular frequency (using the same f as before)
betap_y = 0;                         % Wave number (can be adjusted for direction)
t0_y = 30e-15;                       % Time offset
st_y = 15e-15;                       % Time window
s_y_y = 0;                           % Initial displacement in the y-direction
sty_y = 1.5 * lambda;                % Spatial extent (related to wavelength)
% call to create plane wave
bc{1}.s(2).paras = {mag_y, phi_y, omega_y, betap_y, t0_y, st_y, s_y_y, sty_y, 's'};

Plot.y0 = round(y0/dx);

% Defining the positive/minus boundaries as absorbing boundaries
bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';    % electric field boundary (doesnt absorb)
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

bc{1}.s(2).xm.type = 'a';           
bc{1}.s(2).xp.type = 'a';            
bc{1}.s(2).ym.type = 'a';           
bc{1}.s(2).yp.type = 'a';            

% Create perfectly matched layer
pml.width = 20 * spatialFactor;
pml.m = 3.5;

% Region of the simulation
Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






