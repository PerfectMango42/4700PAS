set(0, 'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2);
set(0, 'Defaultaxeslinewidth', 2)

set(0, 'DefaultFigureWindowStyle', 'docked')

c_c = 299792458;    % m/s TWM speed of light
c_eps_0 = 8.8542149e-12;    % F/m vacuum permittivity
c_eps_0_cm = c_eps_0/100;   % F/cm
c_mu_0 = 1/c_eps_0/c_c^2;   % s^2/Fm
c_q = 1.60217653e-19;       % charge of electron
c_hb = 1.05457266913e-34;   % Dirac constant
c_h = c_hb*2*pi;            % 2 pi Dirac constant

% creating structure InputParasL and assigning values in the structure
% But InputParasR is just a regular scalar value
InputParasL.E0 = 1e5;   % Amplitude of electric field
InputParasL.we = 0;     % Frequency offset
InputParasL.t0 = 2e-12; % Time offset of Gaussian wave
InputParasL.wg = 5e-13; % Standard deviation of the wave
InputParasL.phi = 0;    % Starting phase of the wave
InputParasR = 0;        % No wave starting from the right

n_g = 3.5;              % index of refraction
vg = c_c/n_g*1e2;       % TWM cm/s group velocity
Lambda = 1550e-9;       % wavelength in nm 

% Might be the number of plots
plotN = 10; 

L = 1000e-6*1e2;        % cm
XL = [0,L];             % X axis range in a matrix
YL = [0,InputParasL.E0];% Y axis range in a matrix

Nz = 500;               % total grid steps in the graph
dz = L/(Nz-1);          % spacial step size along the length (L)
dt = dz/vg;             % time step with the corresponding spacial step size
fsync = dt*vg/dz;       % Always equals 1, syncronizing normalized factor

Nt = floor(2*Nz);       % number of time steps (discrete number time points in simulation)
tmax = Nt*dt;           % total simulation time
t_L = dt*Nz;            % time to travel length

% nan(not a number) to be filled in later
z = linspace(0,L,Nz).'; % Nz points, Nz-1 segments (column vector of length Nz, each entry from 0 - L)
time = nan(1, Nt);      % Row vector of size 1 - Nt
InputL = nan(1,Nt);     % Row vector of size 1 - Nt
InputR = nan(1,Nt);     % Row vector of size 1 - Nt
OutputL = nan(1,Nt);    % Row vector of size 1 - Nt
OutputR = nan(1,Nt);    % Row vector of size 1 - Nt

% Initializing Electric fields
Ef = zeros(size(z));    % z has Nz elements 
Er = zeros(size(z));    % z has Nz elements

% SourceFct creates a function handle (allow you to pass functions as
% arguments to other functions, store them...)
Ef1 = @SourceFct;
ErN = @SourceFct;

t = 0;                  % Starting time for the simulation
time(1) = t;            % already initialized as nan matrix, assigns first value to t

% nan values that get assigned values from the SourceFct with the InputParas structures, and the current time
InputL(1) = Ef1(t, InputParasL); 
InputR(1) = ErN(t, InputParasR);

% nan values that get assigned the boundaries of forward and reverse electric fields 
OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);

% Assigning the electric field input values at the boundaries opposite to
% the outputs
Ef(1) = InputL(1);
Er(Nz) = InputR(1);

% Create a figure that contains three sublots, each display different data
% about the electric fields inputs and outputs
% Forward propagating electric field at the left side
figure('name', 'Fields')
subplot(3,1,1)
plot(z*10000,real(Ef),'r');
hold off
xlabel('z(\mum)')
ylabel('E_f')
% Reverse propagating electric field at the right side
subplot(3,1,2)
plot(z*10000,real(Er),'b');
xlabel('z(\mum)')
ylabel('E_r')
hold off
% Inputs and Outputs over time in picoseconds
subplot(3,1,3)
plot(time*1e12, real(InputL), 'r'); hold on
plot(time*1e12, real(OutputR), 'r--'); 
plot(time*1e12, real(InputR), 'b'); hold on
plot(time*1e12, real(OutputL), 'b--');
xlabel('time(ps)')
ylabel('E')
hold off

for i = 2:Nt        % Iterate from 2 to the number of time steps
    t = dt*(i-1);   % Determine next time according to spacial step size and current iteration
    time(i) = t;    % Increment time
    
    % nan values that get assigned values from the SourceFct with the InputParas structures, and the current time
    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t,0);

    Ef(1) = InputL(i);
    Er(Nz) = InputR(i);

    % Updates the current Ef and Er over the spatial grid, ensuring to
    % normalize 
    Ef(2:Nz) = fsync*Ef(1:Nz-1);
    Er(1:Nz-1) = fsync*Er(2:Nz);
    % nan values that get assigned the boundaries of forward and reverse electric fields
    OutputR(i) = Ef(Nz);
    OutputL(i) = Er(1);

    % Create the plots that visualize the forward and reverse propagating
    % electric fields
    if mod(i,plotN) == 0            % updates every plotN iterations
        % Real and imaginary parts of forward propagating wave
        subplot(3,1,1)
        plot(z*10000,real(Ef),'r'); hold on
        plot(z*10000,imag(Ef),'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re','\Im')
        hold off
        % Real and imaginary parts of the reverse propagating wave
        subplot(3,1,2)
        plot(z*10000,real(Er),'b'); hold on
        plot(z*10000,imag(Er),'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re','\Im')
        hold off
        % Input and Output signals over time
        subplot(3,1,3);
        plot(time*1e12, real(InputL), 'r'); hold on
        plot(time*1e12, real(OutputR), 'g'); 
        plot(time*1e12, real(InputR), 'b');
        plot(time*1e12, real(OutputL), 'm');
        xlim([0,Nt*dt*1e12])            % Sets total simulation time (number of time steps * length of each step)
        ylim(YL)                        % Ensures plot is big enough for electric field amplitude
        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', 'Location', 'east')
        hold off
        pause(0.01)                     % Short delay in iterations of for loop (sleep())
    end
end














