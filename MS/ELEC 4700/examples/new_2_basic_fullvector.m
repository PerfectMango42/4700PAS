% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.305;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;         % number of modes to compute


% Define the start and end points, and the number of steps
endValue = 3.44;
numSteps = 10;

% Generate the values using linspace
Changingn2 = linspace(n2, endValue, numSteps);

% Place to keep the Neff values
storedNeff = zeros(1, numSteps);

% Loop through the values
for i = 1:numSteps
    % Replace any rw with values(i)
    
    [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,Changingn2(i),n3],[h1,h2,h3], ...
                                                rh,rw,side,dx,dy); 
    
    % Creating the Meshgrid for surface plots
    [Y, X] = meshgrid(y, x);
    
    % First consider the fundamental TE mode:
    for mode = 1:nmodes
        [Hx,Hy,neff] = wgmodes(lambda,Changingn2(i),nmodes,dx,dy,eps,'000A');

        storedNeff(i) = neff;
        
        fprintf(1,'neff = %.6f\n',neff);
        
        figure();
        subplot(2, 2, 1);
        contourmode(x,y,real(Hx(:,:,mode)));
        title(['Hx (TE mode ', num2str(mode), ')']); xlabel('x'); ylabel('y'); 
        for v = edges, line(v{:}); end
    
        % Add a 3D surface plot of Hx
        subplot(2, 2, 2);
        surf(Y, X, real(Hx(:,:,mode)));  % Surface plot of Hx
        shading interp;  % Smooth out the surface for better visual appearance
        colormap jet;    % Optional: change the color map
        xlabel('y');
        ylabel('x');
        
        
        subplot(2, 2, 3);
        contourmode(x,y,real(Hy(:,:,mode)));
        title(['Hy (TE mode ', num2str(mode), ')']); xlabel('x'); ylabel('y'); 
        for v = edges, line(v{:}); end
    
        subplot(2, 2, 4);
        surf(Y, X, real(Hy(:,:,mode)));  % Surface plot of Hx
        shading interp;  % Smooth out the surface for better visual appearance
        colormap jet;    % Optional: change the color map
        xlabel('y');
        ylabel('x');
    end
end

% plot the saved Neff vs the values() <- the varying ridge widths
figure(2);
plot(Changingn2, storedNeff, 'o-', 'LineWidth', 2);  % Plot stored values vs. the linspace values
xlabel('Ridge Width');
ylabel('Neff');
title('Plot of Ridge Width vs. Neff');
grid on;