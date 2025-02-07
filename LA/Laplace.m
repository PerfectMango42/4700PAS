
nx = 100;
ny = 100;
V = zeros(nx, ny); % Inititalize matrix
num_iterations = 10000;

for n = 1:num_iterations

    for i = 1:nx % loop through each x value

        for j = 1:ny % for each x loop through the y values

            % Updating Boundary Conditions
            %V(1, j) = 1; % left
            %V(nx, j) = 0; % right
            %V(i, 1) = V(i, 2); % bottom
            %V(i, ny) = V(i, ny-1); % top

            % Other Boundary Conditions
            V(1, j) = 1; % left
            V(nx, j) = 1; % right
            V(i, 1) = 0; % bottom
            V(i, ny) = 0; % top

            % Ensuring its doesnt change boundary values
            if (((i > 1) && (j > 1)) && ((i < nx) && (j < ny)))
                % Calculate average the potential
                V(i,j) = (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1))/4;
            end

        end
    end
    
    % filter potential
    V = imboxfilt(V,3);

    % Plot potential
    if (mod(n, 50)) == 0
        surf(V')
        %imagesc(V')
        pause(0.05);
    end
end

% Plot electric fields
[Ex, Ey] = gradient(V);

figure

subplot(1,3,1);
surf(Ex)
title('Ex')

subplot(1,3,2)
surf(Ey)
title('Ey')

subplot(1,3,3);
quiver(-Ey', -Ex', 10)
title('quiver of Ey and Ex')