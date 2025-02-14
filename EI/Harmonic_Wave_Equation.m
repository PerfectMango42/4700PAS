
nx = 50;
ny = 50;
V = zeros(nx,ny);
G = sparse(nx*ny, ny*nx);

% loop through all of G and set the appropriate boundary conditions

for i = 1:nx
    for j = 1:ny

        % mapping to n from i and j
        n = j + (i - 1) * ny;

        % set the boundary conditions into G from known values of n
        if (i == 1 || i == nx || j == 1 || j == ny)
            G(n, :) = 0;
            G(n, n) = 1;
            
        else
            % find the adjacent values to each n using nxm, nxp, nym, and nyp
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = (j - 1) + (i - 1) * ny;
            nyp = (j + 1) + (i - 1) * ny;
    
            if ( i > 10 && i < 20 && j > 10 && j < 20)
                G(n, n) = -2;
            else
                % set the diagonal to -4
                G(n, n) = -4;
            end

            %G(n, n) = -4;
    
            % in the nth row of G matrix set 1 to value indexed at nxm, nxp, nym, and nyp
            G(n, nxm) = 1;
            G(n, nxp) = 1;
            G(n, nym) = 1;
            G(n, nyp) = 1;
        end

    end
end

figure('name', 'Matrix')
spy(G)

nmodes = 20;
[E, D] = eigs(G, nmodes, 'SM');

figure ('name', 'Eigenvalues')
plot(diag(D), '*');

np = ceil(sqrt(nmodes));
figure('name', 'Modes')
for k = 1:nmodes
    M = E(:,k);
    for i = 1:nx
        for j = 1:ny
            n = i + (j-1)*nx;
            V(i,j) = M(n);
        end
        subplot(np,np,k), surf(V,'linestyle','none')
        title(['EV = ' num2str(D(k,k))])
    end
end


