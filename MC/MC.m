% initialize variables
x = 0;
v = 0;
scattering_prob = 0.05;

force = 1;
dt = 1;

num_steps = 100;
num_particles = 10;

position = zeros(num_particles, num_steps);
velocity = zeros(num_particles, num_steps);

for i = 1:num_particles
    v = 0;
    x = 0;


    for t = 1: num_steps

        v = v + force * dt;

        x = x + v * dt;

        if rand() < scattering_prob
            v = 0;
        end

        position(i, t) = x;
        velocity(i, t) = v;
    
        subplot(2, 1, 1)
        plot(1:t, velocity(i, 1:t), 'r'); hold on;
        plot(1:t, position(i, 1:t), 'g');
        xlabel('t');
        ylabel('v and x');
        %legend('velocity', 'position');
     
    
        subplot(2, 1, 2)
        plot(position(i, 1:t), velocity(i, 1:t), 'g'); hold on;
        xlabel('x');
        ylabel('v');
    
        pause(0.01);
    
    end

end

V_drift = mean(velocity(:));
title(["drift velocity = " num2str(V_drift)]);

