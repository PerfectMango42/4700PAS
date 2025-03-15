% Plot Properties
set(0, 'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2)
set(0, 'Defaultaxeslinewidth', 2)
set(0, 'DefaultFigureWindowStyle', 'docked')

Is = 0.01e-12;
Ib = 0.1e-12;
Vb = 1.3;
Gp = 0.1;
n = 200;

V = linspace(-1.97, 0.7, n);

variation = 0.2 * rand(size(V));

I = Is*(exp(1.2 * V/0.025) - 1) + Gp*V - Ib*(exp(-(1.2/0.025)*(V + Vb)) - 1);

I_noise = I + variation; 

% 4th order polynomials
poly_4 = polyfit(V, I, 4);
I_val_4 = polyval(poly_4, V);
% 4th order polynomials noise
poly_4_noise = polyfit(V, I_noise, 4);
I_val_4_noise = polyval(poly_4_noise, V);
% 8th order polynomials
poly_8 = polyfit(V, I, 8);
I_val_8 = polyval(poly_8, V);
% 8th order polynomials noise
poly_8_noise = polyfit(V, I_noise, 8);
I_val_8_noise = polyval(poly_8_noise, V);

% non-linear fit AC
fo_AC = fittype('A*(exp(1.2 * x/0.025) - 1) + Gp*x - C*(exp(-(1.2/0.025)*(x + Vb)) - 1)');
ff_AC = fit(V(:), I(:), fo_AC);
If_AC = ff_AC(V);
% % non-linear fit ABC
fo_ABC = fittype('A*(exp(1.2 * x/0.025) - 1) + B*x - C*(exp(-(1.2/0.025)*(x + Vb)) - 1)');
ff_ABC = fit(V(:), I(:), fo_ABC);
If_ABC = ff_ABC(V);
% % non-linear fit ABCD
fo_ABCD = fittype('A*(exp(1.2 * x/0.025) - 1) + B*x - C*(exp(-(1.2/0.025)*(x + D)) - 1)');
ff_ABCD = fit(V(:), I(:), fo_ABCD);
If_ABCD = ff_ABCD(V);
% non-linear fit AC noise
fo_AC_noise = fittype('A*(exp(1.2 * V/0.025) - 1) + Gp*V - C*(exp(-(1.2/0.025)*(V + Vb)) - 1)', 'independent', 'V');
ff_AC_noise = fit(V(:), I_noise(:), fo_AC_noise);
If_AC_noise = ff_AC_noise(V);
% % non-linear fit ABC noise
fo_ABC_noise = fittype('A*(exp(1.2 * V/0.025) - 1) + B*V - C*(exp(-(1.2/0.025)*(V + Vb)) - 1)', 'independent', 'V');
ff_ABC_noise = fit(V(:), I_noise(:), fo_ABC_noise);
If_ABC_noise = ff_ABC_noise(V);
% % non-linear fit ABCD noise
fo_ABCD_noise = fittype('A*(exp(1.2 * V/0.025) - 1) + B*V - C*(exp(-(1.2/0.025)*(V + D)) - 1)', 'independent', 'V');
ff_ABCD_noise = fit(V(:), I_noise(:), fo_ABCD_noise);
If_ABCD_noise = ff_ABCD_noise(V);

inputs = V.';
targets = I.';
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
view(net)
Inn = outputs;

xp = linspace(-2, 1, 200).';
rgp_I = I.';
m = fitrgp(xp, rgp_I);
[ypred1, yint1] = predict(m, xp);


figure 
subplot(4, 2, 1)
plot (V, I, 'g'); hold on
plot (V, I_noise, 'g--'); 
plot (V, I_val_4, 'r'); 
plot (V, I_val_4_noise, 'r--'); 
plot (V, I_val_8, 'b'); 
plot (V, I_val_8_noise, 'b--'); 
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('I','Inoise', 'I4', 'I4noise', 'I8', 'I8noise')
hold off

subplot(4, 2, 2)
semilogy(V, abs(I), 'g'); hold on
semilogy(V, abs(I_noise), 'g--');
semilogy(V, abs(I_val_4), 'r'); 
semilogy(V, abs(I_val_4_noise), 'r--');
semilogy(V, abs(I_val_8), 'b');
semilogy(V, abs(I_val_8_noise), 'b--');
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('I','Inoise', 'I4', 'I4noise', 'I8', 'I8noise')
hold off

subplot(4, 2, 3)
plot(V, If_AC, 'r'); hold on
plot(V, If_ABC, 'b'); 
plot(V, If_ABCD, 'g'); 
plot(V, If_AC_noise, 'r--'); 
plot(V, If_ABC_noise, 'b--'); 
plot(V, If_ABCD_noise, 'g--'); 
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('AC ','ABC', 'ABCD', 'ACnoise ','ABCnoise', 'ABCDnoise')
hold off;

subplot(4, 2, 4)
semilogy(V, abs(If_AC), 'r'); hold on
semilogy(V, abs(If_ABC), 'b'); 
semilogy(V, abs(If_ABCD), 'g'); 
semilogy(V, abs(If_AC_noise), 'r--'); 
semilogy(V, abs(If_ABC_noise), 'b--'); 
semilogy(V, abs(If_ABCD_noise), 'g--'); 
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('AC ','ABC', 'ABCD', 'ACnoise ','ABCnoise', 'ABCDnoise')
hold off;

subplot(4, 2, 5)
plot(targets, 'b'); hold on
plot(Inn, 'r');
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('data', 'nn')
hold off

subplot(4, 2, 6)
semilogy(abs(targets), 'b'); hold on
semilogy(abs(Inn), 'r');
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('data', 'nn')
hold off

subplot(4, 2, 7)
plot(xp, ypred1, 'g'); 
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('RGP')
hold off;

subplot(4, 2, 8)
semilogy(xp, abs(ypred1), 'g'); 
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('RGP')
hold off;

