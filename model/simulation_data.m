%
% Explain model
%
%

clear;
clc;

%% Open model
fn = 'symmetric_circular_coil.mph';

model = mphopen(fn);
disp("Model loaded succesfully.")


%% Plot geometry and mesh
t = tiledlayout(1,2,'TileSpacing','Tight');
title(t,'Geometry (left) and Mesh (right)')

ax1 = nexttile;
mphgeom(model);
set(ax1,'Visible','off');
ax2 = nexttile;
mphmesh(model);
set(ax2,'Visible','off');


%% Base variables
% The plate thickness in the model equals 5 mm. The baseline variables have
% been chosen such that the skin thickness equals approximately the plate
% thickness. 

t_plate = model.param.evaluate('t_plate'); % [m]
f0 = 10E3;      % Excitation frequency in [Hz]
I0 = 1.0;       % Excitation current in [A]
sigma0 = 1E6;   % Plate conductivity in [S/m]
mu_r0 = 1;      % Relative magnetic permeability [-]
mu0_const = 1.25663706212E-6;

% Function to calculate skin depth
skin_depth = @(f,sigma,mu_r) 1/(sqrt(f*mu0_const*mu_r*pi*sigma));
delta0 = skin_depth(f0, sigma0, mu_r0);

% Run model to obtain base values for induced power and magnetic flux
% density.
disp('Set parameters for baseline simulation.')
model.param.set('I', [num2str(I0), ' [A]']);
model.param.set('mu_r', num2str(mu_r0));
model.param.set('f', [num2str(f0), ' [Hz]']);
model.param.set('sigma', [num2str(sigma0), ' [S/m]']);
disp('Run baseline simulation.')
model.study('std1').run;
power0 = model.result.numerical('int1').getReal;
flux0 = model.result.numerical('pev1').getReal;

% Save data for later use.
% disp('Store data.')
% save('base_values.mat', 'f0', 'sigma0', 'mu_r0', 'I0', 'mu0_const', ...
%     'delta0', 'power0', 'flux0');
disp('Done.')

%% Plot MF and current distribution for small and large skin depth
% Two extreme cases are simulated: one with a large skin depth and one with
% a small skin depth.

% Reset parameters to baseline values
model.param.set('I', [num2str(I0), ' [A]']);
model.param.set('mu_r', num2str(mu_r0));

% Small skin depth
f = 100E3;
sigma = 10E6;
model.param.set('f', [num2str(f), ' [Hz]']);
model.param.set('sigma', [num2str(sigma), ' [S/m]']);
delta_small = skin_depth(f, sigma, mu_r0);

model.study('std1').run;
mf_small = mphplot(model, 'pg1');
J_small = mphplot(model, 'pg2');

% Large skin depth
f = 10E3;
sigma = 10E3;
model.param.set('f', [num2str(f), ' [Hz]']);
model.param.set('sigma', [num2str(sigma), ' [S/m]']);
delta_large = skin_depth(f, sigma, mu_r0);

model.study('std1').run;
mf_large = mphplot(model, 'pg1');
J_large = mphplot(model, 'pg2');
close all;

% Plot data
t = tiledlayout(2,2,'TileSpacing','Tight');
title(t,'Magnetic Field (top) and Current Distribution (bottom)')

ax11 = nexttile;
mphplot(mf_small)
ax11.Title.String = ['\delta/t = ', num2str(delta_small/t_plate, 2)];
set(ax11,'Visible','off');
set(findall(ax11, 'type', 'text'), 'visible', 'on');
ax12 = nexttile;
mphplot(mf_large)
ax12.Title.String = ['\delta/t = ', num2str(delta_large/t_plate, 2)];
set(ax12,'Visible','off');
set(findall(ax12, 'type', 'text'), 'visible', 'on');

ax21 = nexttile;
mphplot(J_small)
ax21.Title.String = ['\delta/t = ', num2str(delta_small/t_plate, 2)];
set(ax21,'Visible','off');
set(findall(ax21, 'type', 'text'), 'visible', 'on');
ax22 = nexttile;
mphplot(J_large)
ax22.Title.String = ['\delta/t = ', num2str(delta_large/t_plate, 2)];
set(ax22,'Visible','off');
set(findall(ax22, 'type', 'text'), 'visible', 'on');

%% Vary conductivity
% Simulations are performed for a large range of conductivities. 

% Reset parameters to baseline values
disp('Initialize.');
model.param.set('f', [num2str(f0), ' [Hz]']);
model.param.set('I', [num2str(I0), ' [A]']);
model.param.set('sigma', [num2str(sigma0), ' [S/m]']);
model.param.set('mu_r', num2str(mu_r0));

% Create range for conductivities and initalize arrays to store data
sigma = [1E3, 2E3, 5E3, 1E4, 2E4, 5E4, 1E5, 2E5, 5E5, 1E6, 2E6, 5E6, 1E7, 2E7, 5E7, 1E8, 2E8, 5E8, 1E9];
induced_power = zeros(size(sigma));
flux_in_coil = zeros(size(sigma));
delta = zeros(size(sigma));

% Run simulations
for i=1:length(sigma)
    disp(['Run simulation: ', num2str(i), '/', num2str(length(sigma)), '.']);
    model.param.set('sigma', [num2str(sigma(i)), ' [S/m]']);
    model.sol().run;
    induced_power(i) = model.result.numerical('int1').getReal;
    flux_in_coil(i) = model.result.numerical('pev1').getReal;
    delta(i) = skin_depth(f0, sigma(i), mu_r0);
end

% Plot data
figure;
loglog(sigma, induced_power);
xlabel('Conductivity [S/m]')
ylabel('Induced power [W]')

figure;
semilogx(sigma, flux_in_coil)
xlabel('Conductivity [S/m]')
ylabel('Magnetic flux in coil [H]')

% Store and save data
% disp('Store data.')
% save('conductivity.mat', 'sigma', 'induced_power', 'flux_in_coil', 'delta')
disp('Done.')

%% Vary frequency
% Simulations are performed for a large range of frequencies. 

% Reset parameters to baseline values
disp('Initialize.');
model.param.set('f', [num2str(f0), ' [Hz]']);
model.param.set('I', [num2str(I0), ' [A]']);
model.param.set('sigma', [num2str(sigma0), ' [S/m]']);
model.param.set('mu_r', num2str(mu_r0));

% Create range for conductivities and initalize arrays to store data
freq = [10, 20, 50, 100, 200, 500, 1E3, 2E3, 5E3, 1E4, 2E4, 5E4, 1E5, 2E5, 5E5, 1E6, 2E6, 5E6, 1E7];
induced_power = zeros(size(freq));
flux_in_coil = zeros(size(freq));
delta = zeros(size(freq));

% Run simulations
for i=1:length(freq)
    disp(['Run simulation: ', num2str(i), '/', num2str(length(freq)), '.']);
    model.param.set('f', [num2str(freq(i)), ' [Hz]']);
    model.sol().run;
    induced_power(i) = model.result.numerical('int1').getReal;
    flux_in_coil(i) = model.result.numerical('pev1').getReal;
    delta(i) = skin_depth(freq(i), sigma0, mu_r0);
end

% Plot data
figure;
loglog(freq, induced_power);
xlabel('frequency [Hz]')
ylabel('Induced power [W]')

figure;
semilogx(freq, flux_in_coil)
xlabel('Frequency [Hz]')
ylabel('Magnetic flux in coil [H]')

% Store and save data
% disp('Store data.')
% save('frequency.mat', 'freq', 'induced_power', 'flux_in_coil', 'delta')
disp('Done.')

%% Vary current
% Simulations are performed for a large range of frequencies. 

% Reset parameters to baseline values
disp('Initialize.');
model.param.set('f', [num2str(f0), ' [Hz]']);
model.param.set('I', [num2str(I0), ' [A]']);
model.param.set('sigma', [num2str(sigma0), ' [S/m]']);
model.param.set('mu_r', num2str(mu_r0));

% Create range for conductivities and initalize arrays to store data
I = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
induced_power = zeros(size(I));
flux_in_coil = zeros(size(I));
delta = zeros(size(I));

% Run simulations
for i=1:length(I)
    disp(['Run simulation: ', num2str(i), '/', num2str(length(I)), '.']);
    model.param.set('I', [num2str(I(i)), ' [A]']);
    model.sol().run;
    induced_power(i) = model.result.numerical('int1').getReal;
    flux_in_coil(i) = model.result.numerical('pev1').getReal;
    delta(i) = skin_depth(f0, sigma0, mu_r0);
end

% Plot data
figure;
loglog(I, induced_power);
xlabel('frequency [Hz]')
ylabel('Induced power [W]')

figure;
semilogx(I, flux_in_coil)
xlabel('Frequency [Hz]')
ylabel('Magnetic flux in coil [H]')

% Store and save data
% disp('Store data.')
% save('current.mat', 'I', 'induced_power', 'flux_in_coil', 'delta')
disp('Done.')