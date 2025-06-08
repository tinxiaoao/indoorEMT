function [SUM_Solid_Angle_mean_R, SUM_Solid_Angle_mean_T, SUM_Solid_Angle_mean_R_T] = ...
    multi_layer_model(frequency, layer_thicknesses, material_id_vector, epsilon_params)
% MULTI_LAYER_MODEL - Calculate solid angle mean reflectance and transmittance for multi-layer materials
%
% Inputs:
%   frequency           - Frequency vector (Hz)
%   layer_thicknesses   - 1x5 real vector for [air, layer1, layer2, layer3, air] thicknesses (m)
%   material_id_vector  - 1x3 integer vector for material IDs corresponding to epsilon_params
%   epsilon_params      - Nx5 matrix with rows [e_inf, e_s, sigma_s, tau, alpha]
%
% Outputs:
%   SUM_Solid_Angle_mean_R   - Mean solid angle reflection loss
%   SUM_Solid_Angle_mean_T   - Mean solid angle transmission
%   SUM_Solid_Angle_mean_R_T - Combined reflection and transmission losses

% Constants
miu0 = 4*pi*1e-7;         % Permeability of free space
epsilon0 = 1/(36*pi)*1e-9; % Permittivity of free space

f = frequency(:);
omega = 2*pi*f;
num_freq = length(f);

% Calculate epsilon (complex permittivity) for each material layer
epsilon_layers = zeros(num_freq, 3);
for idx = 1:3
    params = epsilon_params(material_id_vector(idx), :);
    e_inf = params(1);
    e_s = params(2);
    sigma_s = params(3);
    tau = params(4);
    alpha = params(5);

    epsilon_layers(:, idx) = e_inf + (e_s - e_inf) ./ (1 + (1i*omega*tau).^(1-alpha)) ...
        + sigma_s ./ (1i*omega*epsilon0);
end

% Include air layers
epsilon_r = [ones(num_freq,1), epsilon_layers, ones(num_freq,1)];
layer_thicknesses = layer_thicknesses(:)'; % Ensure row vector

% Initialize parameters
num_layers = length(layer_thicknesses);
gama = 1i * omega .* sqrt(epsilon0 * miu0 * epsilon_r);
phi = gama .* layer_thicknesses;

% Transmission matrix calculations
A = ones(num_freq, num_layers);
B = zeros(num_freq, num_layers);
C = ones(num_freq, num_layers);
D = zeros(num_freq, num_layers);

% Recurrence relations
for n = num_layers-1:-1:1
    Y = sqrt(epsilon_r(:, n+1) ./ epsilon_r(:, n));
    Z = sqrt(epsilon_r(:, n) ./ epsilon_r(:, n+1));

    exp_phi = exp(phi(:, n));
    exp_minus_phi = exp(-phi(:, n));

    A(:, n) = (exp_phi .* (A(:, n+1) .* (1+Y) + B(:, n+1) .* (1-Y)))/2;
    B(:, n) = (exp_minus_phi .* (A(:, n+1) .* (1-Y) + B(:, n+1) .* (1+Y)))/2;

    C(:, n) = (exp_phi .* (C(:, n+1) .* (1+Z) + D(:, n+1) .* (1-Z)))/2;
    D(:, n) = (exp_minus_phi .* (C(:, n+1) .* (1-Z) + D(:, n+1) .* (1+Z)))/2;
end

% Reflection and Transmission coefficients at the first interface
RV = B(:,1)./A(:,1);
RH = D(:,1)./C(:,1);
TV = 1./A(:,1);
TH = 1./C(:,1);

% Numerical integration over solid angle
N_integral = 1000;
n_theta = (pi/2)/N_integral;
theta_i = linspace(0, pi/2, N_integral+1);
weight = cos(theta_i).*sin(theta_i);

reflect_power = 0.5*(abs(RV).^2 + abs(RH).^2);
transmit_power = 0.5*(abs(TV).^2 + abs(TH).^2);

integral_factor = sum(weight)*n_theta;

SUM_Solid_Angle_mean_R = (1 - reflect_power) * integral_factor;
SUM_Solid_Angle_mean_T = transmit_power * integral_factor;
SUM_Solid_Angle_mean_R_T = (reflect_power + transmit_power) * integral_factor;
end
