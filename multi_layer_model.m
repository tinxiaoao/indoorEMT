function result = multi_layer_model(frequency, wall_type, wall_thickness_m, epsilon_params)
% MULTI_LAYER_MODEL - Calculates mean solid angle reflection/transmission properties for layered walls
%
% Inputs:
%   frequency        - Frequency vector (Hz)
%   wall_type        - Integer indicating wall type (1=external wall, 2=internal wall, 3=ceiling/floor, 4=window)
%   wall_thickness_m - Main layer thickness (m)
%   epsilon_params   - Material parameters matrix (N x 5), Cole-Cole parameters [e_inf, e_s, sigma_s, tau, alpha]
%
% Outputs:
%   result - Struct containing:
%       SUM_Solid_Angle_mean_R
%       SUM_Solid_Angle_mean_T
%       SUM_Solid_Angle_mean_R_T
%       RV, RH, TV, TH

mu0 = 4*pi*1e-7;         % Vacuum permeability
epsilon0 = 1/(36*pi)*1e-9; % Vacuum permittivity

f = frequency(:);
omega = 2*pi*f;
num_freq = length(f);

switch wall_type
    case 1  % External wall
        material_id_vector = [15, 19, 15];
        thicknesses = [0.005, wall_thickness_m, 0.005];
    case 2  % Internal wall
        material_id_vector = [15, 16, 15];
        thicknesses = [0.005, wall_thickness_m, 0.005];
    case 3  % Ceiling/floor
        material_id_vector = [15, 19, 15];
        thicknesses = [0.005, 0.1, 0.005];
    case 4  % Window (glass-air-glass)
        material_id_vector = [13, 21, 13];
        air_gap = 0.012;
        glass_thickness = (wall_thickness_m - air_gap)/2;
        thicknesses = [glass_thickness, air_gap, glass_thickness];
    otherwise
        error('Invalid wall_type specified.');
end

N_integral = 1000;
n_theta = (pi/2)/N_integral;
theta_i = linspace(0, pi/2, N_integral);

SUM_R = zeros(num_freq,1);
SUM_T = zeros(num_freq,1);
RV = zeros(num_freq, N_integral);
RH = zeros(num_freq, N_integral);
TV = zeros(num_freq, N_integral);
TH = zeros(num_freq, N_integral);

for theta_idx = 1:N_integral
    theta = theta_i(theta_idx);
    epsilon_layers = zeros(num_freq, 3);
    for idx = 1:3
        params = epsilon_params(material_id_vector(idx), :);
        epsilon_layers(:, idx) = params(1) + (params(2)-params(1))./(1+(1i*omega*params(4)).^(1-params(5))) + params(3)./(1i*omega*epsilon0);
    end

    epsilon_r = [ones(num_freq,1), epsilon_layers, ones(num_freq,1)];
    layer_thicknesses = [0, thicknesses, 0];

    gamma = 1i * omega .* sqrt(epsilon0 * mu0 .* epsilon_r - (sin(theta))^2);
    phi = gamma .* layer_thicknesses;

    A = ones(num_freq, 5);
    B = zeros(num_freq, 5);
    C = ones(num_freq, 5);
    D = zeros(num_freq, 5);

    for n = 4:-1:1
        Y = sqrt(epsilon_r(:, n+1) ./ epsilon_r(:, n));
        Z = sqrt(epsilon_r(:, n) ./ epsilon_r(:, n+1));

        exp_phi = exp(phi(:, n));
        exp_minus_phi = exp(-phi(:, n));

        A(:, n) = (exp_phi .* (A(:, n+1) .* (1+Y) + B(:, n+1) .* (1-Y))) / 2;
        B(:, n) = (exp_minus_phi .* (A(:, n+1) .* (1-Y) + B(:, n+1) .* (1+Y))) / 2;

        C(:, n) = (exp_phi .* (C(:, n+1) .* (1+Z) + D(:, n+1) .* (1-Z))) / 2;
        D(:, n) = (exp_minus_phi .* (C(:, n+1) .* (1-Z) + D(:, n+1) .* (1+Z))) / 2;
    end

    RV(:, theta_idx) = B(:,1)./A(:,1);
    RH(:, theta_idx) = D(:,1)./C(:,1);
    TV(:, theta_idx) = 1./A(:,1);
    TH(:, theta_idx) = 1./C(:,1);
end

for nn = 1:num_freq
    reflect_integrand = (1 - 1/2*(abs(RV(nn,:)).^2 + abs(RH(nn,:)).^2)).*cos(theta_i).*sin(theta_i);
    transmit_integrand = (1/2*(abs(TV(nn,:)).^2 + abs(TH(nn,:)).^2)).*cos(theta_i).*sin(theta_i);

    SUM_R(nn) = sum(reflect_integrand)*n_theta;
    SUM_T(nn) = sum(transmit_integrand)*n_theta;
end

SUM_R_T = SUM_R - SUM_T;

result = struct('SUM_Solid_Angle_mean_R', SUM_R, ...
                'SUM_Solid_Angle_mean_T', SUM_T, ...
                'SUM_Solid_Angle_mean_R_T', SUM_R_T, ...
                'RV', RV(:,1), 'RH', RH(:,1), 'TV', TV(:,1), 'TH', TH(:,1));
end
