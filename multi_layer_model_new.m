% 优化并检查后的 MATLAB 代码
clear; close all; clc;

% 设置频点和墙体参数
f = (1e9:1e8:6e9).'; % 设置频点为列向量
Length_frequency = length(f);

% 定义常数
miu0 = 4 * pi * 1e-7;
e0 = 1 / (36 * pi) * 1e-9;
omega = 2 * pi * f;

num_type_wall = 0;
Reflection_vertical = zeros(Length_frequency, 4);
Reflection_parallel = zeros(Length_frequency, 4);
Transmission_vertical = zeros(Length_frequency, 4);
Transmission_parallel = zeros(Length_frequency, 4);

for type_wall = [1, 2, 3, 4] 
    num_type_wall = num_type_wall + 1;

    if type_wall == 3
        ep = load('epsilon_paper_tabel.txt'); % 材料参数

        % 墙体材料索引
        k_values = [15, 19, 15];
        N_layers = numel(k_values);

        epsilon_r = zeros(Length_frequency, N_layers);

        % 计算epsilon_r
        for t = 1:N_layers
            k = k_values(t);
            e_INF = ep(k, 1);
            e_s = ep(k, 2);
            sigma_s = ep(k, 3);
            r_t = ep(k, 4);
            arf = ep(k, 5);

            e_c2 = e_INF + (e_s - e_INF) ./ (1 + (1i * omega * r_t).^(1 - arf)) + sigma_s ./ (1i * omega * e0);
            epsilon_r(:, t) = e_c2;
        end

        % 空气层
        epsilon_r = [ones(Length_frequency,1), epsilon_r, ones(Length_frequency,1)]; %#ok<AGROW>

        epsilon = e0 .* epsilon_r;
        gama = 1i * omega .* sqrt(epsilon .* miu0);

        % 墙体厚度定义
        d_layers = [0, 0.004, 0.012, 0.004, 0];
        d = repmat(d_layers, Length_frequency, 1);

        % 预分配
        theta = zeros(Length_frequency, N_layers + 2);
        phi = zeros(Length_frequency, N_layers + 2);
        A = zeros(Length_frequency, N_layers + 2);
        B = zeros(Length_frequency, N_layers + 2);
        C = zeros(Length_frequency, N_layers + 2);
        D = zeros(Length_frequency, N_layers + 2);

        % 边界初始化
        A(:, end) = 1;
        C(:, end) = 1;

        % 主循环
        for m = 1:Length_frequency
            theta(m, 1) = 0;
            for n = 2:N_layers + 2
                theta(m, n) = asin(gama(m, n - 1) .* sin(theta(m, n - 1)) ./ gama(m, n));
            end

            phi(m, :) = d(m, :) .* gama(m, :) .* cos(theta(m, :));

            for n = N_layers + 1:-1:1
                Y = (cos(theta(m, n + 1)) / cos(theta(m, n))) * sqrt(epsilon(m, n + 1) / epsilon(m, n));
                Z = (cos(theta(m, n + 1)) / cos(theta(m, n))) / sqrt(epsilon(m, n + 1) / epsilon(m, n));

                exp_phi = exp(phi(m, n));
                exp_minus_phi = exp(-phi(m, n));

                A_next = A(m, n + 1);
                B_next = B(m, n + 1);
                C_next = C(m, n + 1);
                D_next = D(m, n + 1);

                A(m, n) = (exp_phi / 2) * (A_next * (1 + Y) + B_next * (1 - Y));
                B(m, n) = (exp_minus_phi / 2) * (A_next * (1 - Y) + B_next * (1 + Y));
                C(m, n) = (exp_phi / 2) * (C_next * (1 + Z) + D_next * (1 - Z));
                D(m, n) = (exp_minus_phi / 2) * (C_next * (1 - Z) + D_next * (1 + Z));
            end
        end

        % 存储结果
        Reflection_vertical(:, num_type_wall) = B(:, 1) ./ A(:, 1);
        Reflection_parallel(:, num_type_wall) = D(:, 1) ./ C(:, 1);
        Transmission_vertical(:, num_type_wall) = 1 ./ A(:, 1);
        Transmission_parallel(:, num_type_wall) = 1 ./ C(:, 1);

        % 积分计算
        N_integral = 1000;
        theta_i = linspace(0, pi/2, N_integral + 1);

        for nn = 1:Length_frequency
            SUM_Solid_Angle_mean_R(nn, num_type_wall) = trapz(theta_i, (1 - 0.5*(abs(Reflection_vertical(nn, num_type_wall)).^2 + abs(Reflection_parallel(nn, num_type_wall)).^2)) .* cos(theta_i) .* sin(theta_i)); %#ok<*SAGROW>
            SUM_Solid_Angle_mean_T(nn, num_type_wall) = trapz(theta_i, 0.5*(abs(Transmission_vertical(nn, num_type_wall)).^2 + abs(Transmission_parallel(nn, num_type_wall)).^2) .* cos(theta_i) .* sin(theta_i));
            SUM_Solid_Angle_mean_R_T(nn, num_type_wall) = trapz(theta_i, (1 - 0.5*(abs(Reflection_vertical(nn, num_type_wall)).^2 + abs(Reflection_parallel(nn, num_type_wall)).^2 + abs(Transmission_vertical(nn, num_type_wall)).^2 + abs(Transmission_parallel(nn, num_type_wall)).^2)) .* cos(theta_i) .* sin(theta_i));
        end

    end
end
