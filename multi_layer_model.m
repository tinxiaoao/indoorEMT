function [SUM_Solid_Angle_mean_R, SUM_Solid_Angle_mean_T, SUM_Solid_Angle_mean_R_T] = multi_layer_model(frequency, total_thickness, material_id_vector, epsilon_params)
% multi_layer_model  计算三层材料结构的反射/透射立体角功率比
%
% 输入参数：
%   frequency (Hz)            - 频率向量，例如 [1e9; 1.1e9; ...]（列向量）
%   total_thickness (m)       - 中间层材料厚度，例如 0.22 表示0.22米
%   material_id_vector        - 长度为3的材料编号向量，例如 [15 19 15] 表示三层材料的ID（对应 epsilon_params 的行索引）
%   epsilon_params (N×5 矩阵) - 包含N种材料 Cole-Cole模型参数的矩阵，每行格式 [e_inf, e_s, sigma_s, tau, alpha]
%
% 输出参数：
%   SUM_Solid_Angle_mean_R   - 仅反射能量损耗（随频率变化的列向量，长度为 length(frequency)）
%   SUM_Solid_Angle_mean_T   - 仅透射能量比例（同上）
%   SUM_Solid_Angle_mean_R_T - 综合反射+透射损耗（同上）

    % 确保频率为列向量
    f = frequency(:);
    Length_frequency = length(f);
    N_layers = length(material_id_vector);  % 三层结构（应为3）
    
    % 真空介电常数和磁导率常数
    miu0 = (4*pi)*1e-7;             % 真空磁导率 μ0
    e0   = 1/(36*pi) * 1e-9;        % 真空介电常数 ε0
    
    % 每种材料使用 Cole-Cole 模型计算复介电常数 (frequency-dependent)
    epsilon_r_internal = zeros(Length_frequency, N_layers);  % 存储每层材料的相对介电常数 (复数)
    w = 2*pi * f;  % 角频率向量
    
    for ii = 1:N_layers
        material_id = material_id_vector(ii);
        % 提取该材料的 Cole-Cole 参数
        e_inf   = epsilon_params(material_id, 1);
        e_s     = epsilon_params(material_id, 2);
        sigma_s = epsilon_params(material_id, 3);
        tau     = epsilon_params(material_id, 4);
        alpha   = epsilon_params(material_id, 5);
        % Cole-Cole 模型计算复介电常数 (相对介电常数)
        epsilon_complex = e_inf + (e_s - e_inf) ./ (1 + (1i * w * tau).^(1 - alpha)) ...
                                + sigma_s ./ (1i * w * e0);
        epsilon_r_internal(:, ii) = epsilon_complex;
    end
    
    % 构建包含空气层的五层结构介电常数矩阵 (相对介电常数):
    % 三层结构夹在两侧空气中 -> 层序：[空气, 材料1, 材料2, 材料3, 空气]
    atmosphere_col = ones(Length_frequency, 1);  % 空气的相对介电常数约为1
    epsilon_r = [atmosphere_col, epsilon_r_internal, atmosphere_col];  % 尺寸: (Length_frequency × (N_layers+2))
    
    % 定义各层厚度 (m): [空气厚度, 材料1厚度, 材料2厚度, 材料3厚度, 空气厚度]
    % 空气层视为半无限厚度，这里用0表示
    d = zeros(1, N_layers + 2);
    d(1) = 0; 
    d(end) = 0;
    % 两侧材料（材料1和材料3）厚度固定为5mm
    if N_layers >= 1, d(2) = 0.005; end        % 第一层材料厚度 5 mm
    if N_layers >= 3, d(end-1) = 0.005; end    % 第三层材料厚度 5 mm
    if N_layers >= 2, d(3) = total_thickness; end  % 中间层材料厚度 (由参数给定)
    
    % 计算各层的绝对介电常数 ε (ε0 * ε_r) 和传播常数 gama
    epsilon_abs = e0 * epsilon_r;   % (Length_frequency × (N_layers+2)) 每层介质的绝对介电常数
    % 计算传播常数 gama = i * ω * sqrt(ε * μ0)  (长度: Length_frequency × (N_layers+2))
    omega_matrix = repmat(2*pi*f, 1, N_layers+2);  % 将角频率向量扩展成矩阵
    gama = 1i * omega_matrix .* sqrt(epsilon_abs * miu0);
    
    % 计算各层的入射角（Snell定律）和相位厚度 phi
    theta = zeros(Length_frequency, N_layers+2);
    theta(:,1) = 0;  % 入射层（空气）入射角 = 0 (法向入射)
    % 由于法向入射，后续层的 theta 皆为 0，但这里保留一般形式
    for n = 2:(N_layers+2)
        % Snell定律计算折射角： sin(theta_n) = (gama_{n-1}/gama_n) * sin(theta_{n-1})
        theta(:, n) = asin((gama(:, n-1) ./ gama(:, n)) .* sin(theta(:, n-1)));
    end
    % 计算相位因子 phi(m,n) = k_n * d_n * cos(theta_n)，其中 k_n = gama_n (传播常数)
    cos_theta = cos(theta);
    % 将厚度行向量扩展为矩阵 (每个频点一行，厚度相同)
    d_matrix = repmat(d, Length_frequency, 1);
    % phi 为 (Length_frequency × (N_layers+2)) 矩阵
    phi = d_matrix .* gama .* cos_theta;
    
    % 使用传输矩阵法计算反射/透射系数（垂直和水平极化）
    % 初始化递推矩阵 A, B, C, D（每行为一个频率点，每列对应一层边界条件）
    A = zeros(Length_frequency, N_layers+2);
    B = zeros(Length_frequency, N_layers+2);
    C = zeros(Length_frequency, N_layers+2);
    D = zeros(Length_frequency, N_layers+2);
    % 边界条件：在最后一层（空气）界面处，传播波向前的幅度为1，反向为0（两种极化分别设置）
    A(:, end) = 1;
    C(:, end) = 1;
    B(:, end) = 0;
    D(:, end) = 0;
    % 从后向前递推计算各层的系数
    % Y 和 Z 为界面处的阻抗比相关系数，对应垂直极化和水平极化
    for n = (N_layers+1):-1:1  % n 从倒数第二层向前遍历到第1层
        % 计算当前界面 (n 与 n+1 层交界) 的系数 Y 和 Z
        % Y = (cosθ_{n+1}/cosθ_n) * sqrt(ε_{n+1}/ε_n)    （垂直极化）
        % Z = (cosθ_{n+1}/cosθ_n) * sqrt(ε_n/ε_{n+1})    （水平极化）
        Y = (cos_theta(:, n+1) ./ cos_theta(:, n)) .* sqrt(epsilon_abs(:, n+1) ./ epsilon_abs(:, n));
        Z = (cos_theta(:, n+1) ./ cos_theta(:, n)) .* sqrt(epsilon_abs(:, n)   ./ epsilon_abs(:, n+1));
        % 利用递推公式更新当前层的系数 A, B, C, D
        A(:, n) = exp(phi(:, n)) / 2 .* ( A(:, n+1) .* (1 + Y) + B(:, n+1) .* (1 - Y) );
        B(:, n) = exp(-phi(:, n)) / 2 .* ( A(:, n+1) .* (1 - Y) + B(:, n+1) .* (1 + Y) );
        C(:, n) = exp(phi(:, n)) / 2 .* ( C(:, n+1) .* (1 + Z) + D(:, n+1) .* (1 - Z) );
        D(:, n) = exp(-phi(:, n)) / 2 .* ( C(:, n+1) .* (1 - Z) + D(:, n+1) .* (1 + Z) );
    end
    
    % 计算入射端口（第1层界面）的反射和透射系数 (垂直/水平极化)
    RV = B(:,1) ./ A(:,1);   % 垂直极化反射系数（电场垂直于入射面）
    RH = D(:,1) ./ C(:,1);   % 水平极化反射系数（电场平行于入射面）
    TV = 1 ./ A(:,1);        % 垂直极化透射系数（透射到另一侧空气中的比例）
    TH = 1 ./ C(:,1);        % 水平极化透射系数
    
    % 计算立体角平均功率比
    % 反射方向：计算“仅反射能量损耗” (即未被反射的能量比例)
    % 透射方向：计算“仅透射能量比例”
    % 综合：计算“反射+透射总损耗” (即离开墙体的总能量比例)
    N_integral = 1000;                     % 积分等分数
    d_theta = (pi/2) / N_integral;         % 每份对应的角度增量
    theta_i = 0 : d_theta : (pi/2);        % 从0到90°的积分角度划分
    weight = cos(theta_i) .* sin(theta_i); % 积分权重因子 (cosθ * sinθ)
    % 计算各频点的平均值（对每个频率点，反射/透射系数视为对各入射角相同）
    reflect_power_fraction = 0.5 * ( abs(RV).^2 + abs(RH).^2 );  % 平均反射功率系数 (两种极化取平均)
    transmit_power_fraction = 0.5 * ( abs(TV).^2 + abs(TH).^2 ); % 平均透射功率系数
    % 执行数值积分（矩形法）: sum(weight) * d_theta = 0.5（半球立体角积分值）
    integral_factor = sum(weight) * d_theta;
    % 仅反射能量损耗：1 减去反射回来的比例 (即进入墙体的比例，包括透射+吸收)
    SUM_Solid_Angle_mean_R   = (1 - reflect_power_fraction) * integral_factor;
    % 仅透射能量比例：直接透射出去的比例
    SUM_Solid_Angle_mean_T   = transmit_power_fraction * integral_factor;
    % 综合反射+透射损耗：反射或透射离开墙体的总比例
    SUM_Solid_Angle_mean_R_T = (reflect_power_fraction + transmit_power_fraction) * integral_factor;
end
