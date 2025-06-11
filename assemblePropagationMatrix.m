function Mf = assemblePropagationMatrix(sigma_loss, sigma_diag, sigma_cpl)
% 将 sigma_loss、sigma_diag、sigma_cpl 合并为传播矩阵 Mf
%
% 输入:
%   sigma_loss [N x F]：房间自身损耗截面积
%   sigma_diag [N x F]：房间对外传输截面积总和
%   sigma_cpl  [N x N x F]：房间之间耦合截面积
%
% 输出:
%   Mf [N x N x F]：每个频点的传播矩阵

[N, F] = size(sigma_loss);
Mf = zeros(N, N, F);

for fIdx = 1:F
    for i = 1:N
        Mf(i,i,fIdx) = sigma_loss(i,fIdx) + sigma_diag(i,fIdx);
        for j = 1:N
            if i ~= j
                Mf(i,j,fIdx) = -sigma_cpl(j,i,fIdx);
            end
        end
    end
end

end
