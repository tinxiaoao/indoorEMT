function px2m = calculate_px2m(real_value, pixel_value, type)
% calculate_px2m 根据实际尺寸和像素尺寸计算像素到米的转换比例
%
% 输入：
% real_value：实际尺寸（面积或边长）
% pixel_value：对应的像素尺寸（面积或边长）
% type：类型 ('area' 或 'length')
%
% 输出：
% px2m：像素到米的转换系数

switch lower(type)
    case 'area'
        px2m = sqrt(real_value / pixel_value);
    case 'length'
        px2m = real_value / pixel_value;
    otherwise
        error('Unknown type. Use ''area'' or ''length''.');
end
end
