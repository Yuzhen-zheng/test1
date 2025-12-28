clc; clear all;
%% 参数设定
% 电偶极子及介质参数
h = 5;            % 电偶极子所在位置
f = 3;            % 固定频率3Hz
w = 2*pi*f;
I = 5*w*1i;       % 电流大小
l = 2;            % 电偶极子长度
mu = 4*pi*1e-7;   % 磁导率
e1 = 1e-9/(36*pi);% 空气介电常数
r1 = 0;           % 空气电导率
e2 = 80*e1;       % 海水介电常数
r2_range = linspace(0.1, 20, 100); % 电导率范围0.1-20S/m(从0.1开始避免除零错误)

% 固定观测点
obs_point = [100, 25, 30]; % 观测点坐标(x,y,z)
x = obs_point(1);
y = obs_point(2);
z_val = obs_point(3);
p_val = sqrt(x^2 + y^2);

% 预存储结果
E2x_sigma = zeros(size(r2_range));
E2y_sigma = zeros(size(r2_range));
E2z_sigma = zeros(size(r2_range));

%% 设置全局字体
set(0, 'DefaultAxesFontName', '宋体');
set(0, 'DefaultTextFontName', '宋体');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

%% 循环计算不同电导率下的场强
for r2_idx = 1:length(r2_range)
    r2 = r2_range(r2_idx);
    
    % 计算复参数 k12 和 k22
    k12 = w^2*mu*e1;
    k22 = -1i*w*mu*(r2+1i*w*e2);
    
    constant = (mu * I * l) / (2*pi); % 表达式前共有的常数项
    
    %% 定义积分中的被积表达式（使用公用参数）
    u1 = @(kc) sqrt(kc^2 - k12);
    u2 = @(kc) sqrt(kc^2 - k22);
    
    D2x = @(kc) (constant/2)*(kc*(u2(kc)-u1(kc)))/(u2(kc)*(u1(kc)+u2(kc)))*exp(-u2(kc)*h);
    D2z = @(kc) constant * (((kc^2 * (80-1)) / ((u1(kc) + u2(kc)) * (80 * u1(kc) + u2(kc))) ))* exp(-u2(kc) * h);
    
    % 定义随z_val变化的函数
    G1 = @(kc) (kc*D2x(kc)+u2(kc)*D2z(kc))*exp(-u2(kc)*z_val);
    if z_val>h
        P1 = @(kc) (constant*kc/(2*u2(kc)))*exp((h-z_val)*u2(kc));
        G2 = @(kc) (constant*kc^2/(2*u2(kc)))*exp((h-z_val)*u2(kc));
        K1 = @(kc) (constant*kc/2)*exp((h-z_val)*u2(kc));
    else
        P1 = @(kc) (constant*kc/(2*u2(kc)))*exp((z_val-h)*u2(kc));
        G2 = @(kc) (constant*kc^2/(2*u2(kc)))*exp((z_val-h)*u2(kc));
        K1 = @(kc) (constant*kc/2)*exp((z_val-h)*u2(kc));
    end
    
    if p_val == 0
        E2x_sigma(r2_idx) = 0;
        E2y_sigma(r2_idx) = 0;
        E2z_sigma(r2_idx) = 0;
    else
        sin1 = y / p_val;
        cos1 = x / p_val;
        
        % 定义Ex积分被积函数
        E11 = @(kc) (G1(kc)+G2(kc))*besselj(0, kc*p_val)*kc;
        E12 = @(kc) (G1(kc)+G2(kc))*besselj(1, kc*p_val);
        E13 = @(kc) (D2x(kc)*exp(-u2(kc)*z_val)+P1(kc))*besselj(0, kc*p_val);
        
        % 定义Ey积分被积函数
        E21 = E12;
        E22 = E11;
        
        % 定义Ez积分被积函数
        E31 = @(kc) D2z(kc)*exp(-u2(kc)*z_val)*besselj(1, kc*p_val);
        E32 = @(kc) (u2(kc)^2)*E31(kc);
        E33 = @(kc) (u2(kc)*D2x(kc)*exp(-u2(kc)*z_val)+K1(kc))*kc*besselj(1, kc*p_val);
        
        % 计算Ex分量
        integrand1 = @(kc) E11(kc);
        integrand2 = @(kc) E12(kc);
        integrand8 = @(kc) E13(kc);
        E2x_sigma(r2_idx) = (-cos1^2*integral(integrand1, 0, Inf, 'ArrayValued', true)+...
            ((cos1^2-sin1^2)/p_val)*integral(integrand2, 0, Inf, 'ArrayValued', true))/(w*mu*e2*1i+mu*r2)-...
            1i*w*integral(integrand8, 0, Inf, 'ArrayValued', true);
        
        % 计算Ey分量
        integrand3 = @(kc) E21(kc);
        integrand4 = @(kc) E22(kc);
        waypoints = logspace(0, 5, 50); 
        E2y_sigma(r2_idx) = ((2*sin1*cos1/p_val)*integral(integrand3, 0, Inf, 'ArrayValued', true,'Waypoints',waypoints)+...
            (-sin1*cos1)*integral(integrand4, 0, Inf, 'ArrayValued', true,'Waypoints',waypoints))/(w*mu*e2*1i+mu*r2);
        
        % 计算Ez分量
        integrand5 = @(kc) E31(kc);
        integrand6 = @(kc) E32(kc);
        integrand7 = @(kc) E33(kc);
        E2z_sigma(r2_idx) = -1i*w*cos1*integral(integrand5, 0, Inf, 'ArrayValued', true)+...
            (cos1/(1i*w*mu*e2+mu*r2))*(integral(integrand6, 0, Inf, 'ArrayValued', true)+...
            integral(integrand7, 0, Inf, 'ArrayValued', true));
    end
end

%% 绘制场强随电导率变化曲线
figure;

% Ex分量
subplot(3,1,1);
plot(r2_range, abs(E2x_sigma), 'b-', 'LineWidth', 1.5);
xlabel('海水电导率 \sigma_2 (S/m)', 'FontName', '宋体', 'FontSize', 10);
ylabel('|E_x| (V/m)', 'FontName', '宋体', 'FontSize', 10);
title(sprintf('x方向电场强度 (观测点: [%d,%d,%d]m, f=%dHz)', x, y, z_val, f),...
    'FontName', '宋体', 'FontSize', 10);
grid on;

% Ey分量
subplot(3,1,2);
plot(r2_range, abs(E2y_sigma), 'r-', 'LineWidth', 1.5);
xlabel('海水电导率 \sigma_2 (S/m)', 'FontName', '宋体', 'FontSize', 10);
ylabel('|E_y| (V/m)', 'FontName', '宋体', 'FontSize', 10);
title(sprintf('y方向电场强度 (观测点: [%d,%d,%d]m, f=%dHz)', x, y, z_val, f),...
    'FontName', '宋体', 'FontSize', 10);
grid on;

% Ez分量
subplot(3,1,3);
plot(r2_range, abs(E2z_sigma), 'g-', 'LineWidth', 1.5);
xlabel('海水电导率 \sigma_2 (S/m)', 'FontName', '宋体', 'FontSize', 10);
ylabel('|E_z| (V/m)', 'FontName', '宋体', 'FontSize', 10);
title(sprintf('z方向电场强度 (观测点: [%d,%d,%d]m, f=%dHz)', x, y, z_val, f),...
    'FontName', '宋体', 'FontSize', 10);
grid on;

%% 输出结果表格
results_table = table(r2_range', abs(E2x_sigma'), abs(E2y_sigma'), abs(E2z_sigma'),...
    'VariableNames', {'电导率(S/m)', '|Ex|(V/m)', '|Ey|(V/m)', '|Ez|(V/m)'});

% 显示部分结果
disp('电场分量随电导率变化结果（部分显示）：');
disp(results_table(1:10:end,:));

% 保存结果到Excel
filename = '电场分量随电导率变化.xlsx';
writetable(results_table, filename);
disp(['结果已保存到文件: ' filename]);