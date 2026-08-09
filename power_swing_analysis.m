%% 发电机转子摇摆曲线微分方程求解与暂态稳定性分析 demo
% 本程序演示如何利用 MATLAB 数值求解器 (ode45) 求解发电机转子运动方程，
% 并计算电网短路故障期间及切除后的转子功角摇摆曲线 (delta - t)，确定极限切除时间。
%
% 文件路径: D:\Matlab_simulink_system_fangzhen\power_swing_analysis.m

clear; clc; close all;

%% 1. 全局变量与系统物理参数定义
global y0 Tj Pt E U X1

y0 = 2 * pi * 50;   % 系统基准角频率 (2*pi*50 = 314.1593 rad/s)
Tj = 11.28;         % 发电机转子惯性时间常数 (s)
Pt = 1.0;           % 原动机输入机械功率 (p.u.)
E  = 1.47;          % 发电机暂态电势 (p.u.)
U  = 1.0;           % 无限大母线电压 (p.u.)

%% 2. 阶段一：计算故障发生期间的摇摆曲线 (0 ~ 0.3 秒)
X1 = 2.82;          % 短路故障期间的系统转移电抗 (p.u.)
tspan1 = [0.0 0.5]; % 故障持续计算时间范围 (s)
y1_0 = [31.54 * pi / 180; 1.0]; % 初始状态 [初始功角(rad); 初始转速(pu)]

% 使用 4阶-5阶龙格库塔法 (ode45) 求解微分方程组
[t1, YY1] = ode45(@power_tra, tspan1, y1_0);

delta1_deg = YY1(:, 1) * 180 / pi; % 转化为角度 (deg)
omega1_pu  = YY1(:, 2);            % 转速 (pu)

%% 3. 阶段二：计算故障切除后的摇摆曲线 (设定极限切除时间 tc = 0.2416s)
tc = 0.2416; % 极限切除时刻 tc = 0.2416s
delta_c_deg = interp1(t1, delta1_deg, tc); % 插入切除瞬间的功角
omega_c_pu  = interp1(t1, omega1_pu, tc);  % 插入切除瞬间的转速

X1 = 1.062;         % 故障切除后系统的恢复转移电抗 (p.u.)
tspan2 = [tc 2.0];  % 切除后计算时间范围 (s)
y2_0 = [delta_c_deg * pi / 180; omega_c_pu]; % 切除瞬间的状态 (连续不可突变)

[t2, YY2] = ode45(@power_tra, tspan2, y2_0);

delta2_deg = YY2(:, 1) * 180 / pi;
omega2_pu  = YY2(:, 2);

%% 4. 结果可视化与分析
figure('Name', '发电机转子摇摆曲线暂态稳定性分析', 'NumberTitle', 'off');

% 子图1：功角摇摆曲线 (delta - t)
subplot(2,1,1);
plot(t1, delta1_deg, 'r--', 'LineWidth', 1.5); hold on;
plot(t2, delta2_deg, 'b-', 'LineWidth', 2.0);
xline(tc, 'k:', sprintf('切除时刻 %.2fs', tc), 'LineWidth', 1.2);
yline(63.61, 'm--', '极限切除角 63.61°', 'LineWidth', 1.2);
title('发电机转子功角摇摆曲线 (\delta - t)');
xlabel('时间 t / s');
ylabel('功角 \delta / deg');
legend('故障期间 (未切除扩展)', sprintf('%.2fs 切除故障后', tc), 'Location', 'northwest');
grid on;

% 子图2：转速变化曲线 (omega - t)
subplot(2,1,2);
plot(t1, omega1_pu, 'r--', 'LineWidth', 1.5); hold on;
plot(t2, omega2_pu, 'b-', 'LineWidth', 2.0);
xline(tc, 'k:', sprintf('切除时刻 %.2fs', tc), 'LineWidth', 1.2);
title('发电机转子转速变化曲线 (\omega - t)');
xlabel('时间 t / s');
ylabel('转速 \omega / p.u.');
legend('故障期间', sprintf('%.2fs 切除故障后', tc), 'Location', 'northwest');
grid on;

% 命令行打印分析总结
fprintf('=====================================================\n');
fprintf('  发电机暂态稳定性计算结果：\n');
fprintf('  1. 故障发生时间: 0.0 s\n');
fprintf('  2. 故障切除时间: %.2f s\n', tc);
if delta_c_deg <= 63.61
    fprintf('  3. 切除瞬间功角: %.2f° (小于等于极限切除角 63.61°)\n', delta_c_deg);
else
    fprintf('  3. 切除瞬间功角: %.2f° (大于极限切除角 63.61°)\n', delta_c_deg);
end

if max(delta2_deg) < 180
    fprintf('  4. 切除后最大摇摆功角: %.2f° (小于 180°，系统稳态恢复)\n', max(delta2_deg));
    fprintf('  结论: 本切除方案下，电力系统暂态稳定！\n');
else
    fprintf('  4. 切除后最大摇摆功角: %.2f° (超过 180°，系统失步)\n', max(delta2_deg));
    fprintf('  结论: 本切除方案下，电力系统暂态不稳定 (失步)！\n');
end
fprintf('=====================================================\n');

%% 5. 内部 M 函数定义：发电机转子摇摆曲线微分方程
function Yd = power_tra(~, YY)
    % 说明：
    %   t : 时间自变量 (标量)
    %   YY: 状态列向量 [YY(1)=delta (rad); YY(2)=omega (pu)]
    %   Yd: 状态一阶导数向量 [d(delta)/dt; d(omega)/dt]
    %
    % 怎么计算的：
    %   Yd(1) = (omega - 1) * omega0  --> 功角变化率
    %   Yd(2) = (Pt - Pe) / Tj        --> 转速加速度 (不平衡转矩)
    %   其中 Pe = (E * U / X1) * sin(delta) 为电磁功率
    
    global y0 Tj Pt E U X1
    
    Yd = [ (YY(2) - 1.0) * y0 ; 
           (Pt - (E * U / X1) * sin(YY(1))) / Tj ];
end
