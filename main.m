clc
close all
clear all

%% 初始化
global m1 m2 l1 l2 g mu1 mu2 Q Fk threshold  % 声明为全局变量，便于在函数内部使用 Объявить как глобальную переменную для удобства использования внутри функции
m1 = 2;
m2 = 6;
l1 = 6;
l2 = 2;
g  = 9.8;
mu1 = 0.1;
mu2 = 0.1;
% m1 = 1;   % 简单参数，用于调试代码
% m2 = 1;
% l1 = 1;
% l2 = 1;
% g  = 9.8;
% mu1 = 0.01;
% mu2 = 0.01;

threshold = 1;                      % 控制信号阶跃点，threshold用于后面生成控制信号 Точка ступени управляющего сигнала, порог, используется для последующей генерации управляющего сигнала

% X - 4x1, Y - 2x1, u - 4x1
% A - 4X4, B - 4x4, C - 2x4
B = [0, 0, 0, 0;...                 % 根据状态表达式得出 Производная от выражения состояния
    0, 1, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 1;];

C = [1,0,0,0;...                    % 根据系统输出得矩阵C В соответствии с выходом системы матрица C
    0,0,1,0];

dt = 0.01;
t = (0:dt:8)';
N_t = length(t);
X0 = [0;0;0;0];                     % 初始条件 Начальные условия

%% rk4_explicit
X_rk4 = lab_ode_rk4_explicit_x(@derivative_x, t, X0);
Y_rk4 = C * X_rk4;

figure(1)
clf;
subplot(2,1,1)
plot(t, Y_rk4)
xlabel('Время, с')
ylabel('\alpha,рад')
legend('\alpha1','\alpha2');
title('Идеальный сигнал')
grid on
grid minor

subplot(2,1,2)
plot(t, X_rk4(2,:), t, X_rk4(4,:))
xlabel('Время, с')
ylabel('\omega,рад/c')
legend('\omega1','\omega2');
grid on
grid minor

%% 解析法求雅可比矩阵解析表达式 Аналитические выражения для матриц Якоби
syms x1 x2 x3 x4;
B1 = cos(x1-x3);
B2 = sin(x1-x3);
B3 = sin(x1);
B4 = sin(x3);

f1 = x2;
f2 = (-m2*l1^2*x2^2*B1*B2 + m2*g*l1*B1*B4 + l1/l2*mu2*x4*B1 ...
    -m2*l1*l2*x4^2*B2 - (m1+m2)*g*l1*B3 - mu1*x2)/((m1+m2)*l1^2 - m2*l1^2*B1^2);
f3 = x4;
f4 = (m2*l1*l2*x2^2*B2 - m2*l2*g*B4 - mu2*x4 - m2*l1*l2*B1*f2)/(m2*l2^2);

J = jacobian([f1,f2,f3,f4],[x1,x2,x3,x4]);  % 雅可比矩阵解析表达式 Аналитические выражения матрицы Якоби
%% EKF
% X - 4x1, Y - 2x1, u - 4x1
% A - 4X4, B - 4x4, C - 2x4
% 输入信号 Входной сигнал
r_k1 = (t>=threshold) .* ones(N_t,1);
r_k2 = (t>=threshold) .* ones(N_t,1);
r_k = [zeros(N_t, 1), r_k1, zeros(N_t, 1),r_k2]';
% 噪声信号 Шумовые сигналы
q = 4;
r = 2;
G = eye(q);                      %噪声驱动矩阵 Матрица, управляемая шумом
H = eye(r);
Q = B * (0.001 * eye(q)) * B';   %噪声方差   u = r + w   B*u = B*r + B*w  噪声信号为B*w--> w  w服从N(0,BQB')
R = 0.001 * eye(r);
% 无噪声情况 Отсутствие шума
% Q = 0 * G * (0.005 * eye(q)) * G';  
% R = 0 * 0.001 * eye(r);

rng(1000);                              % 指定 MATLAB随机数生成器的种子 Укажите затравку для генератора случайных чисел MATLAB
pdX = makedist('Normal', 'mu', 0, 'sigma', sqrt( Q(2,2)));
pdY = makedist('Normal', 'mu', 0, 'sigma', sqrt( R(1,1)));

% EKF初始化
% P0 = Q;                                 % 初始先验误差协方差矩阵
P0 = zeros(length(X0),length(X0));      % 初始先验误差协方差矩阵 Начальная априорная ковариационная матрица ошибок
I = eye(4);
X_withnoise = zeros(q, length(t));
Y_withnoise = zeros(r, length(t));      % 测量信号 Измерительные сигналы
X_est = zeros(length(X0), length(t));   % 保存滤波结果 Сохранить результаты фильтрации
X_est(:,1) = X0;
P_nrm = zeros(1,length(t));             % 保存误差协方差矩阵的迹 Сохранение следа ковариационной матрицы ошибки
P_pr_nrm = zeros(1,length(t));          % 保存先验误差协方差矩阵的迹 Сохранение следа априорной ковариационной матрицы ошибки

x_kk = X0;
P_kk = P0;
x_k = X0;
h = dt/10;      %龙格库塔法求解步长
for k = 2:1:length(t)      
    %%%%%%%%%%%%预测%%%%%%%%%%%%
    x_k1k1 = x_kk;
    T = (t(k-1): h: t(k));
    x_kk1_rk4 = lab_ode_rk4_explicit_x(@derivative_x, T, x_k1k1);   % 在初值为x_k1k1，[t(k-1),t(k)]区间用rk4解微分方程
    x_kk1 = x_kk1_rk4(:,length(x_kk1_rk4));                         % 取rk4计算结果的最后一个 即x_kk1(tk)
    
    P_k1k1 = P_kk;
%     %разложение Холецкого LL分解
%     L = my_chol(P_kk);
%     P_k1k1 = L'*L;

%     % 解析法得雅可比矩阵
%     Jac = subs(J, [x1,x2,x3,x4], x_k1k1');                          % 带入数值求解析表达式
%     Fk = eval(Jac);                                                % 将表达式结果转换为 普通数值
    % 数值方法得雅可比矩阵
    Fk = linearization_digital(x_k1k1, 0.1);                     
    
    P_kk1 = lab_ode_rk4_explicit_p(@derivative_P, T, P_k1k1);       % 龙格库塔法解微分方程 Метод Лонгакурта для решения дифференциальных уравнений
    P_pr_nrm(:, k) = trace(P_kk1);
   
    %%%%%%%%%%%构造测量信号 y_k %%%%%%%%%%%
    w_k = random(pdX, q, 1);
%     x_k = lab_ode_rk4_explicit_x_withnoise(@derivative_x_withnoise, T, x_k1k1, w_k);
    x_k = lab_ode_rk4_explicit_x_withnoise(@derivative_x_withnoise, T, x_k, w_k);
    
    v_k = random(pdY, r, 1);
%     y_k = C*X_rk4(:,k) + H*v_k;
    y_k = C*x_k + H*v_k;

    X_withnoise(:,k) = x_k;
    Y_withnoise(:,k) = y_k;

    %%%%%%%%%%%%校正%%%%%%%%%%%%  
    S_k = C*P_kk1*C' + R;
    z_k = y_k - C*x_kk1;
    % S_k 直接求逆
    K_k = P_kk1*C'*S_k^(-1);
%     % QR分解求S_k的逆
%     [Q_de,R_de] = qr(S_k);
%     K_k = P_kk1*C' * (Q_de'*R_de);
    x_kk = x_kk1 + K_k * z_k;
    M    = I - K_k * C;
    P_kk = M*P_kk1*M' + K_k*R*K_k';    
    
    X_est(:,k) = x_kk;
    P_nrm(:, k) = trace(P_kk);
end
y_kk = C * X_est;
save('t.mat','t');
save('X_rk4.mat','X_rk4');
save('X_est_digital_Fk.mat','X_est');
save('P_nrm_digital_Fk.mat','P_nrm');

% save('X_est_analy_Fk.mat','X_est');
% save('P_nrm_analy_Fk.mat','P_nrm');

% figure(2)
% clf;
% subplot(2,1,1)
% plot(t, Y_withnoise)
% xlabel('Время, с')
% ylabel('\alpha,рад')
% legend('$\alpha_1$','$\alpha_2$','Interpreter','latex');
% title('Измерение \alpha с шумом')
% grid on
% grid minor
% subplot(2,1,2)
% plot(t, X_withnoise(2,:), t, X_withnoise(4,:))
% xlabel('Время, с')
% ylabel('\omega,рад/c')
% legend('$\omega_1$','$\omega_2$','Interpreter','latex');
% title('Измерение \omega с шумом')
% grid on
% grid minor

figure(3)
clf;
subplot(2,1,1)
plot(t, Y_rk4(1,:), t, y_kk(1,:))
xlabel('Время, с')
ylabel('\alpha,рад')
legend('$\alpha_1$','$\hat{\alpha}_1$','Interpreter','latex');
title('Оценка \alpha')
grid on
grid minor

subplot(2,1,2)
plot(t, Y_rk4(2,:), t, y_kk(2,:))
xlabel('Время, с')
ylabel('\alpha,рад/c')
legend('$\alpha_2$','$\hat{\alpha}_2$','Interpreter','latex');
grid on
grid minor

figure(4)
clf;
subplot(2,1,1)
plot(t, X_rk4(2,:), t, X_est(2,:))
xlabel('Время, с')
ylabel('\omega,рад/c')
legend('$\omega_1$','$\hat{\omega}_1$','Interpreter','latex');
title('Оценка \omega')
grid on
grid minor

subplot(2,1,2)
plot(t, X_rk4(4,:), t, X_est(4,:))
xlabel('Время, с')
ylabel('\omega,рад/c')
legend('$\omega_2$','$\hat{\omega}_2$','Interpreter','latex');
grid on
grid minor

figure(5)
clf;
subplot(2,1,1)
plot(t, P_pr_nrm);
xlabel('Время, с')
ylabel('tr(P_{kk1})')
grid on
grid minor
subplot(2,1,2)
plot(t, P_nrm);
xlabel('Время, с')
ylabel('tr(P_{kk})')
grid on
grid minor

compare_Jacobian

function dx = derivative_x(x, t)
    global m1 m2 l1 l2 g mu1 mu2 threshold;
    u1 = 0; u2 = 0;             % 模拟系统输入信号 Аналоговой входные сигналы системы
    if t >= threshold
        u1 = 1;
        u2 = 1;
    else
        u1 = 0;
        u2 = 0;
    end

    B1 = cos(x(1)-x(3));
    B2 = sin(x(1)-x(3));
    B3 = sin(x(1));
    B4 = sin(x(3));
    
    dx=zeros(length(x),1);      % 4x1向量
    dx(1) = x(2);
    dx(2) = (-m2*l1^2*x(2)^2*B1*B2 + m2*g*l1*B1*B4 + l1/l2*mu2*x(4)*B1 ...
        -m2*l1*l2*x(4)^2*B2 - (m1+m2)*g*l1*B3 - mu1*x(2))/((m1+m2)*l1^2 - m2*l1^2*B1^2) + u1;
    dx(3) = x(4);
    dx(4) = (m2*l1*l2*x(2)^2*B2 - m2*l2*g*B4 - mu2*x(4) - m2*l1*l2*B1*dx(2))/(m2*l2^2) + u2;

%     dx(1) = x(2);     % 简单动力学系统，用于调试代码
%     dx(2) = u1;
%     dx(3) = x(4);
%     dx(4) = u2;
end

function dx = derivative_x_withnoise(x, t, w)
    global m1 m2 l1 l2 g mu1 mu2 threshold;
    u1 = 0; u2 = 0;             % 模拟系统输入信号
    if t >= threshold
        u1 = 1;
        u2 = 1;
    else
        u1 = 0;
        u2 = 0;
    end

    B1 = cos(x(1)-x(3));
    B2 = sin(x(1)-x(3));
    B3 = sin(x(1));
    B4 = sin(x(3));
    
    dx=zeros(length(x),1);      % 4x1向量
    dx(1) = x(2) + w(1);
    dx(2) = (-m2*l1^2*x(2)^2*B1*B2 + m2*g*l1*B1*B4 + l1/l2*mu2*x(4)*B1 ...
        -m2*l1*l2*x(4)^2*B2 - (m1+m2)*g*l1*B3 - mu1*x(2))/((m1+m2)*l1^2 - m2*l1^2*B1^2) + u1 + w(2);
    dx(3) = x(4) + w(3);
    dx(4) = (m2*l1*l2*x(2)^2*B2 - m2*l2*g*B4 - mu2*x(4) - m2*l1*l2*B1*dx(2))/(m2*l2^2) + u2 + w(4);

%     dx(1) = x(2);     % 简单动力学系统，用于调试代码
%     dx(2) = u1;
%     dx(3) = x(4);
%     dx(4) = u2;
end

function dP = derivative_P(P, t)
    global Fk Q
    dP = Fk*P + P*Fk' + Q;
end

function J  = linearization_digital(x,eps)     
    L = length(x);
    J = zeros(L, L);
    for i = 1:L
        m_r = x;
        m_r(i) = m_r(i) + eps;
        m_l = x;
        m_l(i) = m_l(i) - eps;
        J(:,i) = (derivative_x(m_r, 0) - derivative_x(m_l, 0))/(2*eps); % 注意用右侧点-左侧点
    end
end

function [l]=my_chol(A)
    [m,n]=size(A);
    l=zeros(n,n);
    for i=1:n
        for j=1:i
            if i==j
                l(i,i)=sqrt(A(i,i)-sum_l2(i)); 
            else
                l(i,j)=1/l(j,j) * (A(i,j)-sum_ll(i,j));
            end
        end
    end
            
    function res=sum_l2(i)
        res=0;
        if i==1
            res=0;
        else
            for k=1:i-1
                res=res+l(i,k)^2;
            end
        end
    end
    
    function res=sum_ll(i,j)
        res=0;
        if j==1
            res=0;
        else
            for k=1:j-1
                res=res+l(i,k)*l(j,k);
            end
        end
    end
end