% clc
% clear all
% close all

t = load('t.mat').t;
X_rk4 = load('X_rk4.mat').X_rk4;
X_est_digital_Fk = load('X_est_digital_Fk.mat').X_est;
X_est_analy_Fk = load('X_est_analy_Fk.mat').X_est;
P_nrm_digital_Fk = load('P_nrm_digital_Fk.mat').P_nrm;
P_nrm_analy_Fk = load('P_nrm_analy_Fk.mat').P_nrm;

figure(6)
subplot(2,1,1)
plot(t, P_nrm_analy_Fk, t, P_nrm_digital_Fk)
xlabel('Время, с')
ylabel('tr(P_{kk})')
legend('J получен аналитически','J получен численно');
grid on 
grid minor
subplot(2,1,2)
bar([0.8,1.8,2.8,3.8], [norm(X_rk4(1,:)-X_est_analy_Fk(1,:)), norm(X_rk4(2,:)-X_est_analy_Fk(2,:))...
    norm(X_rk4(3,:)-X_est_analy_Fk(3,:)), norm(X_rk4(4,:)-X_est_analy_Fk(4,:))]/length(t),0.3,'b')
hold on
bar([1.2,2.2,3.2,4.2], [norm(X_rk4(1,:)-X_est_digital_Fk(1,:)), norm(X_rk4(2,:)-X_est_digital_Fk(2,:))...
    norm(X_rk4(3,:)-X_est_digital_Fk(3,:)), norm(X_rk4(4,:)-X_est_digital_Fk(4,:))]/length(t),0.3,'r')
xlabel('Номер элемента вектора состояния')
ylabel('Нормированная ошибка')
legend('J получен аналитически','J получен численно','location','best');
grid on 
grid minor
