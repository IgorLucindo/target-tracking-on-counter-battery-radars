clear; close all; clc;

% Parameters of the system
A = [1 1; 0 1];    % State transition matrix
H = [1 0];         % Observation matrix
Q = [0.1 0; 0 0.1]; % Process noise covariance
R = 0.1;           % Measurement noise covariance

% Number of time steps
T = 500;

% True initial state
x_true = [0; 1];

% Initialize arrays to store results
x_true_store = zeros(2, T);
z_store = zeros(1, T);
x_est = zeros(2, T);
P_est = zeros(2, 2, T);

% Initial estimates
x_est(:,1) = [0; 0];
P_est(:,:,1) = eye(2);

% Generate true states and noisy measurements
for k = 2:T
    x_true = A * x_true + mvnrnd([0; 0], Q)';  % True state evolution
    z = H * x_true + normrnd(0, sqrt(R));      % Measurement with noise
    
    % Kalman filter prediction step
    x_pred = A * x_est(:,k-1);
    P_pred = A * P_est(:,:,k-1) * A' + Q;
    
    % Kalman filter update step
    K = P_pred * H' / (H * P_pred * H' + R);
    x_est(:,k) = x_pred + K * (z - H * x_pred);
    P_est(:,:,k) = (eye(2) - K * H) * P_pred;
    
    % Store true state and measurement
    x_true_store(:,k) = x_true;
    z_store(k) = z;
end

% Apply the RTS smoother
x_smooth = x_est;
P_smooth = P_est;

for k = T-1:-1:1
    % RTS smoother gain
    Ck = P_est(:,:,k) * A' / P_pred;
    
    % Smoothed state estimate
    x_smooth(:,k) = x_est(:,k) + Ck * (x_smooth(:,k+1) - A * x_est(:,k));
    
    % Smoothed error covariance
    P_smooth(:,:,k) = P_est(:,:,k) + Ck * (P_smooth(:,:,k+1) - P_pred) * Ck';
end

% Plotting results
error_est = abs(sum(x_est - x_true_store));
error_smoth = abs(sum(x_smooth - x_true_store));
figure;
plot(1:T, error_est(1, :), 'r', 'LineWidth', 2);
hold on;
plot(1:T, error_smoth(1, :), 'g', 'LineWidth', 2);
legend('Kalman Filter Estimate', 'RTS Smoothed Estimate');
xlabel('Time Step');
ylabel('State');
title('Rauch-Tung-Striebel Smoothing');
grid on;