% pegar matrizes de espaco de estado do para Kalman Filter
function [A, B, C, A_rev, B_rev] = getParamsKf(T)
    A = eye(6);
    A(1:3, 4:end) = T*eye(3);
    B = [T^2/2*eye(3); T*eye(3)];
    C = [eye(3) zeros(3)];
    A_rev = inv(A);
    B_rev = -A_rev*B;
end