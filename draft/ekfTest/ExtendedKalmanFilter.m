classdef ExtendedKalmanFilter
    properties
        f; h; F; H;
        Q; R;
        P;
        x;
    end
    


    methods
        % construtor
        function obj = ExtendedKalmanFilter(f, h, F, H, Q, R, P0)
            if nargin > 0
                obj.f = f;
                obj.h = h;
                obj.F = F;
                obj.H = H;
                obj.Q = Q;
                obj.R = R;
                obj.P = P0;
                obj.x = zeros(length(P0), 1);
            end
        end


        % prever
        function obj = predict(obj, u)
            F_ = obj.F(obj.x, u);

            obj.x = obj.f(obj.x, u);
            obj.P = F_*obj.P*F_' + obj.Q;
        end


        % atualizar
        function [obj, y_est, x_est] = update(obj, y, u)
            H_ = obj.H(obj.x);

            K = obj.P*H_'/(H_*obj.P*H_' + obj.R);
            obj.x = obj.x + K*(y - obj.h(obj.x, u));
            obj.P = (eye(length(obj.P)) - K*H_)*obj.P;
            
            y_est = obj.h(obj.x, u);

            x_est = obj.x;
        end


        % definir estado inicial
        function [obj, y_est, x_est] = setInitialState(obj, x0)
            obj.x = x0;
            y_est = x0(1:3);
            x_est = x0;
        end


        % rodar iteracao de kalman filter Extendido
        function [obj, y_est, x_est] = run(obj, y, u)
            obj = obj.predict(u);
            [obj, y_est, x_est] = obj.update(y, u);
        end
    end
end