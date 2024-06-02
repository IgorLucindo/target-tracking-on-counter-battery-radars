classdef KalmanFilter
    properties
        A; B; C;
        Q; R;
        P;
        x;
    end
    


    methods
        % construtor
        function obj = KalmanFilter(A, B, C, Q, R, P0)
            if nargin > 0
                obj.A = A;
                obj.B = B;
                obj.C = C;
                obj.Q = Q;
                obj.R = R;
                obj.P = P0;
                obj.x = zeros(length(A), 1);
            end
        end


        % prever
        function obj = predict(obj, u)
            obj.x = obj.A*obj.x + obj.B*u;
            obj.P = obj.A*obj.P*obj.A' + obj.B*obj.Q*obj.B';
        end


        % atualizar
        function [obj, y_est, x_est] = update(obj, y)
            K = obj.P*obj.C'/(obj.C*obj.P*obj.C' + obj.R);
            obj.x = obj.x + K*(y - obj.C*obj.x);
            obj.P = (eye(length(obj.A)) - K*obj.C)*obj.P;
            
            y_est = obj.C*obj.x;

            x_est = obj.x;
        end


        % definir estado inicial
        function [obj, y_est, x_est] = setX0(obj, y)
            obj.x = [y; 0; 0; 200];
            y_est = y;
            x_est = obj.x;
        end


        % rodar iteracao de kalman filter
        function [obj, y_est, x_est] = run(obj, y, u)
            obj = obj.predict(u);
            [obj, y_est, x_est] = obj.update(y);
        end
    end
end