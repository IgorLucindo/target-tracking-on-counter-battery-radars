classdef Plant
    properties
        A; B; C;
        x;
    end
    


    methods
        % construtor
        function obj = Plant(A, B, C, x0)
            if nargin > 0
                obj.A = A;
                obj.B = B;
                obj.C = C;
                obj.x = x0;
            end
        end


        % rodar iteracao
        function [obj, y_true] = run(obj, u)
            y_true = obj.C*obj.x;
            obj.x = obj.A*obj.x + obj.B*u;
        end


        % definir estado inicial
        function obj = setX0(obj, x)
            obj.x = x;
        end
    end
end