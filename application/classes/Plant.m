classdef Plant
    properties
        f; h;
        x;
    end
    


    methods
        % construtor
        function obj = Plant(f, h, x)
            if nargin > 0
                obj.f = f;
                obj.h = h;
                obj.x = x;
            end
        end


        % rodar iteracao
        function [obj, y_true] = run(obj, u)
            y_true = obj.h(obj.x, u);
            obj.x = obj.f(obj.x, u);
        end

        % definir estado inicial
        function obj = setInitialState(obj, x)
            obj.x = x;
        end
    end
end