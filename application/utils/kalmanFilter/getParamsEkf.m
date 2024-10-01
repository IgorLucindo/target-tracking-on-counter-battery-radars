% pegar funcoes do modelo nao linear para Kalman Filter Extendido
function [f, h, F, H, f_rev, F_rev] = getParamsEkf(T)
    f = @innerf;
    h = @innerh;
    F = @innerF;
    H = @innerH;
    f_rev = @innerf_rev;
    F_rev = @innerF_rev;

    T_rev = -T;
    h_ = 1.01*1e-4;

    % funcao f
    function f = innerf(x, u)
        px = x(1); py = x(2); pz = x(3); vx = x(4); vy = x(5); vz = x(6); gama = x(7);
        g = -u(3);
        v = sqrt(vx^2 + vy^2 + vz^2);

        f = [px + vx*T - gama*v*vx*T^2/4;
             py + vy*T - gama*v*vy*T^2/4;
             pz + vz*T - gama*v*vz*T^2/4 - g*T^2/2;
             vx + (T^2*((gama^2*vx*vy^2)/4 + (gama*vx*((gama*v)/2 + (gama*vx^2)/(2*v))*v)/2 + (gama*vx*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vx*vz*v)/2))/2 - (T*gama*vx*v)/2;
             vy + (T^2*((gama^2*vx^2*vy)/4 + (gama*vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama*vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vy*vz*v)/2))/2 - (T*gama*vy*v)/2;
             vz + (T^2*(((gama*v)/2 + (gama*vz^2)/(2*v))*(g + (gama*vz*v)/2) + (gama^2*vx^2*vz)/4 + (gama^2*vy^2*vz)/4 + (gama*h_*vz^2*v)/2))/2 - T*(g + (gama*vz*v)/2);
             gama + (T^2*(gama*h_^2*vz^2 + gama*h_*(g + (gama*vz*v)/2)))/2 - T*gama*h_*vz];
    end
    

    % funcao h
    function h = innerh(x, u)
        h = x(1:3);
    end

    % jacobiana de f
    function F = innerF(x, u)
        vx = x(4); vy = x(5); vz = x(6); gama = x(7);
        g = -u(3);

        v = sqrt(vx^2 + vy^2 + vz^2);

        F = [1, 0, 0, T - (T^2*gama*v)/4 - (T^2*gama*vx^2)/(4*v), -(T^2*gama*vx*vy)/(4*v), -(T^2*gama*vx*vz)/(4*v), -(T^2*vx*v)/4;
             0, 1, 0, -(T^2*gama*vx*vy)/(4*v), T - (T^2*gama*v)/4 - (T^2*gama*vy^2)/(4*v), -(T^2*gama*vy*vz)/(4*v), -(T^2*vy*v)/4;
             0, 0, 1, -(T^2*gama*vx*vz)/(4*v), -(T^2*gama*vy*vz)/(4*v), T - (T^2*gama*v)/4 - (T^2*gama*vz^2)/(4*v), -(T^2*vz*v)/4;
             0, 0, 0, (T^2*((gama^2*vx*vy)/2 + (gama*vy*((gama*vx)/(2*v) - (gama*vx*vy^2)/(2*v^3))*v)/2 + (gama^2*vx*vy*vz^2)/(4*v^2) + (gama*vx*vy*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) - (gama*vx*vy*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vx*vy*vz)/(2*v)))/2 - (T*gama*vx*vy)/(2*v), (T^2*((gama^2*vx^2)/4 + (gama*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama^2*vy^2*vz^2)/(4*v^2) + (gama*vy^2*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*((3*gama*vy)/(2*v) - (gama*vy^3)/(2*v^3))*v)/2 + (gama*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vz*v)/2 - (gama*vy^2*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy^2*vz)/(2*v)))/2 - (T*gama*v)/2 - (T*gama*vy^2)/(2*v) + 1, (T^2*((gama*vy*((gama*vz)/(2*v) - (gama*vy^2*vz)/(2*v^3))*v)/2 + (gama*vy*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vy*v)/2 + (gama*vy*vz*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) - (gama*vy*vz^2*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy*vz^2)/(2*v)))/2 - (T*gama*vy*vz)/(2*v), (T^2*((gama*vx^2*vy)/2 + (gama*vy*vz^2)/4 + (vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*vy*(v/2 + vy^2/(2*v))*v)/2 + (h_*vy*vz*v)/2))/2 - (T*vy*v)/2;
             0, 0, 0, (T^2*((gama^2*vx*vy)/2 + (gama*vy*((gama*vx)/(2*v) - (gama*vx*vy^2)/(2*v^3))*v)/2 + (gama^2*vx*vy*vz^2)/(4*v^2) + (gama*vx*vy*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) - (gama*vx*vy*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vx*vy*vz)/(2*v)))/2 - (T*gama*vx*vy)/(2*v), (T^2*((gama^2*vx^2)/4 + (gama*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama^2*vy^2*vz^2)/(4*v^2) + (gama*vy^2*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*((3*gama*vy)/(2*v) - (gama*vy^3)/(2*v^3))*v)/2 + (gama*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vz*v)/2 - (gama*vy^2*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy^2*vz)/(2*v)))/2 - (T*gama*v)/2 - (T*gama*vy^2)/(2*v) + 1, (T^2*((gama*vy*((gama*vz)/(2*v) - (gama*vy^2*vz)/(2*v^3))*v)/2 + (gama*vy*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vy*v)/2 + (gama*vy*vz*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) - (gama*vy*vz^2*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy*vz^2)/(2*v)))/2 - (T*gama*vy*vz)/(2*v), (T^2*((gama*vx^2*vy)/2 + (gama*vy*vz^2)/4 + (vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*vy*(v/2 + vy^2/(2*v))*v)/2 + (h_*vy*vz*v)/2))/2 - (T*vy*v)/2;
             0, 0, 0, (T^2*(((gama*vx)/(2*v) - (gama*vx*vz^2)/(2*v^3))*(g + (gama*vz*v)/2) + (gama^2*vx*vz)/2 + (gama*vx*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) + (gama*h_*vx*vz^2)/(2*v)))/2 - (T*gama*vx*vz)/(2*v), (T^2*(((gama*vy)/(2*v) - (gama*vy*vz^2)/(2*v^3))*(g + (gama*vz*v)/2) + (gama^2*vy*vz)/2 + (gama*vy*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) + (gama*h_*vy*vz^2)/(2*v)))/2 - (T*gama*vy*vz)/(2*v), (T^2*((g + (gama*vz*v)/2)*((3*gama*vz)/(2*v) - (gama*vz^3)/(2*v^3)) + ((gama*v)/2 + (gama*vz^2)/(2*v))^2 + (gama^2*vx^2)/4 + (gama^2*vy^2)/4 + (gama*h_*vz^3)/(2*v) + gama*h_*vz*v))/2 - T*((gama*v)/2 + (gama*vz^2)/(2*v)) + 1, (T^2*((g + (gama*vz*v)/2)*(v/2 + vz^2/(2*v)) + (gama*vx^2*vz)/2 + (gama*vy^2*vz)/2 + (vz*((gama*v)/2 + (gama*vz^2)/(2*v))*v)/2 + (h_*vz^2*v)/2))/2 - (T*vz*v)/2;
             0, 0, 0, (T^2*gama^2*h_*vx*vz)/(4*v), (T^2*gama^2*h_*vy*vz)/(4*v), (T^2*(2*gama*vz*h_^2 + gama*((gama*v)/2 + (gama*vz^2)/(2*v))*h_))/2 - T*gama*h_, (T^2*(h_*((gama*vz*v)/2 + 981/100) + h_^2*vz^2 + (gama*h_*vz*v)/2))/2 - T*h_*vz + 1];
    end
    
    % jacobiana de h
    function H = innerH(x)
        H = [eye(3) zeros(3, 4)];
    end

    % funcao f inversa
    function f_rev = innerf_rev(x, u)
        px = x(1); py = x(2); pz = x(3); vx = x(4); vy = x(5); vz = x(6); gama = x(7);
        g = -u(3);
        v = sqrt(vx^2 + vy^2 + vz^2);

        f_rev = [px + vx*T_rev - gama*v*vx*T_rev^2/4;
                 py + vy*T_rev - gama*v*vy*T_rev^2/4;
                 pz + vz*T_rev - gama*v*vz*T_rev^2/4 - g*T_rev^2/2;
                 vx + (T_rev^2*((gama^2*vx*vy^2)/4 + (gama*vx*((gama*v)/2 + (gama*vx^2)/(2*v))*v)/2 + (gama*vx*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vx*vz*v)/2))/2 - (T_rev*gama*vx*v)/2;
                 vy + (T_rev^2*((gama^2*vx^2*vy)/4 + (gama*vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama*vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vy*vz*v)/2))/2 - (T_rev*gama*vy*v)/2;
                 vz + (T_rev^2*(((gama*v)/2 + (gama*vz^2)/(2*v))*(g + (gama*vz*v)/2) + (gama^2*vx^2*vz)/4 + (gama^2*vy^2*vz)/4 + (gama*h_*vz^2*v)/2))/2 - T_rev*(g + (gama*vz*v)/2);
                 gama + (T_rev^2*(gama*h_^2*vz^2 + gama*h_*(g + (gama*vz*v)/2)))/2 - T_rev*gama*h_*vz];
    end

    % jacobiana de f inversa
    function F_rev = innerF_rev(x, u)
        vx = x(4); vy = x(5); vz = x(6); gama = x(7);
        g = -u(3);

        v = sqrt(vx^2 + vy^2 + vz^2);

        F_rev = [1, 0, 0, T_rev - (T_rev^2*gama*v)/4 - (T_rev^2*gama*vx^2)/(4*v), -(T_rev^2*gama*vx*vy)/(4*v), -(T_rev^2*gama*vx*vz)/(4*v), -(T_rev^2*vx*v)/4;
                 0, 1, 0, -(T_rev^2*gama*vx*vy)/(4*v), T_rev - (T_rev^2*gama*v)/4 - (T_rev^2*gama*vy^2)/(4*v), -(T_rev^2*gama*vy*vz)/(4*v), -(T_rev^2*vy*v)/4;
                 0, 0, 1, -(T_rev^2*gama*vx*vz)/(4*v), -(T_rev^2*gama*vy*vz)/(4*v), T_rev - (T_rev^2*gama*v)/4 - (T_rev^2*gama*vz^2)/(4*v), -(T_rev^2*vz*v)/4;
                 0, 0, 0, (T_rev^2*((gama^2*vx*vy)/2 + (gama*vy*((gama*vx)/(2*v) - (gama*vx*vy^2)/(2*v^3))*v)/2 + (gama^2*vx*vy*vz^2)/(4*v^2) + (gama*vx*vy*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) - (gama*vx*vy*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vx*vy*vz)/(2*v)))/2 - (T_rev*gama*vx*vy)/(2*v), (T_rev^2*((gama^2*vx^2)/4 + (gama*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama^2*vy^2*vz^2)/(4*v^2) + (gama*vy^2*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*((3*gama*vy)/(2*v) - (gama*vy^3)/(2*v^3))*v)/2 + (gama*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vz*v)/2 - (gama*vy^2*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy^2*vz)/(2*v)))/2 - (T_rev*gama*v)/2 - (T_rev*gama*vy^2)/(2*v) + 1, (T_rev^2*((gama*vy*((gama*vz)/(2*v) - (gama*vy^2*vz)/(2*v^3))*v)/2 + (gama*vy*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vy*v)/2 + (gama*vy*vz*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) - (gama*vy*vz^2*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy*vz^2)/(2*v)))/2 - (T_rev*gama*vy*vz)/(2*v), (T_rev^2*((gama*vx^2*vy)/2 + (gama*vy*vz^2)/4 + (vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*vy*(v/2 + vy^2/(2*v))*v)/2 + (h_*vy*vz*v)/2))/2 - (T_rev*vy*v)/2;
                 0, 0, 0, (T_rev^2*((gama^2*vx*vy)/2 + (gama*vy*((gama*vx)/(2*v) - (gama*vx*vy^2)/(2*v^3))*v)/2 + (gama^2*vx*vy*vz^2)/(4*v^2) + (gama*vx*vy*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) - (gama*vx*vy*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vx*vy*vz)/(2*v)))/2 - (T_rev*gama*vx*vy)/(2*v), (T_rev^2*((gama^2*vx^2)/4 + (gama*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (gama^2*vy^2*vz^2)/(4*v^2) + (gama*vy^2*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*((3*gama*vy)/(2*v) - (gama*vy^3)/(2*v^3))*v)/2 + (gama*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vz*v)/2 - (gama*vy^2*vz*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy^2*vz)/(2*v)))/2 - (T_rev*gama*v)/2 - (T_rev*gama*vy^2)/(2*v) + 1, (T_rev^2*((gama*vy*((gama*vz)/(2*v) - (gama*vy^2*vz)/(2*v^3))*v)/2 + (gama*vy*(g + (gama*vz*v)/2))/(2*v) + (gama*h_*vy*v)/2 + (gama*vy*vz*((gama*v)/2 + (gama*vy^2)/(2*v)))/(2*v) + (gama*vy*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) - (gama*vy*vz^2*(g + (gama*vz*v)/2))/(2*v^3) + (gama*h_*vy*vz^2)/(2*v)))/2 - (T_rev*gama*vy*vz)/(2*v), (T_rev^2*((gama*vx^2*vy)/2 + (gama*vy*vz^2)/4 + (vy*((gama*v)/2 + (gama*vy^2)/(2*v))*v)/2 + (vy*vz*(g + (gama*vz*v)/2))/(2*v) + (gama*vy*(v/2 + vy^2/(2*v))*v)/2 + (h_*vy*vz*v)/2))/2 - (T_rev*vy*v)/2;
                 0, 0, 0, (T_rev^2*(((gama*vx)/(2*v) - (gama*vx*vz^2)/(2*v^3))*(g + (gama*vz*v)/2) + (gama^2*vx*vz)/2 + (gama*vx*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) + (gama*h_*vx*vz^2)/(2*v)))/2 - (T_rev*gama*vx*vz)/(2*v), (T_rev^2*(((gama*vy)/(2*v) - (gama*vy*vz^2)/(2*v^3))*(g + (gama*vz*v)/2) + (gama^2*vy*vz)/2 + (gama*vy*vz*((gama*v)/2 + (gama*vz^2)/(2*v)))/(2*v) + (gama*h_*vy*vz^2)/(2*v)))/2 - (T_rev*gama*vy*vz)/(2*v), (T_rev^2*((g + (gama*vz*v)/2)*((3*gama*vz)/(2*v) - (gama*vz^3)/(2*v^3)) + ((gama*v)/2 + (gama*vz^2)/(2*v))^2 + (gama^2*vx^2)/4 + (gama^2*vy^2)/4 + (gama*h_*vz^3)/(2*v) + gama*h_*vz*v))/2 - T_rev*((gama*v)/2 + (gama*vz^2)/(2*v)) + 1, (T_rev^2*((g + (gama*vz*v)/2)*(v/2 + vz^2/(2*v)) + (gama*vx^2*vz)/2 + (gama*vy^2*vz)/2 + (vz*((gama*v)/2 + (gama*vz^2)/(2*v))*v)/2 + (h_*vz^2*v)/2))/2 - (T_rev*vz*v)/2;
                 0, 0, 0, (T_rev^2*gama^2*h_*vx*vz)/(4*v), (T_rev^2*gama^2*h_*vy*vz)/(4*v), (T_rev^2*(2*gama*vz*h_^2 + gama*((gama*v)/2 + (gama*vz^2)/(2*v))*h_))/2 - T_rev*gama*h_, (T_rev^2*(h_*(g + (gama*vz*v)/2) + h_^2*vz^2 + (gama*h_*vz*v)/2))/2 - T_rev*h_*vz + 1];
    end
end