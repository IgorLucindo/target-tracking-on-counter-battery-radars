% definir parametros de simulacao a partir da comparacao escolhida
function [simFlags, resultsFilePath, model, T_ref, numOfT, numOfP] = setSimParams(compare)
    % mapa das flags
    flags = ["P_kf" "P_ekf" "T_kf" "T_ekf" "Ti_kf" "Ti_ekf" "Rad_kf" "Rad_ekf"];

    % definir parametros de simulacao
    switch compare
        case "models"
            simFlags = dictionary(flags, [0 0 1 1 0 0 0 0]);
            resultsFilePath = "results\weights_model.mat";
            model = "nonLinear";
            T_ref = 0.05;
            numOfT = 1;
            numOfP = 1;
            return;

        case "P_kf"
            simFlags = dictionary(flags, [1 0 0 0 0 0 0 0]);
            resultsFilePath = "results\weights_P_kf.mat";
            model = "linear";
            T_ref = 0.05;
            numOfT = 1;
            numOfP = 5;
            return;

        case "P_ekf"
            simFlags = dictionary(flags, [0 1 0 0 0 0 0 0]);
            resultsFilePath = "results\weights_P_ekf.mat";
            model = "nonLinear";
            T_ref = 0.05;
            numOfT = 1;
            numOfP = 5;
            return;

        case "T_kf"
            simFlags = dictionary(flags, [0 0 1 0 0 0 0 0]);
            resultsFilePath = "results\weights_T_kf.mat";
            model = "linear";
            T_ref = 0.005;
            numOfT = 4;
            numOfP = 1;
            return;

        case "T_ekf"
            simFlags = dictionary(flags, [0 0 0 1 0 0 0 0]);
            resultsFilePath = "results\weights_T_ekf.mat";
            model = "nonLinear";
            T_ref = 0.005;
            numOfT = 4;
            numOfP = 1;
            return;

        case "interp_kf"
            simFlags = dictionary(flags, [0 0 1 0 1 0 0 0]);
            resultsFilePath = "results\weights_interp_kf.mat";
            model = "linear";
            T_ref = 0.005;
            numOfT = 4;
            numOfP = 1;
            return;

        case "interp_ekf"
            simFlags = dictionary(flags, [0 0 0 1 0 1 0 0]);
            resultsFilePath = "results\weights_interp_kf.mat";
            model = "linear";
            T_ref = 0.005;
            numOfT = 4;
            numOfP = 1;
            return;

        otherwise
            fprintf("comparison unrecognized\n");
            error("");
    end
end