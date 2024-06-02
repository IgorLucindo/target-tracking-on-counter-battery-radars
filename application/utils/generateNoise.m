% gerar ruido
function v = generateNoise(sigma2_v, time, Ts)
    arrayLength = ceil(time/Ts) + 50;
    v = randn(3, arrayLength);
    v = v - mean(v, 2)*ones(1, arrayLength);
    v = v.*sqrt(sigma2_v)./std(v, 0, 2);
end