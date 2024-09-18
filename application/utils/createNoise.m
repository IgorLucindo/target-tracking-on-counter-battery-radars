% gerar ruido
function n = createNoise(sigma2_n, time, T)
    arrayLength = ceil(time/T) + 50;
    n = randn(3, arrayLength);
    n = n - mean(n, 2)*ones(1, arrayLength);
    n = n*sqrt(sigma2_n)./std(n, 0, 2);
end