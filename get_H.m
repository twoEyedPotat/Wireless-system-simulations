% Get the channel co-efficient rayleigh distributed
% E[|h|^2] = 1 ==> complex gaussian circular symmetric with variance 1
function [H] = get_H(Nr,Nt)
    H = sqrt(0.5)*(randn(1, Nr*Nt) + j*randn(1, Nr*Nt));
    H = reshape(H, Nr, Nt);
end