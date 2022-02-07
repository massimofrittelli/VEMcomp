clearvars

L2errors = zeros(1,4);
meshsizes = zeros(1,4);
L2rates = zeros(1,3);

for it=1:4
    tau = 2^(-2-2*it);
    switch it
        case 1
            load('sphere5.mat')
        case 2
            load('sphere9.mat')
        case 3
            load('sphere17.mat')
        case 4
            load('sphere33.mat')
    end
    PARABOLIC_BS_SPHERE
    L2errors(it) = L2err_prod;
    meshsizes(it) = h;
end

for it=1:3
   L2rates(it) = log(L2errors(it+1)/L2errors(it)) / log(meshsizes(it+1)/meshsizes(it));
end