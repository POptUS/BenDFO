load('../data/dfo.dat');
probtypes = ["absnormal", "absuniform", "abswild", "noisy3", "nondiff", "relnormal", "reluniform", "relwild", "smooth", "wild3"];

Results = cell(length(probtypes), size(dfo, 1));
for p = 1:length(probtypes)
    for row = 1:size(dfo, 1)
        nprob = dfo(row, 1);
        n = dfo(row, 2);
        m = dfo(row, 3);
        factor_power = dfo(row, 4);

        X0 = dfoxs(n, nprob, 10^factor_power);

        BenDFO.nprob = nprob;
        BenDFO.m = m;
        BenDFO.n = n;

        [y, F, G, J] = calfun(X0, BenDFO, probtypes(p));

        Results{p, row}.X0 = X0;
        Results{p, row}.y = y;
        Results{p, row}.F = F;
        Results{p, row}.G = G;
        Results{p, row}.J = J;
    end
end
save("fvec_and_gradients_at_starting_values_matlab.mat", "Results");
