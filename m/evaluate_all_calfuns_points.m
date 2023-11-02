load('../data/dfo.dat');
probtypes = ["absnormal", "absuniform", "abswild", "noisy3", "nondiff", "relnormal", "reluniform", "relwild", "smooth", "wild3"];

num_pts = 3;

Results = cell(length(probtypes), size(dfo, 1), num_pts);
for p = 1:length(probtypes)
    for row = 1:size(dfo, 1)

            nprob = dfo(row, 1);
            n = dfo(row, 2);
            m = dfo(row, 3);
            factor_power = dfo(row, 4);

            BenDFO.nprob = nprob;
            BenDFO.m = m;
            BenDFO.n = n;

        for pt = 1:num_pts
            if pt == 1
                X0 = dfoxs(n, nprob, 10^factor_power);
            elseif pt == 2
                X0 = 0.1 * ones(1, n);
            elseif pt == 3
                X0 = 0.1 * [1:n];
            end

            [y, F, G, J] = calfun(X0, BenDFO, probtypes(p));

            Results{p, row, pt}.X0 = X0;
            Results{p, row, pt}.y = y;
            Results{p, row, pt}.F = F;
            Results{p, row, pt}.G = G;
            Results{p, row, pt}.J = J;
        end
    end
end
save("fvec_and_gradients_at_starting_values_matlab.mat", "Results");
