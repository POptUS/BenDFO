% Compares the matlab and python implementations of a method.

M = load("m/fvec_and_gradients_at_starting_values_matlab.mat");
M = M.Results;

P1 = load("py/fvec_and_gradients_at_starting_values_python.mat");
P = cell(size(M));
probtypes = ["absnormal", "absuniform", "abswild", "noisy3", "nondiff", "relnormal", "reluniform", "relwild", "smooth", "wild3"];
num_pts = 3;

for i = 1:53
    for p = 1:length(probtypes)
        for pt = 1:num_pts
            P{p, i, pt} = P1.(['prob_' int2str(p) '_' int2str(i) '_' int2str(pt)]);
        end
    end
end

msg = 'Different %s for prob %d, %s for pt: %s ';

for i = 1:53
    assert(norm(P{1, i}.X0 - M{1, i}.X0') == 0, "different starting point");
    xs = mat2str(P{1, i}.X0);

    for p = [5, 9, 10]
        ps = probtypes(p);
        for pt = 1:num_pts
            assert((P{p, i, pt}.y - M{p, i, pt}.y) / norm(P{p, i, pt}.y) <= 1e-15, sprintf(msg, 'y', 1, ps, xs));
            assert(all(P{p, i, pt}.F - M{p, i, pt}.F' <= 1e-14), sprintf(msg, 'fvec', 1, ps, xs));
            assert(norm(P{p, i, pt}.J - M{p, i, pt}.J) / norm(P{p, i, pt}.J) <= 1e-14, sprintf(msg, 'J', 1, ps, xs));
            assert(norm(P{p, i, pt}.G - M{p, i, pt}.G') / norm(P{p, i, pt}.G) <= 1e-14, sprintf(msg, 'G', 1, ps, xs));
        end
    end
end
