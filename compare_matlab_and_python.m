% Compares the matlab and python implementations of a method.

M = load("m/fvec_and_gradients_at_starting_values_matlab.mat");
M = M.Results;

P1 = load("py/fvec_and_gradients_at_starting_values_python.mat");
P = cell(size(M));
probtypes = ["absnormal", "absuniform", "abswild", "noisy3", "nondiff", "relnormal", "reluniform", "relwild", "smooth", "wild3"];

for probnum = 1:53
    for p = 1:length(probtypes)
        P{p, probnum} = P1.(['prob_' int2str(p) '_' int2str(probnum)]);
    end
end

for probnum = 1:53
    assert(norm(P{1, probnum}.X0 - M{1, probnum}.X0') == 0, "different starting point");
    for p = [5, 9, 10]
        assert((P{p, probnum}.y - M{p, probnum}.y) / norm(P{p, probnum}.y) <= 1e-15, ['different y for X0 for prob: ' int2str(probnum) ', ' char(probtypes(p))]);
        assert(all(P{p, probnum}.F - M{p, probnum}.F' <= 1e-14), ['different fvec for X0 for prob: ' int2str(probnum) ', ' char(probtypes(p))]);
        assert(norm(P{p, probnum}.J - M{p, probnum}.J) <= 1e-14, ['different J for X0 for prob: ' int2str(probnum) ', ' char(probtypes(p))]);
    end
end
