% This is a test to compare the outputs form matlab and python BenDFO versions.
%
% First run
% `python evaluate_all_calfuns_points.py`
% in the `py/` directory. This evaluates `calfun.py` returning 4 outputs
% `(y, fvec, G, J)` for 3 different problem types and at 3 different points,
% and then saving output to a .mat file
%
% Do the same in matlab by running `evaluate_all_calfuns_points` in the `m/`
% directory.
%
% Finally, run this script to compare the differences.

M = load("../m/fvec_and_gradients_at_starting_values_matlab.mat");
M = M.Results;

P1 = load("../py/fvec_and_gradients_at_starting_values_python.mat");
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
            out1(i, p, pt) = abs(P{p, i, pt}.y - M{p, i, pt}.y) / abs(P{p, i, pt}.y);
            out2(i, p, pt) = max(abs(P{p, i, pt}.F - M{p, i, pt}.F') ./ abs(P{p, i, pt}.F));
            out3(i, p, pt) = norm(P{p, i, pt}.J - M{p, i, pt}.J) / norm(P{p, i, pt}.J);
            out4(i, p, pt) = norm(P{p, i, pt}.G - M{p, i, pt}.G') / norm(P{p, i, pt}.G);
            assert(out1(i, p, pt) <= 1e-13, sprintf(msg, 'y', 1, ps, xs));
            assert(out2(i, p, pt) <= 1e-12, sprintf(msg, 'fvec', 1, ps, xs));
            assert(out3(i, p, pt) / norm(P{p, i, pt}.J) <= 1e-14, sprintf(msg, 'J', 1, ps, xs));
            assert(out4(i, p, pt) / norm(P{p, i, pt}.G) <= 1e-14, sprintf(msg, 'G', 1, ps, xs));
        end
    end
end

% Print summary to screen
format short g;
disp(max(out1(:, [5 9 10], :)));
disp(max(out2(:, [5 9 10], :)));
disp(max(out3(:, [5 9 10], :)));
disp(max(out4(:, [5 9 10], :)));
