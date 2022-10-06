% This file covers the BenDFO m files
%   calfun.m dfoxs.m jacobian.m dfovec.m g_dfovec_1d.m
% Updated versions may be found at github.com/POptUS/BenDFO
%
% This file has 100% code coverage for m/*.m with the exception of the
% internal adimat function  g_dfovec_1d>adimat_g_pow_left.

addpath('../data'); % location of dfo.dat
load dfo.dat;
numprobs = size(dfo, 1);

probnames = {'smooth', 'absnormal', 'absuniform', 'reluniform', ...
    'relnormal', 'abswild', 'nondiff', 'relwild',  'wild3', 'noisy3'};

% Global variables needed
global BenDFO

% Initialize the results cell and fvals array
Results = cell(length(probnames), numprobs + 2);
BenDFO.fvals = zeros(1, 4 * numprobs + 2);
BenDFO.np = 0;

for np = 1:numprobs
    BenDFO.nprob = dfo(np, 1);
    BenDFO.n = dfo(np, 2);
    BenDFO.m = dfo(np, 3);
    BenDFO.factor_power = dfo(np, 4);
    BenDFO.sigma = 1e-2;

    % Obtain starting vector
    Xs = dfoxs(BenDFO.n, BenDFO.nprob, 10^BenDFO.factor_power);

    % Loop over the 10 problem types
    for p = 1:10
        BenDFO.probtype = probnames{p};
        BenDFO.nfev = 0;
        BenDFO.np = BenDFO.np + 1;
        switch BenDFO.probtype
            case 'smooth'
                [f, fv, G] = calfun(Xs, BenDFO, 'smooth');
                Results{p, np}.f = f;
                Results{p, np}.Xs = Xs;
                Results{p, np}.Fv = fv;
                Results{p, np}.G = G;
            otherwise % skip gradients for the other problems
                rand('state', 1); % For reproducibility of 'noisy3'
                [f, fv] = calfun(Xs);
                Results{p, np}.f = f;
                Results{p, np}.Fv = fv;
                Results{p, np}.Xs = Xs;
        end
    end
end

% Test problem 5 with x(1)>0 and x(1)=0
np = 9;
BenDFO.nprob = dfo(np, 1);
BenDFO.n = dfo(np, 2);
BenDFO.m = dfo(np, 3);
BenDFO.factor_power = dfo(np, 4);
p = 1;
BenDFO.probtype = 'smooth';
% Obtain original starting vector
Xs = dfoxs(BenDFO.n, BenDFO.nprob, 10^BenDFO.factor_power);
np = numprobs;
for xv = [1 0]
    Xs(1) = xv;
    Xs(2) = 1; % Avoiding a NaN from adimat
    BenDFO.nfev = 0;
    BenDFO.np = BenDFO.np + 1;
    np = np + 1;
    [f, fv, G] = calfun(Xs);
    Results{p, np}.f = f;
    Results{p, np}.Xs = Xs;
    Results{p, np}.Fv = fv;
    Results{p, np}.G = G;
end

% The following file can be compared with ../data/testout.dat
%   diff testout.dat ../data/testout.dat

fid = fopen('testout.dat', 'w');
p = 1;
for np = 1:numprobs + 2
    fprintf(fid, '%3i  %8s  %3i  %3i  %6.5e  %6.5e  %6.5e  %6.5e \n', np, ...
            probnames{p}, size(Results{p, np}.Xs, 1), ...
            size(Results{p, np}.Fv, 1), Results{p, np}.f, ...
            abs(sum(sin(Results{p, np}.Fv))), norm(Results{p, np}.G, 2), ...
            Results{p, np}.G' * Results{p, np}.Xs);
end
for p = 2:4
    for np = 1:numprobs
        fprintf(fid, '%3i  %8s  %3i  %3i  %6.5e  %6.5e \n', np, ...
                probnames{p}, size(Results{p, np}.Xs, 1), ...
                size(Results{p, np}.Fv, 1), ...
                Results{p, np}.f, abs(sum(sin(Results{p, np}.Fv))));
    end
end
fprintf(fid, '\n');
fclose(fid);

% Confirm that problem dimension transposition is caught
try
    [f, fv] = calfun(Xs');
catch
    disp('Caught: Input x is not of size n by 1');
end
