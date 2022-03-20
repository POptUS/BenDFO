% This file covers the BenDFO m files
%   calfun.m dfoxs.m jacobian.m dfovec.m g_dfovec_1d.m
% Updated versions may be found at github.com/POptUS/BenDFO

addpath('../data') % location of dfo.dat
load dfo.dat
numprobs = size(dfo,1);

probnames = {'smooth','nondiff','wild3','noisy3'};

% Global variables needed
global m np nprob probtype fvals nfev

% Initialize the results cell and fvals array
Results = cell(length(probnames),numprobs);
fvals = zeros(1,numprobs);

for np = 1:numprobs
    nprob = dfo(np,1);
    n = dfo(np,2);
    m = dfo(np,3);
    factor_power = dfo(np,4);
    
    % Obtain starting vector
    Xs = dfoxs(n,nprob,10^factor_power);
    
    % Loop over the 4 problem types
    for p = 1:4
        probtype = probnames{p};
        nfev = 0;
        switch probtype
            case 'smooth'
                [f,fv,G] = calfun(Xs);
                Results{p,np}.f = f;
                Results{p,np}.fv = fv;
                Results{p,np}.G = G;
            otherwise % skip gradients for the other problems
                rand('state',1) % For reproducibility of 'noisy3'
                [f] = calfun(Xs);
                Results{p,np}.f = f;
        end
    end
end
%    save('testoutput', 'Results','fvals');
