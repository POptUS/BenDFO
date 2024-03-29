function [J, fvec] = jacobian(m, n, x, nprob)
%     This subroutine computes the Jacobian of the nonlinear equations
%     defining the benchmark problems in
%
%     Benchmarking Derivative-Free Optimization Algorithms
%     Jorge J. More' and Stefan M. Wild
%     SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.
%
%     The latest version of this subroutine is always available at
%          https://github.com/POptUS/BenDFO/
%
%     The dependencies of this function are based on executing the AD
%     software adimat on a modified form of the dfovec function from
%     http://www.mcs.anl.gov/~more/dfo/
%     See the instructions of dfovec.m for additional details on these
%     nonlinear benchmark problems (and appropriate values of m and n).
%
%   J = jacobian(m,n,x,nprob)
%
%       J is an output array of size m-by-n, with J(i,j) denoting the
%         derivative (evaluated at x) of the ith equation with respect to
%         the jth variable.
%       fvec returns the usual dfovec
%       m and n are positive integer input variables. n must not
%         exceed m.
%       x is an input array of length n.
%       nprob is a positive integer input variable which defines the
%         number of the problem. nprob must not exceed 22.
%
%     Argonne National Laboratory
%     Stefan Wild. July 2014.

% addpath('adimat_out');
% This directory contains the dependency g_dfovec_1d.m

% Initialization for adimat objects
t = 0;
g_t = 1;
g_x = zeros(size(x));

J = zeros(m, n);
for ind = 1:n % Do one coordinate direction at a time:
    [g_fvec, fvec] = g_dfovec_1d(g_t, t, ind, m, n, g_x, x, nprob);
    J(:, ind) = g_fvec;
end

end
