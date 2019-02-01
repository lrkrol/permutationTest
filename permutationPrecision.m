% numPermutations = permutationPrecision(precision, alpha, confidenceInterval)
%
%       Returns an estimate of the number of permutations necessary for a
%       permutation test to obtain a p-value with a given precision.
%
%       Based on steffen's answer at
%       https://stats.stackexchange.com/questions/80025/#80879
%
% In:
%       precision - the required precision, i.e., with the actual p being
%                   within p +/- 3*precision
%       alpha - the desired significance level (0-1)
%       confidenceInterval - 1, 2, or 3 corresponding to the size of the
%                            confidence interval, where
%                            | 1 | ~68 % |
%                            | 2 | ~95 % |
%                            | 3 | ~99 % |
%
% Out:
%       numPermutations - estimated minimum number of required permutations

% 2019-02-01 First version

% No rights reserved.


function numPermutations = permutationPrecision(precision, alpha, confidenceInterval)

numPermutations = round(confidenceInterval^2 * (alpha * (1-alpha)) / precision^2);

end