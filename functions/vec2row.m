function [Xrow] = vec2row(X)

% Makes sure vector X will become a row vector
%
% Author: Yui Takeshita
% Carnegie Institution for Science
% Created: 11/16/2015
% Last modified: 11/16/2015

Xrow = X(:)';
