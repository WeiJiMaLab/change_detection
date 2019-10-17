function [output] = myBessel(X,spacing,lookupY)
% Modified bessel function of the first kind, nu == 0
% Reshape X into one column with columns concatenated for qinterp1
tempX = reshape(X,[],1);

% use quick lookup
L5 = qinterp1(spacing,lookupY,tempX);

% reshape for output
output = reshape(L5,size(X));
