function Yi = qinterp1(spacing,Y,xi)

% Performs fast function lookup
%
% S. Keshvari, 2012

xi = spacing*xi+1;
xi = round(xi);
Yi = Y(xi);
