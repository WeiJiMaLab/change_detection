function Yi = qinterp1Lin(spacing,Y,xi)

% Performs fast linear interpolation compared to interp1
%
% S. Keshvari, 2012

xi = spacing*xi+1;
fl = floor(xi);
Yi = (Y(ceil(xi))-Y(fl)).*(xi-fl)+Y(fl);
