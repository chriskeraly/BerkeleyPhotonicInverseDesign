

function[phi] = reinit(phi)

flags=(-1).*(phi<0)+1.*(phi>=0);
%phi=abs(phi);
phi=signedDistold(phi,flags);

