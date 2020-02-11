function psi = heading(pos,targ)
% https://en.wikipedia.org/wiki/Great-circle_navigation
% http://mathworld.wolfram.com/SphericalCoordinates.html
% phi is latitude
% lam is longitude
phi1 = (pi/2)-acos(pos(3));
lam1 = atan2(pos(2),pos(1));
phi2 = (pi/2)-acos(targ(3));
lam2 = atan2(targ(2),targ(1));
lam12 = lam2 - lam1;
psi = (pi/2)-atan2(cos(phi1)*sin(lam12),-cos(phi2)*sin(phi1)+sin(phi2)*cos(phi1)*cos(lam12));
end