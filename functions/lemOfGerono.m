function gndPos = lemOfGerono(pathVar,geomParams)
%LEMOFGERONO Mercator projection of the lemniscate of Gerono onto a sphere
%   INPUTS
%   pathVar = Normalized path variabl (0-1) describing position along the
%   path
%   geomParams(1) = A0, total azimuth sweep angle in degrees
%   geomParams(2) = Z0, total elevation sweep angle in degrees
%   geomParams(3) = A1, mean course azimuth angle in degrees
%   geomParams(4) = Z1, mean course elevation angle in degrees
%   geomParams(5) = R, radius of sphere
%   OUTPUTS
%   gndPos = Nx3 matrix, where N is the number of elements of pathVar
%   Variable names here were chosen to match the notebook
%   lemOfGeronoTanVec.nb

A0 = geomParams(1);
E0 = geomParams(2);
A1 = geomParams(3);
E1 = geomParams(4);
R  = geomParams(5);
freqMult = geomParams(6);

pathVar = pathVar(:);

% Calculate azimuth, elevation, and zenith
a =  (A0/2)*sin(          2*pi*pathVar) + A1;
if freqMult == 1 % Ellipse
    e = -(E0/2)*cos(2*pi*pathVar) + E1;
else % Fig 8
    e = -(E0/2)*sin(4*pi*pathVar) + E1;
end
z =  (pi/2) - e;

% plot(a,e)
% Convert sphereical to cartesian
% http://mathworld.wolfram.com/SphericalCoordinates.html
gndPos = nan(numel(pathVar),3);

gndPos(:,1) = R.*cos(a).*sin(z);
gndPos(:,2) = R.*sin(a).*sin(z);
gndPos(:,3) = R.*cos(z);

end

