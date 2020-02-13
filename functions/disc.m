function [Ad,Bd] = disc(Ac,Bc,step)
Ad = Ac;
Bd = Bc;

CMat = eye(size(Ac.getdatasamples(1)));
DMat = zeros(size(Bc.getdatasamples(1)));
for ii = 1:numel(Ac.Time)
    sysTmp = c2d(ss(Ac.Data(:,:,ii),Bc.Data(:,:,ii),CMat,DMat),step);
    Ad.Data(:,:,ii) = sysTmp.A;
    Bd.Data(:,:,ii) = sysTmp.B;
end

end