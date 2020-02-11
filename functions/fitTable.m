function [CLCoeffs,CDCoeffs] = fitTable(table,alphaRange)

data = [table.alpha(:) table.CL(:) table.CD(:)];

rows = and(table.alpha<=alphaRange(2),table.alpha>=alphaRange(1));

data = data(rows,:,:);

CLCoeffs = polyfit(data(:,1),data(:,2),1);
CDCoeffs = polyfit(data(:,1),data(:,3),2);
% Note coefficients are descending powers (last coefficient is constant
% offset)
end