function idx = cnvrtPathVar2Indx(vars,pathSteps)
% Function to convert path variables (0-1) to indices indicating which path
% step,  function rounds to nearest in range.

% Round to nearest
vars = interp1(pathSteps,pathSteps,vars,'nearest','extrap');
[~,idx] = ismember(vars,pathSteps);
idx = idx(idx~=0);
idx = sort(idx,'ascend');
end