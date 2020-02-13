function Q = wyptSelectionMatrix(stateIdx,nStates,wyptIdx,nSteps)
wyptIdx = wyptIdx(:)'; % Make this a row
% Calculate the number of waypoints
nWyPts = numel(wyptIdx);
% Create zeros sub-matrix
q = zeros(1,nStates);
% Create cell array for block matrix
Q = cell(nSteps,nSteps);
% Set all block elements to zero sub matrics
Q(:) = {q};
% Set specified element of sub matrix to 1
q(stateIdx) = 1;
% Set block diagonal elements of block matrix to new sub matrix
Q(sub2ind(size(Q),wyptIdx,wyptIdx)) = {q};
% Convert cell block matrix to matrix
Q = cell2mat(Q);
% Delete all empty rows
Q = Q(sum(Q,2)~=0,:);
end