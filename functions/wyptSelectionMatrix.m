function Q = wyptSelectionMatrix(stateIdx,wyptIdx,nStates,nSteps)
wyptIdx = wyptIdx(:)'; % Make this a row
% Calculate the number of waypoints
% nWyPts = numel(wyptIdx);
% Create zeros sub-matrix
q = zeros(nStates,nStates);
% Create cell array for block matrix
Q = cell(nSteps,nSteps);
% Set all block elements to zero sub matrics
Q(:) = {q};
% Set specified element of sub matrix to 1
qDiag = zeros(nStates,1);
qDiag(stateIdx) = 1;
% q(stateIdx) = diag(qDiag);
% Set block diagonal elements of block matrix to new sub matrix
Q(sub2ind(size(Q),wyptIdx,wyptIdx)) = {diag(qDiag)};
% Convert cell block matrix to matrix
Q = cell2mat(Q);
% Eliminate empty rows
Q = Q(sum(Q)~=0,:);
end