function Gamma = stateSelectionMatrix(indx,nStates,nSteps)
% Matrix to pick off specified state value at every time step
% Create zeros sub-matrix
gamma = zeros(1,nStates);
% Create cell array for block matrix
Gamma = cell(nSteps,nSteps);
% Set all block elements to zero sub matrics
Gamma(:) = {gamma};
% Set specified element of sub matrix to 1
gamma(indx)= 1;
% Set block diagonal elements of block matrix to new sub matrix
Gamma(1:nSteps+1:nSteps^2) = {gamma};
% Convert cell block matrix to matrix
Gamma = cell2mat(Gamma);
end