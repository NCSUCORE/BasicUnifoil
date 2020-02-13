function tscPath = reparameterize(tscTime)
% Method to overwrite "time" field in all signals with the path variable
tscPath = tscTime;
props = properties(tscTime);

% Remove samples corresponding to non monotonically increasing path steps
nonMonInd = find(diff(tscTime.pathVar.Data)<=0);
if any(nonMonInd)
    for ii = 1:numel(props)
        if isa(tscPath.(props{ii}),'timesignal')
            tscPath.(props{ii}) = tscPath.(props{ii}).delsample('Index',nonMonInd);
        end
    end
end


% Check for monotonicity of path variable
if any(diff(tscTime.pathVar.Data)==0)
   error('Path variable is momentarily constant') 
end
if any(diff(tscTime.pathVar.Data)<0)
   error('Path variable is momentarily decreasing') 
end

for ii = 1:numel(props)
    if isa(tscPath.(props{ii}),'timesignal')
        tscPath.(props{ii}).Time = tscTime.pathVar.Data;
    end
end
end