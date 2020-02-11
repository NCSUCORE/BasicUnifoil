function aeroTable =  buildAirfoilTable(folder,AR,varargin)
p = inputParser;
addRequired(p,'folder');
addRequired(p,'AR');
addOptional(p,'OE',0.8)
parse(p,folder,AR,varargin{:}); 

oswaldEfficiency = p.Results.OE;

folder = lower(folder);
files = dir(fullfile(pwd,folder));
if ispc
    [~,~,data]=xlsread(fullfile(pwd,folder,files(3).name));
    aeroTable.fileName = files(3).name;
    aeroTable.Re = data{4,2};
    jj=1;
    while ~strcmpi(data{jj,1},'alpha')
        jj=jj+1;
    end
    data = cell2mat(data(jj+1:end,1:3));
else
    files = files(3:end);
    for ii = 1:length(files)
        if ~strcmpi(files(ii).name(1),'.')
            data=csvread(fullfile(pwd,folder,files(ii).name),11,0);
            fid = fopen(fullfile(pwd,folder,files(ii).name));
            Re = textscan(fid, '%s','delimiter', '\n');
            fCLose(fid);
            Re=strsplit(Re{1}{4},',');
            aeroTable.Re = str2double(Re{2});
        end
    end
end
aeroTable.alpha = data(:,1)*(pi/180);
aeroTable.CL    = data(:,2);
aeroTable.CL0   = aeroTable.CL(data(:,3)==min(data(:,3)));
aeroTable.CD    = min(data(:,3))+((aeroTable.CL-aeroTable.CL0).^2)./(pi*oswaldEfficiency*AR);
end