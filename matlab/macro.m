
%% simulation params
params.DcostConstr = 1; % cost of decentralized construction;
params.DcostUtil = 1; % cost of decentralized unit per tic; 
params.DoutputUtil = 1; % output of decentralized unit per pixel; 
params.Dradius = 2; % how many tiles away a decentralized center can serve;
params.CcostConstr = 2; % cost of centralized construction;
params.CcostUtil = 1; % cost of centralized unit per tic; 
params.CoutputUtil = 3; % output of centralized unit per pixel; 
params.Cradius = 10; % how many tiles away a centralized center can serve;
params.UtilPerCapita = 3; % how many MW each person needs per tic
params.FinicialConservatism = .5; % how much spending money for new units is disliked. 
params.Head = nan; % head and tail of linked list;
params.Tail = nan;

%% initial conditions

N = 100;
landscapecost = randfield(N,N,5);
population = randfield(N,N,5);
entities = {};
X=1:N;
[params.Y,params.X] = ndgrid(X,X);

%% Run

[entities,landscape,cost] = update(entities,population,landscapecost,params);

%% Simulation Functions
function [entities,landscape,moneyspent] = update(entities,population,landscape,params)

% compute how much utility is needed per tile
deficit = computeDeficit(entities,population,params);
moneyspent = 0;

% for equity, update cells in random order
for i = 1:randperm(numel(entities))
    
    e = entities{i};
    % do something that updates entity
    
end

% look to add new centers
deficit(deficit<0)=0; % get rid of surpluses

% calculate cost of adding a decentralized center
F = double(fspecial('disk',params.Dradius-1)>0);
benefit = imfilter(max( params.DoutputUtil*population./sum(population)-deficit ,0),F);
cost = landscape+params.FinancialConservatism*params.DcostConstr;
[bD,iD] = max(benefit-cost,'all');
cD = cost(iD);

% calculate cost of adding a centralized center
F = double(fspecial('disk',params.Cradius-1)>0);
benefit = imfilter(min(population*params.UtilPerCapita,params.CoutputUtil),F);
cost = landscape+params.FinancialConservatism*params.CcostConstr;
[bC,iC] = max(benefit-cost,'all');
cC = cost(iC);

if bC||bD>0
    
    e = struct;
    
    if bC>bD
        
        % create centralized unit
        e.x = params.X(iC);
        e.y = params.Y(iC);
        e.distribution = double(fspecial('disk',params.Cradius-1)>0);
        moneyspent = moneyspent + cC;
        
    else
        
        % create centralized unit
        e.x = params.X(iD);
        e.y = params.Y(iD);
        moneyspent = moneyspent + cD;
        
    end
    
    % add entity to list
    n = numel(entities);
    if n==0
        params.Head = 1;
    else
        entities{params.Tail}.next = e;
    end
    e.former = params.Tail;
    params.Tail = n+1;
    entities{n+1} = e;
    
end

end

function d = computeDeficit(entities,population,params)

x = params.X;
y = params.Y;
d = population*params.UtilPerCapita;

% construct grids of what decentralized and centralized units are supplying
% to the domain
p = params.pointer;
while ~isnan(p)
    e = entities{p};
    switch entity.type
        case 'D'
            
            isServiced = sqrt((x-e.x).^2+(y-e.y).^2)<params.Dradius;
            d(isServiced) = d(isServiced) - params.DoutputUtil*numel(isServiced)*population(isServiced)./sum(population(isServiced));

        case 'C'
            
            isServiced = e.distribution>0;
            d(isServiced) = d(isServiced) - params.CoutputUtil*numel(isServiced)*population(isServiced)./sum(population(isServiced));
            
%             dtree = 0*d;
%             isServiced = sqrt((x-e.x).^2+(y-e.y).^2)<params.Cradius;
%             dtree(isServiced) = dtree(isServiced) + params.CoutputUtil;
%             
%             ch = e.children;
%             n = numel(ch);
%             p = entity.nextEntity;
%             
%             for c = 1:n
%                 e = entities{ch(c)};
%                 isServiced = sqrt((x-e.x).^2+(y-e.y).^2)<params.Cradius;
%                 dtree(isServiced) = dtree(isServiced) + params.CoutputUtil;
%             end
%             
%             isServiced = dtree>0;
%             d(isServiced) = dtree(isServiced)/sum(isServiced);
            
        otherwise
            error('entity type not recognized');
    end

    p = entity.next;
    
end

end

%% helper functions

function f = randfield(Nx,Ny,cutoff)
    
    % return a smooth, randomly varying field between zero and 1. As cutoff
    % increases, the structure of the field becomes smaller wavelength. 
    f = fftshift(fft2(rand(Ny,Nx)));
    [Nx,Ny] = ndgrid(((1:Ny)-Ny/2),((1:Nx)-Nx/2));
    f(sqrt(Nx.^2+Ny.^2)>cutoff)=0;
    f = ifft2(ifftshift(f),'symmetric');
    f = f-min(f,[],'all');
    f = f/max(f,[],'all');
    
end
