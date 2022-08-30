
%% simulation params
params.Dcost = 2; % cost of decentralized construction;
params.Ccost = 1; % cost of decentralized construction;
params.annualvar = 0; % cost of decentralized unit per tic; 
params.cap = 1; % output of decentralized unit per pixel; 
params.neighborbonus = 0.1;
params.neighborsum = [1 1 1; 1 0 1; 1 1 1];

%% initial conditions

N = 100;
M = 50;
demand = randfield(N,N,M,5);
utilities = zeros(N,N);
[~,ii] = max(demand(:,:,1),[],'all','linear');
utilities(ii) = 2;

%% Run

for i = 1:M
[utilities,cost] = update(utilities,demand(:,:,i),params);
imagesc(utilities);
set(gca,'Ydir','normal');
hold on
contour(demand(:,:,i),'edgecolor','r');
pause(1);
end

%% TODO 
% build plants (?build cost?)
% power plant capacity
% units can downgrade
% cost of construction (landscape?)

%% Simulation Functions
function [e,cost] = update(e,d,params)

% figure out how many centralized neighbors each tile has.
num_neighbors = imfilter(double(e>0),params.neighborsum);

% compute if tiles should upgrade
% dont upgrade if youre a power plant AND
% only can upgrade if you have a neighbor thats a power plant or a
% centralized boi AND
% its cost effective
shouldUpgrade = ( e < 2 ) & ( num_neighbors > 0 ) & (d> ( params.Ccost*params.cap ) / (  params.Dcost-params.annualvar ) );

% compute centralized cost per tile
e(shouldUpgrade) = 1;

cost = sum(params.Dcost*d(e==0),'all')+sum(params.annualvar*d(e>1),'all');

end

%% helper functions

function f = randfield(Nx,Ny,Nt,cutoff)
    
    % return a smooth, randomly varying field between zero and 1. As cutoff
    % increases, the structure of the field becomes smaller wavelength. 
    f = fftshift(fftn(rand(Ny,Nx,Nt)));
    [Nx,Ny,Nt] = ndgrid(((1:Ny)-Ny/2),((1:Nx)-Nx/2),((1:Nt)-Nt/2));
    f(sqrt(Nx.^2+Ny.^2+Nt.^2)>cutoff)=0;
    f = ifftn(ifftshift(f),'symmetric');
    f = f-min(f,[],'all');
    f = f/max(f,[],'all');
    
end
