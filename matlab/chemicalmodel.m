
%% initial conditions

N = 512;
L = 60;
developmentcost = randfield(N,N,5);
u = randfield(N,N,25);

%% simulation params
dx = N/L;
dt = 0.000002;
% L2 = [0,1,0; 1,-4,1; 0,1,0]/dx/dx;
% L4 = [0,0,1,0,0; 0,2,-8,2,0; 1,-8,20,-8,1; 0,2,-8,2,0; 0,0,1,0,0]/dx/dx/dx/dx;
[L2,L4] = diffmat(N,N,L,L);
 
%% Run
i=0;
while true
    u = update(u,L2,L4,dt);
    if mod(i,1000)
        imagesc(u);
        set(gca,'Ydir','normal');
        colormap(bone);
        drawnow
    end
    i=i+1;
end

%% Simulation Functions
function u = update(u,L2,L4,dt)

u1 = RHS(u,             L2,L4);
u2 = RHS(u + dt/2 * u1, L2,L4);
u3 = RHS(u + dt/2 * u2, L2,L4);
u4 = RHS(u + dt * u3,    L2,L4); 
u = u + (u1 + 2*u2 + 2*u3 + u4) * dt/6;

end

%% helper functions

function du = RHS(u,L2,L4)

% du = 2*u-imfilter(u,L4,'circular')-2*imfilter(u,L2,'circular')-u.^3;
du = 2*u-ifft2((L2-2*L4).*fft2(u),'symmetric')-u.^3;
du = fft2(du);
du(abs(L2)>500)=0;
du = ifft2(du,'symmetric');

end

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

function [D2,D4] = diffmat(Nx,Ny,Lx,Ly)

% generate grid data
x = (0:Nx-1)./Nx*Lx-Lx/2; % define x vector corresponding to physical coordinates
y = (0:Ny-1)./Ny*Ly-Ly/2; % define y vector corresponding to physical coordinates

% define frequency grid
if(mod(Nx,2)==0)
   qx = [0:(Nx/2)-1 (Nx/2) -(Nx/2)+1:-1]'*(2*pi/Lx);
else
    qx = [0:ceil(Nx/2)-1 -ceil(Nx/2)+1:-1]'*(2*pi/Lx);
end
if(mod(Ny,2)==0)
    qy = [0:(Ny/2)-1 (Ny/2) -(Ny/2)+1:-1]'*(2*pi/Ly);
else
    qy = [0:ceil(Ny/2)-1 -ceil(Ny/2)+1:-1]'*(2*pi/Ly);
end

% Stored Parameters (Precalculated for speedup. These will be used often)
Dx = 1i*repmat(qx',Ny,1);
Dy = 1i*repmat(qy,1,Nx);
Dx = fftshift(Dx);
Dy = fftshift(Dy);
Dy(1,:) = 0;
Dx(:,1) = 0;
Dx = ifftshift(Dx);
Dy = ifftshift(Dy);

D2 = Dx.^2 + Dy.^2;	% Laplacian
D4 = D2.^2;

end
