close all

N = 32+1;

% four probabilities make up the dynamics
%p_join      = randfield(N,N,1,5)/10; % probability that a de-central node becomes central if its neighbor is central. Compounds linearly with how many neighbors are central
%p_join(1:10,1:10)=0;
p_join      = zeros(N,N)+1;

p_outtage   = 0.005; % probability that a centralized nodes undergoes an outtage
p_recover   = 1.00; % probability that an outtage is fixed. this leaves the tile decentral

%%%%%%%%%%%% create initial state %%%%%%%%%%%%
%%
% make state of zeros
s = zeros(N,N);
s(round(N/2),round(N/2)) = 1;
e = zeros(N^2,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Population and landscape infrastructure %%%%%%%%%%%%

population = zeros(N,N);
population(round(N/2),round(N/2)) = .5;
population(round(N/3),round(N/3)) = 0.3;
D = 0.1;
dt = 0.001;

lands = randfield(N,N,1,5);
[landgradsx, landgradsy] = gradient(lands)
steepness = sqrt(landgradsy.^2+landgradsx.^2);

% Wishlist: find a way to get the gradient of the landscape one
%           so I can use that as my metric for increasing or 
%           decreasing probs. 
% Current approach: use just the plain value as the metric. Less
%           good because you can have a plateau high up that's 
%           still easy to build on. It's the slopes that really 
%           make things hard. 
%           Also, I will be assuming that landscape < 0 is water
%           and thus unbuildable. Idk how to get around that 
%           with only one plant, but that's what I'll try first.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Init for histogram data
binEdges = linspace(0,ceil(sqrt(2)*N),21);
histData = zeros(1, 20);

%%
figure;
resizefigure(2,2);

for k = 1:1000*N^2
    
    %population = rk4(population,D,dt);

    %p_join = population;

    [s,e] = update(s,e,p_join,p_outtage,p_recover);
    
    avgDistance = averageDistances(s, e);
    %fprintf('Average distance from end pixels to origin: %.2f\n', avgDistance);

    % Update the histogram data
    histData = histData + histcounts(avgDistance, binEdges);

    if mod(k,1)==0
        plotstate(s,e,p_join,histEdges,histData);
    end
    
end

function v =mod1(v,N)
v = mod(v-1,N)+1;
end

function [s,e,outtagesize] = update(s,e,p1,p2,p3)

[Ny,Nx] = size(s);
Nn = size(e,2);
neighborsum = conv2(s,[0 1 0; 1 0 1; 0 1 0],'same'); % USE FOR SOLID BOUNDARY
% neighborsum = round(conv_fft2(s,[0 1 0; 1 0 1; 0 1 0],'wrap')); % USE FOR PERIODIC BOUNDARIES
outtagesize = 0;

k = randi(Ny*Nx);
[i,j] = ind2sub([Ny,Nx],k);
sij = s(i,j);

% if current pixel is OUT
if sij == -1
    
    % restore to decentralization with probability p4
    if rand < p3
        s(i,j) = 0;
    end
    
% if current pixel is DECENTRAL
elseif sij == 0
    
    % become central with probability n*p2, where n is the number of neighboring central nodes minus the number of neighboring out nodes
    if rand < neighborsum(i,j)*p1(i,j)
        
        % set pixel to be central
        s(i,j) = 1;
        % we now need to assign one of the neighboring centralized nodes to be the parent. Lets do this in a random fashion
        for n = randperm(Nn)
            if n==1
                ny = mod1(i+1,Ny);
                if s(ny,j)==1
                    kn = sub2ind([Ny,Nx],ny,j);
                    e(kn,n) = k;
                    break % if one of these conditions is set, break, since we only want to set ONE neighbor as the parent
                end
            elseif n==2
                nx = mod1(j-1,Nx);
                if s(i,nx)==1
                    kn = sub2ind([Ny,Nx],i,nx);
                    e(kn,n) = k;
                    break
                end
            elseif n==3
                ny = mod1(i-1,Ny);
                if s(ny,j)==1
                    kn = sub2ind([Ny,Nx],ny,j);
                    e(kn,n) = k;
                    break
                end
            elseif n==4
                nx = mod1(j+1,Nx);
                if s(i,nx)==1
                    kn = sub2ind([Ny,Nx],i,nx);
                    e(kn,n) = k;
                    break
                end
            end
        end
    end
        
% if current pixel is CENTRAL
elseif sij == 1
    % experience an outtage with probability p3
    if rand < p2
        n = sum(s(:)==-1);
        e(e==k) = 0;
        [s,e] = propogateouttage(s,e,k);
         outtagesize = sum(s(:)==-1)-n;
    end
end

if s(round(end/2),round(end/2)) == 0
    s(round(end/2),round(end/2)) = 1;
end

end

function [state,edges] = propogateouttage(state,edges,k)

    % pixel k is experiencing an outtage. find coordinates of k in the NxN grid
    [i,j] = ind2sub(size(state),k);
    
    % figure out how many "children" each pixel can have
    num_edges = size(edges,2);

    % for each edge, if it has a child (edge is non-zero). When non-zero, the edges array stores the k-index of the child
    for e = 1:num_edges
        
        if edges(k,e)>0
            
            % check if that child has children and run the outtage code for the child first, recursively 
            [state,edges] = propogateouttage(state,edges,edges(k,e));
            % set this edge to zero so we know this child is no longer a child anymore
            edges(k,e) = 0;
            
        end

        % set the pixel value here to -1 since it experienced an outtage
        state(i,j) = -1;
        
    end
    
end

% Field maker
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

function P = rk4(P,D,h)
 
% P : population
% D: diffusion coefficient
% h : timestep
 
halfh = h / 2;
h6 = h/6;
 
p1 = dPdt(P,D);
p2 = dPdt(P + halfh * p1,D); 
p3 = dPdt(P + halfh * p2,D); 
p4 = dPdt(P + h * p3,D); 
P = P + (p1 + 2*p2 + 2*p3 + p4) * h6;
 
end
function dP = dPdt(P,D)
 
lapP = del2(P);
 
dP = P.*(1-P)+D*lapP;
dP([1 end],:)=0;
dP(:,[1 end])=0;

end

function plotstate(s,e,p_join,histEdges,histData)

% imagesc(s)
% hold on
% % caxis([-1 1])
% 
% % colormap(flip(gray));
% caxis([-1 1]);
% drawnow
% hold off

[Ny,Nx] = size(s);

subplot(2,2,1)
imagesc(p_join);
% set(gca,'YDir','normal')
hold on
[yyy,xxx] = ind2sub(size(s),find(s==-1));
scatter(xxx,yyy,150,'rs','filled')
reccursiveplot(s,e,round(numel(s)/2));
hold off
axis equal
daspect([1 1 1]);
caxis([-1 1]);
% colormap(flip(rwb))
ylim([1 Ny]);
xlim([1 Nx]);

subplot(2,2,2)
p(1) = sum(s(:)== 0);
p(2) = sum(s(:)== 1);
p(3) = sum(s(:)==-1);
P = p/sum(p);
bar(P);
xticks(1:3);
xticklabels({'decentral','central','out'});
text((1:numel(P))-.2,P+.02,num2cell(P))
title('Average (in time) percent of grid that is...')
text((1:numel(P))-.2,P+.02,num2cell(P))
ylim([0 1]);

subplot(2,2,3)
si = s;
si(si==-1)=0;
si = bwmorph(si,'remove');
% ii = sort(unique(e(:)));
% ii(1)=[];
% ii(sum(e(ii,:),2)>0)=[];
% si = s*0;
% si(ii)=1;
[L,n] = bwlabel(si);
h = hist(L(:),0:n);
[~,ii] = max(h(2:end));
imagesc(L==ii)

subplot(2,2,4)
histogram('BinEdges', histEdges, 'BinCounts', histData); %, 'DisplayStyle', 'stairs'
title('Average Distances Histogram');
xlabel('Average Distance');
ylabel('Frequency');

drawnow

end
function reccursiveplot(s,e,k)

[y1,x1] = ind2sub(size(s),k);

for i = 1:size(e,2)
   if e(k,i)
       [y2,x2] = ind2sub(size(s),e(k,i));
       plot([x1 x2],[y1 y2],'k');
       reccursiveplot(s,e,e(k,i))
   end
end

end

function newmap = rwb(m)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG,
%   COLORMAP, RGBPLOT.
if nargin < 1
    m = size(get(gcf,'colormap'),1);
end

bottom = [0 0 0.5];
botmiddle = [0 0.5 1];
middle = [1 1 1];
topmiddle = [1 0 0];
top = [0.5 0 0];

% bottom = [0 0 0.25];
% botmiddle = [0 0 0.5];
% middle = [1 1 1];
% topmiddle = [0.5 0 0];
% top = [0.25 0 0];

% Find middle
lims = get(gca, 'CLim');
% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG,
%   COLORMAP, RGBPLOT.
if nargin < 1
    m = size(get(gcf,'colormap'),1);
end

bottom = [0 0 0.5];
botmiddle = [0 0.5 1];
middle = [1 1 1];
topmiddle = [1 0 0];
top = [0.5 0 0];

% bottom = [0 0 0.25];
% botmiddle = [0 0 0.5];
% middle = [1 1 1];
% topmiddle = [0.5 0 0];
% top = [0.25 0 0];

% Find middle
lims = get(gca, 'CLim');
% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end



% m = round(m/4);
% newmap = newmap(m+1:end-m,:);
%
% m = 64;
% new = [bottom; botmiddle; middle; topmiddle; top];
% % x = 1:m;
%
% oldsteps = linspace(0, 1, 5);
% newsteps = linspace(0, 1, m);
% newmap = zeros(m, 3);
%
% for i=1:3
%     % Interpolate over RGB spaces of colormap
%     newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
% end
%
% % set(gcf, 'colormap', newmap), colorbar
end

function avgDist = averageDistances(s, e)

    % SUPREME CODE 

    [Ny, Nx] = size(s);
    
    % Find the origin pixel
    origin = [round(Ny/2), round(Nx/2)];
    
    % Initialize variables to store the total distance and count
    totalDist = 0;
    count = 0;
    
    % Iterate over centralized pixels to find the end of branches
    for k = 1:Ny*Nx
        [i, j] = ind2sub([Ny, Nx], k);
        
        if s(i, j) == 1
            % This is a centralized pixel
            isEndOfBranch = true;
            
            % Check if it has children
            for n = 1:size(e, 2)
                childIndex = e(k, n);
                
                if childIndex > 0
                    isEndOfBranch = false;
                    break;
                end
            end
            
            if isEndOfBranch
                % Calculate the Euclidean distance from the end pixel to the origin
                dist = sqrt((i - origin(1))^2 + (j - origin(2))^2);
                totalDist = totalDist + dist;
                count = count + 1;
            end
        end
    end
    
    if count > 0
        % Calculate the average distance
        avgDist = totalDist / count;
    else
        avgDist = NaN; % No end pixels found
    end
    fprintf('Average distance from end pixels to origin: %.8f\n', avgDist);
end