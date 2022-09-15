close all

N = 32+1;

% four probabilities make up the dynamics
p_join      = 0.25; % probability that a de-central node becomes central if its neighbor is central. Compounds linearly with how many neighbors are central
p_outtage   = 0.01; % probability that a centralized nodes undergoes an outtage
p_recover   = 1.00; % probability that an outtage is fixed. this leaves the tile decentral

%%%%%%%%%%%% create initial state %%%%%%%%%%%%

% make state of zeros
s = zeros(N,N);
s(round(N/2),round(N/2)) = 1;
e = zeros(N^2,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
resizefigure(2,1);
p = zeros(3,1);
for k = 1:1000*N^2
    
    [s,e] = update(s,e,p_join,p_outtage,p_recover);
    
    if mod(k,N^2)==1
        subplot(1,2,1)
        imagesc(s);
        axis equal
        %    xlim([1 N]);
        %    ylim([1 N]);
        daspect([1 1 1]);
        caxis([-1 1]);
        colormap(flip(rwb))
        ylim([1 N]);
        xlim([1 N]);
        subplot(1,2,2)
        p(1) = p(1) + sum(s(:)== 0);
        p(2) = p(2) + sum(s(:)== 1);
        p(3) = p(3) + sum(s(:)==-1);
        P = p/sum(p);
        bar(P);
        text((1:numel(P))-.2,P+.02,num2cell(P))
        ylim([0 1]);
        drawnow
    end
    
end

function v =mod1(v,N)
v = mod(v-1,N)+1;
end

function [s,e] = update(s,e,p1,p2,p3)

[Ny,Nx] = size(s);
Nn = size(e,2);
% neighborsum = imfilter(s,[0 1 0; 1 0 1; 0 1 0],'circular');
neighborsum = conv2(s,[0 1 0; 1 0 1; 0 1 0],'same');

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
    if rand < neighborsum(i,j)*p1
        
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
        [s,e] = propogateouttage(s,e,k);
    end
end

s(round(Ny/2),round(Nx/2)) = 1;

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