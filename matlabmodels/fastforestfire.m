close all

N = 128;

p1 = 0.5;
theta = 1; % theta = p1/p2, where p1 is the probability that a tree grows from soil and p2 is the probability a tree burns. must be an integer
assert(mod(theta,1)==0);
p2 = p1/theta;

%%%%%%%%%%%% compute expected stationary state %%%%%%%%%%%%
     
P_exp = [0; 0; 0;];
% coeff = [
%     (16*p1^4 + 16*p1^5);
%     (-32*p1^3 - 64*p1^4 - 48*p1^5); 
%     (24*p1^2 + 72*p1^3 + 104*p1^4 + 56*p1^5);
%     (-8*p1 - 32*p1^2 - 72*p1^3 - 80*p1^4 - 32*p1^5);
%     (-p1 + 4*p1^2 + 18*p1^3 + 20*p1^4 + 7*p1^5 - p2 - 5*p1*p2 - 10*p1^2*p2 - 10*p1^3*p2 - 5*p1^4*p2 - p1^5*p2)
%     p1*p2 + 4*p1^2*p2 + 6*p1^3*p2 + 4*p1^4*p2 + p1^5*p2;
% ];
% r = roots(coeff);
% P_exp = [];
% for k = 1:numel(r)
%     Pk = [r(k)/p1; 1-r(k)*(1+1/p1); r(k)];
%     if all(Pk<=1) && all(Pk>=0)
%         P_exp(:,end+1) = Pk;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% create initial state %%%%%%%%%%%%

% make state of zeros
m = zeros(N,N);

% % make state corresponding to expected equilibrium
% k = 2;
% m = rand(N,N);
% m(m>(1-P_exp(2,k)))=1;
% m(m<P_exp(3,k))=-1;
% m(abs(m)<.9)=0;

% % make state of ones
% m = ones(N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
resizefigure(2,1);
p = zeros(3,1);
for k = 1:1000*N^2
    
    %% Grow stage. plant a tree theta times, at random
    for t = 1:theta
        i = randi(N);
        j = randi(N);
        if m(i,j)==0
            m(i,j)=1;
        end
    end
    
    %% Burn stage. pick a site and if its a tree burn the related cluster.
    i = randi(N);
    j = randi(N);
    lb = bwlabel( m );
    if m(i,j)==1
        m(lb == lb(i,j))=0;
    end

    %% plotting stage
    if mod(k,N)==0
        subplot(1,2,1)
        imagesc(m);
        axis equal
        %    xlim([1 N]);
        %    ylim([1 N]);
        daspect([1 1 1]);
        caxis([-1 1]);
        colormap(flip(rwb))
        subplot(1,2,2)
        p(1) = p(1) + sum(m(:)== 0);
        p(2) = p(2) + sum(m(:)== 1);
        p(3) = p(3) + sum(m(:)==-1);
        P = p/sum(p);
%         P(1,1) = sum(m(:)== 0)/numel(m);
%         P(2,1) = sum(m(:)== 1)/numel(m);
%         P(3,1) = sum(m(:)==-1)/numel(m);
        bar(cat(2,P,P_exp));
        text((1:numel(P))-.2,P+.02,num2cell(P))
        err = min(sqrt(sum((P-P_exp).^2,1)),[],2)/norm(P);
        title(['$\epsilon = ',sprintf('%.4f',err),'$'],'Interpreter','latex');
        ylim([0 1]);
        drawnow
    end
    
end

function v =mod1(v,N)
v = mod(v-1,N)+1;
end

function m = update(m,p1,p2)

r = rand(size(m));

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