close all

if isfile('data\utilityfailurescan.mat')
    
    load('data/utilityfailurescan','plabel','p','clustersize_distribution','outtage_distribution','iii','jjj','p_join','p_outtage','p_recover','N');
    istart = iii;
    jstart = jjj;
    Njoin = numel(p_join);
    Nout = numel(p_outtage);
    
else
    
    N = 33;

    % four probabilities make up the dynamics
    p_join      = 1.00; % probability that a de-central node becomes central if its neighbor is central. Compounds linearly with how many neighbors are central
    p_outtage   = linspace(0,0.1,101); % probability that a centralized nodes undergoes an outtage
    p_outtage(1)=[];
    p_recover   = 1.00; % probability that an outtage is fixed. this leaves the tile decentral

    Njoin = numel(p_join);
    Nout = numel(p_outtage);

    % create variables that store the data we're gonna save
    outtage_distribution = zeros(Njoin,Nout,N^2); % distribution of pixels impacted by an outtage
    % clusternumber_distribution = zeros(Njoin,Nout,N^2); % distribution of connected components
    clustersize_distribution = zeros(Njoin,Nout,N^2); % distribution of connected component size
    load_distribution = zeros(Njoin,Nout,N^2); % distribution of connected component size
    p = zeros(3,1); % probability a state is in empty, filled, or out

    istart = 1;
    jstart = 1;
    
end

for iii = istart:Njoin
    for jjj = jstart:Nout
        
        fprintf('[%g,%g]\tp_join : %f\tp_out : %f\n',iii,jjj,p_join(iii),p_outtage(jjj));

        %% create initial state
        s = zeros(N,N);
        s(round(end/2),round(end/2)) = 1;
        e = zeros(N^2,4);

        %% integrate for some time to make sure our initial condition doesnt matter
        for k = 1:1000*N^2
            [s,e] = update(s,e,p_join(iii),p_outtage(jjj),p_recover);
%             if mod(k,100)==0
%                 plotstate(s,e);
%             end
        end
        
        %% integrate and collect data
        for k = 1:10000*N^2
            
            [s,e,outtagesize] = update(s,e,p_join(iii),p_outtage(jjj),p_recover);
            
            % collect statistics every N^2 steps
            if mod(k,1)==0
                
                if outtagesize>0
                    outtage_distribution(iii,jjj,outtagesize) = outtage_distribution(iii,jjj,outtagesize) + 1;
                end
                cs = sum(s(:)==1);
                if cs > 0
                    clustersize_distribution(iii,jjj,cs) = clustersize_distribution(iii,jjj,cs) + 1;
                end
%                 plotstate(s,e);
                load_distribution(iii,jjj,:) = determineload(squeeze(load_distribution(iii,jjj,:)),e,round(N^2/2));
                
                p(1) = p(1) + sum(s(:)== 0);
                p(2) = p(2) + cs;
                p(3) = p(3) + sum(s(:)==-1);
                
%                 cc = bwlabel( s==1 );
%                 cn = sum(unique(cc)>0);
%                 if cn>0
%                     clusternumber_distribution(iii,jjj,cn) = clusternumber_distribution(iii,jjj,cn) + 1;
%                     for i = 1:cn
%                         cs = sum(cc(:)==i);
%                         clustersize_distribution(iii,jjj,cs) = clustersize_distribution(iii,jjj,cs) + 1;
%                     end
%                 end

            end
            
        end
        
        %% Save data
        plabel = {'empty','filled','outtage'};
        save('data/utilityfailurescan','plabel','p','clustersize_distribution','load_distribution','outtage_distribution','iii','jjj','p_join','p_outtage','p_recover','N');
        
    end
end


function plotstate(s,e)

imagesc(s)
hold on
% caxis([-1 1])
reccursiveplot(s,e,round(numel(s)/2));
% colormap(flip(gray));
caxis([-1 1]);
drawnow
hold off

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

function v =mod1(v,N)
v = mod(v-1,N)+1;
end

function [dist,load] = determineload(dist,e,k)

% [y,x] = ind2sub([sqrt(numel(dist)) sqrt(numel(dist))],k)
load = 1;
for i = 1:size(e,2)
   if e(k,i)
       [dist,childload] = determineload(dist,e,e(k,i));
%        e(k,:)
       load = load + childload;
   end
end
dist(load) = dist(load)+1;

% dist(1:8)
end

function [s,e,outtagesize] = update(s,e,p1,p2,p3)

[Ny,Nx] = size(s);
Nn = size(e,2);
% neighborsum = imfilter(s,[0 1 0; 1 0 1; 0 1 0],'circular');
% neighborsum = conv2(s,[0 1 0; 1 0 1; 0 1 0],'same');
neighborsum = round(conv_fft2(s,[0 1 0; 1 0 1; 0 1 0],'wrap'));
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