close all

N = 33;

% four probabilities make up the dynamics
p_join      = linspace(0.5,1,11); % probability that a de-central node becomes central if its neighbor is central. Compounds linearly with how many neighbors are central
p_outtage   = linspace(0.1,1,11); % probability that a centralized nodes undergoes an outtage
p_recover   = 1.00; % probability that an outtage is fixed. this leaves the tile decentral

Njoin = numel(p_join);
Nout = numel(p_outtage);

% create variables that store the data we're gonna save
outtage_distribution = zeros(Njoin,Nout,N^2); % distribution of pixels impacted by an outtage
% clusternumber_distribution = zeros(Njoin,Nout,N^2); % distribution of connected components
clustersize_distribution = zeros(Njoin,Nout,N^2); % distribution of connected component size
load_distribution = zeros(Njoin,Nout,N^2); % distribution of connected component size
p = zeros(3,1); % probability a state is in empty, filled, or out

for iii = 1:Njoin
    for jjj = 1:Nout
        
        fprintf('[%g,%g]\tp_join : %f\tp_out : %f\n',iii,jjj,p_join(iii),p_outtage(jjj));

        %% create initial state
        graph = zeros(N,N);
        graph = num2cell(graph);
        graph{round(N/2),round(N/2)} = Node(round(N/2),round(N/2));

        %% integrate for some time to make sure our initial condition doesnt matter
        for k = 1:1000*N^2
            graph = update(graph,p_join(iii),p_outtage(jjj),p_recover);
            if mod(k,1)==0
                plotstate(graph);
            end
        end
        
        %% integrate and collect data
        for k = 1:10000*N^2
            
            [graph,outtagesize] = update(graph,p_join(iii),p_outtage(jjj),p_recover);
            
            % collect statistics every N^2 steps
            if mod(k,1)==0
                
                
                s = cellfun(@state, graph);
                cs = sum(s(:)== 1);
                p(1) = p(1) + cs;
                p(2) = p(2) + sum(s(:)== 0);
                p(3) = p(3) + sum(s(:)==-1);
                
                if outtagesize>0
                    outtage_distribution(iii,jjj,outtagesize) = outtage_distribution(iii,jjj,outtagesize) + 1;
                end
                if cs > 0
                    clustersize_distribution(iii,jjj,cs) = clustersize_distribution(iii,jjj,cs) + 1;
                end
                load_distribution(iii,jjj,:) = determineload(squeeze(load_distribution(iii,jjj,:)),graph{round(N/2),round(N/2)});
                
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
        save('data/utilityfailurescan','plabel','p','clustersize_distribution','outtage_distribution','iii','jjj','p_join','p_outtage','p_recover','N');
        
    end
end

function [dist,nodesize] = determineload(dist,node)

nodesize = 1;
for i = 1:numel(node.children)
   [dist,childsize] = determineload(dist,node.children{i});
   nodesize = nodesize + childsize;
end
dist(nodesize) = dist(nodesize)+1;

end

function s = state(s)

if isstruct(s)
    s = 1;
end

end

function v =mod1(v,N)
v = mod(v-1,N)+1;
end

function n = Node(i,j)

n.x = j;
n.y = i;
n.children = {};

end

function [g,outtagesize] = update(g,p1,p2,p3)

[Ny,Nx] = size(g);
Nn = 4;
% neighborsum = imfilter(s,[0 1 0; 1 0 1; 0 1 0],'circular');
s = cellfun(@state, g);
neighborsum = conv2(s,[0 1 0; 1 0 1; 0 1 0],'same');
outtagesize = 0;

k = randi(Ny*Nx);
[i,j] = ind2sub([Ny,Nx],k);
sij = state(g{i,j});

% if current pixel is OUT
if sij == -1
    
    % restore to decentralization with probability p4
    if rand < p3
        g{i,j} = 0;
    end
    
% if current pixel is DECENTRAL
elseif sij == 0
    
    % become central with probability p1 if you have a central neighbor
    if neighborsum(i,j)>0 && rand < p1
        
        % we need to assign one of the neighboring centralized nodes to be the parent. Lets do this in a random fashion
        for n = randperm(Nn)
            if n==1
                ny = mod1(i+1,Ny);
                if isstruct(g{ny,j})
                    g{i,j} = Node(i,j);
                    g{ny,j}.children{end+1} = g{i,j};
                    break % if one of these conditions is set, break, since we only want to set ONE neighbor as the parent
                end
            elseif n==2
                nx = mod1(j-1,Nx);
                if isstruct(g{i,nx})
                    g{i,j} = Node(i,j);
                    g{i,nx}.children{end+1} = g{i,j};
                    break
                end
            elseif n==3
                ny = mod1(i-1,Ny);
                if isstruct(g{ny,j})
                    g{i,j} = Node(i,j);
                    g{ny,j}.children{end+1} = g{i,j};
                    break
                end
            elseif n==4
                nx = mod1(j+1,Nx);
                if isstruct(g{i,nx})
                    g{i,j} = Node(i,j);
                    g{i,nx}.children{end+1} = g{i,j};
                    break
                end
            end
        end
    end
        
% if current pixel is CENTRAL
elseif sij == 1
    % experience an outtage with probability p3
    if rand < p2
        ny = mod1(i+1,Ny);
        g{ny,j}.children{cellfun(@(x) isequal(x,g{i,j}))} = [];
        nx = mod1(j-1,Nx);
        g{i,nx}.children{cellfun(@(x) isequal(x,g{i,j}))} = [];        
        ny = mod1(i-1,Ny);
        g{ny,j}.children{cellfun(@(x) isequal(x,g{i,j}))} = [];       
        nx = mod1(j+1,Nx);
        g{i,nx}.children{cellfun(@(x) isequal(x,g{i,j}))} = [];   
        
        [g,outtagesize] = delete(g,g{i,j});
    end
end


if state(g{round(Ny/2),round(Nx/2)}) == 0
    g{round(Ny/2),round(Nx/2)} = Node(round(Ny/2),round(Nx/2));
end

end

function [g,outtagesize] = delete(g,node)
    
outtagesize = 1;
g{node.x,node.y} = -1;
for i = 1:numel(node.children)
   [g,childouttagesize] = delete(g,node.children{i});
   outtagesize = outtagesize + childouttagesize;
end

end

function plotstate(graph)

s = cellfun(@state, graph);
imagesc(s<1)
hold on
% caxis([-1 1])
reccursiveplot(graph{round(end/2),round(end/2)});
% colormap(flip(gray));
drawnow
hold off

end
function reccursiveplot(node)

for i = 1:numel(node.children)
    reccursiveplot(node.children{i});
    plot([node.x node.children{i}.x],[node.y node.children{i}.y],'k');
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