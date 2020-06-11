% Create periodic MLU
% Author: Ryan Mahutga
% Barocas Research Group
% Date: 11-01-18


% close all
% clear all
% clc

rng('shuffle','twister')

netset = 8; % netset for saving

plot_net=1;

%fiber radii distribution
mu = 10e-9;
sigma = 2e-9;

min_len = 0.05;

numEnodes = 250; % number of elastic lamina nodes
numAnodes = 100; % number of actin nodes

keep_sets = [0.5,1] ; % keep nodes for connecting actin (1:keep_set(i):end) are nodes kept

island_scale_factor = 0.95 ; % if this is less than one nodes come from a smaller region (we use this to prevent things being pulled across boundary)
 
nodesE = island_scale_factor.*(rand(numEnodes,3) - 0.5) ;
nodesA = island_scale_factor.*(rand(numAnodes,2) - 0.5) ;
nodesA(:,3) = 0.35.*(rand(numAnodes,1)) + 0.1 ;
nodesA(round(0.5*numAnodes):end,3) = -nodesA(round(0.5*numAnodes):end,3);

% select fixed nodes to stabilize networks (closest to [0,0,0])
% elastin fixed node
diststabE=[];
% find a node closest to (0,0,0)
for m=1:length(nodesE)
    diststabE(m) = norm(nodesE(m,1:2));
end
[~,indexstabE] = min(diststabE);

%% Generating ECM Networks

num_nodes = length(nodesE);

% percolating elastin network
[pE_fibers] = periodicDelaunay(nodesE);
[pE_fibers] = removeDupes( pE_fibers );
%[pE_fibers,deg_pE] = NetworkPare(pE_fibers, num_nodes,12 );
pE_fibers(:,5) = 0; % remove intersections across z-boundary

% flatten nodes into plane
nodesE(:,3) = 0;
nodes_All = [nodesE;nodesA]; % all nodes;

% percolating collagen network
[pC_fibers] = periodicDelaunay2D(nodesE);
[pC_fibers] = removeDupes(pC_fibers);
[~, init_lens_C ] = networkFeatures( nodesE, pC_fibers, mu, sigma);
  [idx] = find(init_lens_C < min_len); % remove short fibers
  pC_fibers(idx,:) = [];
% [pC_fibers,deg_p] = NetworkPare(pC_fibers, num_nodes,4);

% non-percolating collagen network
% [npC_fibers] = interiorDelaunay2D(nodesE);
% [~, init_lens_C ] = networkFeatures( nodesE, npC_fibers, mu, sigma);
%   [idx] = find(init_lens_C < min_len); % remove short fibers
%   npC_fibers(idx,:) = [];
%[npC_fibers,deg_np] = NetworkPare(npC_fibers, num_nodes,4);

%% Generating Actin Networks

% unconnected actin
[uA_fibers] = periodicDelaunay(nodesA);
[uA_fibers] = removeDupes( uA_fibers );
uA_fibers(:,1) = uA_fibers(:,1) + numEnodes;
uA_fibers(:,2) = uA_fibers(:,2) + numEnodes;

j=1;
for p = 1:length(uA_fibers)
   node1 = uA_fibers(p,1);
   node2 = uA_fibers(p,2);
   real_coord2 = nodes_All(node2,:) + uA_fibers(p,3:5);
   uA_init_lens(p) = norm(real_coord2-nodes_All(node1,:));
   
   if nodes_All(node1,3)*nodes_All(node2,3) < 0 && uA_fibers(p,5) == 0 % z-dim have opposite signs (i.e. crosses x,y plane at z=0
       cross_idx(j) = p; 
       [ x,y ] = XYintercept( nodes_All(node1,:), real_coord2) ;
       
       cross_coord(j,:) = [x,y,0];
       
       short_lens(j,2) = norm(real_coord2-cross_coord(j,:)); % node2 to center
       short_lens(j,1) = norm(cross_coord(j,:)-nodes_All(node1,:)); % center to node1
       
       fibcross(2*j-1,1)=0; fibcross(2*j-1,2)=0;
       fibcross(2*j,1)=0; fibcross(2*j,2)=0;
       
       tagx = 0; tagy=0;
       while x<-0.5 || x>0.5
           if x<0
              x = x+1; 
           else
               x=x-1;
           end
           fibcross(2*j-1,1) = uA_fibers(p,3);
           tagx=1;
       end
       while y<-0.5 || y>0.5
           if y<0
               y = y+1;
           else
               y=y-1;
           end
           fibcross(2*j-1,2) = uA_fibers(p,4);
           tagy=1;
       end
       
       if tagx==0
           fibcross(2*j,1) = uA_fibers(p,3);
       end
       if tagy==0
           fibcross(2*j,2) = uA_fibers(p,4);
       end
       
       fibcross(:,3) = 0;
       
       cross_coord(j,:) = [x,y,0];
       
       j=j+1; % increment
   end
end

%TRI = delaunayn(nodesE(:,1:2));
% close_idx = dsearchn(nodesE(:,1:2),TRI,cross_coord(:,1:2));

for i = 1:length(cross_coord)
    I = knnsearch(nodesE(:,1:2),cross_coord(i,1:2),'K',2);
    close_idx1(i) = I(1);
    close_idx2(i) = I(2);
end

% figure
% plot3(cross_coord(:,1),cross_coord(:,2),cross_coord(:,3),'rx')
% hold on
% plot3(nodesE(:,1),nodesE(:,2),nodesE(:,3),'bo')
% plot3(nodesE(close_idx,1),nodesE(close_idx,2),nodesE(close_idx,3),'ko')

    Afibers=[];
    idx_keep=[];
    Alens=[];
    Afibers = uA_fibers;
    Alens = uA_init_lens;
    
       idx_keep=1:length(cross_idx);
    
%     num_keep = round(length(cross_coord)*keep_sets(n));
    
%     idx_keep=datasample(1:length(cross_idx),num_keep,'Replace',false); % randomly select nodes to connect to cell
    
    for p = 1:length(idx_keep)
        fib_idx = cross_idx(idx_keep(p));
        
        % tacking split fibers onto the end
        Afibers = [Afibers; Afibers(fib_idx,1), close_idx1(idx_keep(p)), fibcross(2*idx_keep(p)-1,:)];
        Alens = [Alens, short_lens(idx_keep(p),1)];
        
        Afibers = [Afibers; close_idx2(idx_keep(p)), Afibers(fib_idx,2), fibcross(2*idx_keep(p),:)];
        Alens = [Alens, short_lens(idx_keep(p),2)];
    end
    Afibers(cross_idx(idx_keep),:) = []; % remove fiber after split
    Alens(cross_idx(idx_keep))=[]; % remove lens after split
    
    fibs = Afibers;
    lens = Alens;

id1 = {'0','10','25','50'};
id2 = {'pE_pC_'}; %,'pE_npC_','npE_npC_'}; % indetifiers for percolating and nonpercolating networks of collagen and elastin

m=1;
for k=1:4
    fibers1=[];
    
    fibers1 = fibs;
    stab_nodes = [indexstabE, -1,-1];
        
    len1=length(fibers1);
    
        fibers2=[];
        if k==1 % both percolating
            fibers2 = [pE_fibers; pC_fibers];
            len2 = length(pE_fibers);
            len3 = length(pC_fibers);
        elseif k==2 % 10% removed elastin 
            numE = length(pE_fibers);
            v = randsample(numE, round(9*numE/10));
            fibers2 = [pE_fibers(v,:) ; pC_fibers];
            len2 = length(pE_fibers(v,:));
            len3 = length(pC_fibers);
        elseif k==3 % 25% removed elastin
            numE = length(pE_fibers);
            v = randsample(numE, round(3*numE/4));
            fibers2 = [pE_fibers(v,:) ; pC_fibers];
            len2 = length(pE_fibers(v,:));
            len3 = length(pC_fibers);
         elseif k==4 % 50% removed elastin
            numE = length(pE_fibers);
            v = randsample(numE, round(numE/2));
            fibers2 = [pE_fibers(v,:) ; pC_fibers];
            len2 = length(pE_fibers(v,:));
            len3 = length(pC_fibers);
        end
        
        % combining actin and elastin/collagen fibers
        fibers = [fibers1;fibers2];
        
        % setting fiber types
        fib_type = ones(length(fibers),1);
        fib_type(len1+1:len1+len2) = 2;
        fib_type(len1+len2+1:len1+len2+len3) = 3;
        
        [ nodes, fibers, stab_nodes ] = MinDegreeReorder( nodes_All, fibers, stab_nodes ) ; % minimum degree reordering
        [ fib_rads, init_lens] = networkFeatures( nodes, fibers, mu, sigma);
        
        %init_lens(1:len1) = lens{k};
        
        dir_id = strcat(id1{k},id2{1});
        
        WriteNet3( stab_nodes, fibers, fib_type, init_lens, fib_rads, nodes, 10e-6, dir_id, netset)
        
        if plot_net==1
            
%             if k<5
%                 subplot(2,3,m+1)
%             else
%                 subplot(2,3,1)
%             end
            rve_stretch = [1 1 1];
            fiber_vol_fract = 0.35;
            
            [ nodes_n, bnd_nodes_n, net_stress, fib_stress, fib_forces, fibers_n ]...
                = solve_periodic_BCs2( nodes, fibers, fib_type, init_lens', fib_rads,...
                fiber_vol_fract, rve_stretch,[] ) ;
            
            ConvPeriodic2FixedNetwork( nodes_n, bnd_nodes_n, fibers_n, fib_type)
            
            plot_net_single_fib_type(nodes_n, bnd_nodes_n, fibers_n, fib_type, rve_stretch)
            
            plot3(nodes(stab_nodes(1),1), nodes(stab_nodes(1),2),nodes(stab_nodes(1),3),'bo','LineWidth',3,'MarkerSize',7)
            
            m=m+1;
        end
end
        