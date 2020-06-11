function plot_net_single_fib_type(scale, nodes, bnd_nodes, fibers, fib_type, fib_rads, rve_stretch)


% plot_net(nodes, fibers)
%
% in:
%
% nodes            N x 3 coordinates for N nodes
% fibers           N x 2 start-end nodes for N fibers
%
% last rev:
%
% tue nov 6 2012 mfh


figure;

nodes = nodes.*scale;
bnd_nodes = bnd_nodes.*scale;
rve_stretch = rve_stretch.*scale;

plot3(nodes(:,1) , nodes(:,2), nodes(:,3),'o', 'LineWidth',0.2, 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',5);
hold on
plot3(bnd_nodes(:,1) , bnd_nodes(:,2), bnd_nodes(:,3),'X', 'LineWidth',0.4, 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',5);

m = 1 ;

max_rad = max(fib_rads);

for n = 1 : length(fibers) % count rows
    
    if fib_rads(n) > 2e-9
        
        lw = 6*(fib_rads(n)/max_rad);
        
        node1 = nodes(fibers(n,1),:) ;
        
        for i = 1:3
            for k = 1:abs(fibers(n,2+i))
                node2 = bnd_nodes(m,:) ;
                
                x1(1) = node1(1); % node 1 x coord
                y1(1) = node1(2); % node 1 y coord
                z1(1) = node1(3); % node 1 z coord
                
                x1(2) = node2(1); % node 2 x coord
                y1(2) = node2(2); % node 2 y coord
                z1(2) = node2(3); % node 2 z coord
                
                if fib_type(n) == 1
                    C = [0.8 0.8 0.4];
                    %lw = 1;
                elseif fib_type(n) == 2
                    C = [0 0.2 0];
                    %lw=4;
                elseif fib_type(n) == 3
                    C = [0.7 0 0];
                    %lw=3;
                end
                
                plot3(x1, y1, z1, 'Color', C,'LineWidth',lw);
                
                m = m+1 ;
                
                node1 = bnd_nodes(m,:) ;
                
                m = m+1 ;
            end
        end
        
        node2 = nodes(fibers(n,2),:) ;
        
        x1(1) = node1(1); % node 1 x coord
        y1(1) = node1(2); % node 1 y coord
        z1(1) = node1(3); % node 1 z coord
        
        x1(2) = node2(1); % node 2 x coord
        y1(2) = node2(2); % node 2 y coord
        z1(2) = node2(3); % node 2 z coord
        
        if fib_type(n) == 1
            C = [0.8 0.8 0.4];
            lw = 1;
        elseif fib_type(n) == 2
            C = [0 0.2 0];
            lw = 4;
        elseif fib_type(n) == 3
            C = [0.7 0 0];
            lw = 3;
        end
        
        plot3(x1, y1, z1, 'Color', C, 'LineWidth',lw);
        
    end
end

row1x = [-0.5,(-0.5+rve_stretch(1))];
row2x = (-0.5+rve_stretch(1)).*ones(length(row1x),1) ;
lowrow = [-0.5,-0.5] ;
row1y = [-0.5,(-0.5+rve_stretch(2))];
row2y = (-0.5+rve_stretch(2)).*ones(length(row1y),1) ;
row1z = [-0.5,(-0.5+rve_stretch(3))];
row2z = (-0.5+rve_stretch(3)).*ones(length(row1z),1) ;

plot3(row1x, lowrow, lowrow,'k--')
plot3(row1x, row2y, lowrow,'k--')
plot3(row1x, lowrow, row2z,'k--')
plot3(row1x, row2y, row2z,'k--')
plot3(lowrow, row1y, lowrow,'k--')
plot3(row2x, row1y, lowrow,'k--')
plot3(lowrow, row1y, row2z,'k--')
plot3(row2x, row1y, row2z,'k--')
plot3(lowrow, lowrow, row1z,'k--')
plot3(lowrow, row2y, row1z,'k--')
plot3(row2x, lowrow, row1z,'k--')
plot3(row2x, row2y, row1z,'k--')

set(gcf, 'color', 'white');
axis equal;

% axis([-0.65 -0.45+rve_stretch(1) -0.65 -0.45+rve_stretch(2) -0.65 -0.45+rve_stretch(3) ]);

xlabel('Circumferential [um]')'; ylabel('Axial [um]'); zlabel('Radial [um]');
set(gca,'FontSize',16);

end
