function [R] = calc_periodic_orient(nodes, fibers)


% calculate the 3D orientation tensor for a network -- in netmat
%
% nodes -- N x 3
% fibers -- N x 2
%
% last rev -- sept 2012 -- mfh


om11=0.0; om12=0.0; om13=0.0;
          om22=0.0; om23=0.0;
                    om33=0.0;
                    
total_length = 0.0;

num_fibers = size(fibers, 1); % rows = fibers

for i = 1 : num_fibers

    a = fibers(i, 1); % node 1 num
    b = fibers(i, 2); % node 2 num
    
    real_node2 = nodes(b,:) + fibers(:,3:5)';
    
    del_x = nodes(a, 1) - real_node2(1); % x1 - x2
    del_y = nodes(a, 2) - real_node2(2); % y1 - y2
    del_z = nodes(a, 3) - real_node2(3); % z1 - z2

    fiber_length = sqrt(del_x^2 + del_y^2 + del_z^2);

    total_length = total_length + fiber_length;
    
    % direction cosines for current fiber  
    
    cosa = del_x / fiber_length;
    cosb = del_y / fiber_length;
    cosg = del_z / fiber_length;
    
    om11 = om11 + fiber_length*cosa*cosa;
    om12 = om12 + fiber_length*cosa*cosb;
    om13 = om13 + fiber_length*cosa*cosg;
    om22 = om22 + fiber_length*cosb*cosb;
    om23 = om23 + fiber_length*cosb*cosg;
    om33 = om33 + fiber_length*cosg*cosg;

end

R = [om11 om12 om13; 
     om12 om22 om23; 
     om13 om23 om33] ./ total_length; 

end

