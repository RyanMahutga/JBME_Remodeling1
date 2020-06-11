function [fibers] = interiorDelaunay2D(nodes)
% interiorDelaunay makes a delaunay tesselation of the interior nodes
%   Detailed explanation goes here
TRI = delaunay(nodes(:,1),nodes(:,2));
fibers = zeros(3*length(TRI),5);

for p=1:length(TRI)
    idx = 3*p;
    fibers(idx-2,1) = TRI(p,1);
    fibers(idx-2,2) = TRI(p,2);
    fibers(idx-1,1) = TRI(p,2);
    fibers(idx-1,2) = TRI(p,3);
    fibers(idx,1) = TRI(p,3);
%     fibers(4*p-3+2,2) = TRI(p,4);
%     fibers(4*p-3+3,1) = TRI(p,4);
    fibers(idx,2) = TRI(p,1);
end

end

