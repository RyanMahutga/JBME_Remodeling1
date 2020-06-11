function WriteNet2( stab_nodes, fibers, fibtype, init_lens, fib_rads, nodes, x_scale, m, netset)
%WriteNet.m writes the periodic network to a file for use in c++ 
%   Detailed explanation goes here

num_nodes = length(nodes);
num_fibers=length(fibers);

if m==1
    dir = strcat('UNP',num2str(netset));
elseif m==2
    dir = strcat('UP',num2str(netset));
elseif m==3
    dir = strcat('CNP',num2str(netset));
elseif m==4
    dir = strcat('CP',num2str(netset));
else
    print('Error!')
end
mkdir(dir);
old=cd(dir);

file='PeriodicNetwork.txt';

fileID = fopen(file,'w');

stab_nodes = stab_nodes-1;

fprintf(fileID,'% d % d % d % d % d % d\r\n',num_fibers, num_nodes, x_scale, stab_nodes);

for  j=1:num_fibers
    
    fibers(j,1) = fibers(j,1)-1;
    fibers(j,2) = fibers(j,2)-1;
    fprintf(fileID,'% d % d % d % d % d % d % d % d\r\n',fibers(j,:), fibtype(j), init_lens(j) ,fib_rads(j));
    
end

for k = 1:num_nodes
    if k<num_nodes
        fprintf(fileID,'% d % d % d\r\n',nodes(k,:));
    else
       fprintf(fileID,'% d % d % d',nodes(k,:)); 
    end
end
fclose(fileID);
cd(old)
end

