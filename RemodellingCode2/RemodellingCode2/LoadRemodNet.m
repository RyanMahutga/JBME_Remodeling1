% comparison script
clear all
clc

dir1 = 'Remodelled';

old = cd(dir1);

k = 1060;  % time at final step

rve_stretch=load(strcat('rve_stretch',num2str(k),'.txt'));
net_stress=load(strcat('net_stress',num2str(k),'.txt'));
fib_vol_fract=load(strcat('fib_vol_fract',num2str(k),'.txt'));
pressure = load(strcat('pressure',num2str(k),'.txt'));
fib_rads = load(strcat('fib_rads',num2str(k),'.txt'));
init_lens = load(strcat('init_lens',num2str(k),'.txt'));
fibers = load(strcat('fibers',num2str(k),'.txt'));
nodes = load(strcat('nodes',num2str(k),'.txt'));
fib_type = load(strcat('fib_type',num2str(k),'.txt'));

num_nodes = length(nodes);
num_fibers = length(fibers);

cd(old)

WriteNet( num_nodes, num_fibers, fibers, fib_type, init_lens, fib_rads, nodes, fib_vol_fract )








