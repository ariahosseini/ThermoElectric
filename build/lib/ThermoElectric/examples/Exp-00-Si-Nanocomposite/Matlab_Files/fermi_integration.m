%% clear workspace

clear
clc
close all

%% read Fermi energy level


Ef_0pct = dlmread("Ef-no-inc.out");
Ef_1pct = dlmread("Ef-inc-1pct.out");
Ef_5pct_dir_up = dlmread("Ef-inc-5pct-dir-up.out");
Ef_5pct_dir_down = dlmread("Ef-inc-5pct-dir-down.out");


for i = 1:size(Ef_0pct,2)

f1_0pct(1,i) = fermi(0.5, Ef_0pct(1,i));
f2_0pct(1,i) = fermi(-0.5, Ef_0pct(1,i));
end

dlmwrite('fermi-0pct.out', [f1_0pct; f2_0pct]);

for i = 1:size(Ef_1pct,2)

f1_1pct(1,i) = fermi(0.5, Ef_1pct(1,i));
f2_1pct(1,i) = fermi(-0.5, Ef_1pct(1,i));
end

dlmwrite('fermi-1pct.out',[f1_1pct; f2_1pct]);

for i = 1:size(Ef_5pct_dir_up,2)

f1_5pct_dir_up(1,i) = fermi(0.5, Ef_5pct_dir_up(1,i));
f2_5pct_dir_up(1,i) = fermi(-0.5, Ef_5pct_dir_up(1,i));

end

dlmwrite('fermi-5pct-dir-up.out', [f1_5pct_dir_up; f2_5pct_dir_up]);

for i = 1:size(Ef_5pct_dir_down,2)
f1_5pct_dir_down(1,i) = fermi(0.5, Ef_5pct_dir_down(1,i));
f2_5pct_dir_down(1,i) = fermi(-0.5, Ef_5pct_dir_down(1,i));
end
dlmwrite('fermi-5pct-dir-down.out',[f1_5pct_dir_down; f2_5pct_dir_down]);

