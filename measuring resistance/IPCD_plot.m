%% Read and analyze excel sheet of IPCD data
load allData


%% growth and delay
for drug = 1:size(allData,2)
    M{drug}=growthanddelay(allData{drug},allT{drug}/12,string(allStrains{drug}(:,1)),drug);
end

%% resistances

for drug = 1:size(allData,2)
    N{drug}=resistances(M{drug},allDrug{drug},string(allStrains{drug}(:,1)));
end

%% plot scatters

names = {'carb';'cipro';'gent';'h2o2';'spec'};

for drug = 1:size(allData,2)
    plot_screen_cat(N{drug},names{drug},min(min(N{drug}(:,1:2))),max(max(N{drug}(:,1:2))),string(allStrains{drug}(:,2)));
end

%%
% H = human infection
% A = animal infection
% E = environmental
% C = CF lung
% U = unknown ("2732" in H2O2)

% strains{:} = carb, cip, gent, h2o2, spec

load allData

%N = {N_carb,N_cipro,N_gent,N_h2o2,N_spec};
allTab = cell(size(allStrains{4},1),12);
allTab(:,1:2) = allStrains{4};

for drug = 1:size(allData,2)
    for strain = 1:size(allStrains{drug},1)
        compare1 = allStrains{drug}(strain);
        idx1 = find(string(compare1)==string(allTab(:,1)));

        % Growth, dyn, ss resist
        allTab(idx1,drug*2+1:drug*2+2) = num2cell(N{drug}(strain,:));
    end
end

allTab = [{'strain' 'category' 'carb dyn' 'carb ss' ...
    'cipro dyn' 'cipro ss' 'gent dyn' 'gent ss' ...
    'h2o2 dyn' 'h2o2 ss' 'spec dyn' 'spec ss'};allTab];

writecell(allTab,'IPCD_ed.xlsx');
