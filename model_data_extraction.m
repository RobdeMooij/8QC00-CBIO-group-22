close all, clear all

% directory of .mat files
dir = '';

cd(dir)
all_stuff = string(ls);
index = 0;
nr_reps = 25;
step_size = 15; %[days]
nr_sit = 9;
nr_steps_1 = 60/step_size;
nr_steps_2 = 120/step_size;
nr_steps_3 = 60/step_size;
total_steps = nr_steps_1+nr_steps_2+nr_steps_3;
end_values_array = NaN(nr_reps*nr_sit,8);
full_values_array = NaN(nr_reps*nr_sit,total_steps,4);
IT = ones(nr_reps*nr_sit,1);

for i = 1:length(all_stuff)
    filename = all_stuff{i};
    if ~isfolder(filename)
        [~,name,ext] = fileparts(all_stuff{i});        
        if strcmp(ext(1:min(end,4)),'.mat')
            load(filename)
            for j = 1:nr_reps
                summary = summaryOne{j};
                for k = 1:size(summaryOne{j},2)
                    full_values_array(index*nr_reps+j,k,1) = summary{k}.TU_Num;
                    full_values_array(index*nr_reps+j,k,2) = summary{k}.IM1_Num;
                    full_values_array(index*nr_reps+j,k,3) = summary{k}.TU_purity;
                    full_values_array(index*nr_reps+j,k,4) = summary{k}.pblock_avg;
                end
                
                summary = summaryTwo{j};
                end_values_array(index*nr_reps+j,1) = summary{end}.TU_Num;
                end_values_array(index*nr_reps+j,2) = summary{end}.IM1_Num;
                end_values_array(index*nr_reps+j,3) = summary{end}.TU_purity;
                end_values_array(index*nr_reps+j,4) = summary{end}.pblock_avg;
                for k = 1:size(summaryTwo{j},2)
                    full_values_array(index*nr_reps+j,k+nr_steps_1,1) = summary{k}.TU_Num;
                    full_values_array(index*nr_reps+j,k+nr_steps_1,2) = summary{k}.IM1_Num;
                    full_values_array(index*nr_reps+j,k+nr_steps_1,3) = summary{k}.TU_purity;
                    full_values_array(index*nr_reps+j,k+nr_steps_1,4) = summary{k}.pblock_avg;
                end
                
                summary = summaryThree{j};
                if summary{end}.stepsDone <= 1
                    IT(index*nr_reps+j) = 0;
                else
                    end_values_array(index*nr_reps+j,5) = summary{end}.TU_Num;
                    end_values_array(index*nr_reps+j,6) = summary{end}.IM1_Num;
                    end_values_array(index*nr_reps+j,7) = summary{end}.TU_purity;
                    end_values_array(index*nr_reps+j,8) = summary{end}.pblock_avg;

                    for k = 1:size(summaryThree{j},2)
                        full_values_array(index*nr_reps+j,k+nr_steps_1+nr_steps_2,1) = summary{k}.TU_Num;
                        full_values_array(index*nr_reps+j,k+nr_steps_1+nr_steps_2,2) = summary{k}.IM1_Num;
                        full_values_array(index*nr_reps+j,k+nr_steps_1+nr_steps_2,3) = summary{k}.TU_purity;
                        full_values_array(index*nr_reps+j,k+nr_steps_1+nr_steps_2,4) = summary{k}.pblock_avg;
                    end
                end
            end
            index = index+1;
        end
    end
end
means = nanmean(end_values_array);

% 1 = TU_Num // 2 = IM1_Num // 3 = TU_purity // 4 = pblock_avg
plot_value = 3;

day_list = linspace(0,(total_steps)*step_size,total_steps);
titles = {'Low pblock_avg and Low IM1_Num';
          'High pblock_avg and Low IM1_Num';
          'Low pblock_avg and High IM1_Num';
          'High pblock_avg and High IM1_Num'};
ylabels = {'Tumor [cells]';
           'CD8+ [cells]';
           'Tumor purity [-]';
           'pblock_avg [-]'};
       
% before immunotherapy
low_IM1_Num = end_values_array(:,2) <= means(2);
low_pblock_avg = end_values_array(:,4) <= means(4);
LL_LH_HH_HL = [boolean(low_IM1_Num.*low_pblock_avg)...
                boolean(low_IM1_Num.*(1-low_pblock_avg))...
                boolean((1-low_IM1_Num).*(1-low_pblock_avg))...
                boolean((1-low_IM1_Num).*low_pblock_avg)];

figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
max_y = max([max(full_values_array(:,:,plot_value)) 1]); 
for i = 1:4
    subplot(2,2,i)
    hold on
    ylabel(ylabels{plot_value})
    plot(transpose(repmat(day_list,sum(LL_LH_HH_HL(:,i)),1)),transpose(full_values_array(LL_LH_HH_HL(:,i),:,plot_value)),'r:','LineWidth',0.5)
    plot(day_list,transpose(nanmean(full_values_array(LL_LH_HH_HL(:,i),:,plot_value),1)),'r','LineWidth',1.5)
    plot(step_size*[nr_steps_1 nr_steps_1],[0 max_y],'k:','LineWidth',1.5)
    plot(step_size*[nr_steps_1+nr_steps_2 nr_steps_1+nr_steps_2],[0 max_y],'k:','LineWidth',1.5)
    xlabel('Time [days]')
    axis([0 (total_steps)*step_size 0 max_y])
    title(titles{i},'Interpreter','none')
end

% after immunotherapy
change_pblock_avg = end_values_array(:,5)./end_values_array(:,1);
PD = (change_pblock_avg >= 1).*IT;
CR = (change_pblock_avg == 0).*IT;
PR = (1-PD).*(1-CR).*IT;

figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
for i = 1:4
    subplot(2,2,i)
    hold on
    ylabel(ylabels{plot_value})
    plot(transpose(repmat(day_list,sum(LL_LH_HH_HL(:,i).*PD),1)),transpose(full_values_array(boolean(LL_LH_HH_HL(:,i).*PD),:,plot_value)),'r:','LineWidth',0.5,'HandleVisibility','off')    
    plot(transpose(repmat(day_list,sum(LL_LH_HH_HL(:,i).*PR),1)),transpose(full_values_array(boolean(LL_LH_HH_HL(:,i).*PR),:,plot_value)),'g:','LineWidth',0.5,'HandleVisibility','off')
    plot(transpose(repmat(day_list,sum(LL_LH_HH_HL(:,i).*CR),1)),transpose(full_values_array(boolean(LL_LH_HH_HL(:,i).*CR),:,plot_value)),'b:','LineWidth',0.5,'HandleVisibility','off')
    plot(day_list,transpose(nanmean(full_values_array(boolean(LL_LH_HH_HL(:,i).*PD),:,plot_value),1)),'r','LineWidth',1.5)
    plot(day_list,transpose(nanmean(full_values_array(boolean(LL_LH_HH_HL(:,i).*PR),:,plot_value),1)),'g','LineWidth',1.5)
    plot(day_list,transpose(nanmean(full_values_array(boolean(LL_LH_HH_HL(:,i).*CR),:,plot_value),1)),'b','LineWidth',1.5)
    plot(step_size*[nr_steps_1+nr_steps_2 nr_steps_1+nr_steps_2],[0 max_y],'k:','LineWidth',1.5,'HandleVisibility','off')
    plot(step_size*[nr_steps_1 nr_steps_1],[0 max_y],'k:','LineWidth',1.5,'HandleVisibility','off')
    legend({ 'Progressive dissease','Partial response','Complete response'},'Location','best')
    xlabel('Time [days]')
    axis([0 (total_steps)*step_size 0 max_y])
    title(titles{i},'Interpreter','none')
end


%make the boxplots op punt step_size*(nr_steps_1+nr_steps_2)
figure('Name','Boxplots of the tumor purities w.r.t the different classes','NumberTitle','off','units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
subplot(2,2,1)
valLL=nonzeros(end_values_array(:,3).*LL_LH_HH_HL(:,1));
boxplot(valLL(~isnan(valLL)))
ylabel(ylabels{3})
title(titles{1})
subplot(2,2,2)
valHL=nonzeros(end_values_array(:,3).*LL_LH_HH_HL(:,4));
boxplot(valHL(~isnan(valHL)))
ylabel(ylabels{3})
title(titles{3})
subplot(2,2,3)
valLH=nonzeros(end_values_array(:,3).*LL_LH_HH_HL(:,2));
boxplot(valLH(~isnan(valLH)))
ylabel(ylabels{3})
title(titles{2})
subplot(2,2,4)
valHH=nonzeros(end_values_array(:,3).*LL_LH_HH_HL(:,3));
boxplot(valHH(~isnan(valHH)))
ylabel(ylabels{3})
title(titles{4})

matrix=[valLL(~isnan(valLL)); valHL(~isnan(valHL)); valLH(~isnan(valLH)); valHH(~isnan(valHH))];
cellf=strings(length(matrix),1);
cellf(1:length(valLL(~isnan(valLL))))=titles{1};
cellf(length(valLL(~isnan(valLL))):length([valLL(~isnan(valLL)); valHL(~isnan(valHL))]))=titles{3};
cellf(length([valLL(~isnan(valLL)); valHL(~isnan(valHL))]):length([valLL(~isnan(valLL)); valHL(~isnan(valHL)); valLH(~isnan(valLH))]))=titles{2};
cellf(length([valLL(~isnan(valLL)); valHL(~isnan(valHL)); valLH(~isnan(valLH))]):length(matrix))=titles{4};
figure
set(gcf,'color','w');
boxplot(matrix,cellf)
set(gca,'XTickLabelRotation',30)
ylabel('Tumor purity')

clear cellf
%make the bar plots
figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
cat = categorical({'High pblock_avg and High IM1_Num';
       'Low pblock_avg and Low IM1_Num';
       'Low pblock_avg and High IM1_Num';
       'High pblock_avg and Low IM1_Num'}); 
%high high
PR1=length(~isnan(nonzeros(LL_LH_HH_HL(:,3).*PR)));
PD1=length(~isnan(nonzeros(LL_LH_HH_HL(:,3).*PD)));
CR1=length(~isnan(nonzeros(LL_LH_HH_HL(:,3).*CR)));
%low low
PR2=length(~isnan(nonzeros(LL_LH_HH_HL(:,1).*PR)));
PD2=length(~isnan(nonzeros(LL_LH_HH_HL(:,1).*PD)));
CR2=length(~isnan(nonzeros(LL_LH_HH_HL(:,1).*CR)));
%high cd8 and low pdl1
PR3=length(~isnan(nonzeros(LL_LH_HH_HL(:,4).*PR)));
PD3=length(~isnan(nonzeros(LL_LH_HH_HL(:,4).*PD)));
CR3=length(~isnan(nonzeros(LL_LH_HH_HL(:,4).*CR)));
%low cd8 and high pdl1
PR4=length(~isnan(nonzeros(LL_LH_HH_HL(:,2).*PR)));
PD4=length(~isnan(nonzeros(LL_LH_HH_HL(:,2).*PD)));
CR4=length(~isnan(nonzeros(LL_LH_HH_HL(:,2).*CR)));
[PR1 PD1 CR1; PR2 PD2 CR2; PR3 PD3 CR3; PR4 PD4 CR4]
bar(cat,[PR1 PD1 CR1; PR2 PD2 CR2; PR3 PD3 CR3; PR4 PD4 CR4],'stacked')
set(gca,'TickLabelInterpreter','none');
ylabel('Number of occurrences')
legend('Partial Response','Progressive Disease','Complete Response','Location','northwest')
posx=[1 1 1 4 4 4 3 3 3 2 2 2];
posy=[PR1/2 PR1+PD1/2 PR1+PD1+CR1/2 PR2/2 PR2+PD2/2 PR2+PD2+CR2/2 PR3/2 PR3+PD3/2 PR3+PD3+CR3/2 PR4/2 PR4+PD4/2 PR4+PD4+CR4/2];
cellf={num2str(PR1/(PR1+PD1+CR1)), num2str(PD1/(PR1+PD1+CR1)),num2str(CR1/(PR1+PD1+CR1)),num2str(PR2/(PR2+PD2+CR2)),num2str(PD2/(PR2+PD2+CR2)),num2str(CR2/(PR2+PD2+CR2)),num2str(PR3/(PR3+PD3+CR3)),num2str(PD3/(PR3+PD3+CR3)),num2str(CR3/(PR3+PD3+CR3)),num2str(PR4/(PR4+PD4+CR4)),num2str(PD4/(PR4+PD4+CR4)),num2str(CR4/(PR4+PD4+CR4))};
idx = strcmp(cellf,"0");
posx(idx) =  [];
posy(idx) = [];
cellf(idx) = [];
text(posx,posy,cellf,'HorizontalAlignment','center');

