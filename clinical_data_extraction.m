clear all
close all
worksheet=1;

[file_location,path]=uigetfile('*.csv');
[numerical_data,textt,all_data]=xlsread([path file_location],worksheet);

PD_id=0;
TP_id=0;
CD8_id=0;

response_id=0;
PDL1=zeros(size(all_data,1)-1,1);
CD8=zeros(size(all_data,1)-1,1);
TP=zeros(size(all_data,1)-1,1);
for k=1:length(all_data)
    line=all_data{k};
    splitted=strsplit(line,',');

    if k==1
        kop=splitted;
        for i=1:length(kop)
            if strcmp(erase(convertCharsToStrings(kop{i}),'"'),"PDL1")
                PD_id=i;


            elseif strcmp(erase(convertCharsToStrings(kop{i}),'"'),"T_cells_CD8")
                CD8_id=i;

            elseif strcmp(erase(convertCharsToStrings(kop{i}),'"'),"TumorPurity")
                TP_id=i;

            elseif strcmp(erase(convertCharsToStrings(kop{i}),'"'),"response")
                response_id=i;
                response=cell(size(all_data,1)-1,1);
                Progr_id=[];
                Part_id=[];
                Comp_id=[];
            end
        end
    end

    if k>1 
        if CD8_id>0
            CD8(k-1)=str2num(splitted{CD8_id});
        end
        if TP_id>0
            TP(k-1)=str2num(splitted{TP_id});
        end
        if PD_id>0
            PDL1(k-1)=str2num(splitted{PD_id});
        end
        if response_id>0
            response{k-1}=splitted{response_id};
            if strcmp(erase(convertCharsToStrings(splitted{response_id}),'"'),"Progressive Disease")
                Progr_id=[Progr_id k-1];
            elseif strcmp(erase(convertCharsToStrings(splitted{response_id}),'"'),"Partial Response")
                Part_id=[Part_id k-1];
            elseif strcmp(erase(convertCharsToStrings(splitted{response_id}),'"'),"Complete Response")
                Comp_id=[Comp_id k-1];
            end
        end
    end  
end

%delete patient with an unknown tumor purity
if TP_id>0
    ind=find(1-isnan(TP));
    TP=TP(ind);
    PDL1=PDL1(ind);
    CD8=CD8(ind);
end
      

mean_PD=mean(PDL1);
mean_CD=mean(CD8);
low_low=[];
lowP_highC=[];
highP_lowC=[];
high_high=[];
clear i
for i=1:length(PDL1)
    if PDL1(i)<mean_PD 
        if CD8(i)<mean_CD
            low_low=[low_low i];
        elseif CD8(i)>mean_CD
            lowP_highC=[lowP_highC i];
        end
    elseif PDL1(i)>mean_PD 
        if CD8(i)<mean_CD
            highP_lowC=[highP_lowC i];
        elseif CD8(i)>mean_CD
            high_high=[high_high i];
        end
    end
end

%for without immunotherapy 
if response_id==0
figure('Name','Boxplots of the tumor purities w.r.t the different classes','NumberTitle','off')
set(gcf,'color','w');
subplot(2,2,1)
boxplot(TP(low_low))
ylabel('Tumor purity')
title('Low PDL1 and Low CD8+')
subplot(2,2,2)
boxplot(TP(lowP_highC))
ylabel('Tumor purity')
title('Low PDL1 and High CD8+')
subplot(2,2,3)
boxplot(TP(highP_lowC))
ylabel('Tumor purity')
title('High PDL1 and Low CD8+')
subplot(2,2,4)
boxplot(TP(high_high))
ylabel('Tumor purity')
title('High PDL1 and High CD8+')

matrix=[TP(low_low); TP(lowP_highC); TP(highP_lowC); TP(high_high)];
cellf=strings(length(matrix),1);
cellf(1:length(TP(low_low)))='Low PDL1 and Low CD8+';
cellf(length(TP(low_low)):length([TP(low_low); TP(lowP_highC)]))='Low PDL1 and High CD8+';
cellf(length([TP(low_low); TP(lowP_highC)]):length([TP(low_low); TP(lowP_highC); TP(highP_lowC)]))='High PDL1 and Low CD8+';
cellf(length([TP(low_low); TP(lowP_highC); TP(highP_lowC)]):length(matrix))='High PDL1 and High CD8+';
figure
set(gcf,'color','w')
boxplot(matrix,cellf)
set(gca,'XTickLabelRotation',30)
ylabel('Tumor purity')
end
%for with immunotherapy

if response_id>0
    PR1=0;
    PD1=0;
    CR1=0;
    PR2=0;
    PD2=0;
    CR2=0;
    PR3=0;
    PD3=0;
    CR3=0;
    PR4=0;
    PD4=0;
    CR4=0;
    clear i
    for i=1:length(high_high)
        re1=response{high_high(i)};
        if strcmp(re1,'"Partial Response"')
            PR1=PR1+1;
        elseif strcmp(re1,'"Progressive Disease"')
            PD1=PD1+1;
        elseif strcmp(re1,'"Complete Response"')
            CR1=CR1+1;
        end
    end
    clear i
    for i=1:length(low_low)
        re2=response{low_low(i)};
        if strcmp(re2,'"Partial Response"')
            PR2=PR2+1;
        elseif strcmp(re2,'"Progressive Disease"')
            PD2=PD2+1;
        elseif strcmp(re2,'"Complete Response"')
            CR2=CR2+1;
        end
    end
    clear i
    for i=1:length(lowP_highC)
        re3=response{lowP_highC(i)};
        if strcmp(re3,'"Partial Response"')
            PR3=PR3+1;
        elseif strcmp(re3,'"Progressive Disease"')
            PD3=PD3+1;
        elseif strcmp(re3,'"Complete Response"')
            CR3=CR3+1;
        end
    end
    clear i
    for i=1:length(highP_lowC)
        re4=response{highP_lowC(i)};
        if strcmp(re4,'"Partial Response"')
            PR4=PR4+1;
        elseif strcmp(re4,'"Progressive Disease"')
            PD4=PD4+1;
        elseif strcmp(re4,'"Complete Response"')
            CR4=CR4+1;
        end
    end
    figure
    set(gcf,'color','w')
    cat=categorical({'High PDL1 and High CD8+','Low PDL1 and Low CD8+','Low PDL1 and High CD8+','High PDL1 and Low CD8+'});
    bar(cat,[PR1 PD1 CR1; PR2 PD2 CR2; PR3 PD3 CR3; PR4 PD4 CR4],'stacked')
    ylabel('Number of occurrences')
    posx=[1 1 1 2 2 3 3 3 4 4 4];
    posy=[0.5 2 4 1 2.5 0.5 2.5 4.5 3 9.5 13.5];
    cellf={num2str(PR1/(PR1+PD1+CR1)), num2str(PD1/(PR1+PD1+CR1)),num2str(CR1/(PR1+PD1+CR1)),num2str(PR4/(PR4+PD4+CR4)),num2str(PD4/(PR4+PD4+CR4)),num2str(PR3/(PR3+PD3+CR3)),num2str(PD3/(PR3+PD3+CR3)),num2str(CR3/(PR3+PD3+CR3)),num2str(PR2/(PR2+PD2+CR2)),num2str(PD2/(PR2+PD2+CR2)),num2str(CR2/(PR2+PD2+CR2))};
    text(posx,posy,cellf,'HorizontalAlignment','center');
    legend('Partial Response','Progressive Disease','Complete Response','Location','northwest')
end




     