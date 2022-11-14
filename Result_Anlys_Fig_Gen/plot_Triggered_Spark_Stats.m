%%% Plots Spark FDHM + Amp (Swarm Plot) & %Sim with Sparks (Bar Chart)
% Simulation Condition (Sim Cond): 1 Dyad 0 Trigger

st_t={'_0'};
IP3R_Names=[0 10 20 40];
SimAnlys_Start_t=950;
SimAnlys_End_t=1070;
nSim_Total=50;
Output_file_remark='';

%%% Pre-allocate result matrix for analyses
%%% Row: Sim Cond, Col: n(IP3R)
% Stores FDHM, Amp, & Time of ALL sparks
Cai_FDHM_Result=cell(numel(st_t),numel(IP3R_Names));
Cai_Amp_Result=cell(numel(st_t),numel(IP3R_Names));
Cai_Amp_Indx_Result=cell(numel(st_t),numel(IP3R_Names));

% Stores n(runs) with sparks
nSparks_in_Sim=cell(numel(st_t),numel(IP3R_Names));

%%% Algorithm to extract and store data in main result matrix
%%% loop each stimulation case, then each n(IP3R)
for TriggerType=1:numel(st_t) % 3 stimulation types
    for IP3R_Indx=1:numel(IP3R_Names)
        file_struct=load(['gryr002_preliminary_tr13_ip3r',num2str(IP3R_Names(IP3R_Indx)),st_t{TriggerType},'_st_one_dyad.mat']);
        
        % removes n(sparks) that occur outside of t(analysis)
        nSpark_remove=cellfun(@(t) sum(t>SimAnlys_End_t|t<SimAnlys_Start_t),file_struct.Cai_Amp_Indx_UD);
        Cai_Num_dummy=file_struct.Cai_Num_UD-nSpark_remove;
        
        % cellfun turns row vectors in cell --> col vectors so cell2mat can work!
        Cai_FDHM_dummy=cell2mat(cellfun(@transpose,file_struct.Cai_FDHM_UD,'UniformOutput',false));
        Cai_Amp_dummy=cell2mat(cellfun(@transpose,file_struct.Cai_Amp_UD,'UniformOutput',false));
        Cai_Amp_Indx_dummy=cell2mat(cellfun(@transpose,file_struct.Cai_Amp_Indx_UD,'UniformOutput',false));

        % chooses in between where we wanna analyse
        Anlys_t_Rng=((Cai_Amp_Indx_dummy>=SimAnlys_Start_t) & (Cai_Amp_Indx_dummy<=SimAnlys_End_t));
        
        % determines count of unique n(sparks) in integers
        [uniq_nSpark_count,uniq_nSpark,~]=histcounts(Cai_Num_dummy,'BinMethod','integers');
        
        % Saves each scenario into main result cell
        nSparks_in_Sim{TriggerType,IP3R_Indx}=Cai_Num_dummy;

        Cai_FDHM_Result{TriggerType,IP3R_Indx}=Cai_FDHM_dummy(Anlys_t_Rng);
        Cai_Amp_Result{TriggerType,IP3R_Indx}=Cai_Amp_dummy(Anlys_t_Rng);
        Cai_Amp_Indx_Result{TriggerType,IP3R_Indx}=Cai_Amp_Indx_dummy(Anlys_t_Rng);

    end
end
clear file_struct;

%%% Plot Generation & its required sub-analyses
Plot_Titles={'1 Dyad 0 Trigger'};
bar_xNames1={'0 IP_3R','5 IP_3R','10 IP_3R','20 IP_3R'};
color_transparency=0.3;
LineWidth=1;
LineStyle='none';
FontName='Sans Serif';
SigstarFontSize=7;
TickFontSize=9;
LabelFontSize=10;
TickLength=[0 0];
PlotBoxAR=[1 1 1];
FigureSize=[0 0 3.98 3.25]*2;
FigNum=1;

% Swarm plot of Cai Amplitude
figure(FigNum);
plotSpread(Cai_Amp_Result,'xNames',bar_xNames1,'showMM',4,'distributionColors',lines(1))
ylim([0 Inf])
xtickangle(30)
set(gca,'FontSize',TickFontSize,'FontName',FontName)
ylabel('Ca^{2+} spark amplitude (\muM)','FontSize',LabelFontSize)
grid off
box off
set(gcf,'color','w','units','centimeters','position',FigureSize);
saveas(gcf,['one_dyad_Cai_Amp_swarm_', Output_file_remark, 'plots.fig'])
FigNum=FigNum+1;

% Swarm plot of Cai FDHM for each sim cond
figure(FigNum);
plotSpread(Cai_FDHM_Result,'xNames',bar_xNames1,'showMM',4,'distributionColors',lines(1))
ylim([0 10])
xtickangle(30)
set(gca,'FontSize',TickFontSize,'FontName',FontName)
ylabel('Ca^{2+} spark FDHM (ms)','FontSize',LabelFontSize)
grid off
box off
set(gcf,'color','w','units','centimeters','position',FigureSize);
saveas(gcf,['one_dyad_Cai_FDHM_swarm_', Output_file_remark, 'plots.fig'])

%%% Statistical Testing
% reshapes all main result data into a column vector of cells
Cai_FDHM_Result_1Col=reshape(Cai_FDHM_Result',[],1);
Cai_Amp_Result_1Col=reshape(Cai_Amp_Result',[],1);
nSpark_Result_1Col=reshape(nSparks_in_Sim',[],1);

% converts all response data (FDHM, Amp) into a row matrix
Cai_FDHM_Response=cell2mat(Cai_FDHM_Result_1Col)';
Cai_Amp_Response=cell2mat(Cai_Amp_Result_1Col)';
nSpark_Response=cell2mat(nSpark_Result_1Col)';

% repeat names of numIP3R for Amp and FDHM
NumIP3R_Rep_SimCond1=num2cell(repmat((IP3R_Names/2),1,numel(st_t))');% cells containing doubles of numIP3R
NumIP3R_Rep_Event1=num2cell(cellfun(@numel,Cai_FDHM_Result_1Col));    % gets numel of each FDHM based on numIP3R & simcond
FDHM_Amp_IP3R_Num=cell2mat(cellfun(@repelem,NumIP3R_Rep_SimCond1,NumIP3R_Rep_Event1,'UniformOutput',false)'); % row matrix of repeated numIP3R

% Using result row matrices above to do 1 way ANOVA and Tukey Post Hoc
[p_FDHM,tbl_FDHM,stats_FDHM] = anova1(Cai_FDHM_Response,FDHM_Amp_IP3R_Num);
figure;
[c_FDHM,m_FDHM,h_FDHM,nms_FDHM] = multcompare(stats_FDHM,'Dimension',1);
ylabel('n(IP_3R)')
saveas(gcf,['one_dyad_FDHM_PostHoc_', Output_file_remark, 'plots.fig'])

[p_Amp,tbl_Amp,stats_Amp] = anova1(Cai_Amp_Response,FDHM_Amp_IP3R_Num);
figure;
[c_Amp,m_Amp,h_Amp,nms_Amp] = multcompare(stats_Amp,'Dimension',1);
ylabel('n(IP_3R)')
saveas(gcf,['one_dyad_Amp_PostHoc_', Output_file_remark, 'plots.fig'])

close all   % close all figures

%%% Plotting significance bars for Amplitude swarm plot
open(['one_dyad_Cai_Amp_swarm_', Output_file_remark, 'plots.fig']);
p_Threshold=[0.05 0.01 0.001 0.0001];
xTick=find(ismember(IP3R_Names/2,cellfun(@str2num,nms_Amp)));
xTick=repelem(xTick,2);
xTick=xTick(2:end-1);
xOffset=0.025;
xSigLine=xTick+(-1).^((1:numel(xTick))+1)*xOffset;
hold on
for i=1:numel(nms_Amp)-1
    RowIndx=find(ismember(c_Amp(:,1:2),[i i+1],'rows'));
    if c_Amp(RowIndx,end) < max(p_Threshold)
        plot(xSigLine([2*i-1 2*i]),[1 1]*max(vertcat(Cai_Amp_Result{xTick(2*i-1:2*i)}))*1.05,'-k');%,mean(sigx([2*i-1 2*i])),max(Cai_Amp_Response)*1.1,'*k')
        plot(xSigLine([2*i-1 2*i-1]),max(vertcat(Cai_Amp_Result{xTick(2*i-1:2*i)}))*[1.025 1.05],'-k')
        plot(xSigLine([2*i 2*i]),max(vertcat(Cai_Amp_Result{xTick(2*i-1:2*i)}))*[1.025 1.05],'-k')

        if c_Amp(RowIndx,end) < min(p_Threshold)
            text(mean(xSigLine(2*i-1:2*i)),max(vertcat(Cai_Amp_Result{xTick(2*i-1:2*i)}))*1.1,['p<',num2str(min(p_Threshold),'%.4f')],...
                HorizontalAlignment="center",FontName='Sans Serif',FontSize=SigstarFontSize)
        else
            text(mean(xSigLine(2*i-1:2*i)),max(vertcat(Cai_Amp_Result{xTick(2*i-1:2*i)}))*1.1,['p=',num2str(c_Amp(RowIndx,end),'%.4f')],...
                HorizontalAlignment="center",FontName='Sans Serif',FontSize=SigstarFontSize)
        end
    end
end
hold off
ylim([0 ceil(max(Cai_Amp_Response)*1.1/10)*10+5])
saveas(gcf,['one_dyad_Cai_Amp_swarm_', Output_file_remark, 'plots.fig'])

close all
disp('done!')
