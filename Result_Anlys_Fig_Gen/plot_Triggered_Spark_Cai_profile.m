%%% CONTROL PANEL - CHANGE THESE FOR DIFFERENT SIMULATIONS
RyR_Total=30;
trigger = 'LTCC';
remark = '';
initiation_settings={'_tr13_try_initiate_LTCC.mat'};...

% Domain Parameters
Sim_t=1070;         % ms
Sim_rec_t=950;      % ms
Sim_dt=0.01;        % ms
Sim_t_elements=(Sim_t-Sim_rec_t)/Sim_dt;
Sim_Total=200;
Sim_Width=8;        % um
Sim_dx=2/50;
ls_length=ceil(Sim_Width/Sim_dx);
Anlys_tCut=35/Sim_dt:Sim_t_elements-1000;
Cai_inset_t=2000;

% Data Processing Parameters
RyR_Spark_Qualifier=5;  % if n(RyR open) greater than ==> run had a spark
RyR_place=1;            % distance relative to central dyad element
IP3R_place=2;
Dyad_place=4;
Dyad_size=9;

for IP3R_Total=[0 5 10 20]*2
    % preallocates for graph plots (i.e. not linescan)
    RyR_open_data=nan(Sim_Total,Sim_t_elements);	% saves n(RyR) = open
    IP3R_open_data=nan(Sim_Total,Sim_t_elements);	% saves n(IP3R) = open
    CaJSR_at_RyR_loc=nan(2*Sim_Total,Sim_t_elements);   % saves [Ca2+]_jSR at each RyR location
    CaJSR_at_IP3R_loc=nan(2*Sim_Total,Sim_t_elements);  % saves [Ca2+]_jSR at each IP3R location
    CaJSR_at_mid_loc=nan(Sim_Total,Sim_t_elements);
    Cai_at_RyR_loc=nan(2*Sim_Total,Sim_t_elements);
    Cai_at_IP3R_loc=nan(2*Sim_Total,Sim_t_elements);
    Cai_at_mid_loc=nan(Sim_Total,Sim_t_elements);

    spark_count=0;    % counts n(sim) which max_val > sparking_run, used to determine mean F/F0

    for Sim_num=1:Sim_Total
        if Sim_num==0
            continue
        else
            file_name=['one_d_Cannell_Cao_ryrs_',num2str(RyR_Total),...
            '_ip3rs_',num2str(IP3R_Total),...
            '_Tsim_',num2str(Sim_t),...
            'smaller_set_test_v',num2str(Sim_num),...
            initiation_settings{1}];
            load(file_name);

            mid_cell=ls_length/2;

            max_RyR_open=max(sum(RyR_open_save(mid_cell-Dyad_place:mid_cell+Dyad_place,:)));

            if max_RyR_open>RyR_Spark_Qualifier
                RyR_open_data(Sim_num,:)=sum(RyR_open_save(mid_cell-Dyad_place:mid_cell+Dyad_place,:));
                IP3R_open_data(Sim_num,:)=sum(IP3R_open_save(mid_cell-Dyad_place:mid_cell+Dyad_place,:));
                
                CaJSR_at_RyR_loc(Sim_num*2-1,:)=Cajsr_save(mid_cell-RyR_place,:);
                CaJSR_at_RyR_loc(Sim_num*2,:)=Cajsr_save(mid_cell+RyR_place,:);
                CaJSR_at_IP3R_loc(Sim_num*2-1,:)=Cajsr_save(mid_cell-IP3R_place,:);
                CaJSR_at_IP3R_loc(Sim_num*2,:)=Cajsr_save(mid_cell+IP3R_place,:);
                CaJSR_at_mid_loc(Sim_num,:)=Cajsr_save(mid_cell,:);

                Cai_at_RyR_loc(Sim_num*2-1,:)=Cai_save(mid_cell-RyR_place,:);
                Cai_at_RyR_loc(Sim_num*2,:)=Cai_save(mid_cell+RyR_place,:);
                Cai_at_IP3R_loc(Sim_num*2-1,:)=Cai_save(mid_cell-IP3R_place,:);
                Cai_at_IP3R_loc(Sim_num*2,:)=Cai_save(mid_cell+IP3R_place,:);
                Cai_at_mid_loc(Sim_num,:)=Cai_save(mid_cell,:);

                spark_count=spark_count+1;
            end
        end
    end
    CaJSR_at_dyad=[CaJSR_at_RyR_loc;CaJSR_at_IP3R_loc;CaJSR_at_mid_loc];
    Cai_at_dyad=[Cai_at_RyR_loc;Cai_at_IP3R_loc;Cai_at_mid_loc];
    
    fprintf('Total Spark Count = %g for n(IP_3R) = %g\n',spark_count,IP3R_Total/2);

    % Figure Control Panel: Figures to be plotted:
    color_transparency=0.3;
    LineWidth=1;
    LineStyle='none';
    FontName='Sans Serif';
    TickFontSize=6;  %get(groot, 'defaultAxesFontSize');
    LabelFontSize=7;
    TickLength=[0 0];
    PlotBoxAR=[1 1 1];
    FigSize=[0 0 3.98 3.25];
    time=time(:,Anlys_tCut);
    FigNum = 1;

    figure(FigNum);
    time2=[time,flip(time)];
    time3=[time(:,end-Cai_inset_t:end),flip(time(:,end-Cai_inset_t:end))];
    [ms1,conf1]=data_shade(Cai_at_mid_loc(:,Anlys_tCut),1);
    hold 'on'
    s1=fill(time2,[conf1(1,:),conf1(2,:)],lines(1),'LineStyle',LineStyle);
    alpha(s1,color_transparency)
    plot(time,ms1,'Color',lines(1),'LineWidth',LineWidth)
    hold 'off'
    xlim([time(1),time(end)])
    ylim([0 400])
    set(gca,'FontSize',TickFontSize,'FontName',FontName)
    xlabel('Time (ms)','FontSize',LabelFontSize)
    ylabel('Dyad Ca^{2+} (\muM)','FontSize',LabelFontSize)
    grid off
    box off
    MainAxPos=get(gca,'Position');
    inset=axes('Position',[MainAxPos(1)+0.4 MainAxPos(2)+0.4 MainAxPos(3)-0.5 MainAxPos(4)-0.5]);
    hold(inset,'on')
    s1=fill(time3,[conf1(1,end-Cai_inset_t:end),conf1(2,1:1+Cai_inset_t)],lines(1),'LineStyle',LineStyle);
    alpha(s1,color_transparency)
    plot(inset,time(:,end-Cai_inset_t:end),ms1(:,end-Cai_inset_t:end),'Color',lines(1),'LineWidth',LineWidth);
    set(inset,'XLim',[Sim_t-20 Sim_t],'YLim',[0 4]);
    set(gca,'FontSize',TickFontSize-1,'FontName',FontName)
    xlabel('Time (ms)','FontSize',LabelFontSize-1)
    ylabel('Dyad Ca^{2+} (\muM)','FontSize',LabelFontSize-1)
    pbaspect(inset,PlotBoxAR);
    hold(inset,'off')
    grid off
    box off
    set(gcf,'color','w','units','centimeters','position',FigSize);
    saveas(gcf,['one_dyad_',trigger,'_',num2str(IP3R_Total),'IP3R_',remark,'plots_Ca_results.fig'])
    FigNum=FigNum+1;

    figure(FigNum);
    time2=[time,flip(time)];
    [ms1,conf1]=data_shade(CaJSR_at_dyad(:,Anlys_tCut),1);
    s1=fill(time2,[conf1(1,:),conf1(2,:)],lines(1),'LineStyle',LineStyle);
    alpha(s1,color_transparency)
    hold on
    plot(time,ms1,'Color',lines(1),'LineWidth',LineWidth)
    hold off
    xlim([time(1) time(end)])
    ylim([0 par.casr_0])
    set(gca,'FontSize',TickFontSize,'FontName',FontName)
    xlabel('Time (ms)','FontSize',LabelFontSize)
    ylabel('JSR Ca^{2+} (\muM)','FontSize',LabelFontSize)
    grid off
    box off
    set(gcf,'color','w','units','centimeters','position',FigSize);
    saveas(gcf,['one_dyad_',trigger,'_',num2str(IP3R_Total),'IP3R_',remark,'plots_JSR_results.fig'])
    FigNum=FigNum+1;

    figure(3);
    [ms1,conf1]=data_shade(RyR_open_data(:,Anlys_tCut),1);
    s1=fill(time2,[conf1(1,:),conf1(2,:)],lines(1),'LineStyle',LineStyle);
    alpha(s1,color_transparency)
    hold on
    plot(time,ms1,'Color',lines(1),'LineWidth',LineWidth)
    hold off
    xlim([time(1) time(end)])
    ylim([0 RyR_Total-10])
    set(gca,'FontSize',TickFontSize,'FontName',FontName)
    xlabel('Time (ms)','FontSize',LabelFontSize)
    ylabel('No. open RyRs','FontSize',LabelFontSize)
    grid off
    box off
    set(gcf,'color','w','units','centimeters','position',FigSize);
    saveas(gcf,['one_dyad_',trigger,'_',num2str(IP3R_Total),'IP3R_',remark,'plots_RyR_results.fig'])
    FigNum=FigNum+1;

    if IP3R_Total>0
        figure(FigNum);
        [ms1,conf1]=data_shade(IP3R_open_data(:,Anlys_tCut),1);
        s1=fill(time2,[conf1(1,:),conf1(2,:)],lines(1),'LineStyle',LineStyle);
        alpha(s1,color_transparency)
        hold on
        plot(time,ms1,'Color',lines(1),'LineWidth',LineWidth)
        hold off
        xlim([time(1) time(end)])
        ylim([0 4])
        set(gca,'FontSize',TickFontSize,'FontName',FontName)
        xlabel('Time (ms)','FontSize',LabelFontSize)
        ylabel('No. open IP_3Rs','FontSize',LabelFontSize)
        grid off
        box off
        set(gcf,'color','w','units','centimeters','position',FigSize);
        saveas(gcf,['one_dyad_',trigger,'_',num2str(IP3R_Total),'IP3R_',remark,'plots_IP3R_results.fig'])
    end
    close all   % closes all figures
end