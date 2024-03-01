%%% Plots the mean trace of [Ca2+] of all 200 simulation files
%%% CONTROL PANEL - CHANGE THESE FOR DIFFERENT SIMULATIONS
Sim_Total=200;
RyR_Total=30;
vary_rslt_name={'Young','Old'};
trigger = 'LTCC';
remark = 'all';
initiation_settings={'_tr13_try_initiate_LTCC.mat'};...

% Domain Parameters
Sim_t=1100;         % ms
Sim_rec_t=950;      % change btw 850 and 950 ms
Sim_dt=0.01;        % ms
Sim_t_elements=(Sim_t-Sim_rec_t)/Sim_dt;
Sim_Width=8;    % um
Sim_dx=2/50;
ls_length=ceil(Sim_Width/Sim_dx);
ls_xlabel=(1:ls_length)*Sim_dx;   % rescales linescan from indices --> um
t_cut=35/Sim_dt:Sim_t_elements-1000;
Cai_inset_t=2000;

% Data Processing Parameters
RyR_Spark_Qualifier=5; % if max_RyR_open value is higher than this, we say run had a spark
Min_Peak_Dist_t=(Sim_t-Sim_rec_t)-1;
Min_Peak_Dist_x=Sim_Width-0.3;
Min_Peak_Prom=0.06;
F_min=1;
F_max=100;
RyR_place=1;  % position of RyR cluster relative to RyR_release_place
IP3R_place=2;
Dyad_place=4;
Dyad_size=9;
Cai_Elem_Avg=IP3R_place;

% PSF Generation
FWHM = 0.41; % um
sigma = FWHM./(2*sqrt(2*log(2))); % getting SD from FWHM
coord=ls_xlabel(:);
mu = mean(coord);   % mean
S=diag(sigma);
PSF_compressed=mvnpdf(coord,mu,S);
PSF_1D=PSF_compressed./sum(PSF_compressed);

for Age_Grp=1:numel(vary_rslt_name)
    % preallocates for graph plots (i.e. not linescan)
    RyR_open_data=nan(Sim_Total,Sim_t_elements);	% saves n(RyR) = open
    CaJSR_at_RyR_loc=nan(2*Sim_Total,Sim_t_elements);   % saves [Ca2+]_jSR at each RyR location
    CaJSR_at_IP3R_loc=nan(2*Sim_Total,Sim_t_elements);  % saves [Ca2+]_jSR at each IP3R location
    CaJSR_at_mid_loc=nan(Sim_Total,Sim_t_elements);
    Cai_at_RyR_loc=nan(2*Sim_Total,Sim_t_elements);
    Cai_at_IP3R_loc=nan(2*Sim_Total,Sim_t_elements);
    Cai_at_mid_loc=nan(Sim_Total,Sim_t_elements);

    spark_count=0;    % counts n(sim) which max_val > sparking_run, used to determine mean F/F0

    for aa=1:Sim_Total
        if aa==0
            continue
        else
            file_name=['one_d_Cannell_Cao_Sheep_',...
                        vary_rslt_name{Age_Grp},...
                        '_ryrs_',num2str(RyR_Total),...
                        '_ip3rs_0_Tsim_',num2str(Sim_t),...
                        'smaller_set_test_v',num2str(aa),...
                        initiation_settings{1}];
            load(file_name);

            mid_cell=ls_length/2;
            max_RyR_open=max(sum(RyR_open_save(mid_cell-Dyad_place:mid_cell+Dyad_place,:)));

            if max_RyR_open>RyR_Spark_Qualifier
                RyR_open_data(aa,:)=sum(RyR_open_save(mid_cell-Dyad_place:mid_cell+Dyad_place,:));

                CaJSR_at_RyR_loc(aa*2-1,:)=Cajsr_save(mid_cell-RyR_place,:);
                CaJSR_at_RyR_loc(aa*2,:)=Cajsr_save(mid_cell+RyR_place,:);
                CaJSR_at_IP3R_loc(aa*2-1,:)=Cajsr_save(mid_cell-IP3R_place,:);
                CaJSR_at_IP3R_loc(aa*2,:)=Cajsr_save(mid_cell+IP3R_place,:);
                CaJSR_at_mid_loc(aa,:)=Cajsr_save(mid_cell,:);

                Cai_at_RyR_loc(aa*2-1,:)=Cai_save(mid_cell-RyR_place,:);
                Cai_at_RyR_loc(aa*2,:)=Cai_save(mid_cell+RyR_place,:);
                Cai_at_IP3R_loc(aa*2-1,:)=Cai_save(mid_cell-IP3R_place,:);
                Cai_at_IP3R_loc(aa*2,:)=Cai_save(mid_cell+IP3R_place,:);
                Cai_at_mid_loc(aa,:)=Cai_save(mid_cell,:);
                
                spark_count=spark_count+1;
            end
        end
    end
    CaJSR_at_dyad=[CaJSR_at_RyR_loc;CaJSR_at_IP3R_loc;CaJSR_at_mid_loc];
    Cai_at_dyad=[Cai_at_RyR_loc;Cai_at_IP3R_loc;Cai_at_mid_loc];
    
    % xlswrite(['Data_Ca_JSR_',vary_rslt_name{Age_Grp},'.csv'],CaJSR_at_dyad)
    % xlswrite(['Data_Ca_Dyad_',vary_rslt_name{Age_Grp},'.csv'],Cai_at_mid_loc)
    % xlswrite(['Data_No_OpenRyR_',vary_rslt_name{Age_Grp},'.csv'],RyR_open_data)
    % xlswrite(['Data_time_',vary_rslt_name{Age_Grp},'.csv'],time)
    % fprintf(['Total Spark Count = %g for cohort = ',vary_rslt_name{Age_Grp},'\n'],spark_count);

    % Figure Control Panel: Figures to be plotted:
    ColorTransparency=0.3;
    LineWidth=2;
    LineStyle='none';
    FontName='Sans Serif';
    TickFontSize=6;
    LabelFontSize=7;
%     ColorMap={'#009DA3','#F98F00'}; % #52986C,#A49336
    ColorOrder=[0 157 163;82 152 108;164 147 54;249 143 0]/255;
    ColorOrderIndx=linspace(1,4,numel(vary_rslt_name));
    ColorOrder=ColorOrder(ColorOrderIndx,:);
    time=time(:,t_cut)/1000;
    
    figure(1);
%     time2=[time flip(time)];
    [ms1,conf1]=data_shade(Cai_at_mid_loc(:,t_cut),1);
    ms1=ms1/1000;
%     conf1=conf1/1000;
%     s1=fill(time2,[conf1(1,:),conf1(2,:)],ColorOrder(Age_Grp,:),'LineStyle',LineStyle);
%     alpha(s1,ColorTransparency)
    hold on
    plot(time,ms1,'Color',ColorOrder(Age_Grp,:),'LineWidth',LineWidth)
    if Age_Grp==numel(vary_rslt_name)
        hold off
        xlim([time(1),time(end)])
        ylim([0 0.5])
        set(gca,'FontSize',TickFontSize,'FontName',FontName)
        xlabel('Time (s)','FontSize',LabelFontSize)
        ylabel('Dyad Ca^{2+} (mM)','FontSize',LabelFontSize)
        grid off
        box off
        set(gcf,'color','w','units','centimeters','position',[0 0 4 4.5]);
        saveas(gcf,['one_dyad_',trigger,'_Sheep_RyR_',remark,'plots_Ca_results.fig'])
    end
    
    figure(2);
    [ms1,conf1]=data_shade(CaJSR_at_dyad(:,t_cut),1);
    ms1=ms1/1000;
%     conf1=conf1/1000;
%     s1=fill(time2,[conf1(1,:),conf1(2,:)],ColorOrder(Age_Grp,:),'LineStyle',LineStyle);
%     alpha(s1,ColorTransparency)
    hold on
    plot(time,ms1,'Color',ColorOrder(Age_Grp,:),'LineWidth',LineWidth)
    if Age_Grp==numel(vary_rslt_name)
        hold off
        xlim([time(1) time(end)])
        ylim([0 1])
        set(gca,'FontSize',TickFontSize,'FontName',FontName)
        xlabel('Time (s)','FontSize',LabelFontSize)
        ylabel('JSR Ca^{2+} (mM)','FontSize',LabelFontSize)
        grid off
        box off
        set(gcf,'color','w','units','centimeters','position',[0 0 4 4.5]);
        saveas(gcf,['one_dyad_',trigger,'_Sheep_RyR_',remark,'plots_JSR_results.fig'])
    end
    
    figure(3);
    [ms1,conf1]=data_shade(RyR_open_data(:,t_cut),1);
%     s1=fill(time2,[conf1(1,:),conf1(2,:)],ColorOrder(Age_Grp,:),'LineStyle',LineStyle);
%     alpha(s1,ColorTransparency)
    hold on
    plot(time,ms1,'Color',ColorOrder(Age_Grp,:),'LineWidth',LineWidth)
    if Age_Grp==numel(vary_rslt_name)
        hold off
        xlim([time(1) time(end)])
        ylim([0 RyR_Total-0])
        set(gca,'FontSize',TickFontSize,'FontName',FontName)
        xlabel('Time (s)','FontSize',LabelFontSize)
        ylabel('N. RyRs open','FontSize',LabelFontSize)
        grid off
        box off
        set(gcf,'color','w','units','centimeters','position',[0 0 4 4.5]);
        saveas(gcf,['one_dyad_',trigger,'_Sheep_RyR_',remark,'plots_RyR_results.fig'])
    end
end