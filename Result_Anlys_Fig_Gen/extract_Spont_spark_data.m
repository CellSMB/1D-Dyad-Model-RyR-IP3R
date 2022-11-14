%%% Generates .mat file containing Amplitude + FDHM data on uncoupled dyad (UD)
% Number of Receptors (RyR+IP3R) in simulation
RyR_Total=30;
IP3R_Total=[0 10 20 40];

% Set n(runs) and time cut to analyse
frst=1;
Sim_Total=200;
Start_tSim=850; % ms
End_tSim=3000;  % ms
Start_tCut=850; % ms
End_tCut=3000;  % ms
Anlys_tRange=1:(End_tCut-Start_tSim)*100;
initiation_settings={'_tr13_no_trigger_at_all.mat'};

% Data processing parameters
RyR_Spark_Qualifier=5;
Min_Peak_Dist=10;       % ms
Cai_Min_Peak_Prom=100;
RyR_place=1;
IP3R_place=2;
Dyad_size=4;

% pre-allocating result matrices for n(event)
Cai_Num_UD=nan(Sim_Total-frst+1,1);         % matrix of n(spark events)
Cai_Amp_UD=cell(Sim_Total-frst+1,1);        % stores matrices of event amp
Cai_Amp_Indx_UD=cell(Sim_Total-frst+1,1);   % stores matrices of event amp index
Cai_FDHM_UD=cell(Sim_Total-frst+1,1);       % stores matrices of event duration

RyR_OpenNoSpark=zeros(Sim_Total-frst+1,1);

for num_IP3R=IP3R_Total
    for Sim_num=frst:Sim_Total
        if ismember(Sim_num,0)
            continue
        else
            % Load file and assign it to a struct
            file_name=['coupl_rep_uncoupl_mult_clusters_Cannell_Cao_ryrs_',num2str(RyR_Total),...
                        '_ip3rs_',num2str(num_IP3R),...
                        '_Tsim_',num2str(End_tSim),...
                        'smaller_set_test_v',num2str(Sim_num),...
                        initiation_settings{1}];
            file_struct=load(file_name);
            
            % finds midpoint of each dyad
            linescan_length=size(file_struct.Caf4_save,1);
            Dyad_1=(linescan_length/4);
            Dyad_2=linescan_length/2;
            Dyad_3=linescan_length/4*3;
            
            % slices data set based on analysis time range
            time=file_struct.time(:,Anlys_tRange);
            RyR_Open_data=sum(file_struct.RyR_open_save([Dyad_2-RyR_place Dyad_2+RyR_place],:),1); % sums open RyRs
            RyR_Open_data=RyR_Open_data(:,Anlys_tRange);
            
            % finds non-sparking RyR openings
            [RyR_Pks1,RyR_Locs1] = findpeaks(RyR_Open_data,time, ...
                'MinPeakDistance',Min_Peak_Dist, ...
                'MinPeakHeight',0);
            [~,RyR_Locs2] = findpeaks(RyR_Open_data,time, ...
                'MinPeakDistance',Min_Peak_Dist, ...
                'MinPeakHeight',RyR_Spark_Qualifier);
            [~,ia] = setdiff(RyR_Locs1,RyR_Locs2,'stable');
            RyR_OpenNoSpark(Sim_num-frst+1)=numel(RyR_Pks1(ia));
            
            % Analyses results for UNCOUPLED DYAD only
            % Gets [Ca2+] Result
            Cai_UD=file_struct.Cai_save(Dyad_2,:);
            Cai_UD=Cai_UD(:,Anlys_tRange);
            
            % findpeaks records amp, indx, and FDHM as row vectors
            % 1. Look for events where Cai peaks have the specified properties and its elemental index
            % 2. Find which of these events have RyR > 5 --> qualified as spark
            % 3. Look for recorded events and its **time** index
            % 4. From (2)&(3), record only events where RyR > 5
            
            [Cai_Amp,Cai_Amp_Indx,Cai_FDHM,~]=findpeaks(Cai_UD,time, ...
                'MinPeakDistance',Min_Peak_Dist, ...
                'MinPeakProminence',Cai_Min_Peak_Prom, ...
                'WidthReference','halfheight');
            [~,time_Indx]=ismember(Cai_Amp_Indx,time);
            Cai_Qualified_Spark_Indx=find(RyR_Open_data(time_Indx)>RyR_Spark_Qualifier);
            
            Cai_Amp=Cai_Amp(Cai_Qualified_Spark_Indx);
            Cai_Amp_Indx=Cai_Amp_Indx(Cai_Qualified_Spark_Indx);
            Cai_FDHM=Cai_FDHM(Cai_Qualified_Spark_Indx);
            time_Indx=time_Indx(Cai_Qualified_Spark_Indx);

            Cai_Num_UD(Sim_num-frst+1)=numel(Cai_Amp);
            Cai_Amp_UD{Sim_num-frst+1}=Cai_Amp;
            Cai_Amp_Indx_UD{Sim_num-frst+1}=Cai_Amp_Indx;
            Cai_FDHM_UD{Sim_num-frst+1}=Cai_FDHM;

            clear file_struct;
        end
    end
    save(['gryr002_preliminary_tr13_ip3r',num2str(num_IP3R),'_0_st_one_dyad.mat']);
end
disp('done')