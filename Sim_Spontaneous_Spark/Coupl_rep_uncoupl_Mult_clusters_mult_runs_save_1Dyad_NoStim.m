%%% Agne Tilunaite 2019
%%% Edited for clarity by Hilary Hunt 2020
%%% Modified by Joshua Chung 2021
numSimulations=200;
simTime=3000;       % ms
triggerTime=1000;   % ms
rec_startTime=850;  % ms
simWidth=8;         % um
c0=0.1;             % uM
casr0=1000*ones(numSimulations,1);  % uM
IP3Conc=0.15;                       % uM
vary_nIP3R=[0 5 10 20];             % n(IP3R)

for run=1:numel(vary_nIP3R)
    nIP3R=vary_nIP3R(run);
    parfor(fop=1:numel(casr0),12)
        %%% Get RyR simulation parameters & overwrites/adds parameter values in par
        par=RxnDiffusion_parameters(c0,casr0(fop));
        par.casr_0=casr0(fop);
        par.sc=1000;
        par.scale_sercas=1;
        par.gryr=2.8;
        par.grefill=0.2;
        par.dt0=1e-4; %ms
        par.dx=0.2/5; %um
        
        %%% Count number of Receptors (IP3R + RyR)
        total_RyRs=30;
        total_IP3Rs=2*nIP3R;
        total_Rs=total_RyRs+total_IP3Rs;
        
        %%% Get/initialises IP3R simulation parameters & overwrites/adds parameter values in ip3r_par
        ip3r_par=IP3R_parameters(c0,IP3Conc);
        ip3r_par.kipr=par.gryr/2.85;    % 2.85 scaling from Foskett et al
        
        %%% Calculates initial values of h_ij and m_ij for IP3R gating at t=0
        h24_inf=ip3r_par.kn24^ip3r_par.nn24./(ip3r_par.kn24^ip3r_par.nn24+(c0).^ip3r_par.nn24);
        h24=h24_inf*ones(total_IP3Rs,1);
        
        m24_inf=(c0).^ip3r_par.n24./(ip3r_par.k24^ip3r_par.n24+(c0).^ip3r_par.n24);
        m24=m24_inf*ones(total_IP3Rs,1);
        
        m42_inf=(c0).^ip3r_par.n42./(ip3r_par.k42.^ip3r_par.n42+(c0).^ip3r_par.n42);
        m42=m42_inf*ones(total_IP3Rs,1);
        
        h42_inf=ip3r_par.kn42^ip3r_par.nn42./(ip3r_par.kn42^ip3r_par.nn42+(c0).^ip3r_par.nn42);
        h42=h42_inf*ones(total_IP3Rs,1);
        
        %%% Set Receptor location in 1D domain
        linescan_length=ceil(simWidth/par.dx);
        quart_cell=ceil(linescan_length/4*1);
        mid_cell=ceil(linescan_length/4*2);
        tquart_cell=ceil(linescan_length/4*3);

        D2_RyR_Cluster1=mid_cell-1;
        D2_RyR_Cluster2=mid_cell+1;
        D2_IP3R_Cluster1=mid_cell-2;
        D2_IP3R_Cluster2=mid_cell+2;

        par.Dyad1=quart_cell-4:quart_cell+4;
        par.Dyad2=mid_cell-4:mid_cell+4;
        par.Dyad3=tquart_cell-4:tquart_cell+4;

        Dyad_region=zeros(linescan_length,1);
        Dyad_region(par.Dyad2)=1;
        
        %%% Used to set TnC location in 1D domain
        par.positions=Dyad_region;

        %%% Set SERCA location in 1D domain
        par.sercas=ones(linescan_length,1).*(~par.positions);

        %%% Explanation: RyR_sparse = Receptor (RyR + IP3R) storage matrix
        % 1st column - Receptor position in 1D line
        % 2nd column - Receptor state (CLOSE = 1, OPEN = 2)
        % 3rd column - IP3R state (C1 C2 C3 C4 O5 O6)***, RyR follows column 2
        % 4th column - cumulative k*dt until Receptor changes state
        % 5th column - random number to determine Receptor transition times
        RyR_sparse=zeros(total_Rs,5);
        
        % Set Receptor location in RyR_sparse
        div=total_RyRs/2;
        RyR_sparse(1:div,1)=D2_RyR_Cluster1;
        RyR_sparse(div+1:2*div,1)=D2_RyR_Cluster2;
        RyR_sparse(total_RyRs+1:total_RyRs+nIP3R,1)=D2_IP3R_Cluster1;
        RyR_sparse(total_RyRs+1+nIP3R:total_Rs,1)=D2_IP3R_Cluster2;
        
        % Initialise Receptor (RyRs & IP3Rs) State in RyR_sparse
        RyRVector = (1:total_RyRs);
        IP3RVector = (1+total_RyRs:total_Rs);
        RyR_sparse(RyRVector,3)=1;
        RyR_sparse(IP3RVector,3)=4; 
        RyR_sparse(:,2)=1;          % simplified R state
        Ch_mx=[-1;1];               % used to change R state

        % Initialise Random Number between 0 and 1
        RyR_sparse(:,5)=rand(total_Rs,1);
        
        %%% Count Receptors (RyRs + IP3Rs) in OPEN state
        RyR_open_sums=zeros(linescan_length,1);
        RyR_open_sums(D2_RyR_Cluster1)=sum(RyR_sparse(1:div,2)==2);
        RyR_open_sums(D2_RyR_Cluster2)=sum(RyR_sparse(div+1:2*div,2)==2);
        
        IP3R_open_sums=zeros(linescan_length,1);
        IP3R_open_sums(D2_IP3R_Cluster1)=sum(RyR_sparse(total_RyRs+1:total_RyRs+nIP3R,2)==2);
        IP3R_open_sums(D2_IP3R_Cluster2)=sum(RyR_sparse(total_RyRs+1+nIP3R:total_Rs,2)==2);

        %%% Pre-allocating vectors to track RD species & save for analysis
        % Variables I will need to track - RyR & SR related regions
        % Y0=[Cai, Cajsr, Cansr, CaCaM, CaATP, F4Ca, CaTnC, CaCSQ];
        %      1     2      3      4      5     6      7      8
        save_t_int=0.01;
        titer_total=(simTime)/save_t_int;
        titer_save=(simTime-rec_startTime)/save_t_int;
        RyR_open_save=zeros(linescan_length,titer_save);
        IP3R_open_save=zeros(linescan_length,titer_save);
        Cai_save=zeros(linescan_length,titer_save);
        CaJSR_save=zeros(linescan_length,titer_save);
        CaF4_save=zeros(linescan_length,titer_save);
        CaATP_save=zeros(linescan_length,titer_save);
        CaCaM_save=zeros(linescan_length,titer_save);
        CaTnC_save=zeros(linescan_length,titer_save);
        CaNSR_save=zeros(linescan_length,titer_save);
        countifs_save=zeros(1,titer_save);
        time=zeros(1,titer_save);

        %%% Initialise [Ca2+] in Dyad, jSR, nSR
        Cai=c0.*ones(linescan_length,1);
        CaJSR=par.casr_0.*par.positions;
        CaNSR=par.casr_0*ones(linescan_length,1);

        %%% Initialise [Buffer] (at *eqbm* with initial free Ca2+)
        CaM0=par.koffCaM.*par.CaMtotal/(par.koffCaM+par.konCaM*c0);
        ATP0=par.koffATP.*par.ATPtotal/(par.koffATP+par.konATP*c0);
        F40=par.koffF4.*par.F4total/(par.koffF4+par.konF4*c0);
        TnC0=par.koffTnC*par.TnCtotal/(par.koffTnC+par.konTnC*c0);
        CSQ0=par.koffCSQ*par.CSQtotal/(par.koffCSQ+par.konCSQ*par.casr_0);

        CaCaM=(par.CaMtotal-CaM0)*ones(linescan_length,1);
        CaATP=(par.ATPtotal-ATP0)*ones(linescan_length,1);
        CaF4=(par.F4total-F40)*ones(linescan_length,1);
        CaTnC=(par.TnCtotal-TnC0)*ones(linescan_length,1).*~par.positions;
        CaCSQ=(par.CSQtotal-CSQ0).*par.positions;
        
        %%% start timer, initialise dt and sim time
        tic
        par.dt=par.dt0;
        simt=0;
        old_indx=0;

        %%% LTCC pulse length
        LTCC_t=2;     % ms

        for iter=1:titer_total
            dumt=0;

            while dumt<=(par.dt0*100)
                %%% System perturbation settings
                if (simt+dumt)==0
                    % open RyRs at t=0 to allow faster system equilibration
                    num_RyR_open_trigger=5;
                    rand_open_D2C1RyR=randsample((1:div),num_RyR_open_trigger);
                    rand_open_D2C2RyR=randsample((div+1:2*div),num_RyR_open_trigger);
                    RyR_sparse(rand_open_D2C1RyR,2)=2; % state of each receptor
                    RyR_sparse(rand_open_D2C2RyR,2)=2; % state of each receptor
                    RyR_sparse(rand_open_D2C1RyR,3)=2;
                    RyR_sparse(rand_open_D2C2RyR,3)=2;
                end
                
                %%% counts open RyR & IP3R for each iteration
                RyR_open_sums(D2_RyR_Cluster1)=sum(RyR_sparse(1:div,2)==2);
                RyR_open_sums(D2_RyR_Cluster2)=sum(RyR_sparse(div+1:2*div,2)==2);
                IP3R_open_sums(D2_IP3R_Cluster1)=sum(RyR_sparse(total_RyRs+1:total_RyRs+nIP3R,2)==2);
                IP3R_open_sums(D2_IP3R_Cluster2)=sum(RyR_sparse(total_RyRs+1+nIP3R:end,2)==2);
                
                %%% Solves RD Species from t=t --> t=t+dt using RK1
                Yn=[Cai, CaJSR, CaNSR, CaCaM, CaATP, CaF4, CaTnC, CaCSQ];
                Yn1=EulerMethod(Yn,par,ip3r_par,RyR_open_sums,IP3R_open_sums);             % values at the time step n+1
                
                %%% Calculates/Updates Receptor (RyR + IP3R) gating kinetics - k_open/k_close
                % Calculates each IP3R transition rates given its state
                q24=ip3r_par.a24+ip3r_par.v24.*(1-m24.*h24);% q42 Park --> Drive
                q42=ip3r_par.a42+ip3r_par.v42.*m42.*h42;    % q24 Drive --> Park
                
                IP3Rs_q_Old=zeros(total_IP3Rs,6);
                for TM_i=1:total_IP3Rs
                    IP3R_TMs_Old=[0 ip3r_par.q12 0 0 0 0;...
                                    ip3r_par.q21 0 ip3r_par.q23 q24(TM_i) 0 ip3r_par.q26;...
                                    0 ip3r_par.q32 0 0 0 0;...
                                    0 q42(TM_i) 0 0 ip3r_par.q45 0;...
                                    0 0 0 ip3r_par.q54 0 0;...
                                    0 ip3r_par.q62 0 0 0 0];
                    IP3Rs_q_Old(TM_i,:)=IP3R_TMs_Old(RyR_sparse(IP3RVector(TM_i),3),:);
                end
                
                % All Receptor (RyR + IP3R) gating kinetics under Yn
                tanVector = [[zeros(total_RyRs,1) ...
                    min(3.17*1e2*((Yn(RyR_sparse(RyRVector,1),1)./par.sc).^2.8),0.700) ...
                    max(0.250*((Yn(RyR_sparse(RyRVector,1),1)./par.sc).^-0.5),0.900) ...
                    zeros(total_RyRs,1)];...
                    [zeros(total_IP3Rs,1)...
                    sum(IP3Rs_q_Old,2)... 
                    sum(IP3Rs_q_Old,2)... 
                    zeros(total_IP3Rs,1)]];
                
                % All Receptor (RyR + IP3R) gating kinetics under initial Yn1
                lambda_h42=(ip3r_par.vh42-ip3r_par.ah42)*(RyR_sparse(IP3RVector,2)==2)+ip3r_par.ah42;
                [m42_New,h42_New,m24_New,h24_New]=Update_IP3R_vars(Yn1(RyR_sparse(IP3RVector,1),1),m42,h42,m24,h24,lambda_h42,ip3r_par,par);
                
                q24=ip3r_par.a24+ip3r_par.v24.*(1-m24_New.*h24_New);% q24 Drive --> Park
                q42=ip3r_par.a42+ip3r_par.v42.*m42_New.*h42_New;    % q42 Park --> Drive
                
                IP3Rs_q_New=zeros(total_IP3Rs,6);
                for TM_i=1:total_IP3Rs
                    IP3R_TMs_New=[0 ip3r_par.q12 0 0 0 0;...
                                    ip3r_par.q21 0 ip3r_par.q23 q24(TM_i) 0 ip3r_par.q26;...
                                    0 ip3r_par.q32 0 0 0 0;...
                                    0 q42(TM_i) 0 0 ip3r_par.q45 0;...
                                    0 0 0 ip3r_par.q54 0 0;...
                                    0 ip3r_par.q62 0 0 0 0];
                    IP3Rs_q_New(TM_i,:)=IP3R_TMs_New(RyR_sparse(IP3RVector(TM_i),3),:);
                end
                
                tanNewVector = [[zeros(total_RyRs,1)...
                    min(3.17*1e2*((Yn1(RyR_sparse(RyRVector,1),1)/par.sc).^2.8),0.700) ...
                    max(0.250*((Yn1(RyR_sparse(RyRVector,1),1)/par.sc).^-0.5),0.900) ...
                    zeros(total_RyRs,1)];...
                    [zeros(total_IP3Rs,1)...
                    sum(IP3Rs_q_New,2)... 
                    sum(IP3Rs_q_New,2)... 
                    zeros(total_IP3Rs,1)]];
                
                %%% Hybrid Gillespie algorithm to determine new dt
                % working to give current receptor gating kinetics
                tanIndices = [(RyR_sparse(:,2)*2)-1 RyR_sparse(:,2)*2];
                tanVectorState = tanVector(((tanIndices-1)*total_Rs)+(1:total_Rs)');
                tanNewVectorState = tanNewVector(((tanIndices-1)*total_Rs)+(1:total_Rs)');

                gOldVector = RyR_sparse(:,4);
                gNewVector = gOldVector + par.dt.*(sum(tanVectorState,2)+sum(tanNewVectorState,2))/2;
                xiVector = log(RyR_sparse(:,5).^-1);
                gNewGreaterThanXiBoolean = gNewVector >= xiVector;

                dt1Vector = par.dt*ones(total_Rs,1);   % Dummy variable to find real minimal time step between the state change of each receptor.
                dt1Vector(gNewGreaterThanXiBoolean) = (xiVector(gNewGreaterThanXiBoolean) - gOldVector(gNewGreaterThanXiBoolean))./...
                    (gNewVector(gNewGreaterThanXiBoolean) - gOldVector(gNewGreaterThanXiBoolean))*par.dt;
                
                [dt_min,indx_dt_min]=min(abs(dt1Vector));
                
                if (dt_min==par.dt) && (indx_dt_min==1)
                    indx_dt_min=0;
                end
                old_indx=indx_dt_min;
                par.dt=dt_min;  % new dt obtained here

                % Recalculate system with new dt (in "par" struct)
                [Yn1,J]=EulerMethod(Yn,par,ip3r_par,RyR_open_sums,IP3R_open_sums);          % values at the time step n+1      
                
                % Receptor gating kinetics under Yn1 with new dt (overwrites previous gating kinetics under Yn1)
                lambda_h42=(ip3r_par.vh42-ip3r_par.ah42)*(RyR_sparse(IP3RVector,2)==2)+ip3r_par.ah42;
                [m42_New,h42_New,m24_New,h24_New]=Update_IP3R_vars(Yn1(RyR_sparse(IP3RVector,1),1),m42,h42,m24,h24,lambda_h42,ip3r_par,par);
                
                q24=ip3r_par.a24+ip3r_par.v24.*(1-m24_New.*h24_New);% q24 Drive --> Park
                q42=ip3r_par.a42+ip3r_par.v42.*m42_New.*h42_New;    % q42 Park --> Drive
                
                IP3Rs_q_New=zeros(total_IP3Rs,6);
                for TM_i=1:total_IP3Rs
                    IP3R_TMs_New=[0 ip3r_par.q12 0 0 0 0;...
                                    ip3r_par.q21 0 ip3r_par.q23 q24(TM_i) 0 ip3r_par.q26;...
                                    0 ip3r_par.q32 0 0 0 0;...
                                    0 q42(TM_i) 0 0 ip3r_par.q45 0;...
                                    0 0 0 ip3r_par.q54 0 0;...
                                    0 ip3r_par.q62 0 0 0 0];
                    IP3Rs_q_New(TM_i,:)=IP3R_TMs_New(RyR_sparse(IP3RVector(TM_i),3),:);
                end
                
                tan1Vector = [[zeros(total_RyRs,1)...
                    min(3.17*1e2*((Yn1(RyR_sparse(RyRVector,1),1)/par.sc).^2.8),0.7) ...
                    max(0.250*((Yn1(RyR_sparse(RyRVector,1),1)/par.sc).^-0.5),0.9) ...
                    zeros(total_RyRs,1)];...
                    [zeros(total_IP3Rs,1)...
                    sum(IP3Rs_q_New,2)...
                    sum(IP3Rs_q_New,2)...
                    zeros(total_IP3Rs,1)]];
                
                tanIndices = [(RyR_sparse(:,2)*2)-1 RyR_sparse(:,2)*2];
                tanVectorState = tanVector(((tanIndices-1)*total_Rs)+(1:total_Rs)');
                tan1VectorState = tan1Vector(((tanIndices-1)*total_Rs)+(1:total_Rs)');

                RyR_sparse(:,4)=RyR_sparse(:,4)+...
                    par.dt.*(sum(tanVectorState,2)+sum(tan1VectorState,2))/2;

                % Y0=[Cai, Cajsr, Cansr, CaCaM, CaATP, CaF4, CaTnC, CaCSQ];
                %      1     2      3      4      5     6      7      8
                Cai=Yn1(:,1);
                CaJSR=Yn1(:,2);
                CaNSR=Yn1(:,3);
                CaCaM=Yn1(:,4);
                CaATP=Yn1(:,5);
                CaF4=Yn1(:,6);
                CaTnC=Yn1(:,7);
                CaCSQ=Yn1(:,8);
                h24=h24_New;
                m24=m24_New;
                h42=h42_New;
                m42=m42_New;

                %%% Determine new Receptor (RyR/IP3R) State
                % 3 cases: no R, RyR, IP3R change state                
                if indx_dt_min==0   % no Receptors change state
                    par.dt=par.dt0;
                    
                elseif indx_dt_min<=total_RyRs  % RyR change state
                    R_position=RyR_sparse(indx_dt_min,1);   % current R position
                    R_state=RyR_sparse(indx_dt_min,2);      % current R state
                    
                    Tan1=[tan1Vector(indx_dt_min,1:2); ...  % [0 k_open;
                            tan1Vector(indx_dt_min,3:4)];   %  k_close 0]
                    Prev_rate=Tan1(R_state,:);
                    New_state=find((cumsum(Prev_rate)/sum(Prev_rate))>=rand,1);
                    
                    RyR_sparse(indx_dt_min,2)=New_state;
                    RyR_sparse(indx_dt_min,3)=New_state;
                    RyR_sparse(indx_dt_min,4)=0;
                    RyR_sparse(indx_dt_min,5)=rand;
                    RyR_open_sums(R_position)=RyR_open_sums(R_position)+Ch_mx(New_state);
                    
                else    % IP3R change state
                    R_position=RyR_sparse(indx_dt_min,1);
                    R_state=RyR_sparse(indx_dt_min,3);  % note: use col 3 instead of 2
                    
                    Prev_rate=IP3Rs_q_New(indx_dt_min-total_RyRs,:);
                    New_state=find((cumsum(Prev_rate,2)/sum(Prev_rate,2))>=rand,1);
                    
                    RyR_sparse(indx_dt_min,2)=1+1.*(heaviside(New_state-4.5));
                    RyR_sparse(indx_dt_min,3)=New_state;
                    RyR_sparse(indx_dt_min,4)=0;
                    RyR_sparse(indx_dt_min,5)=rand;
                    IP3R_open_sums(D2_IP3R_Cluster1)=sum(RyR_sparse(total_RyRs+1:total_RyRs+nIP3R,2)==2);
                    IP3R_open_sums(D2_IP3R_Cluster2)=sum(RyR_sparse(total_RyRs+1+nIP3R:total_Rs,2)==2);
                end
                dumt=dumt+par.dt;
            end

            % save results
            if iter > rec_startTime/save_t_int
                elem_num_save=iter-rec_startTime/save_t_int;
                RyR_open_save(:,elem_num_save)=RyR_open_sums;
                IP3R_open_save(:,elem_num_save)=IP3R_open_sums;
                Cai_save(:,elem_num_save)=Cai;
                CaJSR_save(:,elem_num_save)=CaJSR;
                CaNSR_save(:,elem_num_save)=CaNSR;
                CaF4_save(:,elem_num_save)=CaF4;
                CaATP_save(:,elem_num_save)=CaATP;
                CaCaM_save(:,elem_num_save)=CaCaM;
                CaTnC_save(:,elem_num_save)=CaTnC;
                time(elem_num_save)=simt+dumt;
            end
            simt=simt+dumt;
        end
        toc
        
        parsave(['coupl_rep_uncoupl_mult_clusters_Cannell_Cao_ryrs_',num2str(total_RyRs),'_ip3rs_',num2str(total_IP3Rs),'_Tsim_',num2str(simTime),'smaller_set_test_v',num2str(fop),'_tr13_no_trigger_at_all.mat'],time,Cai_save,CaJSR_save,CaNSR_save,CaF4_save,RyR_open_save,IP3R_open_save,par,ip3r_par)

    end
    disp('finally!')
end