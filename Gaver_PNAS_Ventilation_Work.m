% Script Name: Gaver_PNAS_Ventilation.m
% Description: This script calculates PV work for the recruitment work from the Upper Inflection Point.
%To do so, it determines the inflection point of the insipiration limb of
%the PV curve and then calculates the area difference between the Ptp vs V loop
%and the secant from the beginning of the inspiration component% 

% Inputs:
% - PNAS_DataFiles.txt: A text file listing trial directories.
% - header.txt: A metadata file in each trial directory containing
% experimental conditions Tlow and Phigh

% Calculations:
% Input work, total work, tissue work, and  recruitment work using the pressure-volume loop.

% Outputs:
% Visualization:
% Creates figures visualizing the PV loop and recruitment region
% Saves figures to designated directories.

% Saves calculated work values to Excel files 
% PNAS_Ventilation_Work.xlsx in trial directories
% PNAS_Ventilation_Work_Combined.xlsx in parent directory

% - Processed data files saved in their respective directories.

% Assumptions:
% - Each directory contains a correctly formatted 'header.txt'.
% - All paths listed in 'PNAS_DataFiles.txt' are valid.
% Author: [Your Name]
% Date: [Insert Date]

%Ventilation_Work_Inflation_Limb_Ptp_secant.m
%This code calculates the recruitment work from the Upper Inflection Point.
%To do so, it determines the inflection point of the insipiration limb of
%the PV curve and then calculates the area difference between the Ptp vs V loop
%and the secant from the beginning of the inspiration component
clear all
%% Change Directories
Input_files = 'PNAS_DataFiles.txt'; %This provides the directories for the data to be analyzed
InputfileID = fopen(Input_files,'rt');
NumTrials= fscanf(InputfileID,'%u');        %Number of different studies (Pigs)
Directory={};

for Trial = 1:NumTrials                     % Read Directory for each pig
    Directory{Trial} = fgetl(InputfileID);
end

for Trial = 1:NumTrials              
cd (Directory{Trial})
mkdir PNAS_Figures_PV_Recruit   %Create figures directory showing recruitment region
mkdir PNAS_Figures_PV                 %Create figures directory showing PV loops

%Read file data
filename = 'header.txt';            %Experimental information Descrip, Tlow, Phigh
fileID_In = fopen(filename,'r');

Descrip = fgetl(fileID_In)
Tlow = fscanf(fileID_In, '%u', 1);
Phigh = fscanf(fileID_In, '%u', 1);

filename = 'files.xls';         %Provides file name and airway resistance
[num,txt] = xlsread(filename);
Res = num';                     %Airway Resistance
datafilenames=txt';

mydata={};
numfiles=length(datafilenames);
for k = 1:numfiles
    myfilename = datafilenames{k}; % name of k th file in files.txt
    mydata{k} = importdata(myfilename); %imports data into array mydata
end

%%
%Loop through the files
%Reset arrays
Descriptor={};
Pig_Descriptor=[];
Directory_Descriptor=[];
Tlow_Descriptor=[];
Phigh_Descriptor=[];
N_filter=[];

Work_Insp_Total=zeros(length(mydata),1);
Work_Total=zeros(length(mydata),1);
Work_tissue=zeros(length(mydata),1);
Work_Recruit=zeros(length(mydata),1);
t_low=zeros(length(mydata),1);
nfilter = 2;                        %low pass filter cutoff

for iloop=1:length(mydata)          %number of timepoints (BL, T0, ...)
    loopdata =mydata{1,iloop};

    Pig_Descriptor{iloop}=Descrip;
    Directory_Descriptor{iloop}=Directory{Trial};
    Tlow_Descriptor{iloop}=Tlow;
    Phigh_Descriptor(iloop)=Phigh;
    N_filter(iloop) = nfilter;
    
    %% Convert to output type
    t = loopdata(:,1);
    P_trach = loopdata(:,2);                    %Trachael Pressure
    Q = loopdata(:,3);                          %Flowrate
    P_esoph = loopdata(:,4);                    %Esophogeal pressure
    delt=t(2)-t(1);  
    P_tp=P_trach - P_esoph;                     %Transpulmonary Pressure   
    P_alv = P_trach - Res(iloop)*Q;             %Interpretted Alveolar Pressure             
    P_tissue=P_alv - P_esoph;                    %Tissue Pressure   
    V=zeros(length(Q),1);                        %Empty array for V
    
%% UNFILTERED ANAlYSIS
    % Estimate Volume and Remove Bias
        V(1)=0; % for initial calculation
    for i= 1:length(t)-1
        V(i+1)=V(i)+((Q(i)+Q(i+1))/2)*(t(i+1)-t(i));
    end
%%############# BIAS DRIFT REMOVAL  #######################################    
for j=1:10                                      %Cycle to remove bias trend in data
    %Find Min and Max values
        %calculate flow-rate bias
        %Find First min/max
        rangeV1 =3:length(t)/2;
        Vsubset1 = V(rangeV1);     
        [~, index1]=max(Vsubset1);
            max_index1= rangeV1(index1);  
        [~, index1]=min(Vsubset1);
            min_index1= rangeV1(index1);

        %Find Second min/max
        start=max(min_index1,max_index1)+10;
        rangeV2 =start:length(t);
        Vsubset2 = V(rangeV2);     
        [~, index2]=max(Vsubset2);
            max_index2= rangeV2(index2);  
        [~, index2]=min(Vsubset2);
            min_index2= rangeV2(index2);

        %determine flow rate bias
        Qbias_min=(V(min_index2)-V(min_index1))/...
            (t(min_index2)-t(min_index1));

        Qbias_max=(V(max_index2)-V(max_index1))/...
            (t(max_index2)-t(max_index1));

        Qbias=(Qbias_min+Qbias_max)/2;

        %Remove Qbias
        for i= 1:length(t)
            Q(i)=Q(i)-Qbias_min;
        end
        
        loop_indexV=min_index1:min_index2;
        
    %% Calculate V(t)
        V(1)=0; % for initial calculation
    for i= 1:length(t)-1
        V(i+1)=V(i)+((Q(i)+Q(i+1))/2)*(t(i+1)-t(i));
    end
end
% Reset Minimum to baseline lung volume
     Vmin = 0; %liters                          
     V=V+(Vmin-min(V));
     Vmax(iloop)= max(V);
%%########################################################################  

%% Calibrated final loop
 %Find Min and Max values
        %Find First min/max
        range1 =1:length(t)/2;   %starts after filtered length
        Vsubset1 = V(range1);     
        [~, index1]=max(Vsubset1);
            Vmax_index1_unfilter= range1(index1);  
        [~, index1]=min(Vsubset1);
            Vmin_index1_unfilter= range1(index1);            
            
        %Find Second min/max
        start=max(min_index1,max_index1)+10;
        range2 =start:length(t);
        Vsubset2 = V(range2);     
        [~, index2]=max(Vsubset2);
            Vmax_index2_unfilter= range2(index2);  
        [~, index2]=min(Vsubset2);
            Vmin_index2_unfilter= range2(index2);
            
 % Determine Loop Min and Max
 % Find Pressure Min value
 % Find First min/max
        range1 = 1:round(length(t)/2);   %first half cycle
        Psubset1 = P_trach(range1);  
        [~, index1]=min(Psubset1);                                    %Pressure minumum
        Pmin_index1_unfilter= range1(index1);                              %Pressure minumum index

        range2 =Pmin_index1_unfilter:Pmin_index1_unfilter+round(length(t)/4);   %starts after min           
        Psubset2 = P_trach(range2);    
        [~, index2]=max(Psubset2);                                    %Pressure maximum
        Pmax_index2_unfilter= range2(index2);                              %Pressure maximum index

    %% ################## CALCULATE WORK ################################################### 
    Work_Total(iloop)=0;                                                    %Total Work
    Work_tissue(iloop)=0;                                                   %Tissue Work
    loop_index=Vmin_index1_unfilter:Vmin_index2_unfilter-1;                   % Volume based loop index
    for i=loop_indexV-1                                                     %Shift by one time point to eliminate index over-run
        P_tp_avg=(P_tp(i+1)+P_tp(i-1))/2;
        P_tissue_avg=(P_tissue(i+1)+P_tissue(i-1))/2;                           %Tissue pressure
        Q_avg=(Q(i+1)+Q(i-1))/2;
        Work_Total(iloop)= Work_Total(iloop)+P_tp_avg*Q_avg*delt;       %Total work dissipation provided by ventilator  
        Work_tissue(iloop)= Work_tissue(iloop)+P_tissue_avg*Q_avg*delt;               %Total work provided due to chest wall and ventiltor
    end
    
%% Inhalation Information
    istart=Pmin_index1_unfilter;                                            %Inflation portion
    iend=Pmax_index2_unfilter;
    
     P_trach_inflate=zeros(iend-istart,1);
     P_tp_inflate=zeros(iend-istart,1);
     Q_inflate=zeros(iend-istart,1);  
     P_tissue_inflate=zeros(iend-istart,1);  
     V_inflate=zeros(iend-istart,1); 
     V_inflate_filter=zeros(iend-istart,1);  
     t_inflate=zeros(iend-istart,1);  

    for i=istart:iend
        index=i-istart+1;
        P_trach_inflate(index)=P_trach(i);  
        P_tp_inflate(index)=P_tp(i);  
        P_tissue_inflate(index)=P_tissue(i);
        t_inflate(index)=delt*(index-1);
        Q_inflate(index)=Q(i);
    end
    
%% Inspiratory Work
    Work_Insp_Total(iloop)=0;
    V_inflate(1)=V(istart); % for initial calculation  
    Insp_loop_index=1:length(Q_inflate)-1;
    for i=Insp_loop_index
        P_trach_avg=(P_trach_inflate(i+1)+P_trach_inflate(i))/2;
        P_tp_avg=(P_tp_inflate(i+1)+P_tp_inflate(i))/2;
        Q_avg=(Q_inflate(i+1)+Q_inflate(i))/2;
        V_inflate(i+1)=V_inflate(i)+Q_avg*delt;
%        Work_Insp_Total(iloop)= Work_Insp_Total(iloop)+P_tp_avg*Q_avg*delt;     %Inspiratory work provided by ventilator    
        Work_Insp_Total(iloop)= Work_Insp_Total(iloop)+P_trach_avg*Q_avg*delt;     %Inspiratory work provided by ventilator    
    end
%% ##############  INFLECTION POINT ANALYSIS #################################
        P_tissue_min = min(P_tissue_inflate);
        P_tissue_max = max(P_tissue_inflate);
        
        low_min_init=1;
            low_min=low_min_init;
            del=1;
            for i=1+del:min(75,round(length(Q_inflate)/2))
                i;
                delP=P_tissue_inflate(i+del)-P_tissue_inflate(i-del);
                if(delP<0) && (V_inflate(i+del)< 0.1*max(V_inflate))
                    low_min = i;
                end
            end

% Truncated inflation leg 
    istart=istart+(low_min-1);                                            %Inflation portion
    iend=Pmax_index2_unfilter;
    
     P_tissue_inflate_2=zeros(iend-istart,1); 
     Q_inflate_2=zeros(iend-istart,1); 
     V_inflate_2=zeros(iend-istart,1); 
     V_inflate_filter=zeros(iend-istart,1);  

    for i=istart:iend
        index=i-istart+1;
        P_tissue_inflate_2(index)=P_tissue(i);
        Q_inflate_2(index)=Q(i);
        V_inflate_2(index)=V(i);
    end
    
%  Filtering of truncated inflation leg
    del_t=t(2)-t(1);
    fs=1/del_t;             %Sampling Frequency
    f_Nyquest = fs/2;       %Nyquest Frequency
    Steep=0.5;

% Use filtered data for determination of inflection point
    P_tissue_inflate_filter=lowpass(P_tissue_inflate_2,nfilter,fs,'ImpulseResponse','iir','Steepness',Steep);
    Q_inflate_filter=lowpass(Q_inflate_2,nfilter,fs,'ImpulseResponse','iir','Steepness',Steep);

    V_inflate_filter(1)=V_inflate_2(1); % for initial calculation    
    for i= 1:length(Q_inflate_filter)-1
        V_inflate_filter(i+1)=V_inflate_filter(i)+((Q_inflate_filter(i)+Q_inflate_filter(i+1))/2)*del_t;
    end

%% Derivatives - for determining curvature
        dP_noflow_filter = gradient(P_tissue_inflate_filter);
        dV_inflow_filter = gradient(V_inflate_filter);

        dVdP_noflow_filter=dV_inflow_filter./dP_noflow_filter;
        d2VdP2_noflow_filter=gradient(dVdP_noflow_filter);

%% Find Inflection Point and Secant line
% Find pressure at inflection point locs(1) 
        range=1:round(length(V_inflate_filter));
        [~,locs]=findpeaks(dVdP_noflow_filter(range));                    %Inflection points in **_inflate_2 framework
        locsort=sort(locs);                                                 %order
        locs=locsort;
        locs=locs+(low_min-1);                                              %locs in framework of **_inflate            
        
        if(length(locs)>=1)                                                 %If inflection point exists continue analysis

%Pick appropriate inflection point
        VT= max(V_inflate);                                                  %Tidal Volume
        
        if(length(locs)>=2) 
            inflect_store=locs(1);
            for i_loc = 2:min(10,length(locs))
                inflect_last=inflect_store; 
                inflect_current=locs(i_loc);             
                if(V_inflate(inflect_last) < VT/5)                          % Discard low inflection point                   
                    inflect_store=locs(i_loc);
                end                
                if(abs(V_inflate(inflect_current)-V_inflate(inflect_last))...
                        <VT/10 ...
                        && V_inflate(inflect_current) > VT/5 ...
                        && V_inflate(inflect_current) <= VT/2)                          %differential                     
                    inflect_store=locs(i_loc);
                end               
                if(P_tissue_inflate(inflect_current)<P_tissue_inflate(inflect_last) ...
                        && V_inflate(inflect_current)<=VT/2)                          %differential                     
                    inflect_store=locs(i_loc);
                end
            end
                locs(1)=inflect_store;                                      %store final inflection point at locs(1)
        else                
                locs(1)=min(locs);                                          %translate location of inflection point
        end                                                                      % to untruncated inflation file

%% Draw Secants
            P_init = P_tissue_inflate(low_min);                           %Initial points 
            V_init = V_inflate(low_min);            
            P_inf = P_tissue_inflate(locs(1));                           %Inflection points 
            V_inf = V_inflate(locs(1));
            t_recruit(iloop)=delt*(locs(1)-low_min);

            %Array for integration            
            range_recruit=low_min:locs(1);                 % Recruitment Range reduce range to prevent overshoot
            V_recruit=zeros(length(range_recruit),1);      % Define arrays adding two for endpoints
            P_recruit=zeros(length(range_recruit),1);
            
            V_recruit(1) = V_init;             % 1st points
            P_recruit(1) = P_init;
            
            for i=range_recruit
                V_recruit(i-low_min+1)=V_inflate(i);
                P_recruit(i-low_min+1)=P_tissue_inflate(i);
            end
            
% Calculate Recruit Work
            Int_Pdv_loop=trapz(V_recruit,P_recruit);                %Loop area in recruit region
            DelV_trunc = (V_inf-V_init);                            %Volume span in truncated region
            P_trunc = (P_inf + P_init)/2;                           %Average pressure
            Int_Pdv_trunc=P_trunc*DelV_trunc;                       %Truncated area
            Work_Recruit(iloop) = Int_Pdv_loop-Int_Pdv_trunc;
            
                        
%% RECRUITMENT FIGURES
           range_inflate=low_min:round(length(V_inflate_filter));           % Range for plot
    cd PNAS_Figures_PV_Recruit %for depositing figures (directory created above)
            figure('visible','off')
            hold on
            plot(P_tissue_inflate(range_inflate),V_inflate(range_inflate));                  % Raw inflation limb
            plot(P_recruit,V_recruit,'.');                                              % Raw data in recruitment range
            plot(P_tissue_inflate(locs),V_inflate(locs),'o');                              % Inflection points on raw data      
            x1 = [P_init, P_inf];
            y1 = [V_init, V_inf];
            plot(x1,y1);                                                                % Tangent line
            hold off
            xlabel({'Pressure_{tp} (cmH2O)'});
            ylabel({'Relative Volume (l)'})%P_trach vs time
            title(Descrip, ...
                ['Tlow = ',num2str(Tlow), ', Phigh=', num2str(Phigh),...
                ',   ', datafilenames{iloop}]);
    FileName=append(Descrip,'_Recruitment_',num2str(iloop),'.png');
    saveas(gcf,FileName,'png')
        
    cd ..  %Moves back to parent directory   
        
%% PV FIGURES
    cd PNAS_Figures_PV                                                                   % for depositing figures (directory created above)
            figure('visible','off')
            hold on

             plot(P_tp(loop_indexV),V(loop_indexV));                                                           % airway pressure
             plot(P_tissue(loop_indexV),V(loop_indexV));                                        % Transpulmonary pressure            
            hold off
            xlabel({'Pressure (cmH2O)'});
            ylabel({'Relative Volume (l)'})%P_trach vs time
            title(Descrip, ...
                ['Tlow = ',num2str(Tlow), ', Phigh=', num2str(Phigh),...
                ',   ', datafilenames{iloop}]);
    FileName=append(Descrip,'_PV_',num2str(iloop),'.png');
    saveas(gcf,FileName,'png')
    
     %Write datafile
        PVTable_local = table(P_trach(loop_indexV), P_tp(loop_indexV), ...
                      P_tissue(loop_indexV), V(loop_indexV), ...
                      'VariableNames', {'P_trach', 'P_tp', 'P_tissue', 'Volume'});
        FileName = append('PVloop_',num2str(iloop),'.xlsx');
        writetable(PVTable_local,FileName,"WriteMode","append");

    cd ..  %Moves back to parent directory  
        end      
        formatSpec = 'File evaluated %s \n';
        fprintf(formatSpec,char(datafilenames(iloop)))       % Display file evaluated        
        Store_Work_Recruit(iloop)=Work_Recruit(iloop);

end                                         %iloop of number of timepoints for each pig

    PVTable = table(datafilenames', Tlow_Descriptor', Phigh_Descriptor', ...
                Work_Insp_Total, Work_Total, Work_tissue, Work_Recruit, ...
                'VariableNames', {'FileName', 'Tlow', 'Phigh', ...
                                  'Input Work', 'Total Work', ...
                                  'Tissue Work', 'Recruit Work'});
    writetable(PVTable,'PNAS_Ventilation_Work.xlsx',"WriteMode","append");

cd ../../;
 T={};
  T = table(Pig_Descriptor', Directory_Descriptor', datafilenames', ...
          Tlow_Descriptor', Phigh_Descriptor', Work_Insp_Total, ...
          Work_Total, Work_tissue, Work_Recruit, ...
          'VariableNames', {'Pig', 'Directory', 'FileNames', ...
                            'Tlow', 'Phigh', 'Insp Work', ...
                            'Total Work', 'WorkPalv', 'WorkRecruit'});

      writetable(T,'PNAS_Ventilation_Work_Combined.xlsx',"WriteMode","append");
end                                         %Loop through pigs (Trial = 1:NumTrials) 
