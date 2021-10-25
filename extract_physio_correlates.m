%%Calculation the correlation of physiological parameters, head motion and DVARS.

%add tapas toolbox to path
addpath(genpath('/data/u_fbeyer_software/spm-fbeyer'))


%Specify variables
subjects_file='/data/gh_gr_agingandobesity_share/life_shared/Results/External_coops/PV240_Dagher_Scholz_headmotion_gwas/Amendment2021/Analysis/all_data_final_30.06.21.txt';

subjID = fopen(subjects_file);
subjects=textscan(subjID,'%s');

all_res=zeros(size(subjects{1},1),13);

size(subjects{1})

for i=1:size(subjects{1},1)

    %if subjects{1}{i}=="LI01273319"%for looking at plots of individual subjects

        %RESTING STATE scan needs to be processed
        if isfile(sprintf('/data/pt_life_restingstate_followup/Data/resting/moco/%s/rest_realigned.nii.gz.par', subjects{1}{i}))
            fileID = fopen(sprintf('/data/pt_life_restingstate_followup/Data/resting/moco/%s/rest_realigned.nii.gz.par', subjects{1}{i}),'r');
            formatSpec = '%f %f %f %f %f %f';
            motion_data = fscanf(fileID,formatSpec);
	    fclose(fileID);
            %reshape and remove first 5 volums
            motion_data=transpose(reshape(motion_data,[6,296]));

            %calculate mean FD
            %parameter_source == 'FSL':
            translations = abs(diff(motion_data(:,4:6)));
            rotations = abs(diff(motion_data(:,1:3)));

            FD_power = sum(translations,2) + sum((50*rotations),2);
            all_res(i,1)=mean([0;sum(translations,2) + sum((50*rotations),2)]);
            all_res(i,2)=max([0;sum(translations,2) + sum((50*rotations),2)]);
            
        
            %DVARS data (standardized BOLD signal in entire brain mask)
            if isfile(sprintf(['/data/pt_life/LIFE_fu/wd_preprocessing/hcp_prep_workflow/resting/'...
                                  'transform_timeseries/_subject_%s/dvars/rest2anat_dvars.tsv'], subjects{1}{i}))

                dvars=tdfread(sprintf(['/data/pt_life/LIFE_fu/wd_preprocessing/hcp_prep_workflow/resting/'...
                                  'transform_timeseries/_subject_%s/dvars/rest2anat_dvars.tsv'], subjects{1}{i}));
                dvars_std=dvars.std_DVARS;   

                %results head motion & DVARS

                all_res(i,3)=mean(dvars_std);
                all_res(i,4)=max(dvars_std);
                all_res(i,5)=corr(dvars_std,FD_power);
            end
        end
        
        
        %PHYSIOLOGICAL data
        if isfile(sprintf("/data/pt_life_restingstate_followup/Data/physio/%s.mat", subjects{1}{i}))
            
            if isfile(sprintf('/data/pt_life_restingstate_followup/Data/resting/moco/%s/rest_realigned.nii.gz.par', subjects{1}{i}))
                phys_data=load(sprintf("/data/pt_life_restingstate_followup/Data/physio/%s.mat", subjects{1}{i}));

                %use respiratory and PPU traces resampled to volume acquisitions.
                resp=phys_data.physio.trace.resp;
                oxy=phys_data.physio.trace.oxy;

                %remove first 5 volumes.
                resp=resp(6:300);
                oxy=oxy(6:300);

                all_res(i,6)=corr(FD_power,resp);
                all_res(i,7)=corr(FD_power,oxy);
                all_res(i,8)=corr(FD_power,phys_data.physio.ons_secs.rvt(6:end));
                all_res(i,9)=corr(FD_power,phys_data.physio.ons_secs.hr(6:end)); 
            end
            if isfile(sprintf(['/data/pt_life/LIFE_fu/wd_preprocessing/hcp_prep_workflow/resting/'...
                                  'transform_timeseries/_subject_%s/dvars/rest2anat_dvars.tsv'], subjects{1}{i}))
                	all_res(i,10)=corr(dvars_std,resp);
                	all_res(i,11)=corr(dvars_std,oxy);
                	all_res(i,12)=corr(dvars_std,phys_data.physio.ons_secs.rvt(6:end));
                	all_res(i,13)=corr(dvars_std,phys_data.physio.ons_secs.hr(6:end)); 
            end

        else 
            continue
        end
    %end
      
end
% 

fclose('all')
T = array2table(all_res, 'VariableNames',{'meanFD','maxFD','meanstdDVARS', 'maxstdDVARS', 'corr_FD_stdDVARS'...
                                            'corr_FD_resp','corr_FD_oxy','corr_FD_RVT','corr_FD_HR',...
                                            'corr_dvars_resp','corr_dvars_oxy','corr_dvars_RVT','corr_dvars_HR'})
  
T_final=[cell2table(subjects{1}),T];
writetable(T,'/data/pt_life_restingstate_followup/Results/Physio/results_rs_motion_physio.csv');
writetable(cell2table(subjects{1}),'/data/pt_life_restingstate_followup/Results/Physio/SIC.csv');




