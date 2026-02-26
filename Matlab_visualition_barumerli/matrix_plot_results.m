%{
This standalone code lauch the barumerli's model evaluations, plots the results as shows on the png joined with this code.
To work, you have to install and launch the amtoolbox and to put the uncompressed sofa files .(Le-roux-1D-subject0.sofa , GT etc.) into a folder designated
with path-norm ( you will have to modify it for all instances of path-norm).
No
%}

%% --- Computation of intersubject errors ---
subject_list=["P0091" "P0092" "P0093" "P0095" "P0096" "P0097" "P0098" "P0099" "P0100"];
subject_eval=["subject_0" "subject_1" "subject_2" "subject_3" "subject_4" "subject_5" "subject_6" "subject_7" "subject_8" ];

Matrix_inter_subject=zeros(9,9,3);
Matrix_mean_model_perf=zeros(9,10,3);
Matrix_mean_model_perf_inv=zeros(9,10,3);
progBar=ProgressBar(81,'Title','Intersub errors');
path_norm="C:\Users\alexa\Documents\MATLAB\sofa_hrtf_estimated\";
for j=1:9
    for i=1:9
        sofa_template=strcat(path_norm,subject_list(j),"_Windowed_44kHz.sofa");
        sofa_target=strcat(path_norm,subject_list(i),"_Windowed_44kHz.sofa");
        evalc('Sim_sofa=SOFAload(sofa_target)');
        evalc('GT_sofa=SOFAload(sofa_template)');
        template_dtf=SOFAhrtf2dtf(GT_sofa);
        target_dtf=SOFAhrtf2dtf(Sim_sofa);
        [template_extract, uwu]=barumerli2023_featureextraction(template_dtf, 'dtf','fs',44100);
        [wuw,target_extract]=barumerli2023_featureextraction(target_dtf,'dtf','fs',44100);
        m = barumerli2023(...
                'template', template_extract, ...
                'target', target_extract, ...
                'num_exp', 50, ...
                'sigma_ild', 0.75, ...
                'sigma_spectral', 4.3, ...
                'sigma_prior', 11.5, ...
                'sigma_motor', 13.45);
        sim = barumerli2023_metrics(m, 'middle_metrics');
        Matrix_inter_subject(j,i,1)=sim.rmsP;   % 1, polar rms 2, lateral rms 3, erreur de quadrant
        Matrix_inter_subject(j,i,2)=sim.rmsL;
        Matrix_inter_subject(j,i,3)=sim.querr;
        progBar([], [], []);

    end
end
progBar.release();


%% --- Plot of intersubject errors ---
subject_list=["P0091" "P0092" "P0093" "P0095" "P0096" "P0097" "P0098" "P0099" "P0100"];
subject_eval=["subject_0" "subject_1" "subject_2" "subject_3" "subject_4" "subject_5" "subject_6" "subject_7" "subject_8" ];
figure('Color','w','Position',[100 100 1100 350])

tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
imagesc(Matrix_inter_subject(:,:,1).')
set(gca,'YDir','normal')
axis image
colormap(gca,"copper")
colorbar
title('Rms local polar error (deg)','FontWeight','bold')
set(gca,'XTick',1:9,'XTickLabel',subject_list)
set(gca,'YTick',1:9,'YTickLabel',subject_list)
xlabel('Reference HRTF')
ylabel('Target HRTF ')


nexttile
imagesc(Matrix_inter_subject(:,:,2).')
set(gca,'YDir','normal')
axis image
colormap(gca,"copper")
colorbar
title('Rms lateral error (deg)','FontWeight','bold')
set(gca,'XTick',1:9,'XTickLabel',subject_list)
set(gca,'YTick',1:9,'YTickLabel',subject_list)
xlabel('Reference HRTF')

nexttile
imagesc(Matrix_inter_subject(:,:,3).')
set(gca,'YDir','normal')
axis image
colormap(gca,"copper")
colorbar
title('Quadrant error (%)','FontWeight','bold')
set(gca,'XTick',1:9,'XTickLabel',subject_list)
set(gca,'YTick',1:9,'YTickLabel',subject_list)
xlabel('Reference HRTF')
save("matrice_erreur_intersujet.mat","Matrix_inter_subject");

%% --- Computation of intermodels errors ---
progBar=ProgressBar(90,'Title','Intermodels errors');
model_type=["GT" "Woo-lee_1d" "Woo-lee_2d" "Woo-lee_3d" "Le-roux_1d" "Le-roux_2d" "Le-roux_3d" "Manlin-Zhao_1d" "Manlin-Zhao_2d" "Manlin-Zhao_3d"];
for j=1:9
    for i=1:10
        sofa_template=strcat(path_norm,subject_list(j),"_Windowed_44kHz.sofa");
        if i==1
            sofa_target=strcat(path_norm,subject_list(j),"_Windowed_44kHz.sofa");
        else
            sofa_target=strcat(path_norm,model_type(i),"_",subject_eval(j),".sofa");
        end
        evalc('Sim_sofa=SOFAload(sofa_target)');
        evalc('GT_sofa=SOFAload(sofa_template)');
        template_dtf=SOFAhrtf2dtf(GT_sofa);
        target_dtf=SOFAhrtf2dtf(Sim_sofa);
        [template_extract, uwu]=barumerli2023_featureextraction(template_dtf, 'dtf','fs',44100);
        [wuw,target_extract]=barumerli2023_featureextraction(target_dtf,'dtf','fs',44100);
        m = barumerli2023(...
                'template', template_extract, ...
                'target', target_extract, ...
                'num_exp', 50, ...
                'sigma_ild', 0.75, ...
                'sigma_spectral', 4.3, ...
                'sigma_prior', 11.5, ...
                'sigma_motor', 13.45);
        sim = barumerli2023_metrics(m, 'middle_metrics');
        Matrix_mean_model_perf(j,i,1)=sim.rmsP;   % 1, polar rms 2, lateral rms 3, erreur de quadrant
        Matrix_mean_model_perf(j,i,2)=sim.rmsL;
        Matrix_mean_model_perf(j,i,3)=sim.querr;
        progBar([], [], []);

    end
end
progBar.release();
%%  --- Computation of Kemar/GT errors ---
subject_list=["P0091" "P0092" "P0093" "P0095" "P0096" "P0097" "P0098" "P0099" "P0100"];
path_norm="C:\Users\alexa\Documents\MATLAB\sofa_hrtf_estimated\";
sofa_target="C:\Users\alexa\Documents\MATLAB\KEMAR_GRAS_EarSim_LargeEars_Windowed_44kHz.sofa";
kemar_errors=zeros(9,3);
progBar=ProgressBar(9,'Title','Kemar computation');
for j=1:9
    sofa_template=strcat(path_norm,subject_list(j),"_Windowed_44kHz.sofa");
    evalc('Sim_sofa=SOFAload(sofa_target)');
    evalc('GT_sofa=SOFAload(sofa_template)');
    template_dtf=SOFAhrtf2dtf(GT_sofa);
    target_dtf=SOFAhrtf2dtf(Sim_sofa);
    [template_extract, uwu]=barumerli2023_featureextraction(template_dtf, 'dtf','fs',44100);
    [wuw,target_extract]=barumerli2023_featureextraction(target_dtf,'dtf','fs',44100);
    m = barumerli2023(...
        'template', template_extract, ...
        'target', target_extract, ...
        'num_exp', 50, ...
        'sigma_ild', 0.75, ...
        'sigma_spectral', 4.3, ...
        'sigma_prior', 11.5, ...
        'sigma_motor', 13.45);
     sim = barumerli2023_metrics(m, 'middle_metrics');
     kemar_errors(j,1)=sim.rmsP;
     kemar_errors(j,2)=sim.rmsL;
     kemar_errors(j,3)=sim.querr;
     progBar([], [], []);
end
progBar.release();
save("kemar_mean.mat","kemar_errors")
%% --- Plot of intermodels errors norm ---


kemar_baseline=mean(kemar_errors,1);


subject_list=["P0091" "P0092" "P0093" "P0095" "P0096" "P0097" "P0098" "P0099" "P0100"];
subject_eval=["subject_0" "subject_1" "subject_2" "subject_3" "subject_4" "subject_5" "subject_6" "subject_7" "subject_8" ];
mat_calc=Matrix_inter_subject;
save("matrice_erreur_intermodel.mat","Matrix_mean_model_perf");
sum_best_P=0;
sum_best_L=0;
sum_best_Qu=0;
for k=1:9
    mat_calc(k,k,:)=0;
    sum_best_P=sum_best_P+Matrix_inter_subject(k,k,1);
    sum_best_L=sum_best_L+Matrix_inter_subject(k,k,2);
    sum_best_Qu=sum_best_Qu+ Matrix_inter_subject(k,k,3);
end

mean_rmsP_top=sum_best_P/9;
mean_rmsL_top=sum_best_L/9;
mean_querr_top=sum_best_Qu/9;

mean_rmsP_baseline=sum(mat_calc(:,:,1),"all")/72;
mean_rmsL_baseline=sum(mat_calc(:,:,2),"all")/72;
mean_querr_baseline=sum(mat_calc(:,:,3),"all")/72;

model_list=["GT" "Woo-Lee-1d" "Woo-Lee-2d" "Woo-Lee-3d" "Le-Roux-1d" "Le-Roux-2d" "Le-Roux-3d" "Manlin-Zhao-1d" "Manlin-Zhao-2d" "Manlin-Zhao-3d"];

figure('Color','w','Position',[100 100 1100 700])

tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile
rmsP_sub_mean=mean(Matrix_mean_model_perf,1);
rmsP_sub_std=std(Matrix_mean_model_perf,0,1);
hold on
plot(0:10,kemar_baseline(1)*ones(1,11),"Color","black",'LineWidth',2)
plot(0:10,mean_rmsP_top*ones(1,11),"Color","blue",'LineWidth', 2)
plot(0:10,mean_rmsP_baseline*ones(1,11),"Color","red",'LineWidth', 2)
plot(1:9,rmsP_sub_mean(1,2:10,1),'+','MarkerSize', 10,'Color',"black",'LineWidth', 2)
errorbar(rmsP_sub_mean(1,2:10,1),rmsP_sub_std(1,2:10,1) ,'.','Color',"black")

y1=kemar_baseline(1);
y2=mean_rmsP_baseline;

fill([0 10 10 0], [y1 y1 y2 y2], ...
     [0.01 0.39 0.25], 'FaceAlpha', 0.3, 'EdgeColor', 'none')



y1 = mean_rmsP_baseline;
y2 = mean_rmsP_top;


fill([0 10 10 0], [y1 y1 y2 y2], ...
     'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

y11 = mean_rmsP_top;
y22 = 30;


fill([0 10 10 0], [y11 y11 y22 y22], ...
     [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

y11 = 50;
y22 = kemar_baseline(1);


fill([0 10 10 0], [y11 y11 y22 y22], ...
     'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off

title('Polar error (°) of models')
set(gca,'XTick',1:9,'XTickLabel',model_list(1,2:10))
xlabel('Models')
legend('kemar baseline','self HRTFs mean performance',"intersubject mean performance","models mean performance",'Location',"northeast")



nexttile
rmsL_sub_mean=mean(Matrix_mean_model_perf,1);
rmsL_sub_std=std(Matrix_mean_model_perf,0,1);
hold on
plot(0:10,kemar_baseline(2)*ones(1,11),"Color","black",'LineWidth',2)
plot(0:10,mean_rmsL_top*ones(1,11),"Color","blue",'LineWidth', 2)
plot(0:10,mean_rmsL_baseline*ones(1,11),"Color","red",'LineWidth', 2)
plot(1:9,rmsL_sub_mean(1,2:10,2),'+','MarkerSize', 10,'Color',"black",'LineWidth', 2)
errorbar(rmsL_sub_mean(1,2:10,2),rmsL_sub_std(1,2:10,2),'.','Color','black');

y1 = mean_rmsL_baseline;
y2 = mean_rmsL_top;


fill([0 10 10 0], [y1 y1 y2 y2], ...
     'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

y11 = mean_rmsL_top;
y22 = 10;


fill([0 10 10 0], [y11 y11 y22 y22], ...
     [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

y11 = 30;
y22 = mean_rmsL_baseline;


fill([0 10 10 0], [y11 y11 y22 y22], ...
     'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

hold off

title('Lateral errors (°) of models')
set(gca,'XTick',1:9,'XTickLabel',model_list(1,2:10))
xlabel('Models')
legend('kemar baseline','self HRTFs mean performance',"intersubject mean performance","models mean performance",'Location',"northeast")



nexttile
rms_QU_sub_mean=mean(Matrix_mean_model_perf,1);
rms_QU_sub_std=std(Matrix_mean_model_perf,0,1);
hold on
plot(0:10,kemar_baseline(3)*ones(1,11),"Color","black",'LineWidth',2)
plot(0:10,mean_querr_top*ones(1,11),"Color","blue",'LineWidth', 2)
plot(0:10,mean_querr_baseline*ones(1,11),"Color","red",'LineWidth', 2)
plot(1:9,rms_QU_sub_mean(1,2:10,3),'+','MarkerSize', 10,'Color',"black",'LineWidth', 2)
errorbar(rms_QU_sub_mean(1,2:10,3),rms_QU_sub_std(1,2:10,3),'.','Color','black');

y1=kemar_baseline(3);
y2=mean_querr_baseline;

fill([0 10 10 0], [y1 y1 y2 y2], ...
     [0.01 0.39 0.25], 'FaceAlpha', 0.3, 'EdgeColor', 'none')



y1 = mean_querr_baseline;
y2 = mean_querr_top;


fill([0 10 10 0], [y1 y1 y2 y2], ...
     'green', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

y11 = mean_querr_top;
y22 = 5;


fill([0 10 10 0], [y11 y11 y22 y22], ...
     [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

y11 = 45;
y22 = kemar_baseline(3);
xl=xlim;

fill([0 10 10 0], [y11 y11 y22 y22], ...
     'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold off

title('Quadrant errors (%) of models')
set(gca,'XTick',1:9,'XTickLabel',model_list(1,2:10))
xlabel('Models')
legend('kemar baseline','self HRTFs mean performance',"intersubject mean performance","models mean  performance",'Location',"northeast")






figure('Color','w','Position',[100 100 1100 700])

tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
imagesc(Matrix_mean_model_perf(:,:,1).')
set(gca,'YDir','normal')
axis image
colormap(gca,"copper")
colorbar
title('Rms local polar error (°)','FontWeight','bold')
set(gca,'XTick',1:9,'XTickLabel',subject_list)
set(gca,'YTick',1:10,'YTickLabel',model_list)
xlabel('Reference HRTF')
ylabel('Target HRTF ')
hold on
for j=1:9
    for i=1:10
        if Matrix_mean_model_perf(j,i,1)<=kemar_baseline(1)
            plot(j,i, 'x', ...
            'Color', 'g', ...
            'LineWidth', 2.5, ...
            'MarkerSize', 12);        
        end
        if Matrix_mean_model_perf(j,i,1)<=mean_rmsP_baseline
           plot(j,i, 'x', ...
            'Color', 'b', ...
            'LineWidth', 2.5, ...
            'MarkerSize', 12);
        end
    end
end
hold off


nexttile
imagesc(Matrix_mean_model_perf(:,:,2).')
set(gca,'YDir','normal')
axis image
colormap(gca,"copper")
colorbar
title('Rms lateral error (deg)','FontWeight','bold')
set(gca,'XTick',1:9,'XTickLabel',subject_list)
set(gca,'YTick',1:10,'YTickLabel',model_list)
xlabel('Reference HRTF')
hold on
for j=1:9
    for i=1:10
        if Matrix_mean_model_perf(j,i,2)<=kemar_baseline(2)
            plot(j,i, 'x', ...
            'Color', 'g', ...
            'LineWidth', 2.5, ...
            'MarkerSize', 12);        
        end
        if Matrix_mean_model_perf(j,i,2)<=mean_rmsL_baseline
           plot(j,i, 'x', ...
            'Color', 'b', ...
            'LineWidth', 2.5, ...
            'MarkerSize', 12);
        end
    end
end
hold off


nexttile
imagesc(Matrix_mean_model_perf(:,:,3).')
set(gca,'YDir','normal')
axis image
colormap(gca,"copper")
colorbar
title('Quadrant error (%)','FontWeight','bold')
set(gca,'XTick',1:9,'XTickLabel',subject_list)
set(gca,'YTick',1:10,'YTickLabel',model_list)
xlabel('Reference HRTF')
hold on
for j=1:9
    for i=1:10
        if Matrix_mean_model_perf(j,i,3)<=kemar_baseline(3)
            plot(j,i, 'x', ...
            'Color', 'g', ...
            'LineWidth', 2.5, ...
            'MarkerSize', 12);        
        end
        if Matrix_mean_model_perf(j,i,3)<=mean_querr_baseline
           plot(j,i, 'x', ...
            'Color', 'b', ...
            'LineWidth', 2.5, ...
            'MarkerSize', 12);
        end
    end
end
hold off



save("matrice_erreur_intermodel.mat","Matrix_mean_model_perf");

