% Survival-analysis
% Binary data
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'Fracture','MI','CAD','IS','Hypertension','T2DM'};

% UK Biobank - SURVIVAL (Birth) ------------------------------------------
names  = {'Fracture','MI','CAD','IS','Hypertension','T2DM'};

oddUKB_B = zeros(length(names),1);
clUKB_B  = zeros(length(names),1);
cuUKB_B  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\MR_UKBiobank\BinaryData\SA_Birth\MR_res.xlsx'],Sheet= names{i});

    oddUKB_B(i) = t.OR(1);
    clUKB_B(i)  = t.CI_LOW_OR(1);
    cuUKB_B(i)  = t.CI_HIGH_OR(1);
end
oddUKB_B = flip(oddUKB_B); clUKB_B = flip(clUKB_B); cuUKB_B = flip(cuUKB_B); 

% UK Biobank - SURVIVAL (Enrolment) ------------------------------------------
oddUKB_E = zeros(length(names),1);
clUKB_E = zeros(length(names),1);
cuUKB_E  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\Enrolment\MR_res.xlsx'],Sheet= names{i});

    oddUKB_E(i) = t.OR(1);
    clUKB_E(i)  = t.CI_LOW_OR(1);
    cuUKB_E(i)  = t.CI_HIGH_OR(1);
end
oddUKB_E = flip(oddUKB_E); clUKB_E = flip(clUKB_E); cuUKB_E = flip(cuUKB_E); 

f = figure(2); hold on; grid on; box on;
f.Position = [680 254 560 844];

errorbar(oddUKB_B,[1:6]+0.05,oddUKB_B-clUKB_B,cuUKB_B-oddUKB_B,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)
errorbar(oddUKB_E,[1:6]-0.05,oddUKB_E-clUKB_E,cuUKB_E-oddUKB_E,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[237 177 32]/255,"Color",[237 177 32]/255)

fill([0.05 4 4 0.05],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.05 4];
ax.XTick = [0.1 0.25 0.5 1 2 4];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('SA since birth date','SA since UKB enrolment','FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
xlabel("Hazard ratio per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)

print(f,"SFig_SAEnrolment.png","-dpng","-r600")
