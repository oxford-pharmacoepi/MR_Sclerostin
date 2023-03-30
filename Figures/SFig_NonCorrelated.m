%% Continuous data
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apolipoprotein A',...
    'Apolipoprotein B','C-Reactive protein','Lipoprotein','HbA1c','Glucose'};

% Meta-analysis ----------------------------------------------------------
names  = {'eBMD'};

betaM = zeros(length(names),1);
seM   = zeros(length(names),1);

for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_gwas.xlsx'],Sheet= names{i});

    betaM(i) = t.b(1);
    seM(i)   = t.SE(1);
end

% UK Biobank -------------------------------------------------------------
names  = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apo-A',...
    'Apo-B','CRP','Lipoprotein','HbA1c','Glucose'};

betaUKB = zeros(length(names),1);
seUKB   = zeros(length(names),1);

for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_continuous.xlsx'],Sheet= names{i});

    betaUKB(i) = t.b(1);
    seUKB(i)   = t.se(1);
end
betaUKB = flip(betaUKB); seUKB = flip(seUKB); names = flip(names);

% ------------------------------------------------------------------------
f = figure(1); hold on; grid on; box on;
f.Position = [680 254 560 844];
errorbar(betaM,11+0.1,betaM - (betaM-1.96*seM), (betaM+1.96*seM) - betaM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]./255,"Color",[86 166 194]./255);
errorbar(betaUKB,[1:11]-0.1,betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)
fill([-0.6 1.6 1.6 -0.6],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)
ax = gca;
ax.YTick = 1:11;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-.5 1.5];
ax.YLim = [0.5 11.5];
ax.FontSize = 11;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
l = legend('GWAS','UK Biobank','FontSize',12);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
text(-0.8, 12.4, 'A)','Interpreter','latex','FontSize',17)

xlabel("MR effect size per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
print(f,"SFig_SingleSnp_1.png","-dpng","-r600")


%% Binary data
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'Fracture*','CAD','MI','IS','Hypertension','T2DM'};

% Meta-analysis ----------------------------------------------------------
names  = {'HF','CAD','MI','IS','Hypertension','T2DM'};

oddM = zeros(length(names),1);
clM  = zeros(length(names),1);
cuM  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_gwas.xlsx'],Sheet= names{i});

    oddM(i) = t.OR(1);
    clM(i)  = t.CI_LOW_OR(1);
    cuM(i)  = t.CI_HIGH_OR(1);
end
oddM = flip(oddM); clM = flip(clM); cuM = flip(cuM);

% UK Biobank - LOGISTIC --------------------------------------------------
names  = {'Fracture','CAD','MI','IS','Hypertension','T2DM'};

oddUKB = zeros(length(names),1);
clUKB  = zeros(length(names),1);
cuUKB  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_BinaryLR.xlsx'],Sheet= names{i});

    oddUKB(i) = exp(t.b(1));
    clUKB(i)  = exp(t.b(1)-1.96*t.se(1));
    cuUKB(i)  = exp(t.b(1)+1.96*t.se(1));
end
oddUKB = flip(oddUKB); clUKB = flip(clUKB); cuUKB = flip(cuUKB); names = flip(names);

% UK Biobank - SURVIVAL (Birth) ------------------------------------------
names  = {'Fracture','CAD','MI','IS','Hypertension','T2DM'};

oddUKB_B = zeros(length(names),1);
clUKB_B  = zeros(length(names),1);
cuUKB_B  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_BinarySA.xlsx'],Sheet= names{i});

    oddUKB_B(i) = exp(t.b(1));
    clUKB_B(i)  = exp(t.b(1)-1.96*t.se(1));
    cuUKB_B(i)  = exp(t.b(1)+1.96*t.se(1));
end
oddUKB_B = flip(oddUKB_B); clUKB_B = flip(clUKB_B); cuUKB_B = flip(cuUKB_B); 


f = figure(2); hold on; grid on; box on;
f.Position = [680 254 560 844];

errorbar(oddM,[1:6]+0.1,oddM-clM,cuM-oddM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(oddUKB,[1:6],oddUKB-clUKB,cuUKB-oddUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
errorbar(oddUKB_B,[1:6]-0.1,oddUKB_B-clUKB_B,cuUKB_B-oddUKB_B,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)

fill([0.05 8 8 0.05],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.05 8];
ax.XTick = [0.1 0.25 0.5 1 2 4 8];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('GWAS [Odds]','UKB-LR [Odds]','UKB-SA(Birth) [Hazard]','FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
xlabel("Ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.07, 7, 'B)','Interpreter','latex','FontSize',17)

print(f,"SFig_SingleSnp_2.png","-dpng","-r600")