%% Bidirectional MR
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'eBMD'};

% Meta-analysis ----------------------------------------------------------
names  = {'ebi-a-GCST006979'};

betaM = zeros(length(names),1);
seM   = zeros(length(names),1);

for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\Bidirectional\gMR_res.xlsx'],Sheet= names{i});

    betaM(i) = t.Estimate(1);
    seM(i)   = t.SE(1);
end

f = figure(1); hold on; grid on; box on;
f.Position = [612 576 560 192];

errorbar(betaM,1,betaM - (betaM-1.96*seM), (betaM+1.96*seM) - betaM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]./255,"Color",[86 166 194]./255)
fill([-.3 .3 .3 -0.3],[0.5 0.5 1.5 1.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-.3 0.3];
ax.XTick = [-0.3:0.15:0.3];
ax.YLim = [0.5 1.5];
ax.FontSize = 11;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
%l = legend('GWAS','FontSize',12);
%l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";

xlabel("MR effect size per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
print(f,"SFig_bidirectional_1.png","-dpng","-r600")


%% Binary data
% clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'Fracture*','MI','CAD','IS','Hypertension','T2DM'};

% Meta-analysis ----------------------------------------------------------
names  = {'HF','ebi-a-GCST011365','ebi-a-GCST005194','IS','ukb-b-14057','ebi-a-GCST006867'};

oddM = zeros(length(names),1);
clM  = zeros(length(names),1);
cuM  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\Bidirectional\gMR_res.xlsx'],Sheet= names{i});

    oddM(i) = t.OR(1);
    clM(i)   = t.C1(1);
    cuM(i)   = t.C2(1);
end
oddM = flip(oddM); clM = flip(clM); cuM = flip(cuM);

f = figure(2); hold on; grid on; box on;
f.Position = [680 254 560 844];

errorbar(oddM,[1:6],oddM-clM,cuM-oddM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
fill([0.1 4 4 0.1],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.1 4];
ax.XTick = [0.1 0.25 0.5 1 2 4];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
xlabel("Odds ratio per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)

print(f,"SFig_bidirectional_2","-dpng","-r600")
