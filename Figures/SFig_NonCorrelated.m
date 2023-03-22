
%% Noncorrelated SNPs
clear all; close all;
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apolipoprotein A',...
    'Apolipoprotein B','C-Reactive protein','Lipoprotein','HbA1c','Glucose'};

f = figure(1); axis off;
f.Position = [3 42 1913 1074];

AX = get_axes([1,2],0.1,0.1,0.01,0.07,.01,0.06);


names1 = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apolipoprotein A',...
    'Apolipoprotein B','C-Reactive protein','Lipoprotein','HbA1c','Glucose'};
% Meta-analysis ----------------------------------------------------------
names  = {'ebi-a-GCST006979'};

betaM = zeros(length(names),1);
seM   = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_metaanalysis.xlsx'],Sheet= names{i});

    betaM(i) = t.Estimate(1);
    seM(i)   = t.SE(1);
end

% UK Biobank -------------------------------------------------------------
names  = {'30690-0.0','30780-0.0','30760-0.0','30870-0.0','30630-0.0',...
    '30640-0.0','30710-0.0','30790-0.0','30750-0.0','30740-0.0'};

betaUKB = zeros(length(names),1);
seUKB   = zeros(length(names),1);

for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_continuous.xlsx'],Sheet= names{i});

    betaUKB(i) = t.b(1);
    seUKB(i)   = t.se(1);
end
betaUKB = flip(betaUKB); seUKB = flip(seUKB); names = flip(names);

f.CurrentAxes = AX(1); hold on; box on; grid on;
% ------------------------------------------------------------------------
errorbar(betaM,11,betaM - (betaM-1.96*seM), (betaM+1.96*seM) - betaM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]./255,"Color",[86 166 194]./255)
errorbar(betaUKB,[1:10],betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)
fill([-1.2 1.2 1.2 -1.2],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:11;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-1.2 1.2];
ax.YLim = [0.5 11.5];
ax.FontSize = 12;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
text(-1.2, 11.5, 'A)','Interpreter','latex','FontSize',17,'VerticalAlignment','bottom')
xlabel("MR effect size per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',14)

l = legend('GWAS','UK Biobank','Interpreter','latex','Position',[1/2+0.01 1-0.05 0.07 0.025],'FontSize',14);
l.NumColumns = 2;

% Binary data--------------------------------------------------------------
names1 = {'Fracture','MI','CAD','IS','Hypertension','T2DM'};

% Meta-analysis ----------------------------------------------------------
names  = {'HF','ebi-a-GCST011365','ebi-a-GCST005194','IS','ukb-b-14057','ebi-a-GCST006867'};

oddM = zeros(length(names),1);
clM  = zeros(length(names),1);
cuM  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_metaanalysis.xlsx'],Sheet= names{i});

    oddM(i) = t.OR(1);
    clM(i)   = t.C1(1);
    cuM(i)   = t.C2(1);
end
oddM = flip(oddM); clM = flip(clM); cuM = flip(cuM);

f.CurrentAxes = AX(2); box on; grid on; hold on;

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
ax.FontSize = 12;
xlabel("Odds ratio per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',14)
text(4, 6.5, 'B)','Interpreter','latex','FontSize',17,'VerticalAlignment','bottom','HorizontalAlignment','right')

print(f,"SFig_SingleSNP.png","-dpng","-r600")


%% Noncorrelated SNPs
clear all; close all;
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apolipoprotein A',...
    'Apolipoprotein B','C-Reactive protein','Lipoprotein','HbA1c','Glucose'};

% Meta-analysis ----------------------------------------------------------
names  = {'ebi-a-GCST006979'};

betaM = zeros(length(names),1);
seM   = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_metaanalysis.xlsx'],Sheet= names{i});

    betaM(i) = t.Estimate(1);
    seM(i)   = t.SE(1);
end

% UK Biobank -------------------------------------------------------------
names  = {'30690-0.0','30780-0.0','30760-0.0','30870-0.0','30630-0.0',...
    '30640-0.0','30710-0.0','30790-0.0','30750-0.0','30740-0.0'};

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

errorbar(betaM,11,betaM - (betaM-1.96*seM), (betaM+1.96*seM) - betaM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]./255,"Color",[86 166 194]./255)
errorbar(betaUKB,[1:10],betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)
fill([-0.6 1.2 1.2 -0.6],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:11;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-.6 1.2];
ax.YLim = [0.5 11.5];
ax.FontSize = 11;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
l = legend('GWAS','UK Biobank','FontSize',12);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
text(-0.8, 12.4, 'A)','Interpreter','latex','FontSize',17)

xlabel("MR effect size per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
print(f,"SFig_SingleSNP_1.png","-dpng","-r600")

%% Binary data
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'Fracture','MI','CAD','IS','Hypertension','T2DM'};

% Meta-analysis ----------------------------------------------------------
names  = {'HF','ebi-a-GCST011365','ebi-a-GCST005194','IS','ukb-b-14057','ebi-a-GCST006867'};

oddM = zeros(length(names),1);
clM  = zeros(length(names),1);
cuM  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\MR_metaanalysis\gMR_res.xlsx'],Sheet= names{i});

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
l = legend('GWAS','FontSize',12);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
text(0.07, 7, 'B)','Interpreter','latex','FontSize',17)
print(f,"SFig_SingleSNP_2.png","-dpng","-r600")
%%
function [AX] = get_axes(N,lx,mx,rx,dy,my,uy)
 
Nx = N(2);
Ny = N(1);
 
DX = (1-lx-(Nx-1)*mx-rx)/Nx;
DY = (1-dy-(Ny-1)*my-uy)/Ny;
 
for i = 1:Nx
    for j = 1:Ny
        if i==1 && j==1
            AX = axes('Units','normalized','Position',[lx+(i-1)*(DX+mx) dy+(j-1)*(DY+my) DX DY]);
        else
            AX = [ AX; axes('Units','normalized','Position',[lx+(i-1)*(DX+mx) dy+(j-1)*(DY+my) DX DY])];
        end
    end
end
 
end