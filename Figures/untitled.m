%% All the results
% Binary data -------------------------------------------------------------
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

n = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apolipoprotein A',...
    'Apolipoprotein B','C-Reactive protein','Lipoprotein(a)','HbA1c','Glucose'};
% Gwas
t     = readtable([path '\MR_gwas\gMR_res.xlsx'], Sheet = 'eBMD');
betaGWAS = t.Estimate; seGWAS = t.SE;

% UK Biobank 
names  = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apo-A','Apo-B','CRP','Lipoprotein','HbA1c','Glucose'};
betaUKB = zeros(length(names),1); seUKB = zeros(length(names),1);

for i = 1:length(names)
    t = readtable([path '\MR_UKBiobank\ContinuousData\gMR_res.xlsx'],Sheet= names{i});
    betaUKB(i) = t.Estimate(1); seUKB(i) = t.SE(1);
end
betaUKB = flip(betaUKB); seUKB = flip(seUKB); names = flip(names);
% -----------------------------------
f = figure(1); hold on; grid on; box on; f.Position = [680 254 560 844];
errorbar(betaGWAS,11+0.1,betaGWAS - (betaGWAS-1.96*seGWAS), (betaGWAS+1.96*seGWAS) - betaGWAS,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]./255,"Color",[86 166 194]./255)
errorbar(betaUKB,[1:11]-0.1,betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)
fill([-0.6 1.2 1.2 -0.6],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:11;
ax.YTickLabel = flip(n);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-.4 1.2];
ax.YLim = [0.5 11.5];
ax.FontSize = 11;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
l = legend('GWAS','UK Biobank','FontSize',12);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
text(-0.8, 12.4, 'A)','Interpreter','latex','FontSize',17)

xlabel("1-SD increase per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
print(f,[path '\Results\FigAll_1.png'],"-dpng","-r800")

%% Continuous data ---------------------------------------------------------
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

n = {'Fracture*','CAD','MI','IS','Hypertension','T2DM'};
% GWAS
names  = {'HF','CAD','MI','IS','Hypertension','T2DM'};
oddGWAS = zeros(length(names),1); clGWAS = zeros(length(names),1); cuGWAS = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\MR_gwas\gMR_res.xlsx'],Sheet= names{i});

    oddGWAS(i) = t.OR(1);
    clGWAS(i)  = t.CI_LOW_OR(1);
    cuGWAS(i)  = t.CI_HIGH_OR(1);
end
oddGWAS = flip(oddGWAS); clGWAS = flip(clGWAS); cuGWAS = flip(cuGWAS);

% UK Biobank - LOGISTIC --------------------------------------------------
names  = {'Fracture','MI','CAD','IS','Hypertension','T2DM'};
oddUKB = zeros(length(names),1); clUKB  = zeros(length(names),1); cuUKB  = zeros(length(names),1);
for i = 1:length(names)
    t = readtable([path '\MR_UKBiobank\BinaryData\Logistic\gMR_res.xlsx'],Sheet= names{i});

    oddUKB(i) = t.OR(1);
    clUKB(i)  = t.CI_LOW_OR(1);
    cuUKB(i)  = t.CI_HIGH_OR(1);
end
oddUKB = flip(oddUKB); clUKB = flip(clUKB); cuUKB = flip(cuUKB); 

f = figure(2); hold on; grid on; box on;
f.Position = [680 254 560 844];

errorbar(oddGWAS,[1:6]+0.1,oddGWAS-clGWAS,cuGWAS-oddGWAS,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(oddUKB,[1:6],oddUKB-clUKB,cuUKB-oddUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)

fill([0.1 4 4 0.1],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(n);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.1 4];
ax.XTick = [0.1 0.25 0.5 1 2 4];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('GWAS','UK Biobank','FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
xlabel("Odds Ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.07, 7, 'B)','Interpreter','latex','FontSize',17)

print(f,[path '\Results\FigAll_3.png'],"-dpng","-r800")

%% Leave one out analysis
clear all; close all;
path = pwd; path = path(1:end-length('\Figures'));

% UK Biobank -------------------------------------------------------------
names  = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apo-A','Apo-B','CRP','Lipoprotein','HbA1c','Glucose'};

f = figure(1); hold on; box on; axis off; grid on;
f.Position = [3 42 1169 1074];
lx = 0.11; mx = 0.02; rx = 0.01; dy = 0.06; my = 0.07; uy = 0.07; Nx = 4; Ny = 3;

[AX]  = get_axes([Nx,Ny],lx,mx,rx,dy,my,uy);
order = [4,8,12,3,7,11,2,6,10,1,5];
for i = 1:length(names)
    f.CurrentAxes = AX(order(i)); hold on; box on;

    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\ContinuousData\loo.xlsx'],Sheet= names{i});

    betaUKB = flip(t.Estimate); seUKB = flip(t.SE);
    if i == 1
        errorbar(betaUKB,[0.1:0.1:0.7]-0.02,betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)
    else
        errorbar(betaUKB,[0.1:0.1:0.7],betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)
    end
    ax = gca;
    ax.YTick = [0.1:0.1:0.7];
    ax.YTickLabel = '';
    ax.TickLabelInterpreter = 'latex';
    ax.YLim = [0 0.8];
    xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
    if i == 1
        t = readtable([path '\SensitivityAnalysis\LeaveOneOut\Gwas\loo.xlsx'],Sheet = 'eBMD');
        betaGwas = t.Estimate; seGwas = t.SE;
        errorbar(betaGwas,[0.1:0.1:0.7]+0.02, betaGwas - (betaGwas - 1.96*seGwas), (betaGwas + 1.96*seGwas) - betaGwas,'o','horizontal','LineWidth',1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255);
        ax.XLim(1) = -1.4; ax.XLim(2) = 1.4;
    else
       ax.XLim(1) = -.5; ax.XLim(2) = 0.5;
       ax.XTick = ax.XLim(1):0.25:ax.XLim(2);
    end

    % horizontal axis
    if order(i) == 1 || order(i) == 5 || order(i) == 10
        xlabel("MR effect size per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',11)
    end
    % vertical axis
    if order(i) == 1 || order(i) == 2 || order(i) == 3 || order(i) == 4
        ax.YTickLabel = flip(t.SNPS);
    end
        grid on; box on;
end

f.CurrentAxes = AX(9); box off; axis off;
annotation('textbox',[0.052 dy 1-uy-dy-my lx/2],'Color','k','String','Continuous outcomes','Interpreter','latex',...
    'FontSize',17,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90,'BackgroundColor',[200 200 200]./255,'EdgeColor',[200 200 200]./255)
Dx = (1-2*mx-lx-rx)/3; Dy = (1-uy-uy-3*my)/4;

% First row
annotation('textbox',[lx 1-0.9*uy Dx 0.03],'Color','k','String','A) Heel BMD','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[Dx+mx+lx 1-0.9*uy Dx 0.03],'Color','k','String','B) Cholesterol','Interpreter','latex',...
     'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[2*(Dx+mx)+lx 1-0.9*uy Dx 0.03],'Color','k','String','C) LDL cholesterol','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)

% Second row
annotation('textbox',[lx 1-uy-1*Dy-0.975*my Dx 0.03],'Color','k','String','D) HDL cholesterol','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[Dx+mx+lx 1-uy-1*Dy-0.975*my Dx 0.03],'Color','k','String','E) Triglycerides','Interpreter','latex',...
     'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[2*(Dx+mx)+lx 1-uy-1*Dy-0.975*my Dx 0.03],'Color','k','String','F) Apolipoprotein A','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)

% Third row
annotation('textbox',[lx 1-uy-2*Dy-1.975*my Dx 0.03],'Color','k','String','G) Apolipoprotein B','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[Dx+mx+lx 1-uy-2*Dy-1.975*my Dx 0.03],'Color','k','String','H) C-Reactive protein','Interpreter','latex',...
     'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[2*(Dx+mx)+lx 1-uy-2*Dy-1.975*my Dx 0.03],'Color','k','String','I) Lipoprotein (a)','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)

% Third row
annotation('textbox',[lx 1-uy-3*Dy-2.975*my Dx 0.03],'Color','k','String','J) HbA1c','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[Dx+mx+lx 1-uy-3*Dy-2.975*my Dx 0.03],'Color','k','String','K) Glucose','Interpreter','latex',...
     'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)

f.CurrentAxes = AX(9); 
errorbar(nan,nan,0.1,0.1,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255); hold on;
errorbar(nan,nan,0.1,0.1,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255);
box off; axis off;
l = legend('GWAS','UK Biobank','Interpreter','latex','FontSize',14);
l.NumColumns = 2; l.Position = [2*(Dx+mx)+1.85*lx Dy+dy/2+0.003 0.1 0.03];

print(f,[path '\Results\Fig_LOO_Continuous_1.png'],"-dpng","-r800")
%%
% Binary data ------------------------------------------------------------
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'Fracture','CAD','MI','IS','Hypertension','T2DM'};
namesM  = {'Hip fracture','CAD','MI','IS','Hypertension','T2DM'};
namesU  = {'Fracture','CAD','MI','IS','Hypertension','T2DM'};


f = figure(1); hold on; box on; axis off; grid on;
f.Position = [3 42 1913 1074];
lx = 0.09; mx = 0.015; rx = 0.01; dy = 0.1; my = 0.07; uy = 0.07; Nx = 2; Ny = 3;

[AX]  = get_axes([Nx,Ny],lx,mx,rx,dy,my,uy);
order = [2,4,6,1,3,5]; 

for i = 1:length(names1)
    f.CurrentAxes = AX(order(i));

    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\Gwas\loo.xlsx'],Sheet= namesM{i});
    betaM = flip(t.Estimate);
    seM   = flip(t.SE);
    if i == 6
        betaM = [betaM(1:2); NaN; betaM(3:6)];
        seM   = [seM(1:2); NaN; seM(3:6)];
    end
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\BinaryData_LR\loo.xlsx'],Sheet= namesU{i});
    betaLR = flip(t.Estimate);
    seLR    = flip(t.SE);
    
    SNPs    = flip(t.SNPS);

    hold on; grid on; box on;
    errorbar(betaM,[0.1:0.1:0.7]+0.01,betaM - (betaM - 1.96*seM),(betaM + 1.96*seM) - betaM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
    errorbar(betaLR,[0.1:0.1:0.7]-0.01,betaLR - (betaLR - 1.96*seLR),(betaLR + 1.96*seLR) - betaLR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
   
    ax = gca;
    ax.YTick = [0.1:0.1:0.7];
    ax.YTickLabel = '';
    ax.TickLabelInterpreter = 'latex';
    ax.YLim = [0 0.8];
    if i == 1
        ax.XLim = [-2.5 2.5];
    else
        ax.XLim = [-1 1];
    end

    xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
    if i == 4 || i == 5 || i == 6
        xlabel("MR effect size per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
    end
    if i == 1 || i == 4
        ax.YTickLabel = SNPs;
    end
end
Dx = (1-2*mx-lx-rx)/3; Dy = (1-uy-uy-my)/2;

annotation('textbox',[0.05 dy 0.465 lx-0.02],'Color','k','String','Categorical outcomes','Interpreter','latex',...
    'FontSize',17,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90,'BackgroundColor',[200 200 200]./255,'EdgeColor',[200 200 200]./255)
% First column
annotation('textbox',[lx 1-0.9*uy Dx 0.03],'Color','k','String','A) Fracture*','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[lx 1-uy-Dy-0.7*my Dx 0.03],'Color','k','String','D) Ischemic stroke','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
% Second column
annotation('textbox',[Dx+mx+lx 1-0.9*uy Dx 0.03],'Color','k','String','B) Coronary artery disease','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[Dx+mx+lx 1-uy-Dy-0.7*my Dx 0.03],'Color','k','String','E) Hypertension','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
% Third column
annotation('textbox',[2*(Dx+mx)+lx 1-0.9*uy Dx 0.03],'Color','k','String','C) Myocardial infarction','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[2*(Dx+mx)+lx 1-uy-Dy-0.7*my Dx 0.03],'Color','k','String','F) Type II diabetes mellitus','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)

l = legend('GWAS','UK Biobank','Interpreter','latex','FontSize',14);
l.NumColumns = 3; l.Position = [1/2+0.03 dy/6 .03 .03];
print(f,[path '\Results\Fig_LOO_Binary_1.png'],"-dpng","-r800")

%% Uncorrelated SNPs
% Continuous data ---------------------------------------------------------
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

    betaM(i) = t.Estimate(1);
    seM(i)   = t.SE(1);
end

% UK Biobank -------------------------------------------------------------
names  = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apo-A',...
    'Apo-B','CRP','Lipoprotein','HbA1c','Glucose'};

betaUKB = zeros(length(names),1);
seUKB   = zeros(length(names),1);

for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\SingleSNP\MR_res_continuous.xlsx'],Sheet= names{i});

    betaUKB(i) = t.Estimate(1);
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
print(f,[path '\Results\Fig_SingleSNP_1.png'],"-dpng","-r800")
%%
% Binary outcome ----------------------------------------------------------
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

    oddUKB(i) = t.OR;
    clUKB(i)  = t.CI_LOW_OR;
    cuUKB(i)  = t.CI_HIGH_OR;
end
oddUKB = flip(oddUKB); clUKB = flip(clUKB); cuUKB = flip(cuUKB); names = flip(names);

f = figure(2); hold on; grid on; box on;
f.Position = [680 254 560 844];

errorbar(oddM,[1:6]+0.1,oddM-clM,cuM-oddM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(oddUKB,[1:6],oddUKB-clUKB,cuUKB-oddUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)

fill([0.0005 24 24 0.0005],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(names1);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.01 6.5];
ax.XScale = 'log';
ax.XLim = [0.05 12];
ax.XTick = [ 0.01 0.1 0.25 0.5 1 2 4 8 12];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('GWAS','UK Biobank','FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside";
xlabel("Odds Ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.0005, 7, 'B)','Interpreter','latex','FontSize',17)

print(f,[path '\Results\Fig_SingleSNP_3.png'],"-dpng","-r800")
%% FUNCTIONs
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



