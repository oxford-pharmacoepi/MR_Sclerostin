% Leave one out analysis
clear all; close all;
path = pwd; path = path(1:end-length('\Figures'));

% UK Biobank -------------------------------------------------------------
names  = {'3148-0.0','30690-0.0','30780-0.0','30760-0.0','30870-0.0','30630-0.0',...
    '30640-0.0','30710-0.0','30790-0.0','30750-0.0','30740-0.0'};

f = figure(1); hold on; box on; axis off; grid on;
f.Position = [3 42 1913 1074];
lx = 0.09; mx = 0.015; rx = 0.01; dy = 0.06; my = 0.07; uy = 0.07; Nx = 3; Ny = 4;

[AX]  = get_axes([Nx,Ny],lx,mx,rx,dy,my,uy);
order = [3,6,9,12,2,5,8,11,1,4,7,10]; 
for i = 1:length(names)
    f.CurrentAxes = AX(order(i)); hold on; box on;
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\ContinuousData\loo.xlsx'],Sheet= names{i});

    betaUKB = flip(t.Estimate); seUKB = flip(t.SE);
    errorbar(betaUKB,[0.1:0.1:0.7],betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)
    ax = gca;
    ax.YTick = [0.1:0.1:0.7];
    ax.YTickLabel = '';
    ax.TickLabelInterpreter = 'latex';
    ax.YLim = [0 0.8];
    xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
    if i == 1
        ax.XLim(1) = -1.4; ax.XLim(2) = 1.4;
    else
       ax.XLim(1) = -.45; ax.XLim(2) = 0.45;
    end

    % horizontal axis
    if i == 11 || i == 9 ||i == 10 || i == 8
        xlabel("MR effect size per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
    end
    % vertical axis
    if i == 1 || i == 5 || i == 9 
        ax.YTickLabel = flip(t.SNPS);
    end
        grid on; box on;
end

f.CurrentAxes = AX(10); box off; axis off;
annotation('textbox',[0.05 dy 0.4875 lx-0.01],'Color','k','String','Continuous data','Interpreter','latex',...
    'FontSize',17,'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90,'BackgroundColor',[200 200 200]./255,'EdgeColor',[200 200 200]./255)
Dx = (1-3*mx-lx-rx)/4; Dy = (1-uy-uy-2*my)/3;
% First column
annotation('textbox',[lx 1-0.9*uy Dx 0.03],'Color','k','String','A) Heel BMD','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[lx 1-uy-Dy-0.95*my Dx 0.03],'Color','k','String','E) HDL cholesterol','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[lx 1-uy-2*Dy-1.975*my Dx 0.03],'Color','k','String','I) Lipoprotein (a)','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
% Second column
annotation('textbox',[Dx+mx+lx 1-0.9*uy Dx 0.03],'Color','k','String','B) Cholesterol','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[Dx+mx+lx 1-uy-Dy-0.95*my Dx 0.03],'Color','k','String','F) Apolipoprotein A','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[Dx+mx+lx 1-uy-2*Dy-1.975*my Dx 0.03],'Color','k','String','J) HbA1c','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
% Third column
annotation('textbox',[2*(Dx+mx)+lx 1-0.9*uy Dx 0.03],'Color','k','String','C) Triglycerides','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[2*(Dx+mx)+lx 1-uy-Dy-0.95*my Dx 0.03],'Color','k','String','G) Apolipoprotein B','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[2*(Dx+mx)+lx 1-uy-2*Dy-1.975*my Dx 0.03],'Color','k','String','K) Glucose','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
% Forth column
annotation('textbox',[3*(Dx+mx)+lx 1-0.9*uy Dx 0.03],'Color','k','String','D) LDL cholesterol','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)
annotation('textbox',[3*(Dx+mx)+lx 1-uy-Dy-0.95*my Dx 0.03],'Color','k','String','H) C-Reactive protein','Interpreter','latex',...
    'FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[220 220 220]./255,'EdgeColor',[220 220 220]./255)


print(f,['./SFig_LOO_Continuous.png'],"-dpng","-r300")

%% Binary data
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'Fracture','CAD','MI','IS','Hypertension','T2DM'};
namesM  = {'HF','ebi-a-GCST005194','ebi-a-GCST011365','IS','ukb-b-14057','ebi-a-GCST006867'};
namesU  = {'Fracture','CAD','MI','IS','Hypertension','T2DM'};


f = figure(1); hold on; box on; axis off; grid on;
f.Position = [3 42 1913 1074];
lx = 0.09; mx = 0.015; rx = 0.01; dy = 0.1; my = 0.07; uy = 0.07; Nx = 2; Ny = 3;

[AX]  = get_axes([Nx,Ny],lx,mx,rx,dy,my,uy);
order = [2,4,6,1,3,5]; 

for i = 1:length(names1)
    f.CurrentAxes = AX(order(i));

    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\Metaanalysis\loo.xlsx'],Sheet= namesM{i});
    betaM = flip(t.Estimate);
    seM   = flip(t.SE);
    if i == 6
        betaM = [betaM(1:2); NaN; betaM(3:6)];
        seM   = [seM(1:2); NaN; seM(3:6)];
    end
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\BinaryData_LR\loo.xlsx'],Sheet= namesU{i});
    betaLR = flip(t.Estimate);
    seLR    = flip(t.SE);
    
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\BinaryData_SA\loo.xlsx'],Sheet= namesU{i});
    betaSA = flip(t.Estimate);
    seSA   = flip(t.SE);
   
    SNPs    = flip(t.SNPS);

    hold on; grid on; box on;
    errorbar(betaM,[0.1:0.1:0.7]+0.02,betaM - (betaM - 1.96*seM),(betaM + 1.96*seM) - betaM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
    errorbar(betaLR,[0.1:0.1:0.7],betaLR - (betaLR - 1.96*seLR),(betaLR + 1.96*seLR) - betaLR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
    errorbar(betaSA,[0.1:0.1:0.7]-0.02,betaSA - (betaSA - 1.96*seSA),(betaSA + 1.96*seSA) - betaSA,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)

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
        xlabel("MR effect size per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
    end
    if i == 1 || i == 4
        ax.YTickLabel = SNPs;
    end
end
Dx = (1-2*mx-lx-rx)/3; Dy = (1-uy-uy-my)/2;

annotation('textbox',[0.05 dy 0.465 lx-0.02],'Color','k','String','Categorical data','Interpreter','latex',...
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

l = legend('GWAS','UKB-LR','UKB-SA','Interpreter','latex','FontSize',14);
l.NumColumns = 3; l.Position = [1/2+0.03 dy/6 .03 .03];
print(f,['./SFig_LOO_Binary.png'],"-dpng","-r600")

%% Continuous data
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'eBMD','Cholesterol','LDL','HDL','Triglycerides','Apolipoprotein A',...
    'Apolipoprotein B','C-Reactive protein','Lipoprotein','HbA1c','Glucose'};

% UK Biobank -------------------------------------------------------------
names  = {'3148-0.0','30690-0.0','30780-0.0','30760-0.0','30870-0.0','30630-0.0',...
    '30640-0.0','30710-0.0','30790-0.0','30750-0.0','30740-0.0'};

betaUKB = zeros(length(names),1);
seUKB   = zeros(length(names),1);

for i = 1:length(names)
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\ContinuousData\loo.xlsx'],Sheet= names{i});

    betaUKB = flip(t.Estimate);
    seUKB   = flip(t.SE);
    SNPs    = flip(t.SNPS);

    f = figure(i); hold on; grid on; box on;

    errorbar(betaUKB,[0.1:0.1:0.7],betaUKB - (betaUKB - 1.96*seUKB),(betaUKB + 1.96*seUKB) - betaUKB,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]./255,"Color",[167 21 49]./255)

    ax = gca;
    ax.YTick = [0.1:0.1:0.7];
    ax.YTickLabel = SNPs;
    ax.TickLabelInterpreter = 'latex';
    ax.YLim = [0 0.8];
    xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
    xlabel("MR effect size per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
    print(f,['./SFig_LOO_Continuous/Fig_' names1{i} '.png'],"-dpng","-r600")
end

%% Binary data
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

names1 = {'Fracture','MI','CAD','IS','Hypertension','T2DM'};
namesM  = {'HF','ebi-a-GCST011365','ebi-a-GCST005194','IS','ukb-b-14057','ebi-a-GCST006867'};
namesU  = {'Fracture','MI','CAD','IS','Hypertension','T2DM'};

for i = 1:length(names1)
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\Metaanalysis\loo.xlsx'],Sheet= namesM{i});
    betaM = flip(t.Estimate);
    seM   = flip(t.SE);
    if i == 6
        betaM = [betaM(1:2); NaN; betaM(3:6)];
        seM   = [seM(1:2); NaN; seM(3:6)];
    end
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\BinaryData_LR\loo.xlsx'],Sheet= namesU{i});
    betaLR = flip(t.Estimate);
    seLR    = flip(t.SE);
    
    t = readtable([path '\SensitivityAnalysis\LeaveOneOut\BinaryData_SA\loo.xlsx'],Sheet= namesU{i});
    betaSA = flip(t.Estimate);
    seSA   = flip(t.SE);
   
    SNPs    = flip(t.SNPS);

    f = figure(i); hold on; grid on; box on;
    errorbar(betaM,[0.1:0.1:0.7]+0.02,betaM - (betaM - 1.96*seM),(betaM + 1.96*seM) - betaM,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
    errorbar(betaLR,[0.1:0.1:0.7],betaLR - (betaLR - 1.96*seLR),(betaLR + 1.96*seLR) - betaLR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
    errorbar(betaSA,[0.1:0.1:0.7]-0.02,betaSA - (betaSA - 1.96*seSA),(betaSA + 1.96*seSA) - betaSA,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)

    ax = gca;
    ax.YTick = [0.1:0.1:0.7];
    ax.YTickLabel = SNPs;
    ax.TickLabelInterpreter = 'latex';
    ax.YLim = [0 0.8];
    xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
    xlabel("MR effect size per SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
    print(f,['./SFig_LOO_Binary/Fig_' names1{i} '.png'],"-dpng","-r600")
end

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
