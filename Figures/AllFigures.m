%% Supplementary Figure 3
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';
n1 = {"eBMD", "hip fracture"};

% Heel bone mineral density
for i = 1
    t = readtable([path '\Study1\InstrumentSelection\Validation\Validation_Fixedstep' num2str(i) '.txt']);
    t = t(strcmp(t.SNP,'rs66838809') | strcmp(t.SNP,'rs7220711'),:);

    t.beta_outcome(t.beta_exposure < 0)  = -t.beta_outcome(t.beta_exposure < 0);
    t.beta_exposure(t.beta_exposure < 0) = -t.beta_exposure(t.beta_exposure < 0);

    betax = t.beta_exposure; sex = t.se_exposure; betay = t.beta_outcome; sey = t.se_outcome;

    f = figure(1); clf; grid on; hold on; box on;
    errorbar(betax,betay,betay-(betay-1.96*sey),(betay+1.96*sey)-betay,betax-(betax-1.96*sex),(betax+1.96*sex)-betax,"o","LineWidth",0.5,"MarkerFaceColor",'k',"Color",'k');
    ax = gca;
    ax.XLim(1) = -0.12; ax.XLim(2) = 0.12;
    ax.YLim(1) = -0.12; ax.YLim(2) = 0.12;

    ax.TickLabelInterpreter = 'latex';
    text(-0.145,0.13,'A)','FontSize',17,'Interpreter','latex')
    text(betax*1.075,betay,t.SNP,'FontSize',11,'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom')

    xline(0); yline(0);

    xlabel("SNP effect size per 1-SD change in sclerostin levels","Interpreter","latex")
    ylabel(['SNP effect size per 1-SD change in ' n1{i}],"Interpreter","latex")

    print(f,[path '\Figures\SFig3_A.png'],"-dpng","-r800")
end

% Hip fracture risk
for i = 2
    t = readtable([path '\Study1\InstrumentSelection\Validation\Validation_Fixedstep' num2str(i) '.txt']);
    t = t(strcmp(t.SNP,'rs66838809') | strcmp(t.SNP,'rs7220711'),:);

    t.beta_outcome(t.beta_exposure < 0)  = -t.beta_outcome(t.beta_exposure < 0);
    t.beta_exposure(t.beta_exposure < 0) = -t.beta_exposure(t.beta_exposure < 0);

    betax = t.beta_exposure; sex = t.se_exposure; betay = t.beta_outcome; sey = t.se_outcome;

    f = figure(2); clf; grid on; hold on; box on;
    errorbar(betax,betay,betay-(betay-1.96*sey),(betay+1.96*sey)-betay,betax-(betax-1.96*sex),(betax+1.96*sex)-betax,"o","LineWidth",0.5,"MarkerFaceColor",'k',"Color",'k');
    ax = gca;
    ax.XLim(1) = -0.12; ax.XLim(2) = 0.12;
    ax.YLim(1) = -0.24; ax.YLim(2) = 0.24;

    ax.TickLabelInterpreter = 'latex';
    text(-0.145,0.26,'B)','FontSize',17,'Interpreter','latex')
    text(betax*0.925,betay,t.SNP,'FontSize',11,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','bottom')

    xline(0); yline(0);

    xlabel("SNP effect size per 1-SD change in sclerostin levels","Interpreter","latex")
    ylabel(['SNP effect size per log(OR) change in ' n1{i}],"Interpreter","latex")
    print(f,[path '\Figures\SFig3_B.png'],"-dpng","-r800")
end


%% Figure 1.A
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';

% GWAS results
n1 = {'bmd','ldl','hdl','glucose', 'hba1c'};

t1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
    t1  = [t1; t0];
end

% UK Biobank results
n1 = {'cholesterol','ldl','hdl','triglycerides', 'apoA', 'apoB', 'crp', 'lipoprotein',...
    'glucose', 'hba1c'};
t2 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logistic_' n1{1} '.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logistic_' n1{i} '.txt']);
    t2  = [t2; t0];
end

n0 = {'Heel bone mineral density', 'Cholesterol','LDL-Cholesterol','HDL-Cholesterol','Triglycerides',...
    'Apolipoprotein-A', 'Apolipoprotein-B','C-Reactive protein','Lipoprotein (a)','Glucose','HbA1c'};
f = figure(1); hold on; grid on; box on;
f.Position = [392 250 560 677];
errorbar(t1.Beta, [11,9,8,2,1]+0.1, t1.Beta-t1.L95CI, t1.U95CI-t1.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(t2.Beta, flip([1:10])-0.1, t2.Beta-t2.L95CI, t2.U95CI-t2.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)

fill([-0.6 1.2 1.2 -0.6],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:11;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-.4 1.2]; ax.XTick = [-0.4:0.2:1.2];
ax.YLim = [0.5 11.5];
ax.FontSize = 11;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
l1 = legend('GWAS','UK Biobank','FontSize',12);
l1.Interpreter = 'latex'; l1.NumColumns = 3; l1.Location = "northoutside";
text(-1, 12.1, 'A)','Interpreter','latex','FontSize',17)
xlabel("Unit increase per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)

print(f,[path '\Figures\Fig1_A.png'],"-dpng","-r800")

%% Figure 1.B
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';

% GWAS results
n1 = {'hf','cad','mi','is','hypertension', 't2dm'};

t1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
    t1  = [t1; t0];
end

% UK Biobank results
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t2 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logaritmic_' n1{1} '.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logaritmic_' n1{i} '.txt']);
    t2  = [t2; t0];
end

% UK Biobank results - survival
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t3 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{1} '.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{i} '.txt']);
    t3 = [t3; t0];
end


n0 = {'Hip fracture', 'Coronary artery disease','Myocardial infarction','Ischaemic stroke',...
    'Hypertension','Type 2 diabetes'};
f = figure(2); hold on; grid on; box on;
f.Position = [680 254 560 644];
errorbar(t1.OR, flip([1:6])+0.1, t1.OR-t1.L95CI, t1.U95CI-t1.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(t2.OR, flip([1:5]), t2.OR-t2.L95CI, t2.U95CI-t2.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
errorbar(t3.OR, flip([1:5])-0.1, t3.OR-t3.L95CI, t3.U95CI-t3.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)

fill([0.075 4 4 0.075],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.075 4];
ax.XTick = [0.1 0.25 0.5 1 2 4];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('GWAS [Odds]','UKB-LR [Odds]','UKB-SA(Birth) [Hazard]','FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside"; l.Position = [0.1050 0.9370 0.8804 0.0373];
xlabel("Ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.015, 6.8, 'B)','Interpreter','latex','FontSize',17)

print(f,[path '\Figures\Fig1_B.png'],"-dpng","-r800")

%% SFig5
clear all; close all;

path = 'D:\Projects\MR_Sclerostin\Results\';

% UK Biobank results
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t2 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{1} '.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{i} '.txt']);
    t2  = [t2; t0];
end

% UK Biobank results - survival
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t3 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Enrol_' n1{1} '.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Enrol_' n1{i} '.txt']);
    t3 = [t3; t0];
end


n0 = {'Coronary artery disease','Myocardial infarction','Ischaemic stroke',...
    'Hypertension','Type 2 diabetes'};
f = figure(2); hold on; grid on; box on;
f.Position = [680 254 560 644];
errorbar(t2.OR, flip([1:5]), t2.OR-t2.L95CI, t2.U95CI-t2.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)
errorbar(t3.OR, flip([1:5])-0.1, t3.OR-t3.L95CI, t3.U95CI-t3.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[237 177 32]/255,"Color",[237 177 32]/255)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.075 4];
ax.XTick = [0.1 0.25 0.5 1 2 4];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('SA since birth date','SA since UKB first assessment','FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside"; l.Position = [0.1050 0.9370 0.8804 0.0373];
xlabel("Hazard ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)

print(f,[path '\Figures\SFig5.png'],"-dpng","-r800")

%% SFig6
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';

outcomes = {'hf','cad','mi','is','hypertension','t2dm','bmd','ldl','hdl','glucose','hba1c'};

outcomes_names = {'Hip fracture','Coronary artery disease','Myocardial infarction',...
    'Ischaemic stroke','Hypertension','Type 2 diabetes','eBMD','LDL-Cholesterol', 'HDL-Cholesterol',...
    'Fasting glucose','HbA1c'};

alt = [7,7,7,7,7,7,100,10,90,7,7];
let = {'A)','B)','C)','D)','E)','F)','G)','H)','I)','J)','K)'};
order = [3, 6, 9, 12, 2, 5, 8, 11, 1,4,7,10];
f = figure(1); hold on; box on; axis off;
f.Position = [9 79  1924 890];

AX = get_axes([3,4],0.02,0.017,0.01,0.07,0.05,0.03);

genstart = 43753738;
genend   = 43758791;
for i = 1:length(let)
    t = readtable([path 'Colocalization\' outcomes{i} '.txt']);
    pos = t.pos_exposure;
    pval.outcome  = t.pval_outcome;
    pval.exposure = t.pval_exposure;

    f.CurrentAxes = AX(order(i));
    scatter(pos,-log10(pval.exposure),20,'k','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
    scatter(pos,-log10(pval.outcome),20,'r','MarkerFaceColor','r')
    ax = gca;
    ax.XTick = [genstart-20*1000, genstart, genend, genend+20*1000];
    ax.XTickLabel = '';
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 8.5;
    ax.XLim = [genstart-21*1000 genend+21*1000];

    fill([genstart, genend, genend, genstart],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[237 177 32]./255,'EdgeColor',[237 177 32]./255,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    fill([genstart-20*1000, genstart, genstart, genstart-20*1000],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[237 177 32]./255,'EdgeColor',[237 177 32]./255,'FaceAlpha',0.05,'EdgeAlpha',0)
    fill([genend, genend+20*1000, genend+20*1000, genend],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[237 177 32]./255,'EdgeColor',[237 177 32]./255,'FaceAlpha',0.05,'EdgeAlpha',0)

    if order(i) == 1 || order(i) == 4 || order(i) == 7 || order(i) == 11
        xlabel('Chromosome 17 position [bp]','Interpreter','latex','FontSize',14)
        ax.XTick = [genstart-20*1000, genstart, genend, genend+20*1000];
        ax.XAxis.Exponent = 0;
        xtickformat('%.0f');
    end

    if order(i) == 1 || order(i) == 2 || order(i) == 3
        ylabel('-$log_{10}$(P-Value)','Interpreter','latex','FontSize',14)
    end

    if order(i) == 1
        ax.YLim(2) = 90;
    end
    text(genstart-20*1000,alt(i),let{i},'Interpreter','latex','FontSize',17,'VerticalAlignment','bottom')
    l = legend('Sclerostin',outcomes_names{i},'Interpreter','latex','FontSize',11,'FontWeight','bold','Location','northeast');
end

f.CurrentAxes = AX(10);  box off; axis off;
print(f,[path 'Figures\SFig6.png'],"-dpng","-r800")


%% Step-wise pruning - continuous
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';
r2 = {'0.001', '0.8'};

for j = 1:length(r2)

f = figure(j); hold on; box on; axis off;
AX = get_axes([1,2],0.1,0.12,0.01,0.07,0.1,0.1);
f.Position = [9 60 1767 920];

    % GWAS results
    n1 = {'bmd','ldl','hdl','glucose', 'hba1c'};

    t1   = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
    t1_1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_' r2{j} '_250000.txt']);

    for i = 2:length(n1)
        t0   = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
        t0_1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_' r2{j} '_250000.txt']);
        t1   = [t1; t0];
        t1_1 = [t1_1; t0_1];
    end
    % UK Biobank results
    n1 = {'cholesterol','ldl','hdl','triglycerides', 'apoA', 'apoB', 'crp', 'lipoprotein',...
        'glucose', 'hba1c'};
    t2   = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logistic_' n1{1} '.txt']);
    t2_1 = readtable([path 'Study2\MendelianRandomisation\Fixed_' r2{j} '_Logistic_' n1{1} '.txt']);
    for i = 2:length(n1)
        t0   = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logistic_' n1{i} '.txt']);
        t0_1 = readtable([path 'Study2\MendelianRandomisation\Fixed_' r2{j} '_Logistic_' n1{i} '.txt']);
        t2   = [t2; t0];
        t2_1 = [t2_1; t0_1];
    end

    n0 = {'Heel bone mineral density', 'Cholesterol','LDL-Cholesterol','HDL-Cholesterol','Triglycerides',...
        'Apolipoprotein-A', 'Apolipoprotein-B','C-Reactive protein','Lipoprotein (a)','Glucose','HbA1c'};
    f.CurrentAxes = AX(1); hold on; grid on; box on;
    errorbar(t1.Beta, [11,9,8,2,1]+0.3, t1.Beta-t1.L95CI, t1.U95CI-t1.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
    errorbar(t1_1.Beta, [11,9,8,2,1]+0.1, t1_1.Beta-t1_1.L95CI, t1_1.U95CI-t1_1.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[90 139 162]/255,"Color",[90 139 162]/255)
    errorbar(t2.Beta, flip([1:10])-0.1, t2.Beta-t2.L95CI, t2.U95CI-t2.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
    errorbar(t2_1.Beta, flip([1:10])-0.3, t2_1.Beta-t2_1.L95CI, t2_1.U95CI-t2_1.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[118 31 42]/255,"Color",[118 31 42]/255)

    fill([-0.6 1.2 1.2 -0.6],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

    ax = gca;
    ax.YTick = 1:11;
    ax.YTickLabel = flip(n0);
    ax.TickLabelInterpreter = 'latex';
    ax.XLim = [-.4 1.2]; ax.XTick = [-0.4:0.2:1.2];
    ax.YLim = [0.5 11.5];
    ax.FontSize = 11;
    xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
    l1 = legend('GWAS (r2$\leq$0.3)', ['GWAS (r2$\leq$' r2{j} ')'], 'UK Biobank (r2$\leq$0.3)', ['UK Biobank (r2$\leq$' r2{j} ')'],'FontSize',12);
    l1.Interpreter = 'latex'; l1.NumColumns = 2; l1.Location = "northoutside";
    l1.Position    = [0.1886    0.9234    0.2078    0.0466];
    text(-0.5, 12.05, 'A)','Interpreter','latex','FontSize',17)
    xlabel("Unit increase per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)

    % Figure B ------------------------------------------------------------
    % GWAS results
    n1 = {'hf','cad','mi','is','hypertension', 't2dm'};

    t1   = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
    t1_1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_' r2{j} '_250000.txt']);
    for i = 2:length(n1)
        t0   = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
        t0_1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_' r2{j} '_250000.txt']);
        t1   = [t1; t0];
        t1_1 = [t1_1; t0_1];
    end

    % UK Biobank results
    n1 = {'CAD','MI','IS','Hypertension','T2DM'};

    t2   = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logaritmic_' n1{1} '.txt']);
    t2_1 = readtable([path 'Study2\MendelianRandomisation\Fixed_' r2{j} '_Logaritmic_' n1{1} '.txt']);

    for i = 2:length(n1)
        t0   = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logaritmic_' n1{i} '.txt']);
        t0_1 = readtable([path 'Study2\MendelianRandomisation\Fixed_' r2{j} '_Logaritmic_' n1{i} '.txt']);
        t2   = [t2; t0];
        t2_1 = [t2_1; t0_1];
    end

    % UK Biobank results - survival
    n1 = {'CAD','MI','IS','Hypertension','T2DM'};

    t3   = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{1} '.txt']);
    t3_1 = readtable([path 'Study2\MendelianRandomisation\Fixed_' r2{j} '_Birth_' n1{1} '.txt']);

    for i = 2:length(n1)
        t0   = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{i} '.txt']);
        t0_1 = readtable([path 'Study2\MendelianRandomisation\Fixed_' r2{j} '_Birth_' n1{i} '.txt']);

        t3   = [t3; t0];
        t3_1 = [t3_1; t0_1];
    end


    n0 = {'Hip fracture', 'Coronary artery disease','Myocardial infarction','Ischaemic stroke',...
        'Hypertension','Type 2 diabetes'};
    f.CurrentAxes = AX(2); hold on; grid on; box on; hold on; grid on; box on;

    errorbar(t1.OR, flip([1:6])+0.25, t1.OR-t1.L95CI, t1.U95CI-t1.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
    errorbar(t1_1.OR, flip([1:6])+0.15, t1_1.OR-t1_1.L95CI, t1_1.U95CI-t1_1.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[90 139 162]/255,"Color",[90 139 162]/255)
    errorbar(t2.OR, flip([1:5])+0.05, t2.OR-t2.L95CI, t2.U95CI-t2.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
    errorbar(t2_1.OR, flip([1:5])-0.05, t2_1.OR-t2_1.L95CI, t2_1.U95CI-t2_1.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[118 31 42]/255,"Color",[118 31 42]/255)
    errorbar(t3.OR, flip([1:5])-0.15, t3.OR-t3.L95CI, t3.U95CI-t3.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)
    errorbar(t3_1.OR, flip([1:5])-0.25, t3_1.OR-t3_1.L95CI, t3_1.U95CI-t3_1.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[128 76 43]/255,"Color",[128 76 43]/255)

    fill([0.075 4 4 0.075],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

    ax = gca;
    ax.YTick = 1:6;
    ax.YTickLabel = flip(n0);
    ax.TickLabelInterpreter = 'latex';
    ax.YLim = [0.5 6.5];
    ax.XScale = 'log';
    ax.XLim = [0.075 4];
    ax.XTick = [0.1 0.25 0.5 1 2 4];
    xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
    ax.FontSize = 11;
    l = legend('GWAS [Odds, r2$\leq$0.3]', ['GWAS [Odds, r2$\leq$' r2{j} ']'], ...
        'UKB-LR [Odds, r2$\leq$0.3]',['UKB-LR [Odds, r2$\leq$' r2{j} ']'], ...
        'UKB-SA(Birth) [Hazard, r2$\leq$0.3]',['UKB-SA(Birth) [Hazard, r2$\leq$' r2{j} ']'],'FontSize',10.5);
    l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside"; l.Position = [0.6052    0.9175    0.3847    0.0514];
    xlabel("Ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
    text(0.05, 6.8, 'B)','Interpreter','latex','FontSize',17)
    print(f,[path 'Figures\SFig7_' r2{j} '.png'],"-dpng","-r800")
end



%% Random-effects
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';

f = figure(1); hold on; box on; axis off;
AX = get_axes([1,2],0.1,0.12,0.01,0.07,0.1,0.1);
f.Position = [9 60 1767 920];

% GWAS results
n1 = {'bmd','ldl','hdl','glucose', 'hba1c'};

t1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
t11 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Random_0.3_250000.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
    t00 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Random_0.3_250000.txt']);
    t1  = [t1; t0];
    t11  = [t11; t00];
end


% UK Biobank results
n1 = {'cholesterol','ldl','hdl','triglycerides', 'apoA', 'apoB', 'crp', 'lipoprotein',...
    'glucose', 'hba1c'};
t2 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logistic_' n1{1} '.txt']);
t22 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Logistic_' n1{1} '.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logistic_' n1{i} '.txt']);
    t00 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Logistic_' n1{i} '.txt']);

    t2  = [t2; t0];
    t22  = [t22; t00];
end

n0 = {'Heel bone mineral density', 'Cholesterol','LDL-Cholesterol','HDL-Cholesterol','Triglycerides',...
    'Apolipoprotein-A', 'Apolipoprotein-B','C-Reactive protein','Lipoprotein (a)','Glucose','HbA1c'};

f.CurrentAxes = AX(1); hold on; grid on; box on;

errorbar(t1.Beta, [11,9,8,2,1]+0.25, t1.Beta-t1.L95CI, t1.U95CI-t1.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(t11.Beta, [11,9,8,2,1]+0.1, t11.Beta-t11.L95CI, t11.U95CI-t11.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[90 139 162]/255,"Color",[90 139 162]/255)
errorbar(t2.Beta, flip([1:10])-0.1, t2.Beta-t2.L95CI, t2.U95CI-t2.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
errorbar(t22.Beta, flip([1:10])-0.25, t22.Beta-t22.L95CI, t22.U95CI-t22.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[118 31 42]/255,"Color",[118 31 42]/255)

fill([-0.6 1.3 1.3 -0.6],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:11;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-.5 1.3]; ax.XTick = [-0.4:0.2:1.2];
ax.YLim = [0.5 11.5];
ax.FontSize = 11;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
l1 = legend('GWAS (Fixed)', 'GWAS (Random)','UK Biobank (Fixed)', 'UK Biobank (Random)','FontSize',12);
l1.Interpreter = 'latex'; l1.NumColumns = 2; l1.Location = "northoutside"; l1.Position = [0.1886    0.9234    0.2078    0.0466];
text(-0.6, 12.1, 'A)','Interpreter','latex','FontSize',17)
xlabel("Unit increase per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)

% Figure B ------------------------------------------------------------
% GWAS results
n1 = {'hf','cad','mi','is','hypertension', 't2dm'};

t1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
t11 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Random_0.3_250000.txt']);

for i = 2:length(n1)
    t0 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
    t00 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Random_0.3_250000.txt']);

    t1  = [t1; t0];
    t11  = [t11; t00];
end

% UK Biobank results
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t2 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logaritmic_' n1{1} '.txt']);
t22 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Logaritmic_' n1{1} '.txt']);

for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Logaritmic_' n1{i} '.txt']);
    t00 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Logaritmic_' n1{i} '.txt']);

    t2  = [t2; t0];
    t22  = [t22; t00];
end

% UK Biobank results - survival
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t3 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{1} '.txt']);
t33 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Birth_' n1{1} '.txt']);

for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{i} '.txt']);
    t00 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Birth_' n1{i} '.txt']);
    t3 = [t3; t0];
    t33 = [t33; t00];
end


n0 = {'Hip fracture', 'Coronary artery disease','Myocardial infarction','Ischaemic stroke',...
    'Hypertension','Type 2 diabetes'};

f.CurrentAxes = AX(2); hold on; grid on; box on; hold on; grid on; box on;

errorbar(t1.OR, flip([1:6])+0.25, t1.OR-t1.L95CI, t1.U95CI-t1.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(t11.OR, flip([1:6])+0.15, t11.OR-t11.L95CI, t11.U95CI-t11.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[90 139 162]/255,"Color",[90 139 162]/255)
errorbar(t2.OR, flip([1:5])+0.05, t2.OR-t2.L95CI, t2.U95CI-t2.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[167 21 49]/255,"Color",[167 21 49]/255)
errorbar(t22.OR, flip([1:5])-0.05, t22.OR-t22.L95CI, t22.U95CI-t22.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[118 31 42]/255,"Color",[118 31 42]/255)
errorbar(t3.OR, flip([1:5])-0.15, t3.OR-t3.L95CI, t3.U95CI-t3.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)
errorbar(t33.OR, flip([1:5])-0.25, t33.OR-t33.L95CI, t33.U95CI-t33.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[128 76 43]/255,"Color",[128 76 43]/255)

fill([0.05 15 15 0.05],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.05 15];
ax.XTick = [0.1 0.25 0.5 1 2 4 8];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('GWAS [Odds, Fixed]','GWAS [Odds, Random]', ...
    'UKB-LR [Odds, Fixed]','UKB-LR [Odds, Random]', ...
    'UKB-SA(Birth) [Hazard, Fixed]','UKB-SA(Birth) [Hazard, Random]','FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside"; l.Position = [0.6052    0.9175    0.3847    0.0514];
xlabel("Ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.03, 6.8, 'B)','Interpreter','latex','FontSize',17)

print(f,[path '\Figures\Fig9.png'],"-dpng","-r800")



%% Principal components analysis

clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';

f = figure(1); hold on; box on; axis off;
AX = get_axes([1,2],0.1,0.12,0.01,0.07,0.1,0.1);
f.Position = [9 60 1767 920];

% GWAS results
n1 = {'bmd','ldl','hdl','glucose', 'hba1c'};

t1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
t11 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_999_Fixed_PCA.txt.txt']);

for i = 2:length(n1)
    t0 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
    t00 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_999_Fixed_PCA.txt.txt']);
    t1  = [t1; t0];
    t11  = [t11; t00];
end

n0 = {'Heel bone mineral density', 'LDL-Cholesterol','HDL-Cholesterol','Glucose','HbA1c'};

f.CurrentAxes = AX(1); hold on; grid on; box on;

errorbar(t1.Beta, flip([1:5])+0.1, t1.Beta-t1.L95CI, t1.U95CI-t1.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(t11.Beta, flip([1:5])-0.1, t11.Beta-t11.L95CI, t11.U95CI-t11.Beta,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[90 139 162]/255,"Color",[90 139 162]/255)

fill([-0.6 1.3 1.3 -0.6],[5.5 5.5 4.5 4.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:5;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.XLim = [-.5 1.3]; ax.XTick = [-0.4:0.2:1.2];
ax.YLim = [0.5 5.5];
ax.FontSize = 11;
xline(0,"LineWidth",1.5,"LineStyle","--","Color",'k')
l1 = legend('GWAS (Pruning)', 'GWAS (PCA)','FontSize',12); 
text(-0.65, 5.75, 'A)','Interpreter','latex','FontSize',17)
l1.Interpreter = 'latex'; l1.NumColumns = 2; l1.Location = "northoutside";l1.Position = [0.1886    0.9234    0.2078    0.035];
xlabel("Unit increase per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)

% Figure B ------------------------------------------------------------
% GWAS results
n1 = {'hf','cad','mi','is','hypertension', 't2dm'};

t1 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_Fixed_0.3_250000.txt']);
t11 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{1} '_999_Fixed_PCA.txt.txt']);

for i = 2:length(n1)
    t0 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_Fixed_0.3_250000.txt']);
    t00 = readtable([path 'Study1\MendelianRandomisation\Results\' n1{i} '_999_Fixed_PCA.txt.txt']);

    t1  = [t1; t0];
    t11  = [t11; t00];
end

n0 = {'Hip fracture', 'Coronary artery disease','Myocardial infarction','Ischaemic stroke',...
    'Hypertension','Type 2 diabetes'};

f.CurrentAxes = AX(2); hold on; grid on; box on; hold on; grid on; box on;
errorbar(t1.OR, flip([1:6])+0.1, t1.OR-t1.L95CI, t1.U95CI-t1.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[86 166 194]/255,"Color",[86 166 194]/255)
errorbar(t11.OR, flip([1:6])-0.1, t11.OR-t11.L95CI, t11.U95CI-t11.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[90 139 162]/255,"Color",[90 139 162]/255)
fill([0.05 15 15 0.05],[5.5 5.5 6.5 6.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)

ax = gca;
ax.YTick = 1:6;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 6.5];
ax.XScale = 'log';
ax.XLim = [0.05 15];
ax.XTick = [0.1 0.25 0.5 1 2 4 8];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('GWAS (Pruning)','GWAS (PCA)','FontSize',12);
l.Interpreter = 'latex'; l.NumColumns = 2; l.Location = "northoutside"; l.Position = [0.6875    0.9234   0.2078    0.035];
xlabel("Odds ratio per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.03, 6.8, 'B)','Interpreter','latex','FontSize',17)

print(f,[path '\Figures\Fig10.png'],"-dpng","-r800")

%% Survival analysis since enrolment - categorical and survival outcomes
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';

% UK Biobank results - survival
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t3 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{1} '.txt']);
t33 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Enrol_' n1{1} '.txt']);

for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{i} '.txt']);
    t00 = readtable([path 'Study2\MendelianRandomisation\Random_0.3_Enrol_' n1{i} '.txt']);
    t3 = [t3; t0];
    t33 = [t33; t00];
end


n0 = {'Coronary artery disease','Myocardial infarction','Ischaemic stroke',...
    'Hypertension','Type 2 diabetes'};
f = figure(2); hold on; grid on; box on;
f.Position = [680   254   799   644];
errorbar(t3.OR, flip([1:5])-0.15, t3.OR-t3.L95CI, t3.U95CI-t3.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)
errorbar(t33.OR, flip([1:5])-0.25, t33.OR-t33.L95CI, t33.U95CI-t33.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[236 200 69]/255,"Color",[236 200 69]/255)
ax = gca;
ax.YTick = 1:5;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 5.5];
ax.XScale = 'log';
ax.XLim = [0.05 15];
ax.XTick = [0.1 0.25 0.5 1 2 4 8];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('SA since birth', "SA since enrolment",'FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside"; l.Position = [0.1050 0.9370 0.8804 0.0373];
xlabel("HR per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.015, 6.8, 'B)','Interpreter','latex','FontSize',17)

print(f,[path '\Figures\Fig11.png'],"-dpng","-r800")

%% Survival analysis new package - categorical and survival outcomes
clear all; close all;
path = 'D:\Projects\MR_Sclerostin\Results\';

% UK Biobank results - survival
n1 = {'CAD','MI','IS','Hypertension','T2DM'};

t3  = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{1} '.txt']);
for i = 2:length(n1)
    t0 = readtable([path 'Study2\MendelianRandomisation\Fixed_0.3_Birth_' n1{i} '.txt']);
    t3 = [t3; t0];
end

t33 = readtable([path 'SensitivityAnalyses\SurvivalPackage_results.txt']);
t33 = t33(2:6,:);

n0 = {'Coronary artery disease','Myocardial infarction','Ischaemic stroke',...
    'Hypertension','Type 2 diabetes'};
f = figure(2); hold on; grid on; box on;
f.Position = [680   254   799   644];
errorbar(t3.OR, flip([1:5])-0.15, t3.OR-t3.L95CI, t3.U95CI-t3.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[206 109 50]/255,"Color",[206 109 50]/255)
errorbar(t33.OR, flip([1:5])-0.25, t33.OR-t33.L95CI, t33.U95CI-t33.OR,"o","horizontal","LineWidth",1.25,"MarkerFaceColor",[236 200 69]/255,"Color",[236 200 69]/255)
ax = gca;
ax.YTick = 1:5;
ax.YTickLabel = flip(n0);
ax.TickLabelInterpreter = 'latex';
ax.YLim = [0.5 5.5];
ax.XScale = 'log';
ax.XLim = [0.05 15];
ax.XTick = [0.1 0.25 0.5 1 2 4 8];
xline(1,"LineWidth",1.5,"LineStyle","--","Color",'k')
ax.FontSize = 11;
l = legend('SA', "SA with gwasurvivr",'FontSize',10.5);
l.Interpreter = 'latex'; l.NumColumns = 3; l.Location = "northoutside"; l.Position = [0.1050 0.9370 0.8804 0.0373];
xlabel("HR per 1-SD decrease in sclerostin levels","Interpreter","latex",'FontSize',12)
text(0.015, 6.8, 'B)','Interpreter','latex','FontSize',17)

print(f,[path '\Figures\Fig12.png'],"-dpng","-r800")


%% Functions
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


