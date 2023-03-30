clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

t = readtable([path '\Pruning\gMR_res.xlsx'],'Sheet','data_eBMD');
t.beta_outcome(t.beta_exposure < 0)  = -t.beta_outcome(t.beta_exposure < 0);
t.beta_exposure(t.beta_exposure < 0) = -t.beta_exposure(t.beta_exposure < 0);

betax = t.beta_exposure; sex = t.se_exposure;
betay = t.beta_outcome;  sey = t.se_outcome;

f = figure(1); grid on; hold on; box on;
errorbar(betax,betay,betay-(betay-1.96*sey),(betay+1.96*sey)-betay,betax-(betax-1.96*sex),(betax+1.96*sex)-betax,"o","LineWidth",0.5,"MarkerFaceColor",'k',"Color",'k');
ax = gca;
ax.XLim(1) = 0;
ax.TickLabelInterpreter = 'latex';
text(-0.01,0.0025,'A)','FontSize',17,'Interpreter','latex')
t = readtable([path '\Pruning\gMR_res.xlsx'],'Sheet','results_eBMD');
plot(linspace(ax.XLim(1),ax.XLim(2),10),t.Estimate(1)*linspace(ax.XLim(1),ax.XLim(2),10),'--',"Color",'k')

xlabel("MR effect size per 1-SD change in sclerostin levels","Interpreter","latex")
ylabel("MR effect size per 1-SD change in eBMD","Interpreter","latex")
print(f,"Fig6.png","-dpng","-r800")
%%
clear all; close all
path = pwd; path = path(1:end-length('\Figures'));

t = readtable([path '\Pruning\gMR_res.xlsx'],'Sheet','data_HF');
t.beta_outcome(t.beta_exposure < 0)  = -t.beta_outcome(t.beta_exposure < 0);
t.beta_exposure(t.beta_exposure < 0) = -t.beta_exposure(t.beta_exposure < 0);

betax = t.beta_exposure; sex = t.se_exposure;
betay = t.beta_outcome;  sey = t.se_outcome;

f = figure(1); grid on; hold on; box on;
errorbar(betax,betay,betay-(betay-1.96*sey),(betay+1.96*sey)-betay,betax-(betax-1.96*sex),(betax+1.96*sex)-betax,"o","LineWidth",0.5,"MarkerFaceColor",'k',"Color",'k');
ax = gca;
ax.XLim(1) = 0;
ax.TickLabelInterpreter = 'latex';
text(-0.0125,0.25,'B)','FontSize',17,'Interpreter','latex')
t = readtable([path '\Pruning\gMR_res.xlsx'],'Sheet','results_HF');
plot(linspace(ax.XLim(1),ax.XLim(2),10),t.Estimate(1)*linspace(ax.XLim(1),ax.XLim(2),10),'--',"Color",'k')
xlabel("MR effect size per 1-SD change in sclerostin levels","Interpreter","latex")
ylabel("MR effect size for hip fracture ","Interpreter","latex")
print(f,"Fig6.1.png","-dpng","-r800")