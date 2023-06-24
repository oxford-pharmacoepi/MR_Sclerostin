clear all; close all;

col = [192 0 0]./255;
f = figure(1); f.Position = [351.5 91 768.5 621.5]; hold on; axis off;

N = [1,2]; lx = 0.07; mx = 0.05; rx = 0.07; dy = 0.1; my = 0; uy = 0.4; 
Fx = (1-mx-rx-lx)/2; AX = get_axes(N,lx,mx,rx,dy,my,uy);

annotation('textbox',[lx 1-0.07 Fx 0.06],'Color','k','String','Linear regression','Interpreter','latex',...
    'FontSize',17,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[225 225 225]./255,'EdgeColor',[225 225 225]./255)
annotation('textbox',[Fx+lx+mx 1-0.07 Fx 0.06],'Color','k','String','Logistic regression','Interpreter','latex',...
    'FontSize',17,'VerticalAlignment','middle','HorizontalAlignment','center','BackgroundColor',[225 225 225]./255,'EdgeColor',[225 225 225]/255)

f.CurrentAxes = AX(1); grid on; box on; hold on; % ------------------------
rng(2); n = 20;
x1 = 0.15*rand(n,1); x2 = 0.75+(1-0.75).*rand(n,1);
t = [x1;x2] + rand(length([x1;x2]),1);
scatter([zeros(n,1); ones(n,1)],t,'filled','o','MarkerFaceColor',col,'MarkerFaceAlpha',0.2)
ax = gca; ax.XLim(1) = -0.05; ax.XLim(2) = 1.05; ax.YLim(1) = -0.05; ax.YLim(2) = 2.05;
plot([0 1], [mean(t([1:n])) mean(t([n+1:2*n]))],'Color',col,'LineWidth',1.5)
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 12;
ax.XTick = [0:0.2:1]; ax.XTickLabel = {'0','','','','','1'}; 
ax.YTick = [-0.06:0.36:2.34]; ax.YTickLabel = {'0','','','','','','1'}; 
ylim([-0.26, 2.26])
xlabel('SNP','FontSize',14,'Interpreter','latex')
ylabel('Outcome','FontSize',12,'Interpreter','latex')

text(-0.16,4.15,'A)','FontSize',15,'Interpreter','latex')
text(1.1,4.15,'B)','FontSize',15,'Interpreter','latex')

text(0.5,3.56,'Outcome = $\beta_0$ + $\beta_1$ $\cdot$ [SNP]','HorizontalAlignment','center', ...
    'FontSize',16,'Interpreter','latex')
text(-0.05,3.1,'Where:','FontSize',13,'Interpreter','latex','HorizontalAlignment','left')
text(0,2.90,'- SNP = 0,1','FontSize',13,'Interpreter','latex','HorizontalAlignment','left')
text(0,2.70,'- Outcome, $\beta_0$, $\beta_1$ $\in R$','FontSize',13,'Interpreter','latex','HorizontalAlignment','left')

f.CurrentAxes = AX(2); box on; hold on; grid on; % ------------------------
x1 = 0.1*rand(25,1); x2 = 0.9+(1-0.9).*rand(25,1);

x11  = [0.05-x1(randi(length(x1),length(x1),1));0.05+x2(randi(length(x2), length(x2),1))];
fact = round(x11); fact(fact == 0) = -0.05; fact(fact == 1) = 0.05; 
scatter([x1+0.05;x2-0.05],x11+fact,'filled','o','MarkerFaceColor',col,'MarkerFaceAlpha',0.2);
x = [x1+0.05; x2-0.05];
y = [x11+fact];

x11  = [x2(randi(length(x2),5,1))+0.05;0.05-x1(randi(length(x1),5,1))];
fact = round(x11); fact(fact == 0) = -0.05; fact(fact == 1) = 0.05; 
scatter([x1(1:5)+0.05;x2(1:5)-0.05],x11+fact,'filled','o','MarkerFaceColor',col,'MarkerFaceAlpha',0.2); 

x = [x; x1(1:5)+0.05;x2(1:5)-0.05];
y = [y; x11+fact];

[logitCoef,dev] = glmfit(x,round(y),'binomial','logit');
logitFit = glmval(logitCoef,linspace(-0.5,1.5,length(x)),'logit');
plot(linspace(-0.5,1.5,length(x)),logitFit,'Color',col,'LineWidth',1.5); hold on;

ax = gca; ax.XLim(1) = -0.1; ax.XLim(2) = 1.1; ax.YLim(1) = -0.1; ax.YLim(2) = 1.1;
ax.XTick = [-0.1:0.2:1]; ax.XTickLabel = {'','0','','','','1'}; 
ax.YTick = [-0.03:0.18:1.17]; ax.YTickLabel = {'0','','','','','','1'}; 
ylim([-0.13 1.13]);
ax.TickLabelInterpreter = 'latex'; ax.FontSize = 12;
xlabel('SNP','FontSize',14,'Interpreter','latex');
ylabel('Prob(outcome)','FontSize',12,'Interpreter','latex')

text(0.2,1.77,'Prob(Outcome) = ', ...
    'HorizontalAlignment','center','Interpreter','latex','FontSize',16)
%$\frac{1}{1+e^{-(\beta_0+\beta_1 \cdot [SNP])}}$
text(0.85,1.77,'$\frac{1}{1+e^{-(\beta_0+\beta_1 \cdot [SNP])}}$', ...
    'HorizontalAlignment','center','Interpreter','latex','FontSize',21)

text(0,1.55,'Where:','HorizontalAlignment','left','Interpreter','latex','FontSize',13)
text(0.05,1.45,'- SNP, Outcome = 0,1','Interpreter','latex','FontSize',13)
text(0.05,1.35,'- $\beta_0, \beta_1 \in R$','Interpreter','latex','FontSize',13)
text(0.05,1.25,'- Prob(Outcome) $\in$ [0,1]','Interpreter','latex','FontSize',13)

print(f,['LinearLogisticRegression1.png'],'-dpng','-r800')


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
