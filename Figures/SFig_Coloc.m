clear all; close all;
path = pwd; path = path(1:end-length('\Figures'));

outcomes = {'HF','CAD','MI','IS','Hypertension','T2DM','eBMD'};
outcomes_names = {'Hip fracture','Coronary artery disease','Myocardial infarction','Ischaemic stroke','Hypertension','Type II diabetes mellitus','eBMD'};

alt = [7,7,7,7,7,7,90];
let = {'A)','B)','C)','D)','E)','F)','G)'};
order = [3,6,9,2,5,8,1];
f = figure(1); hold on; box on; axis off;
f.Position = [3 42 1913 1074];

AX = get_axes([3,3],0.03,0.02,0.02,0.07,0.05,0.03);

for i = 1:length(let)
    t = readtable([path '\Colocalization\data.xlsx'],'Sheet',outcomes{i});
    pos = t.Pos;
    pval.outcome  = t.pval_outcome;
    pval.exposure = t.pval_exposure;

    f.CurrentAxes = AX(order(i));
    scatter(pos,-log10(pval.exposure),20,'k','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
    scatter(pos,-log10(pval.outcome),20,'r','MarkerFaceColor','r')
    ax = gca;
    ax.XTick = [41831099-20*1000, 41831099, 41836156, 41836156+20*1000];
    ax.XTickLabel = '';
    ax.TickLabelInterpreter = 'latex';
    ax.FontSize = 10.5;
    ax.XLim = [41831099-21*1000 41836156+21*1000];
  
    fill([41831099, 41836156, 41836156, 41831099],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[237 177 32]./255,'EdgeColor',[237 177 32]./255,'FaceAlpha',0.2,'EdgeAlpha',0.2)
    fill([41831099-20*1000, 41831099, 41831099, 41831099-20*1000],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[237 177 32]./255,'EdgeColor',[237 177 32]./255,'FaceAlpha',0.05,'EdgeAlpha',0)
    fill([41836156, 41836156+20*1000, 41836156+20*1000, 41836156],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],[237 177 32]./255,'EdgeColor',[237 177 32]./255,'FaceAlpha',0.05,'EdgeAlpha',0)

    if i == 7 || i == 5 || i == 6
        xlabel('Chromosome 17 position [bp]','Interpreter','latex','FontSize',14)
        ax.XTick = [41831099-20*1000, 41831099, 41836156, 41836156+20*1000];
        ax.XAxis.Exponent = 0;
        xtickformat('%.0f');
    end

    if i == 1 || i == 4 ||i == 7
        ylabel('-$log_{10}$(P-Value)','Interpreter','latex','FontSize',14)
    end
    if i == 7
        ax.YLim(2) = 90;
    end
    text(41831099-20*1000,alt(i),let{i},'Interpreter','latex','FontSize',17,'VerticalAlignment','bottom')
    l = legend('Sclerostin',outcomes_names{i},'Interpreter','latex','FontSize',11,'FontWeight','bold','Location','northeast');
    
end

f.CurrentAxes = AX(7);  box off; axis off;
f.CurrentAxes = AX(4);  box off; axis off;

print(f,'SFig_coloc.png','-dpng','-r700')


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
