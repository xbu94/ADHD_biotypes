%% PLS analysis for gene_brain
data_dir = 'D:\Helab\ADHDsubtyping\AHBA\pls_gene\CTsubtypes\SAspin\weights_disorders\';
cortical_t = load('D:\Helab\ADHDsubtyping\AHBA\surrogatemap\t_CTsubtype1vs0_unth.txt');
gene = csvread('D:\Helab\ADHDsubtyping\AHBA\gene_expression_on_brain\left_DS01_expression.csv',1,0);
X = zscore(gene);
Y = zscore(cortical_t);
dim = 5;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim,'CV',dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);

%align PLS components with desired direction%
R1 = corr([XS(:,1),XS(:,2),XS(:,3)],Y(:,1));
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end
%% spatial autocorrelation corrected permutation test to assess the significance of PLS component variance explained ratios
surrogate_tmaps_CTsubtype1vs0_5000 = csvread('D:/Helab/ADHDsubtyping/AHBA/nullspin_tmap5000_subtype1_unth.csv',1,0);
PCTVARrand = zeros(5000,dim);
Rsq = zeros(5000,1);
for j=1:5000
    disp(j);
    %Yp=zscore(surrogate_tmaps_CTsubtype1vs0_5000(j,:)');
    Yp=zscore(surrogate_tmaps_CTsubtype1vs0_5000(:,j));
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(zscore(gene),Yp,dim);
    PCTVARrand(j,:)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim); 
end 
p_single = zeros(1,dim);
for k=1:dim
    p_single(k)=length(find(PCTVARrand(:,k)>=PCTVAR(2,k)))/5000;
end
p_cum = length(find(Rsq>=Rsquared))/5000;
myStats=[PCTVAR; p_single];
csvwrite(['D:\Helab\ADHDsubtyping\AHBA\pls_gene\CTsubtypes\SAspin\PLS_stats_CTsubtypes1_unth.csv'],myStats);
save([data_dir,'PLS_results_ID.mat']);

%% Draw variance explanation
%load('D:\Helab\ADHDsubtyping\AHBA\pls_gene\PLS_stats_sa_t.csv');
%py = plot(sort(PLS_stats_sa_t(2,:)','descend'),'.-');
py = plot(sort(PCTVAR(2,:)','descend'),'.-');
py.Color = [115 130 184]/255;
py.MarkerSize = 7;
%hold on
%plot(0:11,0.1*ones(1,12),'--','Color',[226,115,134]./255,'LineWidth',0.5);
%hold off

xlabel('PLS Component');
ylabel('Explained variance');
set(gca,'XLim',[0,6]);
set(gca,'YLim',[0,0.04],'YTick',[0:0.01:0.04]);
%t1 = text(1,0.31,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',7);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',7);
box off

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 9 5.6]);
print(gcf,[data_dir,'PLS\PLS_Variance_MDD.tif'],'-dtiff','-r1000')

%% save score,weights, and calculate correlation
csvwrite([data_dir,'SAspin\ADHDgenes\PLS_scores_CTsubtype1vs0.csv'],XS(:,1));
csvwrite([data_dir,'SAspin\ADHDgenes\PLS_weights_CTsubtype1vs0.csv'],stats.W(:,1));

[r_corr,p_corr] = corr(XS(:,1),Y);
corr_real = r_corr;
% permutation 5000times for correlation
corr_surr = zeros(1,5000);
for j = 1:5000
    disp(j);
    corr_surr(1,j) = corr(XS(:,1),surrogate_tmaps_CTsubtype1vs0_5000(:,j));
end
p = zeros(1,1);
p = length(find(corr_surr(1,:)>r_corr))/5000;
save([data_dir,'PLS\PLS_corr_ADHD1.mat'],'corr_surr','corr_real','p_corr');

%[PLS1w,x1] = sort(stats.W(:,1),'descend');
%gene_names = importdata([data_dir,'gene_namesAHBA.txt']);
%PLS1ids=gene_names(x1);
%genes_pos_W = PLS1ids(1:33,1);
%genes_neg_W = PLS1ids(34:end,1);
%writecell(genes_pos_W,['D:\Helab\ADHDsubtyping\AHBA\pls_gene\CTsubtypes\genes_pos_W.csv']);
%writecell(genes_neg_W,['D:\Helab\ADHDsubtyping\AHBA\pls_gene\CTsubtypes\genes_neg_W.csv']);

%% Draw figures for correlations
close all
dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = {'{\itt}-statistic of CTsubtype1';'group difference'};
xlable1 = 'PLS scores';
[xData, yData] = prepareCurveData(XS(:,1),Y);
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );
h=plot( fitresult, xData, yData);
set(h(1),'Marker','.','MarkerSize',6,'Color',dotcolor)
set(h(2),'LineWidth',0.5,'Color',linecolor)
hold on
xFit = linspace(min(xData),max(xData),100);
yPredict = predint(fitresult,xFit,0.95,'functional','off');
fy = cat(2,yPredict(:,2)',flip(yPredict(:,1),1)')';
fx = cat(2,xFit,flip(xFit',1)')';
fill(fx,fy,[0.5 0.5 0.5],'EdgeAlpha',0,'FaceAlpha',0.3);
hold off
legend off
ylabel(ylable1);
xlabel(xlable1);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'XTick',[-0.22:0.1:0.26]);
set(gca,'YTick',[-2.7:0.9:2.7]);
t1 = text(-0.18,2.6,{'{\itr} = 0.55';'{\itP} < 0.001'},'FontName','Arial','FontSize',4);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.95 3.8]);
grid off
box off
print(gcf,[data_dir,'PLS\corr_PLS_SCZ.tif'],'-dtiff','-r1000')

%% calculate corrected weight
gene_name = importdata([data_dir,'gene_namesADHDYL.txt']);
geneindex=1:size(gene,2);
bootnum=5000;
%X = zscore(gene);
%Y = zscore(cortical_t);
dim=1;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

[R1,p1]=corr(XS(:,1),cortical_t);
if R1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=gene_name(x1);
geneindex1=geneindex(x1);
%PLS1_ROIscores=XS(:,1);
%save([data_dir,'GeneAssociation_Main\PLS1_ROIscore.mat'],'PLS1_ROIscores');
%csvwrite([data_dir,'GeneAssociation_Main\PLS1_ROIscores.csv'],XS(:,1));
%PLS1_score=XS(:,1);

%[R2,p2]=corr(XS(:,2),MRIdata);
%if R2(1,1)<0
%    stats.W(:,2)=-1*stats.W(:,2);
%    XS(:,2)=-1*XS(:,2);
%end
%[PLS2w,x2] = sort(stats.W(:,2),'descend');
%PLS2ids=genes(x2);
%geneindex2=geneindex(x2);
%PLS2_ROIscores_280=XS(:,2);
%save([data_dir,'GeneAssociation_Main\PLS2_ROIscore.mat'],'PLS2_ROIscores_280');
%csvwrite([data_dir,'GeneAssociation_Main\PLS2_ROIscores.csv'],XS(:,2));
%PLS2_score=XS(:,2);

PLS1weights = zeros(length(gene_name),bootnum);
%PLS2weights = zeros(10027,10000);

parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights(:,i) = newW;%store (ordered) weights from this bootstrap run
    
    %temp=stats.W(:,2);%extract PLS2 weights
    %newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    %if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
    %    newW=-1*newW;
    %end
    %PLS2weights(:,i) = newW; %store (ordered) weights from this bootstrap run    
end

PLS1sw = std(PLS1weights');
temp1=PLS1w./PLS1sw';
[Z1,ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
fid1 = fopen([data_dir,'SAspin\ADHDgenes\PLS_geneWeights_corrected_subtype1.csv'],'w');
for i=1:length(gene_name)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
end
fclose(fid1);
