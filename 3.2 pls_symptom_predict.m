%% PLS analysis for brain_symptoms
data_dir = 'D:\Helab\ADHDsubtyping\Stats\pls\CTsubtypes\';
brain = csvread([data_dir,'symptoms_all\z_ct1_clinic.csv'],1,1);
clinic1 = readtable([data_dir,'symptoms_all\CTsubtype1_clinic.csv']);
clinic2 = readtable([data_dir,'symptoms_all\CTsubtype2_clinic.csv']);

symptom1 = clinic1(:,["inattention","hyperactive_impulsive","ADHD_index"]);
symptom1 = table2array(symptom1);
symptom2 = clinic2(:,["inattention","hyperactive_impulsive","ADHD_index"]);
symptom2 = table2array(symptom2);

ADHDindex1 = clinic1(:,"ADHD_index");
ADHDindex1 = table2array(ADHDindex1);
ADHDindex2 = clinic2(:,"ADHD_index");
ADHDindex2 = table2array(ADHDindex2);

IA1 = clinic1(:,"inattention");
IA1 = table2array(IA1);
IA2 = clinic2(:,"inattention");
IA2 = table2array(IA2);

HI1 = clinic1(:,"hyperactive_impulsive");
HI1 = table2array(HI1);
HI2 = clinic2(:,"hyperactive_impulsive");
HI2 = table2array(HI2);

X = zscore(brain);
Y = zscore(symptom1);
dim = 10;
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
%% permutation test to assess the significance of PLS component variance explained ratios
nperm = 5000;
perm_index = zeros(nperm,size(symptom2,1));
for j=1:nperm
    perm_index(j,:) = randperm(size(symptom2,1));
end
PCTVARrand = zeros(nperm,dim);
Rsq = zeros(nperm,1);
for j=1:nperm
    disp(j);
    Yp=zscore(symptom2(perm_index(j,:),:));
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(zscore(brain),Yp,dim);
    PCTVARrand(j,:)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim); 
end 
clear Yp XLr YLr XSr YSr BETAr PCTVARr MSEr statsr
p_single = zeros(1,dim);
for k=1:dim
    p_single(k)=length(find(PCTVARrand(:,k)>=PCTVAR(2,k)))/nperm;
end
myStats=[PCTVAR; p_single];
csvwrite([data_dir,'symptoms_all\stats_subtype2.csv'],myStats);

%% Draw variance explanation
py = plot(PCTVAR(2,:)','.-');
py.Color = [115 130 184]/255;
py.MarkerSize = 7;
hold on
plot(0:11,0.1*ones(1,12),'--','Color',[226,115,134]./255,'LineWidth',0.5);
hold off

xlabel('PLS Component');
ylabel('Explained variance');
set(gca,'XLim',[0,11]);
set(gca,'YLim',[0,0.40],'YTick',0:0.10:0.40);
t1 = text(1,0.38,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',7);
%t2 = text(4,0.2,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',7);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',7);
box off

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 9 5.6]);
print(gcf,[data_dir,'symptoms_all\PLS_Variance1.tif'],'-dtiff','-r1000')

%% save score,weights, and calculate correlation
csvwrite([data_dir,'symptoms_all\brain_scores1.csv'],XS(:,1));
csvwrite([data_dir,'symptoms_all\brain_weights1.csv'],stats.W(:,1));
csvwrite([data_dir,'symptoms_all\brain_loadings1.csv'],XL(:,1));
csvwrite([data_dir,'symptoms_all\symptoms_scores1.csv'],YS(:,1));
csvwrite([data_dir,'symptoms_all\symptoms_loadings1.csv'],YL(:,1));

[r_corr,p_corr_self] = corr(XS(:,1),YS(:,1));
% permutation 5000times for correlation
corr_perm = zeros(1,nperm);
for j = 1:nperm
    disp(j);
    Yp=YS(perm_index(j,:),1);
    corr_perm(1,j) = corr(XS(:,1),Yp);
end
p_corr = length(find(corr_perm(1,:)>r_corr))/nperm;
save([data_dir,'symptoms_all\pls_corr_symptom1.mat'],'corr_perm','r_corr','p_corr');
clear Yp
%% Draw figures for correlations
close all
dotcolor = [115 130 184]/255;
linecolor = [0 0 0];
ylable1 = 'Symptoms PLS score';
xlable1 = 'Brain PLS score';
[xData, yData] = prepareCurveData(XS(:,1),YS(:,1));
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
set(gca,'XLim',[-0.3,0.3],'XTick',[-0.30:0.10:0.30]);
set(gca,'YLim',[-50,50],'YTick',[-50:10:50]);
t1 = text(-0.2,47,{'{\itr} = 0.68';'{\itP} < 0.001'},'FontName','Arial','FontSize',4);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.95 3.8]);
grid off
box off
print(gcf,[data_dir,'symptoms_all\pls_corr1.tif'],'-dtiff','-r1000')

%% bootstrapping (5000)
bootnum=5000;
X = zscore(brain);
Y = zscore(symptom1);
dim=1;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
[PLS1w,x1] = sort(XL(:,1),'descend');
[clinic1w,y1] = sort(YL(:,1),'descend');

PLS1loadings = zeros(219,5000);
Clinicloadings = zeros(3,5000);
parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=XL(:,1);%extract brain loadings
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1loadings(:,i) = newW;%store (ordered) weights from this bootstrap run 

    temp=YL(:,1);%extract clinic loadings
    newW=temp(y1); %order the newly obtained weights the same way as initial PLS 
    if corr(clinic1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    Clinicloadings(:,i) = newW;%store (ordered) weights from this bootstrap run
end

writematrix(PLS1loadings,[data_dir,'symptoms_all\brain_boot.csv']);
writematrix(x1,[data_dir,'symptoms_all\brain_boot_index.csv']);
writematrix(Clinicloadings,[data_dir,'symptoms_all\Clinic_boot.csv']);
writematrix(y1,[data_dir,'symptoms_all\Clinic_boot_index.csv']);
