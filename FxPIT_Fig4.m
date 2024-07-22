%%
load PDTchemo_SFdata.mat
load DoseResponse_PDTChemo_params.mat
xrange1 = 1:60; xrange2 = 1:300;
% fun = @(x,xdata)x(1)./(1+exp(x(2).*(xdata-x(3))));
% fun2 = @(x,xdata)x(1).*exp((-xdata.*x(2)));
fun = @(x,xdata)(x(2)-x(1))./(1+exp(x(3).*log(xdata)-log(x(4))))+x(1)


figure(2);clf
clf
B = panel();
B.pack(1,3)

B(1,1).select()
errorbar(light_dose,nanmean(Tumor_BPD_RAW'),nanstd(Tumor_BPD_RAW')./sqrt(sum(~isnan(Tumor_BPD_RAW'))),'r.','linewidth',1.5,'markersize',20)
hold on
plot(xrange1,fun(p_Tumor_BPD,xrange1),'r--','linewidth',1.5)

errorbar(light_dose,nanmean(Tcell_BPD_RAW'),nanstd(Tcell_BPD_RAW')./sqrt(sum(~isnan(Tcell_BPD_RAW'))),'k.','linewidth',1.5,'markersize',20)
plot(xrange1,fun(p_Tcell_BPD,xrange1),'k--','linewidth',1.5)

h1 = hline(1,'k-')

title('Verteporfin + PDT'); xlabel('Dose (J/cm^2)'); ylabel('Survival Fraction')

set(gca,'tickdir','out','linewidth',1.5,'fontsize',14,'xtick',[1 3 10 30 60],...
    'ytick',[0 1 2],'xscale','log','xlim',[0.589,590],'ylim',[-0.05,2.37]); box off; 
axis square
axis([1-0.15 60 0-0.04 2.1])

B(1,2).select()
errorbar(light_dose,nanmean(Tumor_CetBPD_RAW'),nanstd(Tumor_CetBPD_RAW')./sqrt(sum(~isnan(Tumor_CetBPD_RAW'))),'r.','linewidth',1.5,'markersize',20)
hold on
plot(xrange1,fun(p_Tumor_CetBPD,xrange1),'r--','linewidth',1.5)

errorbar(light_dose,nanmean(Tcell_CetBPD_RAW'),nanstd(Tcell_CetBPD_RAW')./sqrt(sum(~isnan(Tcell_CetBPD_RAW'))),'k.','linewidth',1.5,'markersize',20)
plot(xrange1,fun(p_Tcell_CetBPD,xrange1),'k--','linewidth',1.5)
hline(1,'k-')

title('Verteporfin-Cetuximab + PDT'); xlabel('Dose (J/cm^2)'); ylabel('Survival Fraction')
%legend('Tumor Cells', 'T-cells'); box off; set(gca,'linewidth',1.5,'fontsize',12)
%legend('location','northeast','box','off')
set(gca,'tickdir','out','linewidth',1.5,'fontsize',14,'xtick',[1 3 10 30 60],'ytick',[0 1 2],'xscale','log'); box off; 
axis([1-0.15 60 0-0.04 2.1])
axis square

B(1,3).select()
errorbar(cisplatin_dose,nanmean(Tumor_Cisplatin_RAW'),nanstd(Tumor_Cisplatin_RAW')./sqrt(sum(~isnan(Tumor_Cisplatin_RAW'))),'r.','linewidth',1.5,'markersize',20)
hold on
plot(xrange2,fun(p_Tumor_Cisplatin,xrange2),'r--','linewidth',1.5)

errorbar(cisplatin_dose,nanmean(Tcell_Cisplatin_RAW'),nanstd(Tcell_Cisplatin_RAW')./sqrt(sum(~isnan(Tcell_Cisplatin_RAW'))),'k.','linewidth',1.5,'markersize',20)
plot(xrange2,fun(p_Tcell_Cisplatin,xrange2),'k--','linewidth',1.5)
hline(1,'k-')

title('Cisplatin'); xlabel('Dose (\muM)'); ylabel('Survival Fraction')
%legend('Tumor Cells', 'T-cells'); box off; set(gca,'linewidth',1.5,'fontsize',12)
%legend('location','northeast','box','off')
set(gca,'tickdir','out','linewidth',1.5,'fontsize',14,'xtick',10.^[0 1 2],'ytick',[0 1 2],'xscale','log'); box off; 
axis([1-0.15 300 0-0.04 2.1])

legend('Tumor Cells','fit', 'T-cells','fit'); box off; set(gca,'linewidth',1.5,'fontsize',12)
legend('location','northeast','box','off')

axis square




B.de.margin = 15
   B.fontsize = 12;
    B.export('Fig4_PDTChemo_DoseResponse.tiff','-w225','-h155', '-rp')
