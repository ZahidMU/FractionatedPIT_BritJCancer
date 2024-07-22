%% Figure 4: chemo vs. PDT
load DoseResponse_PDTChemo_params.mat
figure(6); clf
A = panel();
fun = @(x,xdata)(x(2)-x(1))./(1+exp(x(3).*log(xdata)-log(x(4))))+x(1);

A.pack([60 40],3)
color_vec = [0.9000    0.7000    0.8000;
             0.9000    0.7000         0;
             0.9000         0    0.8000;
             0.7000         0    0.7000;
                  0    0.9000    0.7000;
             0.7000    0.9000    0.7000;
             0.7000    0.9000    0.7000;
             0.7000    0.9000    0.7000;
             0.7000    0.9000    0.7000;
             0         0.8000         0];
         
         q_map = colormap(copper(61));
         q_map = flipud(q_map);

         q_map2 = colormap(copper(101));
         q_map2 = flipud(q_map2);

I0 = [0.5 450];

Dose_vec1 = 0:100;
Dose_vec1 = Dose_vec1';

Dose_vec2 = 0:60;
Dose_vec2 = Dose_vec2';

Dose_coeff_chemo = [1 fun(p_Tcell_Cisplatin,Dose_vec1(2:end))';1 fun(p_Tumor_Cisplatin,Dose_vec1(2:end))'];
Dose_coeff_chemo = Dose_coeff_chemo';

Dose_coeff_PDT = [1 fun(p_Tcell_BPD,Dose_vec2(2:end))';1 fun(p_Tumor_BPD,Dose_vec2(2:end))'];
Dose_coeff_PDT = Dose_coeff_PDT';

Dose_coeff_PDTAb = [1 fun(p_Tcell_CetBPD,Dose_vec2(2:end))';1 fun(p_Tumor_CetBPD,Dose_vec2(2:end))'];
Dose_coeff_PDTAb = Dose_coeff_PDTAb';


A(1,1).select()
title('Chemotherapy')
ABC = A(2,1);
ABC.pack(2,1)

A(1,1).select()

r = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1)+100);

set(r,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on
    
   g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))];
b = [0                                  0 30  30];

q = patch(a,b,'b');
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on


clear h
for k = 1:length(Dose_vec2)
A(1,1).select()
    
initCond = I0.*Dose_coeff_PDT(k,:);

   sols = solve(initCond);
   y_vec = sols.y(2,:)/max(y)*(Npoints-1);
   x_vec = sols.y(1,:)/max(x)*(Npoints-1);
   t_vec = sols.x;


    hold on

if k>1
    x_seg = [initCond(1)/max(x)*(Npoints-1) I0(1)/max(x)*(Npoints-1)];
    y_seg = [initCond(2)/max(y)*(Npoints-1) I0(2)/max(y)*(Npoints-1)];
    h10 = plot(x_seg,y_seg,'.','linewidth',1.5,'markersize',13);
    
    set(h10,'color',q_map(k,:))
end

hold on
ABC(2,1).select()
m = plot(t_vec,x_vec,'-','linewidth',1.5);
set(m,'color',q_map(k,:))
hold on

hold on
ABC(1,1).select()
m = plot(t_vec,y_vec,'-','linewidth',1.5);
set(m,'color',q_map(k,:))

hold on





end





A(1,1).select()
    x_seg = [I0(1)/max(x)*(Npoints-1) I0(1)/max(x)*(Npoints-1)];
    y_seg = [I0(2)/max(y)*(Npoints-1) I0(2)/max(y)*(Npoints-1)];
    h11 = plot(x_seg,y_seg,'.','linewidth',1.5,'markersize',13);
    set(h11,'color',[0 0.75 0.75])

g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'b--','linewidth',1.2);
set(g,'color',[0 0 0.75])

axis([0 14 0 30])

set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xtick',[0 7 14],'ytick',[0 15 30])

     box off
 
        axis square


ABC(2,1).select()

set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xlim',[0 12],'ylim',[0 6],'xtick',[0 6 12],'xticklabel',[0 2 4],'ytick',[0 3 6])
    h1 = plot(0,x_seg,'.','linewidth',1.5,'markersize',13);
    set(h1,'color',[0 0.75 0.75])
        h3 = hline(x_seg,'--');
        set(h3,'color',[0 0.75 0.75],'linewidth',1)

ABC(1,1).select()
    
set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xlim',[0 12],'ylim',[0 33],'xtick',[0 6 12],'xticklabel',[0 2 4],'ytick',[0 15 30])
    h2 = plot(0,y_seg,'.','linewidth',1.5,'markersize',13);
    set(h2,'color',[0 0.75 0.75])    
        %h4 = hline(y_seg,'--');
        %set(h4,'color',[0 0.75 0.75],'linewidth',1)

A(1,3).select()
title('PIT')

DEF = A(2,3);
DEF.pack(2,1)


r = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1)+100);

set(r,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on
    
   g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))];
b = [0                                  0 30  30];

q = patch(a,b,'b');
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on


clear h
for k = 1:length(Dose_vec2)
   A(1,3).select()
 
initCond = I0.*Dose_coeff_PDTAb(k,:);

   sols = solve(initCond);
   y_vec = sols.y(2,:)/max(y)*(Npoints-1);
   x_vec = sols.y(1,:)/max(x)*(Npoints-1);
   t_vec = sols.x;

if k>1
    x_seg = [initCond(1)/max(x)*(Npoints-1) I0(1)/max(x)*(Npoints-1)];
    y_seg = [initCond(2)/max(y)*(Npoints-1) I0(2)/max(y)*(Npoints-1)];
    h(k-1) = plot(x_seg,y_seg,'.','linewidth',1.5,'markersize',13);
    
    set(h(k-1),'color',q_map(k,:))
end

set(gca,'tickdir','out','linewidth',1,'fontsize',14)


hold on
DEF(2,1).select()
m = plot(t_vec,x_vec,'-','linewidth',1.5);
set(m,'color',q_map(k,:))

hold on

hold on
DEF(1,1).select()
m = plot(t_vec,y_vec,'-','linewidth',1.5);
set(m,'color',q_map(k,:))

hold on

end
A(1,3).select()
    x_seg = [I0(1)/max(x)*(Npoints-1) I0(1)/max(x)*(Npoints-1)];
    y_seg = [I0(2)/max(y)*(Npoints-1) I0(2)/max(y)*(Npoints-1)];
    h = plot(x_seg,y_seg,'.','linewidth',1.5,'markersize',13);
    
    set(h,'color',[0 0.75 0.75])
g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
set(g,'color',[0 0 0.75])
    axis([0 14 0 30]); box off
 set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xtick',[0 7 14],'ytick',[0 15 30])
 
    axis square
    
DEF(2,1).select()

set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xlim',[0 12],'ylim',[0 6],'xtick',[0 6 12],'xticklabel',[0 2 4],'ytick',[0 3 6])
    h1 = plot(0,x_seg,'.','linewidth',1.5,'markersize',13);
    set(h1,'color',[0 0.75 0.75])
        h3 = hline(x_seg,'--');
        set(h3,'color',[0 0.75 0.75],'linewidth',1)


DEF(1,1).select()
    h2 = plot(0,y_seg,'.','linewidth',1.5,'markersize',13);
    set(h2,'color',[0 0.75 0.75])    
        %h4 = hline(y_seg,'--');
        %set(h4,'color',[0 0.75 0.75],'linewidth',1)

set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xlim',[0 12],'ylim',[0 33],'xtick',[0 6 12],'xticklabel',[0 2 4],'ytick',[0 15 30])

GHI = A(2,2);
GHI.pack(2,1) 

A(1,2).select()
title('PDT')

r = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1)+100);

set(r,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on
    
   g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))];
b = [0                                  0 30  30];

q = patch(a,b,'b');
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on


clear h
for k = 1:length(Dose_vec1)
    
initCond = I0.*Dose_coeff_chemo(k,:);

   sols = solve(initCond);
   y_vec = sols.y(2,:)/max(y)*(Npoints-1);
   x_vec = sols.y(1,:)/max(x)*(Npoints-1);
   t_vec = sols.x;


%pause

A(1,2).select()

if k>1
    x_seg = [initCond(1)/max(x)*(Npoints-1) I0(1)/max(x)*(Npoints-1)];
    y_seg = [initCond(2)/max(y)*(Npoints-1) I0(2)/max(y)*(Npoints-1)];
    h(k-1) = plot(x_seg,y_seg,'.','linewidth',1.5,'markersize',13);
    
    set(h(k-1),'color',q_map2(k,:))
end

set(gca,'tickdir','out','linewidth',1,'fontsize',14)

hold on
GHI(2,1).select()
m = plot(t_vec,x_vec,'-','linewidth',1.5);
set(m,'color',q_map2(k,:))

hold on

hold on
GHI(1,1).select()
m = plot(t_vec,y_vec,'-','linewidth',1.5);
set(m,'color',q_map2(k,:))

hold on


end

A(1,2).select()
    x_seg = [I0(1)/max(x)*(Npoints-1) I0(1)/max(x)*(Npoints-1)];
    y_seg = [I0(2)/max(y)*(Npoints-1) I0(2)/max(y)*(Npoints-1)];
    h = plot(x_seg,y_seg,'.','linewidth',1.5,'markersize',13);
    
    set(h,'color',[0 0.75 0.75])
g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
set(g,'color',[0 0 0.75])
    axis([0 14 0 30]); box off
 set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xtick',[0 7 14],'ytick',[0 15 30])

        axis square

GHI(2,1).select()
axis([0 12 0 5])
set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xlim',[0 12],'ylim',[0 6],'xtick',[0 6 12],'xticklabel',[0 2 4],'ytick',[0 3 6])
    h1 = plot(0,x_seg,'.','linewidth',1.5,'markersize',13);
    set(h1,'color',[0 0.75 0.75])
        h3 = hline(x_seg,'--');
        set(h3,'color',[0 0.75 0.75],'linewidth',1)

GHI(1,1).select()
axis([0 12 0 33])
set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xtick',[0 6 12],'ytick',[0 15 30],'xticklabel',[0 2 4])
    h2 = plot(0,y_seg,'.','linewidth',1.5,'markersize',13);
    set(h2,'color',[0 0.75 0.75])
        %h4 = hline(y_seg,'--');
        %set(h4,'color',[0 0.75 0.75],'linewidth',1)

   A.de.margin = 10; 
    A.export('Fig5_PhasePlaneDoseResponse_PDTChemo.tiff','-w150','-h150', '-rp')