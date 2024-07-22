%close all
I0 = [0.5 450];
I0_orig = I0;

%num_fx = 3;
fx_dt = 2;

% Dose_vec_str = {'0 J/cm^2' '1 J/cm^2' '3 J/cm^2' '10 J/cm^2' '30 J/cm^2' '60 J/cm^2'};
% Dose_vec = Dose';

%Dose_vec = [1 3 10 30 60]

Dose_vec1 = round(linspace(0,100,5));
Dose_vec1 = Dose_vec1';

Dose_vec2 = round(linspace(0,60,5));
%Dose_vec2 = [0  .1]
Dose_vec2 = Dose_vec2';

Dose_coeff_chemo = [1 fun(p_Tcell_Cisplatin,Dose_vec1(2:end))';1 fun(p_Tumor_Cisplatin,Dose_vec1(2:end))']
Dose_coeff_chemo = Dose_coeff_chemo'

Dose_coeff_PDT = [1 fun(p_Tcell_BPD,Dose_vec2(2:end))';1 fun(p_Tumor_BPD,Dose_vec2(2:end))']
Dose_coeff_PDT = Dose_coeff_PDT'

Dose_coeff_PDTAb = [1 fun(p_Tcell_CetBPD,Dose_vec2(2:end))';1 fun(p_Tumor_CetBPD,Dose_vec2(2:end))']
Dose_coeff_PDTAb = Dose_coeff_PDTAb'


    
       %Dose_vec = [40 10 77 40];
       dose_vec = [15 3 60 6];

       num_fx_vec = [1 5 1 10];
       %dose_vec = [40 20 5]
       
       
       %dose_vec = fliplr(dose_vec)
       %num_fx_vec = round(Dose_vec./dose_vec)


%for j = 1:length(num_fx_vec)

num_sims = length(dose_vec);


figure(21); clf
%hold on
B = panel();
B.pack(3,2)

pos_array = [1,1;
             1,2;
             2,1;
             2,2];

for j = 1:length(dose_vec)
 
   I0 = I0_orig;
start_time = 0;

   dose = dose_vec(j);
   num_fx = num_fx_vec(j);

   
   subB = B(pos_array(j,1), pos_array(j,2));
   subB.pack(1,2)
   subB(1,1).select();
   
   r = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1)+100)

set(r,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on
    
   g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1))

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1))

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  33 33 max(curve.x./max(x)*(Npoints-1))]
b = [0                                  0 33  33]
%patch([xu(:)' fliplr(xu(:)')], [ymin fliplr(ymax)], [0.1 0.9 0.1], 'FaceAlpha',0.3) % Plot Filled Background
q = patch(a,b,'b')
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on


        plot(I0_orig(1)/max(x)*(Npoints-1),I0_orig(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',12)
        hold on
   
   subC=subB(1,2);
   subC.pack(2,1)
   subC(2,1).select();
        plot(start_time,I0_orig(1)/max(x)*(Npoints-1),'b.','linewidth',1.5,'markersize',12)
        hold on
        h3 = hline(I0_orig(1)/max(x)*(Npoints-1),'--');
        set(h3,'color',[0 0.75 0.75],'linewidth',1)

   subC(1,1).select();
        plot(start_time,I0_orig(2)/max(y)*(Npoints-1),'b.','linewidth',1.5,'markersize',12)
        hold on

   for k = 1:num_fx
       %start_time = (k-1)*fx_dt+1;
       
       initCond = I0.*[fun(p_Tcell_CetBPD,dose) fun(p_Tumor_CetBPD,dose)];

       sols = solve(initCond);
       y_vec = sols.y(2,:)/max(y)*(Npoints-1);
       x_vec = sols.y(1,:)/max(x)*(Npoints-1);
       
       if k < num_fx
           
   subB(1,1).select();
            plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
            plot(x_vec(1:fx_dt),y_vec(1:fx_dt),'k','linewidth',1.5); 

   subC(2,1).select();
            plot(start_time,initCond(1)/max(x)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
            plot(start_time:fx_dt:start_time+fx_dt,x_vec(1:fx_dt),'k','linewidth',1.5); 

   subC(1,1).select();
            plot(start_time,initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
            plot(start_time:fx_dt:start_time+fx_dt,y_vec(1:fx_dt),'k','linewidth',1.5); 
   start_time = start_time + fx_dt;

           
       elseif k == num_fx
   subB(1,1).select();
              plot(x_vec,y_vec,'k','linewidth',1.5); 
              plot(initCond(1)/max(x)*(Npoints-1),initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)

   subC(2,1).select();
              plot(start_time,initCond(1)/max(x)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
              plot(start_time+(0:length(x_vec)-1),x_vec,'k','linewidth',1.5); 

   subC(1,1).select();
              plot(start_time,initCond(2)/max(y)*(Npoints-1),'r.','linewidth',1.5,'markersize',12)
              plot(start_time+(0:length(y_vec)-1),y_vec,'k','linewidth',1.5); 
   start_time = start_time + fx_dt;

       end


        I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);
   
       
   end

   
   subB(1,1).select();
          set(gca,'tickdir','out','linewidth',1.5,'xtick',[0 15 30],'ytick',[0 15 30])
           g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0 0.75])
           axis([0 33 0 31]); box off; axis square
           title([num2str(dose) ' J/cm^2 x ' num2str(num_fx) ])

   subC(2,1).select();
          set(gca,'tickdir','out','linewidth',1.5)
          ylabel('E'); %axis square
        h3 = hline(I0_orig(1)/max(x)*(Npoints-1),'--');
        set(h3,'color',[0 0.75 0.75],'linewidth',1)
        
          if j == 1
             set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 5],'ytick',[0 2.5 5])
          elseif j == 2
              set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 30],'ytick',[0 15 30])
          elseif j == 3
              set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 10],'ytick',[0 5 10])              
          elseif j == 4
              set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 33],'ytick',[0 15 30])
          end
% 
   subC(1,1).select();
          set(gca,'tickdir','out','linewidth',1.5)
          ylabel('T'); %axis square
          
          if j == 1
             set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 30],'ytick',[0 15 30])
          elseif j == 2
             set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 30],'ytick',[0 15 30])              
          elseif j == 3
              set(gca,'xlim',[0 20],'xtick',[0 10 20],'ylim',[0 30],'ytick',[0 15 30])
          elseif j == 4
              set(gca,'xlim',[0 60],'xtick',[0 30 60],'ylim',[0 30],'ytick',[0 15 30])
          end
%    
end

load OutcomeFields.mat

% Fractionated PIT, fixed fraction size

num_ICs = sum(sum(outcome_mat == 0));

IC_vec = zeros(num_ICs,2);


count = 1;
for i = 1:length(x)
    for j = 1:length(y)
        
        if outcome_mat(i,j) == 0
            
            IC_vec(count,:) = initCond_mat{i,j};
            count = count+1;
            
        end
        
    end
end


   E = B(3,1);
   
   E.pack(1,2);
   E(1,1).select();
   
   dose = 1;
   
   [NumFx_map,NumFx_vec] = MinPITFxNum(dose,IC_vec,num_ICs,fun,p_Tcell_CetBPD,p_Tumor_CetBPD);
   NumFx_map(NumFx_map==0)=1000; 

    image(NumFx_map); axis xy

    a = autumn;
    %colorbar; 
    colormap(flipud(a(end-25:end,:)))
    colormap(flipud(autumn))

    hold on

    set(gca,'linewidth',1.5,'tickdir','out','xtick',[1 7 14],'ytick',[1 15 30],'xticklabel',[0 7 14],'yticklabel',[0 15 30])
    %xlabel('Immune Effector Cells'); ylabel('Tumor Cells')
    title([num2str(dose) ' J/cm^2/fx'])
       g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0.75 0.75])
                  axis([1 14 1 30]); box off; axis square

              g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))];
b = [0                                  0 30  30];

q = patch(a,b,'b');
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on

g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0 0.75])   
   E(1,2).select();
   dose = 3
   [NumFx_map,NumFx_vec] = MinPITFxNum(dose,IC_vec,num_ICs,fun,p_Tcell_CetBPD,p_Tumor_CetBPD);
   NumFx_map(NumFx_map==0)=1000; 

    image(NumFx_map); axis xy

    a = autumn;
    %colorbar; 
    colormap(flipud(a(end-25:end,:)))
    colormap(flipud(autumn))

    hold on

    set(gca,'linewidth',1.5,'tickdir','out','xtick',[1 7 14],'ytick',[1 15 30],'xticklabel',[0 7 14],'yticklabel',[0 15 30])
    %xlabel('Immune Effector Cells'); ylabel('Tumor Cells')
    title([num2str(dose) ' J/cm^2/fx'])
       g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0.75 0.75])
                  axis([1 14 1 30]); box off; axis square

              g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))];
b = [0                                  0 30  30];

q = patch(a,b,'b');
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on

g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0 0.75])       
   F = B(3,2);
   F.pack(1,2);

   F(1,1).select();
   dose = 5;
   [NumFx_map,NumFx_vec] = MinPITFxNum(dose,IC_vec,num_ICs,fun,p_Tcell_CetBPD,p_Tumor_CetBPD);
   NumFx_map(NumFx_map==0)=1000; 

    image(NumFx_map); axis xy

    a = autumn;
    %colorbar; 
    colormap(flipud(a(end-25:end,:)))
    colormap(flipud(autumn))


    hold on

    set(gca,'linewidth',1.5,'tickdir','out','xtick',[1 7 14],'ytick',[1 15 30],'xticklabel',[0 7 14],'yticklabel',[0 15 30])
    %xlabel('Immune Effector Cells'); ylabel('Tumor Cells')
    title([num2str(dose) ' J/cm^2/fx'])
       g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0.75 0.75])
                  axis([1 14 1 30]); box off; axis square

              g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))];
b = [0                                  0 30  30];

q = patch(a,b,'b');
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on

g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0 0.75])       
       F(1,2).select();
       dose = 7;
   [NumFx_map,NumFx_vec] = MinPITFxNum(dose,IC_vec,num_ICs,fun,p_Tcell_CetBPD,p_Tumor_CetBPD);

   NumFx_map(NumFx_map==0)=1e5; 
   image(NumFx_map); axis xy

    a = flipud(copper);
    %colorbar; 

    colormap(a)

    %colormap(flipud(copper))

%colormap hsv
    hold on

    set(gca,'linewidth',1.5,'tickdir','out','xtick',[1 7 14],'ytick',[1 15 30],'xticklabel',[0 7 14],'yticklabel',[0 15 30])
    %xlabel('Immune Effector Cells'); ylabel('Tumor Cells')
    title([num2str(dose) ' J/cm^2/fx'])
                  axis([1 14 1 30]); box off; axis square

              g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1));

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))];
b = [0                                  0 30  30];

q = patch(a,b,'b');
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on

g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',1);
           set(g,'color',[0 0 0.75])

       
B.fontsize = 14;

   B.de.margin = 15 ;
    B.export('Fig6_FxPIT.tiff','-w200','-h150', '-rp')
    