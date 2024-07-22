%% DEFINING MODEL PARAMETERS
% sigma = 0.118; rho = 1.131; eta = 20.19; mu = 0.00311;
% delta = 0.374; alpha = 1.636; beta = 0.002;

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1; 
 
%% DEFINING THE MODEL (INLINE FUNCTION)
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-x(1,:).*x(2,:)]);
 
%% FUNCTION RETURNING MODEL SOLUTION ON [0,100] FOR GIVEN INITIAL CONDITION
 options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));
%% EVALUATING THE VECTOR FIELD IN [0, 3.5]x[0, 450]
Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 
dx = x(2)-x(1);
 
dy = y(2)-y(1);
[X, Y] = meshgrid(x,y);
G = rhs([],[reshape(X,1,[]); reshape(Y,1,[])]);
U = reshape(G(1,:),Npoints,Npoints);
V = reshape(G(2,:),Npoints,Npoints)*dx/dy;
N = sqrt(U.^2+V.^2);
U = U./N; V = V./N;

%% EVALUATING MODEL SOLUTIONS FOR DIFFERENT INITIAL CONDITIONS

sigma = 0.118; rho = 0.95;    eta = 20.19;  mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002; gamma = 1;   
    [curve] = Kuznetsov_SeparatrixCalc(sigma,rho,eta,mu,delta,alpha,beta,gamma);

Npoints = 30
initCond = [0.15, 100;
            0.8,  400;
            2.5,    375;
            0.95,    110;]; 
    
sols = cell(1,size(initCond,1));
 options = odeset('Refine',100);


for i = 1:size(initCond,1)
sols{i} = solve(initCond(i,:));
end


figure(4)
clf
A = panel();
A.pack(1,2)
clf
AB = A(1,1);
AB.pack(2,2)

A(1,2).select()


r = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1)+100)

set(r,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on
    
   g1 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1))

set(g1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])

   g2 = area(curve.x./max(x)*(Npoints-1),curve.y./max(y)*(Npoints-1))

set(g2,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
a = [max(curve.x./max(x)*(Npoints-1))  30 30 max(curve.x./max(x)*(Npoints-1))]
b = [0                                  0 30  30]
%patch([xu(:)' fliplr(xu(:)')], [ymin fliplr(ymax)], [0.1 0.9 0.1], 'FaceAlpha',0.3) % Plot Filled Background
q = patch(a,b,'b')
set(q,'FaceColor',[0 1 0],'EdgeColor',[0 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.1,'linewidth',0.01)
hold on
[X1, Y1] = meshgrid(0:Npoints-1,0:Npoints-1);

 q = streamslice(X1,Y1,U,V); %plotting vector field
set(q,'color',[0.5 0.5 0.5])
 %q.AutoScaleFactor = 0.5;
hold on
plot(initCond(:,1)/max(x)*(Npoints-1),initCond(:,2)/max(y)*(Npoints-1),'k.','markersize',15)

color_vec = [0.5  0    0.75;
             0.75 0    0;
             0    0.75  0;
             0    0.5 0.75]

for i = 1:length(sols) %plotting each solution

    x_start = initCond(i,1)/max(x)*(Npoints-1);
    y_start = initCond(i,2)/max(y)*(Npoints-1);
t_vec = sols{i}.x;

    
A(1,2).select()

h1 = plot(sols{i}.y(1,:)/max(x)*(Npoints-1),sols{i}.y(2,:)/max(y)*(Npoints-1),'-','linewidth',1.25)
%h1 = streamline(X1,Y1,U,V,x_start,y_start)
set(h1,'color',color_vec(i,:),'linewidth',1.5)
axis([0 25 0 29])

hold on

    if i == 1
    
        AB(1,1).select()
        h2=plot(t_vec,sols{i}.y(2,:)/max(y)*(Npoints-1),'-','linewidth',1.5)
        hold on
        %set(gca,'ylim',[0 30],'xlim',[0 40])
        set(h2,'color',color_vec(i,:))
        axis([0 12 0 30])
    
        hold on
        %AB(2,1).select()
        h3=plot(t_vec,sols{i}.y(1,:)/max(x)*(Npoints-1),'-.','linewidth',1.5)
        hold on
        set(h3,'color',color_vec(i,:))
        axis([0 12 0 30])

        
    elseif i == 2
        AB(1,2).select()
        h4=plot(t_vec,sols{i}.y(2,:)/max(y)*(Npoints-1),'-','linewidth',1.5)
        hold on
        %set(gca,'ylim',[0 30],'xlim',[0 40])
        set(h4,'color',color_vec(i,:))
        axis([0 12 0 30])
    
        hold on
        %AB(2,2).select()
        h5=plot(t_vec,sols{i}.y(1,:)/max(x)*(Npoints-1),'-.','linewidth',1.5)
        hold on
        set(h5,'color',color_vec(i,:))
        axis([0 12 0 30])
    
    elseif i == 3
        AB(2,1).select()
        h4=plot(t_vec,sols{i}.y(2,:)/max(y)*(Npoints-1),'-','linewidth',1.5)
        hold on
        %set(gca,'ylim',[0 30],'xlim',[0 40])
        set(h4,'color',color_vec(i,:))
        axis([0 60 0 30])
    
        hold on
        %AB(2,2).select()
        h5=plot(t_vec,sols{i}.y(1,:)/max(x)*(Npoints-1),'-.','linewidth',1.5)
        hold on
        set(h5,'color',color_vec(i,:))
        axis([0 60 0 30])

    elseif i ==4
        AB(2,2).select()
        h4=plot(t_vec,sols{i}.y(2,:)/max(y)*(Npoints-1),'-','linewidth',1.5)
        hold on
        %set(gca,'ylim',[0 30],'xlim',[0 40])
        set(h4,'color',color_vec(i,:))
        axis([0 60 0 30])
    
        hold on
        %AB(2,2).select()
        h5=plot(t_vec,sols{i}.y(1,:)/max(x)*(Npoints-1),'-.','linewidth',1.5)
        hold on
        set(h5,'color',color_vec(i,:))
        axis([0 60 0 30])       
    end
        plot(0,initCond(i,2)/max(y)*(Npoints-1),'k.','markersize',15)
        plot(0,initCond(i,1)/max(x)*(Npoints-1),'k.','markersize',15)

    
end


A(1,2).select()


 set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xtick',[0 12 24],'ytick',[0 14 28])
% hold on
g = plot((curve.x)./max(x)*(Npoints-1),(curve.y)./max(y)*(Npoints-1),'g--','linewidth',2);
set(g,'color',[0 0 0.75])
set(gca,'xtick',[0 12 24],'ytick',[0 14 28])%,'xticklabel',[],'yticklabel',[])
axis square



AB(1,1).select()
set(gca,'tickdir','out','linewidth',1,'fontsize',14)
set(gca,'xtick',[0 6 12],'ytick',[0 15 30])
set(gca,'xticklabel',round([0 6 12]*9.91/30,1))

axis square

AB(2,1).select()
set(gca,'tickdir','out','linewidth',1,'fontsize',14)
set(gca,'xtick',[0 30 60],'ytick',[0 15 30])
set(gca,'xticklabel',round([0 30 60]*9.91/30))

axis square

AB(1,2).select()
set(gca,'tickdir','out','linewidth',1,'fontsize',14)
set(gca,'xtick',[0 6 12],'ytick',[0 15 30])
set(gca,'xticklabel',round([0 6 12]*9.91/30,1))

axis square

AB(2,2).select()
set(gca,'tickdir','out','linewidth',1,'fontsize',14)
set(gca,'xtick',[0 30 60],'ytick',[0 15 30])
set(gca,'xticklabel',round([0 30 60]*9.91/30))
%plot(0,initCond(3:4,1)/max(x)*(Npoints-1),'k.','markersize',15)

axis square


A.de.margin = 10
   A.fontsize = 12;
    A.export('Fig3_CancerImmunePhasePortrait+TimeTraces_alt.tiff','-w150','-h75', '-rp')
    
