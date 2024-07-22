function [NumFx_map,NumFx_vec] = MinPITFxNum(dose,IC_vec,num_ICs,fun,p_Tcell_CetBPD,p_Tumor_CetBPD)

%% DEFINING MODEL PARAMETERS
sigma = 0.118; rho = 1.131; eta = 20.19; mu = 0.00311;
delta = 0.374; alpha = 1.636; beta = 0.002;
 
%% DEFINING THE MODEL (INLINE FUNCTION)
rhs = @(t,x)([sigma+rho*x(1,:).*x(2,:)./(eta+x(2,:))-mu*x(1,:).*x(2,:)-delta*x(1,:);...
              alpha*x(2,:).*(1-beta*x(2,:))-x(1,:).*x(2,:)]);
 
%% FUNCTION RETURNING MODEL SOLUTION ON [0,100] FOR GIVEN INITIAL CONDITION

 options = odeset('Refine',100);

solve = @(init)(ode45(rhs,[0 100],init,options));

Npoints = 30;
x = linspace(0,3.5,Npoints);
y = linspace(0,450,Npoints);
 

fx_dt = 2;
NumFx_vec = zeros(1,num_ICs);
for  j = 1:num_ICs
    j;
    I0 = IC_vec(j,:);
    I0_orig = I0;


       I0 = I0_orig;

    start_time = 1;
    k = 1;
    y_final =11;
        while y_final > 10

                   initCond = I0.*[fun(p_Tcell_CetBPD,dose) fun(p_Tumor_CetBPD,dose)];

                   sols = solve(initCond);
                   y_vec = sols.y(2,:)/max(y)*(Npoints-1);
                   x_vec = sols.y(1,:)/max(x)*(Npoints-1);

                    I0 = [x_vec(fx_dt)*max(x) y_vec(fx_dt)*max(y)]./(Npoints-1);

                   sols = solve(I0);
                   y_final = sols.y(2,end)/max(y)*(Npoints-1);

                    start_time = start_time + fx_dt;
                    k = k+1;

            %end
        end
                           NumFx_vec(j) = k;

end


NumFx_map = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        
        indx = find((IC_vec(:,1) == x(i)) .* (IC_vec(:,2) == y(j)));
        
        if isempty(indx) == 0
        NumFx_map(j,i) = NumFx_vec(indx);
        end
        
    end

end


end

