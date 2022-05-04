classdef adaptive_filter
    %Performs sample-wise adaptive filtering
    %   input x and e, outputs r
    
    properties
        length;  %number of filter coefficients
        type;    %filter type
        adaptation_constant;
        w=[];           %filter coefficients
        x_delayed = []; %delayed input signals
        w_history = []; %keep track of filter coefficients as they evolve
        r;
        
        %add your own variables here
        %add estimate of the input signal power
        deltax;
         %add forgetting factor ,g(k+1),Rxverse(k+1),rex(k+1)for RLS
        gamma;
        g=[];
        Rxverse=[];
        rex=[];
        e_delayed=[];
        %%FDAF
        F=[];
        P=[];
        
    end
    
    methods
        function obj = adaptive_filter(length,type,adaptation_constant)
            if(nargin>0)
            obj.length=length;
            obj.type=type;
            obj.adaptation_constant=adaptation_constant;
            obj.w=zeros(length,1);
            obj.x_delayed=zeros(length,1);
            obj.w_history=zeros(0,length);
            
            %initialize your variables here
            %add initialize of the input signal power
            obj.deltax=0;
            %add initialize of the factors for RLS
            obj.gamma=1-10^(-4);
            obj.g=zeros(length,1);
            rdel=[0.01,0.01];
            obj.Rxverse=diag(rdel);
            obj.rex=zeros(length,1);
            obj.e_delayed=zeros(length,1);
            %add initialize of FDAF
            obj.F = dftmtx(length);
            obj.P= eye(length);
            %obj.x_hist=zeros(length,1);
            %obj.X_freq=zeros(length,1)
            %obj.w_old=zeros(1,length);
            %obj.PI=inv(diag(x_hist));

            end
        end
          function obj = filter(obj,x,e)
            obj.x_delayed=[x ; obj.x_delayed(1:obj.length-1)];
            %%e[k+1]
            e_hat=obj.w.' * obj.x_delayed;
            obj.r=e-e_hat;
            obj.w_history=[obj.w_history ; obj.w.']; %you may want to remove this line to gain speed
            %add update of deltax
            obj.deltax=(obj.x_delayed).' * obj.x_delayed/200 +0.05;
            %add update of RLS
            obj.gamma=1-10^(-1);
            obj.g= obj.Rxverse * obj.x_delayed/((obj.gamma)^2 + (obj.x_delayed).' * obj.Rxverse * obj.x_delayed);
            obj.Rxverse=(obj.Rxverse - obj.g * (obj.x_delayed).' * obj.Rxverse)/(obj.gamma)^2;
            obj.rex=(obj.gamma)^2 * obj.rex + obj.x_delayed * e;
             %%FDAF
            pn=diag(abs(obj.F * obj.x_delayed).^2/obj.length);
            obj.P=0.3* obj.P + 0.7*pn;
            obj= update_filter(obj,e); 
           
            %obj.x_hist = [obj.x_hist ; x];
            %obj.PI=0.3* obj.PI + 0.7*(obj.PI*inv(X_freq).')^(2)/200;

            %obj.x_hist = [obj.x_hist ; x.']; % Something to store history in a column vector
            %nv= length(obj.x_hist);
            %F = dftmtx(nv);
            %X_freq = F*obj.x_hist; % Frequency domain representation of x  (X[k])
            %w_new = obj.w_old + 2*obj.adaptation_constant/200*obj.PI*conj(X_freq)*obj.r; % Implement FDAF update f
            %obj.PI=0.3*obj.PI + 0.7*(inv(X_freq)*inv(X_freq).')^(2)/200;
            %obj.w_old = w_new; % Store new values as old values for later
            %obj.w = w_new*F; % Make sure that w*x => w*F*x = w*X
      
        end
    end
    
end

