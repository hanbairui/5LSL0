function [ filter ] = update_filter( filter,e )

x=filter.x_delayed;
w_old = filter.w;
% r=filter.r;
filter_type=filter.type;

filter.w=w_old; %default output

%% A1 scenario 1:f
if strcmpi(filter_type,'SGD')
%     implement the SGD update rule here
    alpha=filter.adaptation_constant;
    Rx=[2 -1;-1 2];
    rex=[0;3];
    filter.w=filter.w+2*alpha*(rex-Rx*filter.w);
end

%% A1 scenario 1:i
if strcmpi(filter_type,'Newton')
    %implement the Newton update rule here
    alpha=filter.adaptation_constant;
    Rx=[2 -1;-1 2];
    rex=[0;3];
    filter.w=filter.w+2*alpha*inv(Rx)*(rex-Rx*filter.w);
end

%% A1 scenario 2:a
%%add deltax for NLMS
deltax=filter.deltax;
%%
filter.w=w_old; %default output
%%2:a
if strcmpi(filter_type,'LMS')
    %implement the LMS update rule here
    alpha=filter.adaptation_constant;
    filter.w = filter.w + 2*alpha*x*r; 
    
end
%% A1 scenario 2:b
if strcmpi(filter_type,'NLMS')
    %implement the NLMS update rule here
    alpha=filter.adaptation_constant;
    filter.w= filter.w + 2*alpha*x*r/deltax;
end
%% A1 scenario 2:d
if strcmpi(filter_type,'RLS')
    %implement the RLS update rule here
   % lambda=filter.adaptation_constant;   
    filter.w= filter.Rxverse*filter.rex;
end
%% A1 scenario 2:e
if strcmpi(filter_type,'FDAF')
    %implement the FDAF update rule here
    alpha=filter.adaptation_constant;
    X_freq = filter.F * filter.x_delayed; % Frequency domain representation of x  (X[k])
    w_old = inv(filter.F)*filter.w; % Implement FDAF update f
    w_new = w_old+ 2*alpha*inv(filter.P)*conj(X_freq)*filter.r; % Store new values as old values for later
    filter.w = filter.F * w_new; % Make sure that w*x => w*F*x = w*X
end
% 
% 
 end

