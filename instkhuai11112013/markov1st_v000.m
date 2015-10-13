%Convert d->c or c->d 1st order markov parameters
%x'=ax+bu
%mode 1:c->d 2:d->c, fr=frequency of discretication
%function [cm, dm]=markov(a,b,ss_std,mode,fr)

function [cm, dm]=markov1st(a,b,ss_std,mode,fr)
dt=1/fr;
if mode==1
    if (~(isempty(a) || isempty(b)))
        if (a~=0)
            a_c=a;
            b_c=b;
            a_d=exp(a_c*dt);
            var_c=b_c^2/(-2*a_c);
            var_d=var_c;
            b_d=sqrt(var_c*(1-exp(2*a_c*dt)));
        else %randon walk
            a_c=0;
            b_c=b;
            var_c=ss_std;
            var_d=ss_std;
            a_d=1;
            b_d=b*sqrt(dt);
        end 
    elseif (~(isempty(a) || isempty(ss_std)))
        a_c=a;
        var_c=ss_std^2;
        var_d=var_c;
        b_c=sqrt(var_c*(-2*a_c));
        a_d=exp(a_c*dt);
        b_d=sqrt(var_c*(1-exp(2*a_c*dt)));
    elseif (~(isempty(b) || isempty(ss_std)))
        b_c=b;
        var_c=ss_std^2;
        var_d=var_c;
        a_c=b_c^2/(-2)/var_c;
        a_d=exp(a_c*dt);
        b_d=sqrt(var_c*(1-exp(2*a_c*dt)));
    end
else
    if (~(isempty(a) || isempty(b)))
        if (a~=1)
            a_d=a;
            b_d=b;
            var_d=b_d^2/(1-a_d^2);
            var_c=var_d;
            a_c=log(a_d)/dt;
            b_c=sqrt(var_c*(-2*a_c));
        else
            a_d=a;
            b_d=b;
            var_d=ss_std;
            var_c=ss_std;
            a_c=0;
            b_c=b_d/sqrt(dt);
        end
    elseif (~(isempty(a) || isempty(ss_std)))
        a_d=a;
        var_d=ss_std^2;
        var_c=var_d;
        b_d=sqrt(var_d*(1-a_d^2));
        a_c=log(a_d)/dt;
        b_c=sqrt(var_c*(-2*a_c));
    elseif (~(isempty(b) || isempty(ss_std)))
        b_d=b;
        var_d=ss_std^2;
        var_c=var_d;
        a_d=sqrt(var_d*(1-a_d^2)); %%note only +ve sol is used. (but -ve sol is also valid!)
        a_c=log(a_d)/dt;    %-ve solution creates problem in this step
        b_c=sqrt(var_c*(-2*a_c));
    end
end
T=-1/a_c;
cm=[a_c;b_c;var_c;T];
dm=[a_d;b_d;var_d;T];
