function [vb,vs,mub,mus] = binodal(f,mu,dmu)
% Adam Optimization used for calculating the binodal and spinodal of a
% 2-component system determined through inline functions.
%
% Input:
% f = @(v)f(v) : free-energy density
% mu = @(v)mu(v) : chemical potential (1st derivative of f)
% dmu = @(v)dmu(v) : 1st derivative of mu
%
% Output:
% vb = [dilute phase vol.frac., dense phase vol.frac.]
% vs = [dilute spinodal vol.frac., dense spinodal vol.frac]
% mub = binodal chemical potential
% mus = mu(vs)

% Adam Optimization parameters
adamparam{1}=0.1;
adamparam{2}=0.9;
adamparam{3}=0.99;
adamparam{4}=1e-8;
adamparam{5}=1e-5;

% Calculate spinodal
vs = calvs(dmu,adamparam);
mus = mu(vs);

% Parameterized variable to bound binodal mu
mubg = @(g)mus(1)*(1./(1+exp(-g))) + mus(2)*(1-1./(1+exp(-g)));

% Calculate binodal
mub = mubg(adamopt(@(x)dmub(f,mu,mubg(x),vs,adamparam),0,adamparam));
vb = calvb(mu,mub,vs,adamparam);
return

function f = adamopt(df, f, param, varargin)
a=param{1};b1=param{2};b2=param{3};ep=param{4};tol=param{5};
if ~isempty(varargin)
    plt = true;
else
    plt = false;
end
x = f + 1;
m1 = zeros(size(x));
m2 = zeros(size(x));
t = 0;
while abs(x-f) > tol
    t = t + 1;
    x = f;
    g = df(x);
    m1 = b1*m1 + (1-b1)*g;
    m2 = b2*m2 + (1-b2)*g^2;
    m1t = m1/(1-b1^t);
    m2t = m2/(1-b2^t);
    f = x - a * m1t/(sqrt(m2t) + ep);
    if plt
        plot(f,df(f),'.k')
        pause(0.1)
    end
end
return

function vs = calvs(dmu,param)
vg = @(g)1./(1+exp(-g));
vs(1) = adamopt(@(x)-dmu(vg(x)),-5,param);
vs(2) = adamopt(@(x)dmu(vg(x)),5,param);
vs = vg(vs);
return

function vb = calvb(mu,mub,vs,param)
vg = @(g)1./(1+exp(-g));
vb(1) = adamopt(@(x)mu(vg(x))-mub,-log(1/vs(1)-1)-5,param);
vb(2) = adamopt(@(x)mu(vg(x))-mub,-log(1/vs(2)-1)+5,param);
vb = vg(vb);
return

function ermu = dmub(f,mu,mub,vs,param)
vb = calvb(mu,mub,vs,param);
ermu = mub - diff(f(vb))/diff(vb);
return
