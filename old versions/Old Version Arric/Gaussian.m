function [ integral ] = Gaussian( fx,~,GaussPoints,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch GaussPoints
    case 1
        gp = 0;
        gw = 2;
        
    case 2
        gp = [-sqrt(1/3) sqrt(1/3)]';
        gw = [1 1]';
        
    case 3
        gp = [-sqrt(3/5) 0 sqrt(3/5)]';
        gw = [5/9 8/9 5/9]';
        
    case 4
        gp = [-sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5))]';
        gw = [(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36]';
        
    case 5
        gp = [-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))]';
        gw = [(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900]';
        
    case GaussPoints
        n=GaussPoints;
        x=sym('x','real');
        dg=(x^2-1)^n;
        for i=1:n
            dg=diff(dg,x);
        end
        P=1/(2^n*factorial(n))*dg;
        gp=solve(P);
        gw=zeros(n,1);
        for i=1:length(gp)
            gw(i)=2/((1-gp(i)^2)*(subs(P,'x',gp(i)))^2);
        end
end

fg=zeros(GaussPoints,1);

if verLessThan('matlab','8.2.0')
    for k=1:GaussPoints
        fg(k,1)=subs(fx,symvar(fx),gp(k));
    end
else
    for k=1:GaussPoints
        fg(k,1)=eval(subs(fx,symvar(fx),gp(k)));
    end
end

integral=gw'*fg;

end
