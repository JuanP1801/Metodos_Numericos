function [n, xn, fm, E] = secante(x0,x1,Tol,niter)
syms x
f=sin(x)-(100*10^-9)
c=0;
fm(c+1)=eval(subs(f,x0));
fe0=fm(c+1) ;
E(c+1) =Tol+2;
xn (c+1)=x0;
N(c+1)=c;
    if fe0==0
        s=X0;
        n=c;
        fprintf ('% es raiz de f(x) ' , x0)
    end 
c=c+1;
fm(c+1) = eval(subs (f, x1)) ;
fe=fm(c+1);
den=fe-fe0;
E(c+1) =Tol+1;
error=E(c+1);
xn(c+1)=x1;
N(c+1)=c;
    while error>Tol && fe~=0 && den~=0 && c<niter
        xn(c+2)=x1-fe*(x1-x0)/den;
        %Este es decimales correctos
        E(c+2) =abs ( (xn (c+2) -x1));
        %este es cifras significativas
        %E(c+2) =abs ( (xn (c+2) -x1))/abs(xn(c+2));
        error=E(c+2);
        x0=x1;
        fe0=fe;
        x1=xn(c+2);
        fm(c+2) =eval(subs (f, x1)) ;
        fe=fm(c+2); den=fe-fe0;
        N (c+2)=c+1; c=c+1;
    end 
format long
disp('              n                 Xn                  Fx               Error  ')
D=[N' xn' fm' E'];
disp (D)
    if fe==0
        n=C;
        fprintf('%f es raiz de f(x) \n', x1)
    elseif error<Tol
        n=c;
        fprintf('%f es una aproximación de una raiz de f(x) con una tolerancia= %f \n', x1)
    elseif den==0
        n=c;
        fprintf('%f es una posible raiz de f(x) \n',x0)
    else
        s=x0;
        n=c;
        fprintf('Fracasó en %f iteraciones \n',niter)
    end
