function [x,histout,costdata] = optim_kelley_levmar(x0,f,tol,maxit)
//
// C. T. Kelley, Dec 14, 1997
//
// This code comes with no guarantee or warranty of any kind.
//
// function [x,histout,costdata] = optim_kelley_levmar(x0,f,tol,maxit)
//
// Levenberg-Marquardt code, trust region control of LM parameter
//
//
// Input: x0 = initial iterate
//        f = r^T r/2 = objective function,
//            the calling sequence for f should be
//            [fout,gout,jac]=f(x) where fout=f(x) is a scalar
//              gout = jac^T r = grad f(x) is a COLUMN vector
//              and jac = r' = Jacobian of r is an M x N matrix
//        tol = termination criterion norm(grad) < tol
//        maxit = maximum iterations (optional) default = 100
//
// Output: x = solution
//         histout = iteration history   
//             Each row of histout is      
//       [norm(grad), f, number of stepsize cuts, iteration count] 
//       costdata = [num f, num grad, num hess] (for levmar, num hess=0)
//
// At this stage all iteration parameters are hardwired in the code.
//
//

[nargout, nargin] = argn();

maxarm = 10;

if (nargin < 4) then
  maxit = 100; 
end

itc = 1;
xc  = x0; 

[fc,gc,jac] = f(xc);
nu0   = .001d0;
mvar  = length(gc);
nvar  = length(xc);
nfun  = 1;
ngrad = 1;
numh  = 0;
numf  = 1;
numg  = 1;
numh  = 0;

ithist(1,1) = norm(gc);
ithist(1,2) = fc;
ithist(1,4) = itc-1;
ithist(1,3) = 0; 
nu          = norm(gc);

while(norm(gc) > tol & itc <= maxit)
  itc = itc+1;
  hc  = (jac'*jac)+ nu*eye(nvar,nvar);
  dc  = hc\gc;
  xt  = xc-dc; 
  [xp,nup,idid] = trtestlm(f,xc,xt,fc,jac,gc,nu,nu0);
  if (idid > 20) then
    error('too many iterations in TR solve')
  end
  xc = xp;
  nu = nup;
  numf = numf+idid;
  if idid > 1
    [fc,gc,jac] = f(xc);
    numf        = numf+1;
    numg        = numg+1;
  end
  ithist(itc,1) = norm(gc);
  ithist(itc,2) = fc; 
  ithist(itc,4) = itc-1;
  ithist(itc,3) = idid;
end
x = xc;
histout  = ithist(1:itc,:);
costdata = [numf, numg, numh];
endfunction

//
//
//

function [xp,nup,idid] = trtestlm(f,xc,xt,fc,jac,gc,nu,nu0)
mu0    = .1;
mulow  = .25;
muhigh = .75;
mup    = 2.d0;
mdown  = .5d0;

z     = xc; 
nvar  = length(z);
numft = 0;
numgt = 0;
itr   = 0;

while ((z==xc) & (itr <= 20))
  ft   = f(xt); 
  itr  = itr+1;
  s    = xt-xc;
  ared = fc-ft;
  pred = -gc'*s;
  _rat = ared/pred;
  if (_rat < mu0) then
    nu = max(nu*mup,nu0);
    hc = (jac'*jac)+ nu*eye(nvar,nvar);
    dc = hc\gc;
    xt = xc-dc;
  elseif (_rat < mulow) then
    z  = xt;
    nu = max(nu*mup,nu0);
	else
    z = xt;
    if (_rat > muhigh) then
      nu = mdown*nu;
      if (nu < nu0) then
        nu = 0.d0;
      end
    end
  end
end
nup  = nu;
xp   = z;
idid = itr;
endfunction
