function [xc,histout,costdata] = optim_kelley_cgtrust(x0,f,parms,resolution)
//
//
// C. T. Kelley, Jan 13, 1998
//
// This code comes with no guarantee or warranty of any kind.
//
// function [x,histout] = optim_kelley_cgtrust(x0,f,parms,resolution)
//
// Steihaug Newton-CG-Trust region algorithm
// 
// Input: x0 = initial iterate
//        f = objective function,
//            the calling sequence for f should be
//            [fout,gout]=f(x) where fout=f(x) is a scalar
//              and gout = grad f(x) is a COLUMN vector
//        parms = [tol, eta, maxitnl, maxitl]
//        tol = termination criterion norm(grad) < tol
//        eta = forcing term in linear iteration (optional) default = .1
//        maxitnl = maximum nonlinear iterations (optional) default = 100
//        maxitl = maximum linear iterations (optional) default = 20
//        resolution = estimated accuracy in functions/gradients (optional)
//                     default = 1.d-12
//                     The finite difference increment in the difference
//                     Hessian is set to sqrt(resolution).
//
//
//
// Output: x = solution
//         histout = iteration history
//             Each row of histout is
//              [norm(grad), f, TR radius, inner iteration count]
//         costdata = [num f, num grad]
//
// Requires: dirdero.m
//
// trust region algorithm parameters
//

[nargout, nargin] = argn();

omegaup   = 2;
omegadown = .5;
mu0       = .25;
mulow     = .25; 
muhigh    = .75;

//
// set maxit, the resolution, and the difference increment
//

_debug = 0;
tol    = parms(1);
np     = length(parms);
itc    = 1;

if (np < 2) then
  eta = .1;
else
  eta = parms(2);
end

if (np < 3) then
  maxit = 100;
else
  maxit = parms(3);
end
if (np < 4) then
  maxitl = 20;
else
  maxitl = parms(4);
end
if (nargin < 5) then
  resolution = 1.d-12;
end

numf  = 0;
numg  = 0;
hdiff = sqrt(resolution);

t   = 100;
itc = 0;
xc  = x0;
n   = length(x0);
[fc,gc] = f(xc);
numf   = 1;
numg   = 1; 
trrad  = norm(x0);
ithist = zeros(maxit,4);
ithist(itc+1,:) = [norm(gc), fc, trrad, 0];

//

paramstr = [eta,maxitl];
trcount  = 1;

while ((norm(gc) > tol & itc <=maxit) & trcount < 30)
  itc = itc+1;
  [s,dirs] = trcg(xc, f, gc, trrad, paramstr);
  vd   = size(dirs);
  itsl = vd(2);
  numf = numf+itsl;
  numg = numg+itsl;
  xt   = xc+s;
  ft   = f(xt);
  numf = numf+1;
  ared = ft-fc;
  w    = dirdero(xc, s, f, gc);
  numg = numg+1;
  pred = gc'*s + .5*(w'*s);
  _rat  = ared/pred;
  if (_rat > mu0) then
    xc = xc+s;
    [fc,gc] = f(xc);
    numf = numf+1;
    if (_rat > muhigh) & (norm(s) > trrad-1.d-8) then
      trrad = omegaup*trrad;
    end
    if (_rat < mulow) then
      trrad = omegadown*trrad;
    end
  else
    for k=1:itsl
      dsums(k) = norm(sum(dirs(:,1:k),2));
    end
    trcount = 1;
    while ((_rat <= mu0) & (trcount < 30))
      trrad = omegadown*min(trrad,norm(s));
      [s, kout] = tradj(trrad, dirs, itsl);
      xt   = xc+s;
      ft   = f(xt);
      numf = numf+1;
      ared = ft-fc; 
      //
      //      Only compute a new pred if ared < 0 and there's some hope
      //      that rat > mu0
      //
      if (ared < 0) then
        w    = dirdero(xc, s, f, gc);
        numg = numg+1;
        pred = gc'*s + .5*(w'*s); 
      end
      _rat     = ared/pred;
      itsl    = kout;
      trcount = trcount+1;
      if (trcount > 30) then
        //           ithist(itc+1,:)=[norm(gc), fc, trrad, itsl];
        //           histout=ithist(1:itc+1,:); costdata=[numf,numg];
        disp(' stagnation in CGTR')
      end
    end
    if (trcount < 30) then
      xc      = xt;
      [fc,gc] = f(xc);
      numf    = numf+1;
      numg    = numg+1;
    end
  end
  ithist(itc+1,:) = [norm(gc), fc, trrad, itsl];
end

histout  = ithist(1:itc+1,:);
costdata = [numf,numg];
endfunction

//
// find the point of intersetion of the TR boundary and the PL path
//

function [st, kout] = tradj(trrad, dirs, itsl)
st     = dirs(:,1);
inside = 1;
if (norm(st) > trrad) | (itsl == 1) then
  st   = st*trrad/norm(st);
  kout = 1;
else
  for k=2:itsl
    if ((norm(st+dirs(:,k)) > trrad) & (inside == 1)) then
      kout   = k;
      p      = dirs(:,k);
      ac     = p'*p;
      bc     = 2*(st'*p);
      cc     = st'*st - trrad*trrad;
      alpha  = (-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac);
      st     = st+alpha*p;
      inside = 0;
    else
      st = st+dirs(:,k);
    end
  end
end
endfunction

//
//
//

function [x, directions]  = trcg(xc, f, gc, delta, params, pcv)
//
// Solve the trust region problem with preconditioned conjugate-gradient
//
// C. T. Kelley, January 13, 1997
//
// This code comes with no guarantee or warranty of any kind.
// function [x, directions]
//                    = trcg(xc, f, gc, delta)
//
//
//
// Input:        xc=current point
//               b=right hand side
//           f = function, the calling sequence is
//                               [fun,grad]=f(x)
//           gc = current gradient
//                gc has usually been computed
//                before the call to dirdero
//           delta = TR radius
//           params = two dimensional vector to control iteration
//                params(1) = relative residual reduction factor
//                params(2) = max number of iterations
//           pcv, a routine to apply the preconditioner
//                if omitted, the identity is used.
//                The format for pcv is 
//                       function px = pcv(x).
//
// Output:   x = trial step
//           directions = array of search directions TR radius reduction
// 
//

//
// initialization
//

[nargout, nargin] = argn();

n = length(xc);
errtol   = params(1);
maxiters = params(2); 
x = zeros(n,1);
b = -gc;
r = b - dirdero(xc, x, f, gc);

if (nargin == 5) then
  z = r;
else
  z = pcv(r);
end

rho        = z'*r;
tst        = norm(r);
terminate  = errtol*norm(b);
it         = 1;
directions = zeros(n,1);
hatdel     = delta*(1-1.d-6);

while((tst > terminate) & (it <= maxiters) & norm(x) <= hatdel)
  if (it==1) then
    p = z;
  else
    _beta = rho/rhoold;
    p = z + _beta*p;
  end
  w     = dirdero(xc, p, f, gc);
  alpha = p'*w;

  //
  // If alpha <=0 head to the TR boundary and return
  //

  ineg = 0;
  if (alpha <= 0) then
    ac    = p'*p;
    bc    = 2*(x'*p);
    cc    = x'*x - delta*delta;
    alpha = (-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac);
    disp(' negative curvature')
  else
    alpha = rho/alpha;
    if (norm(x+alpha*p) > delta) then
      ac    = p'*p;
      bc    = 2*(x'*p);
      cc    = x'*x - delta*delta;
      alpha = (-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac);
    end
  end

  x = x+alpha*p;
  directions(:,it) = alpha*p;
  r      = r - alpha*w;
  tst    = norm(r);
  rhoold = rho;

  if (nargin < 6) then
    z=r;
  else
    z = pcv(r);
  end

  rho = z'*r;
  it  = it+1;
end
endfunction
