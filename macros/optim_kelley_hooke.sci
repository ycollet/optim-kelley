function [x, histout] = optim_kelley_hooke(x0, f, budget, scales, tol, v)
//
// Nelder-Mead optimizer, No tie-breaking rule other than MATLAB's sort
//
// C. T. Kelley, July 10, 1998
//
//
// This code comes with no guarantee or warranty of any kind.
//
// function [x, histout] = optim_kelley_hooke(x0, f, budget, scales, tol, v)
//
// inputs:
//
//       x0= initial iterate 
//       f = objective function 
//       budget = max f evals (default=50*number of variables)
//                 The iteration will terminate after the iteration that
//                 exhausts the budget
//       scales = the decreasing sequence of stencil sizes
//                 This is an optional argument and the default is
//                 1, 1/2, ... 1/128
//       tol = termination tolerance 
//            difference is function values at successive iterations
//            (default = 1.d-6)
//       v = matrix of search directions (default $v = I$)
//       target = value of f at which the iteration should be terminated
//             This is an optional argument, which you SHOULD set to
//             something reasonable for your problem. The default is
//             as close to no limit as we can get, -1.d8
//
//
// outputs:
//
//        final results = x
//        iteration histor = histout itout x 5
//        histout = iteration history, updated after each nonlinear iteration
//                 = lhist x 2 array, the rows are
//                   [fcount, fval, sval]
//                   fcount = cumulative function evals
//                   fval = current best function value
//                   sval = current scale
//
// set debug = 1 to print iteration stats
//

[nargout, nargin] = argn();

_debug = 0;

//
//

n      = length(x0);
fcount = 0;
bud    = 50*n;
scal   = -(0:8)';
scal   = 2.^scal;
_dir   = eye(n,n);
tolf   = 1.d-6;
targ   = -1.d8;

//
// cache control paramters
//

global cache_size cache cache_ptr cache_fvals

cache_size  = 4*n;
cache       = zeros(n,cache_size);
cache_ptr   = 0;
cache_fvals = zeros(cache_size,1);

//
// count the incoming arguments
//

if (nargin > 2) then bud=budget;  end
if (nargin > 3) then scal=scales; end
if (nargin > 4) then tolf=tol;    end
if (nargin > 5) then targ=target; end
if (nargin > 6) then targ=target; end
if (nargin > 7) then _dir = v;    end

//
// main loop to sweep the scales
//

fcount  = 1;
h       = scal(1);
fv      = geval(f,x0); 
histout = [fcount,fv,h];
fbest   = fv;
x       = x0;
nscal   = length(scal);
ns      = 0; 

while (ns < nscal & fcount <= bud)
  ns = ns+1; 
  [x, sf, lhist, fcount] = hjsearch(x, f, scal(ns), _dir, fcount,bud);
  histout = [histout',lhist']'; 
end
endfunction

//
// Search with a single scale
//

function [x, sf, lhist, finc] = hjsearch(xb, f, h, _dir, fcount,bud)
x  = xb;
xc = x;
sf = 0;
finc  = fcount;
lhist = [];

[x, fv, sf, numf] = hjexplore(xb, xc, f, h, _dir);
finc = finc+numf;

if (sf==1) then
  thist = [finc, fv, h];
  lhist = [lhist',thist']';
end

while ((sf==1) & (finc < bud))
  //
  // pattern move
  //
  d  = x-xb;
  xb = x;
  xc = x+d;
  fb = fv;
  
  [x, fv, sf, numf] = hjexplore(xb, xc, f, h, _dir,fb);
  finc = finc+numf; 
  
  if (sf == 0) then // pattern move fails!
    [x, fv, sf, numf] = hjexplore(xb, xb, f, h, _dir, fb);
    finc = finc+numf; 
  end
  if (sf==1) then
    thist = [finc, fv, h];
    lhist = [lhist',thist']';
  end
end
endfunction

//
// exploratory move
//

function [x, fv, sf, numf] = hjexplore(xb, xc, f, h, _dir, fbold);
global cache_size cache cache_ptr cache_fvals

[nargout, nargin] = argn();

n    = length(xb);
x    = xb;
numf = 0;

if (nargin == 5) then
  [fb,ctr] = geval(f,x); 
  numf     = numf+ctr; 
else
  fb = fbold; 
end

fv    = fb; 
xt    = xc;
sf    = 0;
dirh  = h*_dir;
fvold = fv;

for k=1:n
  p    = xt+dirh(:,k);
  ft   = f(p);
  numf = numf+1;
  if (ft >= fb) then
    p        = xt-dirh(:,k);
    [ft,ctr] = geval(f,p);
    numf     = numf+ctr;
  end
  if (ft < fb) then
    sf = 1;
    xt = p;
    fb = ft;
  end
end

if (sf==1) then
  x  = xt;
  fv = fb;
end
endfunction

//
//

function [fs,ctr] = geval(fh,xh)
global cache_size cache cache_ptr cache_fvals

for i=1:cache_size
  nz(i) = norm(xh-cache(:,i));
end

[vz,iz] = min(nz);

if ((vz == 0) & (cache_ptr ~=0)) then
  fs  = cache_fvals(iz);
  ctr = 0;
else
  fs  = fh(xh);
  ctr = 1;
  cache_ptr              = modulo(cache_ptr,cache_size)+1;
  cache(:,cache_ptr)     = xh;
  cache_fvals(cache_ptr) = fs;
end

