// TestKelley.sce

Title     = 'rosenbrock';
execstr('functoopt = ' + Title);
execstr('get_min   = get_min_bound_' + Title);
execstr('get_max   = get_max_bound_' + Title);
MaxIt     = 1000;
Tol       = 1e-6;
GradTol   = 1e-12;

Min = get_min();
Max = get_max();

lines(0);

//
// Display the cost function
//

clf();

x1 = Min(1):(Max(1)-Min(1))/100:Max(1);
x2 = Min(2):(Max(2)-Min(2))/100:Max(2);
for i=1:length(x1)
  for j=1:length(x2)
    z(i,j) = functoopt([x1(i) x2(j)]);
  end
end
xset('fpf',' ');
contour(x1, x2, z, 20);

function [f, g] = f_optim(x)
f = functoopt(x);
g = derivative(functoopt, x); g = g'; // ??
endfunction

function [f, g, jac] = f_optim_2(x)
f = functoopt(x);
[g, jac] = derivative(functoopt, x); g = g'; jac = jac';// ??
endfunction

x0 = (Max - Min).*rand(2,1) + Min; //x0 = x0';

xtitle(Title + ' function - optim_kelley_bfgswopt','x1','x2');

printf('optim_kelley_bfgswopt:\n');
printf('steepest descent/bfgs with polynomial line search\n');
printf(' and Steve Wright storage of H^-1\n');

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_bfgswopt(x0, f_optim, Tol, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_bfgswopt: Starting point x0:'); disp(x0);
printf('optim_kelley_bfgswopt: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_bfgswopt: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - cgtrust','x1','x2');

printf('optim_kelley_cgtrust:\n');
printf('Steihaug Newton-CG-Trust region algorithm\n')

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_cgtrust(x0, f_optim, [Tol, 0.1, MaxIt, ceil(0.1*MaxIt)], GradTol);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_cgtrust: Starting point x0:'); disp(x0);
printf('optim_kelley_cgtrust: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_cgtrust: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_gaussn','x1','x2');

printf('optim_kelley_gaussn:\n');
printf('Damped Gauss-Newton with Armijo rule\n');
printf('simple divide by 2 stepsize reduction\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_gaussn(x0, f_optim_2, Tol, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_gaussn: Starting point x0:'); disp(x0);
printf('optim_kelley_gaussn: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_gaussn: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_gradproj','x1','x2');

printf('optim_kelley_gradproj:\n');
printf('gradient projection with Armijo rule, simple linesearch\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_gradproj(x0, f_optim, Max, Min, Tol, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_gradproj: Starting point x0:'); disp(x0);
printf('optim_kelley_gradproj: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_gradproj: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_hooke','x1','x2');

printf('optim_kelley_hooke:\n');
printf('Nelder-Mead optimizer, No tie-breaking rule other than MATLAB''s sort\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout] = optim_kelley_hooke(x0, f_optim, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_hooke: Starting point x0:'); disp(x0);
printf('optim_kelley_hooke: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_hooke: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_imfil','x1','x2');

printf('optim_kelley_imfil:\n');
printf('Unconstrained implicit filtering code\n');
printf('IMPLICIT FILTERING with SR1 and BFGS quasi-Newton methods\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_imfil(x0, f_optim, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_imfil: Starting point x0:'); disp(x0);
printf('optim_kelley_imfil: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_imfil: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_levmar','x1','x2');

printf('optim_kelley_levmar:\n');
printf('Levenberg-Marquardt code, trust region control of LM parameter\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_levmar(x0, f_optim_2, Tol, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_levmar: Starting point x0:'); disp(x0);
printf('optim_kelley_levmar: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_levmar: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_nelder','x1','x2');

printf('optim_kelley_nelder:\n');
printf('Nelder-Mead optimizer, No tie-breaking rule other than MATLAB''s sort\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

x0_nelder = [];
x0_nelder(:,1) = (Max - Min).*rand(2,1) + Min;
x0_nelder(:,2) = (Max - Min).*rand(2,1) + Min;
x0_nelder(:,3) = (Max - Min).*rand(2,1) + Min;
//x0_nelder = x0_nelder';

plot(x0_nelder(1,:), x0_nelder(2,:), 'ro');

[x_opt, histout, costdata] = optim_kelley_nelder(x0_nelder, f_optim, Tol, MaxIt);

plot(x_opt(1,:), x_opt(2,:), 'go');

printf('optim_kelley_nelder: Starting simpled x0_nelder:'); disp(x0_nelder);
printf('optim_kelley_nelder: Final simplex x_opt:'); disp(x_opt);

xclick();
clf();

xtitle(Title + ' function - optim_kelley_mds','x1','x2');

printf('optim_kelley_mds:\n');
printf('Multidirectional search\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0_nelder(:,1), x0_nelder(:,2), 'ro');

[x_opt, histout, costdata] = optim_kelley_mds(x0_nelder, f_optim, Tol, MaxIt);

plot(x_opt(:,1), x_opt(:,2), 'go');

printf('optim_kelley_mds: Starting simplex x0:'); disp(x0_nelder);
printf('optim_kelley_mds: Final simplex x_opt:'); disp(x_opt);
printf('optim_kelley_mds: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_ntrust','x1','x2');

printf('optim_kelley_ntrust:\n');
printf('Dogleg trust region, Newton model, dense algorithm \n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_ntrust(x0, f_optim, MaxIt, GradTol);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_ntrust: Starting point x0:'); disp(x0);
printf('optim_kelley_ntrust: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_ntrust: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_projbfgs','x1','x2');

printf('optim_kelley_projbfgs:\n');
printf('optim_kelley_projected BFGS with Armijo rule, simple linesearch \n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_projbfgs(x0, f_optim, Max, Min, Tol, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_projbfgs: Starting point x0:'); disp(x0);
printf('optim_kelley_projbfgs: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_projbfgs: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));

xclick();
clf();

xtitle(Title + ' function - optim_kelley_steep','x1','x2');

printf('optim_kelley_steep:\n');
printf('steepest descent with Armijo rule, polynomial linesearch\n');

xset('fpf',' ');
contour(x1, x2, z, 20);

plot(x0(1), x0(2), 'ro');

[x_opt, histout, costdata] = optim_kelley_steep(x0, f_optim, Tol, MaxIt);

plot(x_opt(1), x_opt(2), 'go');

printf('optim_kelley_steep: Starting point x0:'); disp(x0);
printf('optim_kelley_steep: Final point x_opt:'); disp(x_opt);
printf('optim_kelley_steep: Initial fobj = %f - Final fobj = %f\n', functoopt(x0), functoopt(x_opt));


