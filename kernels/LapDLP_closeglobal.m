function [u ux uy info] = LapDLP_closeglobal(t, s, tau, side)
% LAPDLP_CLOSEGLOBAL - Laplace DLP potential & deriv w/ global close-eval quad
%
% u = LapDLP_closeglobal(t,s,dens,side) returns potentials at targets t.x
%  due to double-layer potential with real-valued density dens sampled on the
%  nodes s.x of a smooth global quadrature rule on the curve s, either inside
%  outside the curve. The global scheme of [hel08] based on barycentric Cauchy
%  close-evaluation (see references in Cau_closeglobal.m) is used.
% 
% Our definition of the DLP on curve \Gamma is, associating R^2 with the complex
% plane, with x = (x_1,x_2) a point in the complex plane,
%
%    u(x) = Re (1/2\pi i) \int_\Gamma \tau(y) / (x-y) dy
%
% Inputs:
%  t = target struct containing
%      t.x = M-by-1 list of targets, as points in complex plane
%      (optionally also t.nx target normals as unit complex numbers)
%  s = curve struct containing N-by-1 vector s.x of source nodes (as complex
%     numbers), and other fields as generated by setupquad.
%  dens = double-layer density values at nodes. Note, must be real-valued to
%     evaluate Laplace layer potentials (note Stokes feeds it complex dens).
%     If dens has multiple columns, the evaluation is done for each.
%     If dens is empty, outputs are matrices mapping density to values, etc.
%  side = 'i','e' to indicate targets are interior or exterior to the curve.
%
% Outputs:
%  u = column vector of M potential values where M = numel(t.x), or if
%      dens has n columns, the n column vector outputs are stacked (M-by-n).
%  [u un] = LapDLPeval_closeglobal(t,s,dens,side) also returns target-normal
%           derivatives (M-by-n) using t.nx
%  [u ux uy] = LapDLPeval_closeglobal(t,s,dens,side) instead returns x- and y-
%           partials of u at the targets (ignores t.nx)
%  [u ux uy info] = LapDLPeval_closeglobal(t,s,dens,side) also gives diagnostic:
%           info.vb = vector of v boundary values (M-by-n)
%           info.imv = imag part of v at targets (M-by-n)
%
% References:
%
%  [hel08] J. Helsing and R. Ojala, On the evaluation of layer potentials close
%          to their sources, J. Comput. Phys., 227 (2008), pp. 2899–292
%
%  [lsc2d] Spectrally-accurate quadratures for evaluation of layer potentials
%          close to the boundary for the 2D Stokes and Laplace equations,
%          A. H. Barnett, B. Wu, and S. Veerapaneni, SIAM J. Sci. Comput.,
%          37(4), B519-B542 (2015)   https://arxiv.org/abs/1410.2187
%
% See also: LAPDLP, SETUPQUAD, CAU_CLOSEGLOBAL, LAPINTDIRBVP
%
% complexity O(N^2) for evaluation of v^+ or v^-, plus O(NM) for globally
% compensated quadrature to targets

% Barnett 2013; multiple-column generalization by Gary Marple, 2014.
% Repackaged Barnett 6/12/16, 6/27/16 col vec outputs

if nargin==0, test_LapDLP_closeglobal; return; end
N = numel(s.x); M = numel(t.x);      % # source, target nodes
if isempty(tau), tau = eye(N); end   % case of filling matrices
n = size(tau,2);                     % # density columns

% Helsing step 1: eval bdry limits at nodes of v = complex DLP(tau)...
% (note to future fortran/C versions: this also needs to be able to handle
% complex input, real output, since the Stokes DLP feeds that in)
vb = zeros(N,n);                    % alloc v^+ or v^- bdry data of holom func
taup = zeros(N,n);                 % alloc tau'
for k=1:n
    taup(:,k) = perispecdiff(tau(:,k));  % numerical deriv of each dens
end
for i=1:N, j = [1:i-1, i+1:N];      % skip pt i. Eqn (4.2) in [lsc2d]
  vb(i,:) = sum((tau(j,:)-ones(N-1,1)*tau(i,:))./((s.x(j)-s.x(i))*ones(1,n)).*(s.cw(j)*ones(1,n)),1) + taup(i,:)*s.w(i)/s.sp(i);   % vectorized over cols
end
vb = vb*(1/(-2i*pi));               % prefactor
if side=='i', vb = vb - tau; end   % JR's add for v^-, cancel for v^+ (Eqn 4.3)
info.vb = vb;                       % diagnostics

% Helsing step 2: compensated close-evaluation of v & v', followed by take Re:
if nargout>1                                % want derivatives
  [v vp] = Cau_closeglobal(t.x,s,vb,side);  % does Sec. 3 of [lsc2d]
  ux = real(vp); uy = -imag(vp);            % col vecs, leave as partials...
  if nargout==2           % or dot w/ targ nor...
    ux = bsxfun(@times,ux,real(t.nx)) + bsxfun(@times,uy,imag(t.nx));
  end
else
  v = Cau_closeglobal(t.x,s,vb,side);        % does Sec. 3 of [lsc2d]
end
u = real(v); info.imv = imag(v);

%%%%%%%%%%%%%%%%%%%
function test_LapDLP_closeglobal  % check far-field matches the native rule
verb = 0;       % to visualize
s = wobblycurve(0.3,5,200); if verb, figure;showsegment(s); end
tau = 0.7+sin(3*s.t);                    % pick smooth density w/ nonzero mean
nt = 100; t.nx = exp(2i*pi*rand(nt,1));  % target normals
for side = 'ie'
  if side=='e', t.x = 1.5+1i+rand(nt,1)+1i*rand(nt,1);         % distant targs
  else, t.x = 0.6*(rand(nt,1)+1i*rand(nt,1)-(0.5+0.5i)); end % targs far inside
  if verb, plot(t.x,'.'); end
  [u un] = LapDLP(t,s,tau);
  [uc unc] = LapDLP_closeglobal(t,s,tau,side);
  [uc uxc uyc] = LapDLP_closeglobal(t,s,tau,side);
  fprintf('density case, side=%s: max abs errors in u, un, and [ux,uy]...\n',side)
  disp([max(abs(u-uc)), max(abs(un-unc)), max(abs(un - (uxc.*real(t.nx)+uyc.*imag(t.nx))))])
  [u un] = LapDLP(t,s);   % matrix cases....
  [uc unc] = LapDLP_closeglobal(t,s,[],side);
  [uc uxc uyc] = LapDLP_closeglobal(t,s,[],side);
  fprintf('matrix fill case, side=%s: max abs errors in u, un, and [ux,uy]...\n',side)
  disp([max(abs(u(:)-uc(:))), max(abs(un(:)-unc(:))), max(max(abs(un - (bsxfun(@times,uxc,real(t.nx))+bsxfun(@times,uyc,imag(t.nx))))))])
end

