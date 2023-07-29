function A = HelmDLP_closepanel(omega,t,s,a,b,side, meth)
% special quadrature for 2d helmholtz
% via kernel split Laplace DLP special quadrature

p = numel(s.x); % GL nodes order
if isfield('be',s), be = s.be; else, be = 1; end  % upsample factor of GL nodes
if be > 1
  [x, w] = gauss(p); w = w';
  [xf, wf] = gauss(ceil(be*p)); wf = wf';
  Imn = interpmat(x, xf);
  hl = mean(s.w./(w.*s.sp));
  s.w = wf; % overwrite source panel
  s.x = Imn*s.x; s.xp = Imn*s.xp;
  s.sp = abs(s.xp); s.tang = s.xp./s.sp; s.nx = -1i*s.tang;    % outward unit normals
  s.wxp = hl*s.w.*s.xp; % complex speed weights (Helsing's wzp)
  s.w = abs(s.wxp); % speed weights
end

% Variants of an explicit kernel-split panel-based Nystrom discretization scheme for Helmholtz boundary value problems.
% Johan Helsing & Anders Holst. p 696, eq(26), pay attention to possible typo in laplace DLP coefficient 1/(2*pi) instead of 1/pi
lapslpmat = -2*pi*LapSLP_closepanel(t,s,a,b,side, meth); % close eval for log kernel (Laplace SLP)
lapdlpmat = -2*pi*LapDLP_closepanel(t,s,a,b,side, meth); % close eval for Cauchy kernel (Laplace DLP)
u = 1i/4*(omega*abs(s.x(:).'-t.x(:))).*besselh(1,omega*abs(s.x(:).'-t.x(:))).*imag(s.wxp(:).'./(t.x(:)-s.x(:).'));
A = ( u + 2/pi*log(abs(s.x(:).'-t.x(:))).*imag(u) + 1/(2*pi)*real(s.nx(:).'./(s.x(:).'-t.x(:))).*s.w' ) ...
         - 2/pi*(lapslpmat./s.w').*imag(u) - 1/(2*pi)*real(lapdlpmat); 
if be > 1, A = A*Imn; end

end