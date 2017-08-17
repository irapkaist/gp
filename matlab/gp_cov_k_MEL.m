function Kij = gp_cov_k_MEL (Xi, Xj, l, s, cfg)
%
%-----------------------------------------------------------------
%    History:
%    Date            Who         What
%    -----------     -------     -----------------------------
%    2013.12.26      ak          created and written

M = diag(1./[cfg.slam.along_dx,cfg.slam.cross_dx,1.0].^2);
d = sqrt((Xi-Xj)'*M*(Xi-Xj));

%d = norm (Xi-Xj);

if (d < l) % if distance is smaller than l
    term1 = ((2+cos(2*pi*d/l))/3)*(1-d/l);
    term2 = 1/(2*pi)*sin(2*pi*d/l);
    Kij = s*(term1 + term2);
else
    Kij = 0;
end
