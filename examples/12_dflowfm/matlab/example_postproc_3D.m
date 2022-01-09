clear all
close all

addpath('d:/pijl/openearthtools_svn/trunk/matlab\');
try
    dflowfm.readMap;
catch
    oetsettings;
end

DIR = 'd:/pijl/unstruc_svn/cases/3D/harlingen_model_3d/output';
MDUnam = 'har';

FNAM=sprintf('%s/%s_map.nc', DIR,MDUnam);

% read data
G = dflowfm.readNet(FNAM);
D = dflowfm.readMap(FNAM,'vel',true,'sal',true);

% specify polygon for cross section
pol.x = [1.561342e+05 1.561360e+05 1.561377e+05 1.561394e+05 1.561410e+05 1.561426e+05 1.561441e+05 1.561457e+05 1.561473e+05 1.561489e+05 1.561506e+05 1.561524e+05 1.561542e+05 1.561560e+05 1.561580e+05 1.561600e+05 1.561621e+05 1.561644e+05 1.561668e+05 1.561718e+05 1.561745e+05 1.561772e+05 1.561800e+05 1.561828e+05 1.561857e+05 1.561966e+05 1.562073e+05 1.562179e+05 1.562283e+05 1.562384e+05 1.562483e+05 1.562672e+05 1.562764e+05 1.562854e+05 1.562944e+05 1.563123e+05 1.563214e+05 1.563306e+05 1.563400e+05 1.563496e+05 1.563594e+05 1.563692e+05 1.563791e+05 1.563890e+05 1.563990e+05 1.564090e+05 1.564191e+05 1.564290e+05 1.564389e+05 1.564486e+05 1.564580e+05 1.564678e+05 1.564780e+05 1.564880e+05 1.565065e+05 1.565151e+05 1.565233e+05 1.565311e+05 1.565386e+05 1.565458e+05 1.565527e+05 1.565592e+05 1.565660e+05 1.565729e+05 1.565799e+05 1.565874e+05 1.565952e+05 1.566029e+05 1.566108e+05 1.566272e+05 1.566359e+05 1.566452e+05 1.566553e+05 1.566665e+05 1.566784e+05 1.567028e+05 1.567154e+05 1.567287e+05 1.567426e+05 1.567564e+05 1.567703e+05 1.567981e+05 1.568128e+05 1.568281e+05 1.568582e+05 1.568870e+05 1.569017e+05 1.569167e+05 1.569319e+05 1.569470e+05 1.569620e+05 1.569768e+05 1.569915e+05 1.570204e+05 1.570343e+05 1.570480e+05 1.570614e+05 1.570746e+05 1.570875e+05 1.571002e+05 1.571126e+05 1.571250e+05 ];
pol.y = [5.766081e+05 5.766018e+05 5.765957e+05 5.765897e+05 5.765838e+05 5.765780e+05 5.765723e+05 5.765665e+05 5.765607e+05 5.765547e+05 5.765487e+05 5.765425e+05 5.765362e+05 5.765298e+05 5.765231e+05 5.765163e+05 5.765091e+05 5.765017e+05 5.764939e+05 5.764771e+05 5.764681e+05 5.764586e+05 5.764487e+05 5.764385e+05 5.764278e+05 5.764307e+05 5.764333e+05 5.764359e+05 5.764384e+05 5.764407e+05 5.764429e+05 5.764468e+05 5.764486e+05 5.764501e+05 5.764516e+05 5.764540e+05 5.764552e+05 5.764562e+05 5.764573e+05 5.764584e+05 5.764595e+05 5.764607e+05 5.764622e+05 5.764638e+05 5.764657e+05 5.764679e+05 5.764705e+05 5.764735e+05 5.764770e+05 5.764809e+05 5.764853e+05 5.764903e+05 5.764966e+05 5.765037e+05 5.765202e+05 5.765294e+05 5.765391e+05 5.765493e+05 5.765598e+05 5.765706e+05 5.765811e+05 5.765908e+05 5.766005e+05 5.766100e+05 5.766193e+05 5.766286e+05 5.766383e+05 5.766477e+05 5.766569e+05 5.766743e+05 5.766826e+05 5.766906e+05 5.766985e+05 5.767064e+05 5.767146e+05 5.767300e+05 5.767370e+05 5.767439e+05 5.767509e+05 5.767578e+05 5.767647e+05 5.767779e+05 5.767844e+05 5.767914e+05 5.768055e+05 5.768183e+05 5.768245e+05 5.768306e+05 5.768369e+05 5.768431e+05 5.768493e+05 5.768554e+05 5.768616e+05 5.768737e+05 5.768798e+05 5.768858e+05 5.768917e+05 5.768974e+05 5.769031e+05 5.769087e+05 5.769142e+05 5.769195e+05 ];

% get polygon
% pol=dflowfm.get_xy_from_figure(1)

% interpolate in cross section
crs = dflowfm.interpolate_on_polygon(G,D,pol);


% plot top layer salt
figure(1)
dflowfm.plotMap(G,D.cen.sal(:,end));
hold on
plot(crs.x,crs.y,'-w')
axis([min(pol.x), max(pol.x), min(pol.y), max(pol.y)]);

% plot contours in cross section
figure(2)
contourf(crs.x,crs.z,crs.cen.sal);
colorbar;

figure(3)
contourf(crs.x,crs.z,sqrt(crs.cen.u.^2+crs.cen.v.^2));
colorbar;