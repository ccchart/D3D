function h=colquiver(hline,cdata);
%COLQUIVER Color quiver plot.
%   Replaces the uncolored quiver plot generated by QUIVER or QUIVER3
%   by a colored quiver plot.
%
%   COLQUIVER(HL,C)
%   where HL is the vector of line handles returned by QUIVER or QUIVER3
%   and C is a matrix of the same size as used to generate the quiver plot.
%
%   HP = COLQUIVER(...) returns a vector of patch handles.
%
%   Example
%      [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%      z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
%      contour(x,y,z), hold on
%      colquiver(quiver(x,y,px,py),sqrt(px.^2+py.^2))
%      hold off, axis image
%
%   See also QUIVER, QUIVER3.

%----- LGPL --------------------------------------------------------------------
%                                                                               
%   Copyright (C) 2011-2013 Stichting Deltares.                                     
%                                                                               
%   This library is free software; you can redistribute it and/or                
%   modify it under the terms of the GNU Lesser General Public                   
%   License as published by the Free Software Foundation version 2.1.                         
%                                                                               
%   This library is distributed in the hope that it will be useful,              
%   but WITHOUT ANY WARRANTY; without even the implied warranty of               
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
%   Lesser General Public License for more details.                              
%                                                                               
%   You should have received a copy of the GNU Lesser General Public             
%   License along with this library; if not, see <http://www.gnu.org/licenses/>. 
%                                                                               
%   contact: delft3d.support@deltares.nl                                         
%   Stichting Deltares                                                           
%   P.O. Box 177                                                                 
%   2600 MH Delft, The Netherlands                                               
%                                                                               
%   All indications and logos of, and references to, "Delft3D" and "Deltares"    
%   are registered trademarks of Stichting Deltares, and remain the property of  
%   Stichting Deltares. All rights reserved.                                     
%                                                                               
%-------------------------------------------------------------------------------
%   http://www.deltaressystems.com
%   $HeadURL$
%   $Id$

x=get(hline,'xdata');
y=get(hline,'ydata');
z=get(hline,'zdata');
quiv3=~isempty(z{1});
Nc=prod(size(cdata));
c=zeros(4,Nc);
c(1:3,:)=repmat(cdata(:)',3,1);
c(4,:)=NaN;
if quiv3
    hpatch(2)=patch(repmat(x{2}',1,2),repmat(y{2}',1,2),repmat(z{2}',1,2),repmat(c(:),1,2), ...
        'parent',get(hline(1),'parent'), ...
        'edgecolor','flat','facecolor','none');
    c=c(2:4,:);
    hpatch(1)=patch(repmat(x{1}',1,2),repmat(y{1}',1,2),repmat(z{1}',1,2),repmat(c(:),1,2), ...
        'parent',get(hline(1),'parent'), ...
        'edgecolor','flat','facecolor','none');
else
    if isequal(length(x{2}),Nc)
        cc=c(1,:);
        hpatch(2)=patch(repmat(x{2}',1,2)',repmat(y{2}',1,2)',repmat(cc(:),1,2)', ...
            'parent',get(hline(1),'parent'), ...
            'edgecolor','flat','facecolor','none', ...
            'linestyle','none','marker',get(hline(2),'marker'));
    else
        hpatch(2)=patch(repmat(x{2}',1,2),repmat(y{2}',1,2),repmat(c(:),1,2), ...
            'parent',get(hline(1),'parent'), ...
            'edgecolor','flat','facecolor','none');
    end
    c=c(2:4,:);
    hpatch(1)=patch(repmat(x{1}',1,2),repmat(y{1}',1,2),repmat(c(:),1,2), ...
        'parent',get(hline(1),'parent'), ...
        'edgecolor','flat','facecolor','none');
end
delete(hline);
if nargout>0
    h=hpatch;
end
