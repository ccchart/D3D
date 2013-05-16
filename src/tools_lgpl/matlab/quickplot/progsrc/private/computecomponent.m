function [data,scalar,vpt]=computecomponent(data,Ops)
%COMPUTECOMPONENT Compute component of vector data set.
%
%   NewData = COMPUTECOMPONENT(Data,Component)
%   where Data is a vector data structure obtained from QPREAD and
%   Component equals one of the following strings: 'magnitude',
%   'magnitude in plane', 'angle', 'x component', 'y component',
%   'z component', 'm component', 'n component', 'k component'

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

if ischar(Ops)
    vpt=Ops;
    scalar=0;
else
    vpt=Ops.vectorcomponent;
    switch lower(vpt)
        case {'edge'}
            N=length(data);
            data=cat(2,data,data);
            for dU = 1:N
                dV = dU+N;
                if ~isempty(Ops.vectorcolour)
                    data(dU).Val=abs(data(dU).XComp(:,2:end));
                    data(dV).Val=abs(data(dV).YComp(2:end,:));
                end
                %
                dX = data(dU).X(:,2:end)-data(dU).X(:,1:end-1);
                dY = data(dU).Y(:,2:end)-data(dU).Y(:,1:end-1);
                Mg = sqrt(dX.^2+dY.^2);
                data(dU).X=(data(dU).X(:,1:end-1)+data(dU).X(:,2:end))/2;
                data(dU).Y=(data(dU).Y(:,1:end-1)+data(dU).Y(:,2:end))/2;
                data(dU).YComp=-data(dU).XComp(:,2:end).*dX./Mg;
                data(dU).XComp=data(dU).XComp(:,2:end).*dY./Mg;
                %
                dX = data(dV).X(2:end,:)-data(dV).X(1:end-1,:);
                dY = data(dV).Y(2:end,:)-data(dV).Y(1:end-1,:);
                Mg = sqrt(dX.^2+dY.^2);
                data(dV).X=(data(dV).X(1:end-1,:)+data(dV).X(2:end,:))/2;
                data(dV).Y=(data(dV).Y(1:end-1,:)+data(dV).Y(2:end,:))/2;
                data(dV).XComp=-data(dV).YComp(2:end,:).*dY./Mg;
                data(dV).YComp=data(dV).YComp(2:end,:).*dX./Mg;
            end
            scalar=0;
            vpt='vector';
            return
        case {'vector (split x,y)'}
            N=length(data);
            data=cat(2,data,data);
            for d = 1:N
                data(d).YComp=zeros(size(data(d).YComp));
                data(N+d).XComp=data(d).YComp;
                if ~isempty(Ops.vectorcolour)
                    data(d).Val=abs(data(d).XComp);
                    data(N+d).Val=abs(data(N+d).YComp);
                end
            end
            scalar=0;
            vpt='vector';
            return
        case {'vector (split m,n)'}
            N=length(data);
            data=cat(2,data,data);
            for d = 1:N
                data(d).XComp=real(data(d).XComp);
                data(d).YComp=real(data(d).YComp);
                data(N+d).XComp=imag(data(N+d).XComp);
                data(N+d).YComp=imag(data(N+d).YComp);
                if ~isempty(Ops.vectorcolour)
                    data(d).Val=sqrt(data(d).XComp.^2+data(d).YComp.^2);
                    data(N+d).Val=abs(data(N+d).YComp.^2+data(N+d).YComp.^2);
                end
            end
            scalar=0;
            vpt='vector';
            return
        case {'vector','patch centred vector'}
            scalar=0;
            if ~isempty(Ops.vectorcolour)
                vpt=Ops.vectorcolour;
            else
                return
            end
        otherwise
            scalar=1;
    end
end
for d=1:length(data)
    switch lower(vpt)
        case {'vector','vector (split x,y)'}
        case 'magnitude'
            data(d).Val=data(d).XComp.^2;
            if isfield(data,'YComp')
                data(d).Val=data(d).Val+data(d).YComp.^2;
            end
            if isfield(data,'ZComp')
                data(d).Val=data(d).Val+data(d).ZComp.^2;
            end
            data(d).Val=sqrt(data(d).Val);
        case 'magnitude in plane'
            if isfield(data,'XComp') && size(data(d).XComp,1)>1
                data(d).Val=data(d).XComp.^2;
                if isfield(data,'YComp') && size(data(d).YComp,2)>1
                    data(d).Val=data(d).Val+data(d).YComp.^2;
                end
            elseif isfield(data,'YComp') && size(data(d).YComp,2)>1
                data(d).Val=data(d).YComp.^2;
            end
            if isfield(data,'ZComp') && size(data(d).ZComp,3)>1
                data(d).Val=data(d).Val+data(d).ZComp.^2;
            end
            data(d).Val=sqrt(data(d).Val);
        case 'angle'
            sf = qp_unitconversion('radians',Ops.units);
            switch Ops.angleconvention
                case {'Nautical','Nautical Positive'}
                    data(d).Val=sf*atan2(data(d).XComp,data(d).YComp); % Nautical convention
                otherwise
                    data(d).Val=sf*atan2(data(d).YComp,data(d).XComp); % Cartesian convention
            end
            if strfind(Ops.angleconvention,'Positive')
                neg = data(d).Val<0;
                data(d).Val(neg) = data(d).Val(neg)+2*pi*sf;
            end
            data(d).Units=Ops.units;
            vpt='angle';
        case 'x component'
            data(d).Val=data(d).XComp;
        case 'y component'
            data(d).Val=data(d).YComp;
        case 'z component'
            data(d).Val=data(d).ZComp;
        case 'm component'
            data(d).Val=data(d).XComp; %MComp
        case 'n component'
            data(d).Val=data(d).YComp; %NComp
        case 'k component'
            data(d).Val=data(d).ZComp; %KComp
        case 'normal component' % only for a vertical slice
            if  size(data(d).Val,1)>1
                data(d).Val=data(d).YComp; %NComp
            else
                data(d).Val=data(d).XComp; %MComp
            end
        otherwise
            ui_message('error','Unexpected colour/plot type encountered: %s.',vpt);
            scalar=0;
    end
    if isfield(data,'Val')
        if isfield(data,'X') && isequal(size(data(d).Val),size(data(d).X))
            data(d).Val(isnan(data(d).X))=NaN;
        end
        if isfield(data,'Y') && isequal(size(data(d).Val),size(data(d).Y))
            data(d).Val(isnan(data(d).Y))=NaN;
        end
        if isfield(data,'Z') && isequal(size(data(d).Val),size(data(d).Z))
            data(d).Val(isnan(data(d).Z))=NaN;
        end
    end
end
