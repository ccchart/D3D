function delwaq2raster(ini_file)
%DELWAQ2RASTER rasterize delwaq MAP results
%   DELWAQ2RASTER(INIFILE) reads the INIFILE and converts selected data
%   from the delwaq MAP results file referenced therein to raster grid
%   files for use with for instance the HABITAT module. The INIFILE
%   should contain the following elements:
%
%   Block [general] containing keywords:
%     raster : name of reference raster file (determining raster dimensions
%              and raster file format). Currently, supports bil/hdr pairs
%              and asc files. Output files will use the same settings.
%     data   : name of delwaq MAP file
%     grid   : reference to Delft3D curvilinear cco/lga files or flexible
%              mesh netCDF file.
%     method : data mapping method, one of:
%              center   : raster value determined by simulation cell in
%                         which the raster cell center is located
%              maxarea  : raster value determined by biggest overlapping
%                         (aggregated) simulation cell
%              weighted : raster value determined by area weighted average
%                         of overlapping simulationn cells; the total
%                         overlapping area fraction must be at least equal
%                         to minarea.
%     minarea: indicates for method 'weighted' the minimum overlap area to
%              determine the average (in case of smaller overlap, the
%              raster cell is assigned the nodata value). Default: 0.
%     outdir : output directory for the raster files.
%     tstart : start index of time steps to include in processing (by
%              default 1, that is the first time step).
%     tstop  : stop index of time steps to include in processing (by
%              default the last time step in the file).
%
%     One or more [action] blocks containing the following keywords:
%     time_op: name of operation in time. The tool supports the following
%              operations:
%              max, mean, min, std, and ident (default)
%              The default time operator "ident" will convert every time
%              step in the specified period (from tstart until tstop). In
%              the other cases, it will take the max, mean, min, std of all
%              values in the selected time period.
%     tstart : start index of time steps to include in processing (default
%              value taken from [general] block).
%     tstop  : stop index of time steps to include in processing (default
%              value taken from [general] block).
%     include: name of a variable to be converted to raster (the include
%              keyword may be repeated to convert multiple variables using
%              the same settings). The names specified should either match
%              the short DELWAQ name, or the expanded names as used by
%              Delft3D-QUICKPLOT.
%
%   See also: DELWAQ, ARCGRID.

%----- LGPL --------------------------------------------------------------------
%
%   Copyright (C) 2011-2015 Stichting Deltares.
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

if isstandalone
   fprintf(1,'--------------------------------------------------------------------------------\n')
   fprintf(1,'Delft3D-DELWAQ2RASTER conversion tool.\n')
   fprintf(1,'Version <VERSION> (<CREATIONDATE>)\n')
   fprintf(1,'--------------------------------------------------------------------------------\n')
end

if nargin==0
    fprintf(1,'Syntax: delwaq2raster <configuration.ini>\n');
    error('No configuration.ini file specified')
end
ini_info       = inifile('open',ini_file);
%
raster_reffile = inifile('get',ini_info,'general','raster');
gridfile       = inifile('get',ini_info,'general','grid');
waqfile        = inifile('get',ini_info,'general','data');
flwfile        = gridfile;
method         = inifile('get',ini_info,'general','method');
outdir         = inifile('get',ini_info,'general','outdir');
%
raster_info    = qpfopen(raster_reffile);
for f = {'qp_filetype','QP_FileType'}
    fs = f{1};
    if isfield(raster_info,fs)
        raster_info    = rmfield(raster_info,fs);
    end
end
%
switch raster_info.FileType
    case 'HDR Raster File - BIP/BIL/BSQ'
        raster_xco     = bil('read',raster_info,'xco');
        raster_xcc     = bil('read',raster_info,'xcc');
        raster_yco     = bil('read',raster_info,'yco');
        raster_ycc     = bil('read',raster_info,'ycc');
        raster_sz      = [raster_info.nRows raster_info.nCols];
        cellarea       = raster_info.Xdim*raster_info.Ydim;
    case 'arcgrid'
        tmp = qpread(raster_info,'grid','grid');
        raster_xco     = tmp.X(:,1)';
        raster_xcc     = raster_xco(1:end-1)+raster_info.CellSize(1)/2;
        raster_yco     = tmp.Y(1,:);
        raster_ycc     = raster_yco(1:end-1)-raster_info.CellSize(2)/2;
        raster_sz      = [raster_info.NRows raster_info.NCols];
        cellarea       = prod(raster_info.CellSize);
end
%
yflip          = raster_yco(2)<raster_yco(1);
if yflip
    raster_yco = raster_yco(end:-1:1);
end
%
waq_info       = qpfopen(waqfile,gridfile);
dwq            = waq_info.DwqBin;
waq_qnt        = qpread(waq_info);
waq_qnt(strcmp('-------',{waq_qnt.Name}'))=[];
%
flw_info       = qpfopen(flwfile);
flw_qnt        = qpread(flw_info);
flw_qnt(strcmp('-------',{flw_qnt.Name}'))=[];
nElm = waq_info.Grid.CoordDims{1};
% filter quantities at cell centers
for i = length(flw_qnt):-1:1
    if ~strcmp(flw_qnt(i).DimName{3},nElm)
        flw_qnt(i) = [];
    end
end
%
tstart         = inifile('get',ini_info,'general','tstart',1);
tstop          = inifile('get',ini_info,'general','tstop' ,inf);
%
if exist(outdir,'dir')~=7
    error('Invalid output directory: %s',outdir)
end
%
ndatasets = 0;
chp = inifile('chapters',ini_info);
for ifld = 1:length(chp)
    switch chp{ifld,1}
        case 'general'
        case 'action'
            chp{ifld,1} = inifile('get',ini_info,ifld,'time_op','ident');
            chp{ifld,2} = inifile('get',ini_info,ifld,'include',{});
            [chp{ifld,3},chp{ifld,2}] = filter_qnt(chp{ifld,2},waq_qnt,flw_qnt);
            info = [];
            info.tstart = inifile('get',ini_info,ifld,'tstart',tstart);
            info.tstop  = inifile('get',ini_info,ifld,'tstop' ,tstop);
            chp{ifld,4} = info;
            %
            if strcmp(chp{ifld,1},'ident')
                ndatasets = ndatasets + length(chp{ifld,2})*(info.tstop-info.tstart+1);
            else
                ndatasets = ndatasets + length(chp{ifld,2});
            end
        otherwise
            error('Unknown chapter encountered: %s',chp{ifld,1})
    end
end
%
qnts = cat(1,chp{:,2});
iqnt = cat(1,chp{:,3});
if any(iqnt==0)
    qnts = unique(qnts(iqnt==0));
    if length(qnts)==1
        err = 'Unable to locate the following quantity in the data file(s)';
    else
        err = 'Unable to locate the following quantities in the data file(s)';
    end
    error('%s\n',err,qnts{:})
end
%
for ifld = 1:length(chp)
    if ~isempty(chp{ifld,3})
        nq = length(chp{ifld,3});
        info = chp{ifld,4};
        %
        AllTimes  = cell(2,0);
        for i = 1:nq
            if chp{ifld,3}(i)>0
                times = delwaq('read',dwq,[],[]);
            else
                times = qpread(flw_info,flw_qnt(-chp{ifld,3}(i)),'times');
            end
            match=0;
            for it = 1:size(AllTimes,2)
                if isequal(AllTimes{1,it},times)
                    % add elem
                    AllTimes{2,it}(end+1,1) = chp{ifld,2}(i);
                    match=1;
                end
            end
            if ~match
                AllTimes{1,end+1} = times;
                AllTimes{2,end} = chp{ifld,2}(i);
            end
        end
        if size(AllTimes,2)>1
            if nq>2
                all_ = 'all ';
            else
                all_ = '';
            end
            err = {sprintf('The output times for the selected variables in block %i don''t %smatch.',ifld,all_)};
            for it = 1:size(AllTimes,2)
                if length(AllTimes{2,it})==1
                    err{end+1,1} = sprintf('%s: %i time step(s)',AllTimes{2,it}{1},length(AllTimes{1,it}));
                else
                    err{end+1,1} = sprintf('%i shared time step(s) for',length(AllTimes{1,it}));
                    for i = 1:length(AllTimes{2,it})
                        err{end+1,1} = sprintf('  %s',AllTimes{2,it}{i});
                    end
                end
            end
            error('%s\n',err{:})
        end
        if ~isfinite(info.tstop)
            info.tstop = length(AllTimes{1});
        elseif info.tstart>=info.tstop
            error('Invalid time range specified: tstart (now %i) should be less than tstop (now %i)',info.tstart,info.tstop)
        elseif info.tstop>length(AllTimes{1})
            error('Time range selected (%i:%i) does not fall inside the time range available (1:%i)',info.tstart,info.tstop,length(AllTimes{1}))
        end
        %
        chp{ifld,4} = info;
    end
end
%
G = qpread(waq_info, 'segment number', 'gridcell');
if size(G.X,2)>1
else
    idx = find(~isnan(G.X));
    id2 = find(diff([-100;idx;-100])~=1);
    %
    ncell = length(id2)-1;
    %
    G.CellX = cell(ncell,1);
    G.CellY = cell(ncell,1);
    G.CellIdx = zeros(ncell,1);
    for m = 1:ncell
        m1 = idx(id2(m));
        m2 = idx(id2(m+1)-1);
        %
        G.CellX{m} = G.X(m1:m2-1)';
        G.CellY{m} = G.Y(m1:m2-1)';
        G.CellIdx(m) = m1;
    end
    G.CellVal = G.Val(G.CellIdx);
end
%
switch method
    case 'center'
        min_area = 0;
        yflip = 0;
        %
        % create mapping
        %
        h = lcwaitbar(0,'Mapping centers ...');
        if size(G.X,2)>1
            wght = cell(size(G.X)-1);
            lidx = cell(size(G.X)-1);
            for m = 1:size(G.X,1)-1
                for n = 1:size(G.X,2)-1
                    x = [G.X(m,n) G.X(m,n+1) G.X(m+1,n+1) G.X(m+1,n)];
                    y = [G.Y(m,n) G.Y(m,n+1) G.Y(m+1,n+1) G.Y(m+1,n)];
                    %
                    if ~any(isnan(x) | isnan(y)) && ~isnan(G.Val(m,n))
                        [wght{m,n},lidx{m,n}] = inside_raster_poly(x,y,raster_xcc,raster_ycc);
                    end
                end
                lcwaitbar(m/(size(G.X,1)-1),h)
            end
        else
            wght = cell(ncell,1);
            lidx = cell(ncell,1);
            for m = 1:ncell
                x = G.CellX{m};
                y = G.CellY{m};
                if ~any(isnan(x) | isnan(y)) && ~isnan(G.CellVal(m))
                    [wght{m},lidx{m}] = inside_raster_poly(x,y,raster_xcc,raster_ycc);
                end
                lcwaitbar(m/ncell,h)
            end
        end
        delete(h)
    case {'weighted','maxarea'}
        min_area = inifile('get',ini_info,'general','minarea',0);
        if ~isnumeric(min_area) || min_area<0 || min_area>1
            error('Invalid value for minarea')
        else
            min_area = min_area*cellarea;
        end
        %
        % create mapping
        %
        h = lcwaitbar(0,'Determining weights ...');
        if size(G.X,2)>1
            wght = cell(size(G.X)-1);
            lidx = cell(size(G.X)-1);
            for m = 1:size(G.X,1)-1
                for n = 1:size(G.X,2)-1
                    x = [G.X(m,n) G.X(m,n+1) G.X(m+1,n+1) G.X(m+1,n)];
                    y = [G.Y(m,n) G.Y(m,n+1) G.Y(m+1,n+1) G.Y(m+1,n)];
                    %
                    if ~any(isnan(x) | isnan(y)) && ~isnan(G.Val(m,n))
                        [wght{m,n},lidx{m,n}] = intersect_raster_poly(x,y,raster_xco,raster_yco);
                    end
                end
                lcwaitbar(m/(size(G.X,1)-1),h)
            end
        else
            wght = cell(ncell,1);
            lidx = cell(ncell,1);
            for m = 1:ncell
                x = G.CellX{m};
                y = G.CellY{m};
                if ~any(isnan(x) | isnan(y)) && ~isnan(G.CellVal(m))
                    [wght{m},lidx{m}] = intersect_raster_poly(x,y,raster_xco,raster_yco);
                end
                lcwaitbar(m/ncell,h)
            end
        end
        delete(h)
        drawnow
    otherwise
        error('Unsupported method: %s',method)
end
%
nwght=cellfun('length',wght);
Wght = zeros(sum(nwght(:)),3);
offset = 0;
if size(G.X,2)>1
    for m = 1:size(G.X,1)-1
        for n = 1:size(G.X,2)-1
            i = offset+(1:nwght(m,n));
            if ~isempty(i)
                Wght(i,1) = wght{m,n};
                Wght(i,2) = lidx{m,n};
                Wght(i,3) = G.Val(m,n);
                offset = i(end);
            end
        end
    end
else
    for m = 1:ncell
        i = offset+(1:nwght(m));
        if ~isempty(i)
            Wght(i,1) = wght{m};
            Wght(i,2) = lidx{m};
            Wght(i,3) = G.CellVal(m);
            offset = i(end);
        end
    end
end
%
switch method
    case 'maxarea'
        % choose the maximum area
        [Ws,I]=sort(Wght(:,2));
        Wght=Wght(I,:);
        dW=find(diff([0;Wght(:,2);0]));
        num_dW=length(dW)-1;
        h = lcwaitbar(0,'Locating maximum weights ...');
        for ix=1:num_dW
            W=Wght(dW(ix):dW(ix+1)-1,:);
            NoData=cellarea-sum(W(:,1));
            %
            if size(W,1)>1
                [c,ia,ic]=unique(W(:,3));
                szc=size(c);
                W=[accumarray(ic,W(:,1),size(c)) repmat(W(1,2),szc) c];
                [mx,im]=max(W(:,1));
                W=W(im,:);
            end
            %
            Wght(dW(ix):dW(ix+1)-1,:)=NaN;
            if W(1)>NoData
                Wght(dW(ix),:)=W;
            end
            lcwaitbar(ix/num_dW,h)
        end
        delete(h)
        drawnow
        Wght(isnan(Wght(:,1)),:)=[];
    case 'center'
        % check if only single value is used
        if any(diff(sort(Wght(:,2)))==0)
            warning('Raster cell centre coincides with simulation cell edge: using average value');
        end
end
%
% process data sets
%
idataset  = 0;
h = lcwaitbar(0,'Processing data sets ...');
nfld = size(chp,1);
for ifld = 1:nfld
    %
    % get data
    %
    tmop = chp{ifld,1}; % time_operator
    qstr = chp{ifld,2};
    if isempty(qstr)
        continue
    end
    iqnt = chp{ifld,3};
    info = chp{ifld,4};
    %
    iwaq = iqnt>0;
    iflw = find(~iwaq);
    iqnt_waq = iqnt(iwaq);
    idwq = [waq_qnt(iqnt_waq).Val1];
    iqnt_flw = -iqnt(~iwaq);
    %
    ntim = info.tstop - info.tstart + 1;
    %
    if strcmp(tmop,'ident')
        it1 = info.tstart;
        it2 = info.tstop;
    else
        it1 = 1;
        it2 = 1;
        switch tmop
            case {'max','min'}
                D = NaN;
            otherwise
                D  = 0;
                D2 = 0;
        end
        for it = info.tstart:info.tstop
            if any(iwaq)
                [t,Dt] = delwaq('read',dwq,idwq,0,it);
                Dt(Dt==-999)=0;
            end
            if any(~iwaq)
                if any(iwaq)
                    % expand Dt for flow fields
                    Dt(iwaq,:) = Dt;
                end
                for i = length(iflw):-1:1
                    X = qpread(flw_info,flw_qnt(iqnt_flw(i)),'data',it,0);
                    Dt(iflw(i),1:size(X.Val,1)) = X.Val';
                end
            end
            switch tmop
                case 'mean'
                    D = D+Dt;
                    if it==info.tstop
                        D = D/ntim;
                    end
                case 'max'
                    D = max(D,Dt);
                case 'min'
                    D = min(D,Dt);
                case 'std'
                    D  = D  + Dt;
                    D2 = D2 + Dt.^2;
                    if it==info.tstop
                        D = sqrt(D2 - (D.^2)/ntim)/max(1,ntim-1);
                    end
            end
        end
    end
    %
    for it = it1:it2
        if strcmp(tmop,'ident')
            [t,D] = delwaq('read',dwq,idwq,0,it);
            D(D==-999)=NaN;
            tmopstr = sprintf('time_step_%i',it);
        else
            tmopstr = sprintf('%s_%i-%i',tmop,info.tstart,info.tstop);
        end
        for i = 1:length(iqnt)
            data = accumarray(Wght(:,2),Wght(:,1).*D(i,Wght(:,3))',[prod(raster_sz) 1]);
            wght = accumarray(Wght(:,2),Wght(:,1),[prod(raster_sz) 1]);
            if min_area>0
                data(wght<min_area)=0;
                wght(wght<min_area)=0;
            end
            data = reshape(data./wght,raster_sz);
            %
            if yflip
                data = flipud(data);
            end
            %
            filename = fullfile(outdir,sprintf('%s.%s',qstr{i},tmopstr));
            switch raster_info.FileType
                case 'HDR Raster File - BIP/BIL/BSQ'
                    bil('write',filename,raster_info,data);
                case 'arcgrid'
                    FD = raster_info;
                    FD.Data = data';
                    arcgrid('write',FD,[filename '.asc']);
            end
            idataset = idataset+1;
            lcwaitbar(idataset/ndatasets,h)
        end
    end
end
delete(h)

function x=lcwaitbar(varargin)
persistent alpha
if nargout==1
    x = waitbar(varargin{:});
    alpha = 0.01;
else
    a = varargin{1};
    if a==1 || a>=alpha
        alpha = alpha+0.01;
        waitbar(a,varargin{2:end})
    end
end

function [iselqnt,selqnt] = filter_qnt(selqnt,waq_qnt,flw_qnt)
if ~iscell(selqnt)
    selqnt = {selqnt};
end
selqnt = unique(selqnt);
[mem,idx1]=ismember(selqnt,{waq_qnt.Name}');
[mem,idx2]=ismember(selqnt,{waq_qnt.ShortName}');
iselqnt=max(idx1,idx2);
%
[mem,idx3]=ismember(selqnt,{flw_qnt.Name}');
iselqnt(idx3>0)=-idx3(idx3>0);

function [area,lidx] = inside_raster_poly(x,y,raster_xcc,raster_ycc)
col = find(min(x)<raster_xcc & max(x)>raster_xcc);
row = find(min(y)<raster_ycc & max(y)>raster_ycc);
nrowt = length(raster_ycc);
%
[xcc,ycc]=ndgrid(raster_xcc(col),raster_ycc(row));
[col,row]=ndgrid(col,row);
area = inpolygon(xcc,ycc,x,y);
lidx = row + (col-1)*nrowt;
chk = area>0;
area = area(chk);
lidx = lidx(chk);

function [area,lidx] = intersect_raster_poly(x,y,raster_xco,raster_yco)
col = find(min(x)<raster_xco(2:end) & max(x)>raster_xco(1:end-1));
ncol = length(col);
col = col(1):col(end)+1;
row = find(min(y)<raster_yco(2:end) & max(y)>raster_yco(1:end-1));
nrow = length(row);
row = row(1):row(end)+1;
nrowt = length(raster_yco)-1;
%
% figure
% plot(x([1:end 1]),y([1:end 1]))
% set(gca,'ytick',raster_yco(row))
% set(gca,'xtick',raster_xco(col))
% grid on
% hold on
%
area = zeros(nrow,ncol);
lidx = zeros(nrow,ncol);
dx = diff(raster_xco);
dy = diff(raster_yco);
for n = row(1:end-1)
    for m = col(1:end-1)
        in = n-row(1)+1;
        im = m-col(1)+1;
        area(in,im) = intersect_cell_poly((x-raster_xco(m))/dx(m),(y-raster_yco(n))/dx(n))*dx(m)*dy(n);
        lidx(in,im) = n + (m-1)*nrowt;
    end
end
chk  = area>1e-3;
area = area(chk);
lidx = lidx(chk);

function area = intersect_cell_poly(x,y)
area = zeros(1,length(x));
for i = 1:length(x)
    x0 = x(i);
    y0 = y(i);
    if i<length(x)
        x1 = x(i+1);
        y1 = y(i+1);
    else
        x1 = x(1);
        y1 = y(1);
    end
    area(i) = intersect_cell_tri(x0,y0,x1,y1);
end
area = abs(sum(area));

function area = intersect_cell_tri(x0,y0,x1,y1)
%   1 ____
%    |    |       tri (0,0)-(x0,y0)-(x1,y1)
%    |____|
%   0      1
%
if x0<=0 && x1<=0
    % whole triangle to the left of the cell
    area = 0;
elseif y0<=0 && y1<=0
    % whole triangle beneath cell
    area = 0;
else
    % triangle may overlap with cell
    % intersect (0,0)-(x0,y0) with (1,0)-(1,1) ==> (1,y2)
    y2 = y0/x0; % 0 + (y0-0)*(1-0)/(x0-0);
    % intersect (0,0)-(x0,y0) with (0,1)-(1,1) ==> (x2,1)
    x2 = x0/y0; % 0 + (x0-0)*(1-0)/(y0-0);
    %
    % intersect (x0,y0)-(x1,y1) with (0,1)-(1,1) ==> (x3,1)
    x3 = x0 + (x1-x0)*(1-y0)/(y1-y0);
    % intersect (x0,y0)-(x1,y1) with (1,0)-(1,1) ==> (1,y3)
    y3 = y0 + (y1-y0)*(1-x0)/(x1-x0);
    % intersect (x0,y0)-(x1,y1) with (0,0)-(1,0) ==> (x3z,0)
    x3z = x0 + (x1-x0)*(0-y0)/(y1-y0);
    % intersect (x0,y0)-(x1,y1) with (0,0)-(0,1) ==> (0,y3z)
    y3z = y0 + (y1-y0)*(0-x0)/(x1-x0);
    %
    % intersect (0,0)-(x1,y1) with (1,0)-(1,1) ==> (1,y4)
    y4 = y1/x1; %0 + (y1-0)*(1-0)/(x1-0);
    % intersect (0,0)-(x1,y1) with (0,1)-(1,1) ==> (x4,1)
    x4 = x1/y1; %0 + (x1-0)*(1-0)/(y1-0);
    %
    xa = [0 NaN 1 NaN 1 NaN];
    ya = [0 NaN 1 NaN 1 NaN];
%     idx = idx+1;
%     if ismember(idx,39)
%         figure; plot([0 0 1 1 0],[0 1 1 0 0],'b:',[0 x0 x1 0],[0 y0 y1 0],'b',x0,y0,'b*');
%     end
    %
    if (x0<0 && y3z<0) || (y0<0 && x3z>=0)
        if x3z>1 && x0>=1 && x1<1
            xa(2)=1;
            ya(2)=0;
            ya(3)=max(0,min(y3,1));
        else
            xa(3)=max(0,min(x3z,1));
            ya(3)=0;
        end
    elseif x0<0 || y0<0
        if y3z>1 && y0>=1 && y1<1
            xa(2)=0;
            ya(2)=1;
            xa(3)=max(0,min(x3,1));
        else
            xa(3)=0;
            ya(3)=max(0,min(y3z,1));
        end
    elseif x2>=1 && x0>1
        xa(2)=1;
        ya(2)=y2;
        if x1<1
            if y3>1 && x3>0 && x3<1
                xa(4)=x3;
                ya(4)=1;
            else
                ya(3)=max(0,min(y3,1));
            end
        end
    elseif y2>=1 && y0>1
        xa(2)=x2;
        ya(2)=1;
        if y1<1
            if x3>1 && y3>0 && y3<1
                xa(4)=1;
                ya(4)=y3;
            else
                xa(3)=max(0,min(x3,1));
            end
        end
    else
        xa(3)=x0;
        ya(3)=y0;
    end
    %
    if (x1<0 && y3z<0) || (y1<0 && x3z>=0)
        if x3z>1 && x1>=1 && x0<1
            xa(6)=1;
            ya(6)=0;
            ya(5)=max(0,min(y3,1));
        else
            xa(5)=max(0,min(x3z,1));
            ya(5)=0;
        end
    elseif x1<0 || y1<0
        if y3z>1 && y1>=1 && y0<1
            xa(6)=0;
            ya(6)=1;
            xa(5)=max(0,min(x3,1));
        else
            xa(5)=0;
            ya(5)=max(0,min(y3z,1));
        end
    elseif x4>=1 && x1>1
        xa(6)=1;
        ya(6)=y4;
        if x0<1
            if y3>1 && x3>0 && x3<1
                xa(4)=x3;
                ya(4)=1;
            else
                ya(5)=max(0,min(y3,1));
            end
        end
    elseif y4>=1 && y1>1
        xa(6)=x4;
        ya(6)=1;
        if y0<1
            if x3>1 && y3>0 && y3<1
                xa(4)=1;
                ya(4)=y3;
            else
                xa(5)=max(0,min(x3,1));
            end
        end
    else
        xa(5)=x1;
        ya(5)=y1;
    end
    %
    eq=0;
    ip=1;
    for i=2:5
        if isnan(xa(i))
            continue
        elseif eq==1
            if xa(i)==xa(ip)
                xa(ip)=NaN;
                ya(ip)=NaN;
                ip=i;
                continue
            end
        elseif eq==2
            if ya(i)==ya(ip)
                xa(ip)=NaN;
                ya(ip)=NaN;
                ip=i;
                continue
            end
        end
        %
        if xa(i)==xa(ip) && ya(i)==ya(ip)
            xa(ip)=NaN;
            ya(ip)=NaN;
        elseif xa(i)==xa(ip)
            eq = 1;
        elseif ya(i)==ya(ip)
            eq = 2;
        else
            eq = 0;
        end
        ip = i;
    end
    xa(isnan(xa))=[];
    ya(isnan(ya))=[];
    %area=polyarea(xa,ya)*clockwise([0 x0 x1],[0 y0 y1]);
    area=-sum(xa(2:end-1).*(ya(3:end)-ya(2:end-1))-ya(2:end-1).*(xa(3:end)-xa(2:end-1)))/2;
    %
%     if ismember(idx,39)
%         figure; plot([0 0 1 1 0],[0 1 1 0 0],'b:',[0 x0 x1 0],[0 y0 y1 0],'b',[xa xa(1)],[ya ya(1)],'r',x0,y0,'b*');
%         title(sprintf('%g',area))
%     end
end
