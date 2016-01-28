function varargout=wldep(cmd,varargin)
%WLDEP Read/write Delft3D field files (e.g. depth files).
%   WLDEP can be used to read and write Delft3D field files for spatial
%   varying quantities. This includes files for depth (most common
%   application), roughness (2 fields in one file!), eddy viscosity,
%   critical bed shear stress (for mud fractions), bed composition files
%   (per layer per fraction), initial conditions (water level, velocities,
%   concentrations), etc.
%
%   DEPTH = WLDEP('read',FILENAME,SIZE) reads the file assuming that the
%   file contains a single data matrix of the indicate size.
%
%   DEPTH = WLDEP('read',FILENAME,GRID) reads the file assuming that the
%   file contains a single data matrix of the size that matches the grid
%   data structure GRID as generated by WLGRID.
%
%   STRUCT = WLDEP(...,'multiple') reads multiple fields from the file. The
%   function returns a structure vector with one field: Data. In case of a
%   roughness file which contains 2 fields, the U data is accessible via
%   STRUCT(1).Data and the V data is accessible via STRUCT(2).Data.
%
%   [FLD1,FLD2,...] = WLDEP(...,'multiple') reads multiple fields from the
%   file. The function returns one field to one output argument. In case of
%   a roughness file you would use this as [URGH,VRGH] = WLDEP(...).
%
%   WLDEP('write',FILENAME,MATRIX) writes a depth file for the specified
%   matrix. It will interactively ask whether the data should be multiplied
%   by -1 (e.g. for depth versus elevation) and whether one data line
%   should be added at the high end of both grid directions (the Delft3D
%   depth files are one grid line larger than the associated grid files).
%   To suppress these interactive questions use the command:
%   WLDEP('write',FILENAME,'',MATRIX)
%
%   WLDEP('write',FILENAME,MATRIX1,MATRIX2,MATRIX3,...) writes multiple
%   fields to the specified file. Optionally you may add keywords between
%   data blocks by inserting Keyword strings between the matrices:
%   WLDEP('write',FILENAME,KEYWORD1,MATRIX1,KEYWORD2,MATRIX2,...)
%   The keywords will be enclosed in single quotes in the file:
%
%   'KEYWORD1'
%   DATA of MATRIX1
%   'KEYWORD2'
%   DATA of MATRIX2
%
%   If the keyword is empty, nothing will be written. The file format with
%   keywords is not supported by Delft3D.
%
%   WLDEP('write',FILENAME,STRUCT) writes the STRUCT(i).Data fields to the
%   specified file. Optionally the STRUCT may contain Keyword fields which
%   will be used in the same manner as the keyword arguments mentioned
%   above.
%
%   The default number format used by WLDEP is a fixed point %15.7f. For
%   small parameters like variable grain size diameters this setting may
%   not be appropriate. You can change the format setting by inserting the
%   pair 'format',FORMAT after the FILENAME and before MATRIX argument that
%   should be influenced. Here FORMAT should contain a single valid format
%   specification like '%15.7f', '%16.7e' or '%g'.
%
%   See also: WLGRID

%----- LGPL --------------------------------------------------------------------
%
%   Copyright (C) 2011-2014 Stichting Deltares.
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

if nargin==0
    if nargout>0
        varargout=cell(1,nargout);
    end
    return
end

switch cmd
    case 'read'
        if nargout>1 % automatically add 'multiple'
            if nargin==3 % cmd, filename, size
                varargin{3}='multiple';
            end
        end
        Dep=Local_depread(varargin{:});
        if isstruct(Dep)
            if length(Dep)<nargout
                error('Too many output arguments.')
            end
            if nargout==1
                varargout{1}=Dep;
            else
                varargout = cell(1,max(1,nargout));
                for i=1:length(varargout)
                    varargout{i}=Dep(i).Data;
                end
            end
        elseif nargout>1
            error('Too many output arguments.')
        else % nargout<=1
            varargout={Dep};
        end
    case 'write'
        Out=Local_depwrite(varargin{:});
        if nargout>0
            varargout{1}=Out;
        end
    otherwise
        error('Unknown command: %s',var2str(cmd))
end


function DP=Local_depread(filename,dimvar,option)
% DEPREAD reads depth information from a given filename
%    DEPTH=DEPREAD('FILENAME.DEP',SIZE)
%    or
%    DEPTH=DEPREAD('FILENAME.DEP',GRID)
%    where GRID was generated by GRDREAD.
%
%    ...,'multiple') to read multiple fields from the file.

DP=[];

if nargin<2
    error('No size or grid specified.')
end

if nargin<3
    multiple=0;
else
    multiple=strcmp(option,'multiple');
    if ~multiple
        error('Unknown 3rd argument: %s.',var2str(option))
    end
end

if strcmp(filename,'?')
    [fname,fpath]=uigetfile('*.*','Select depth file');
    if ~ischar(fname)
        return
    end
    filename=[fpath,fname];
end

fid=fopen(filename);
if fid<0
    error('Cannot open %s.',filename)
end
try
    if isstruct(dimvar) % new grid format G.X, G.Y, G.Enclosure
        dim=size(dimvar.X)+1;
    elseif iscell(dimvar) % old grid format {X Y Enclosure}
        dim=size(dimvar{1})+1;
    else
        dim=dimvar;
    end
    i=1;
    while 1
        %
        % Skip lines starting with *
        %
        line='*';
        cl=0;
        while ~isempty(line) && line(1)=='*'
            if cl>0
                DP(i).Comment{cl,1}=line;
            end
            cl=cl+1;
            currentpoint=ftell(fid);
            line=fgetl(fid);
        end
        fseek(fid,currentpoint,-1);
        %
        [DP(i).Data,NRead]=fscanf(fid,'%f',dim);
        if NRead==0 % accept dredging input files containing 'keyword' on the first line ...
            str=fscanf(fid,['''%[^' char(10) char(13) ''']%['']']);
            if ~isempty(str) && isequal(str(end),'''')
                [DP(i).Data,NRead]=fscanf(fid,'%f',dim);
                DP(i).Keyword=str(1:end-1);
            end
        end
        DP(i).Data=DP(i).Data;
        %
        % Read remainder of last line
        %
        Rem=fgetl(fid);
        if ~ischar(Rem)
            Rem='';
        else
            Rem=deblank(Rem);
        end
        if NRead<prod(dim)
            if isempty(Rem)
                Str=sprintf('Not enough data in the file for complete field %i (only %i out of %i values).',i,NRead,prod(dim));
                if i==1 % most probably wrong file format
                    error(Str)
                else
                    warning(Str)
                end
            else
                error('Invalid string while reading data: %s',Rem)
            end
        end
        pos=ftell(fid);
        if isempty(fscanf(fid,'%f',1))
            break % no more data (at least not readable)
        elseif ~multiple
            fprintf('More data in the file. Use ''multiple'' option to read all fields.\n');
            break % don't read data although there seems to be more ...
        end
        fseek(fid,pos,-1);
        i=i+1;
    end
    fclose(fid);
catch
    fclose(fid);
    rethrow(lasterror)
end
if ~multiple
    DP=DP.Data;
end


function OK=Local_depwrite(filename,varargin)
% DEPWRITE writes depth information to a given filename
%
% Usage: depwrite('filename',Matrix)
%
%    or: depwrite('filename',Struct)
%        where Struct is a structure vector with one field: Data

fid=fopen(filename,'w');
Keyword='';
interactive = length(varargin)==1;
format = '%15.8f';
idp = 0;
while idp<length(varargin)
    idp = idp+1;
    DP = varargin{idp};
    %
    if ischar(DP)
        if idp==length(varargin)
            % ignore this last argument
        else
            DP2 = varargin{idp+1};
            if strcmpi(DP,'format') && ischar(DP2)
                format = DP2;
                idp = idp+1;
            else
                Keyword = DP;
            end
        end
    elseif isstruct(DP)
        % DP(1:N).Data=Matrix;
        
        for i=1:length(DP)
            if isfield(DP,'Keyword')
                Keyword = DP(i).Keyword;
            else
                Keyword = '';
            end
            writeblock(fid,DP(i).Data,Keyword,format);
        end
        Keyword = '';
    else
        % DP=Matrix;
        if DP(end,end)~=-999 && interactive
            switch input('Negate data points? (Y/N) ','s')
                case {'Y','y'}
                    DP=-DP;
                otherwise
            end
            switch input('Grid extension: 9 (-999 values)/B (boundary values) /N (Don''t extend) ','s')
                case {'9'}
                    DP=[DP -999*ones(size(DP,1),1); ...
                        -999*ones(1,size(DP,2)+1)];
                case {'B','b'}
                    DP=[DP DP(:,end); ...
                        DP(end,:) DP(end,end)];
                otherwise
            end
        end
        writeblock(fid,DP,Keyword,format);
        Keyword = '';
    end
end
fclose(fid);
OK=1;


function writeblock(fid,DP,Keyword,format)
if ~isempty(Keyword)
    fprintf(fid,'''%s''\n',Keyword);
end

DP(isnan(DP))=-999;

lformat = length(format)+2;
Frmt=repmat([format '  '],[1 size(DP,1)]);
k=lformat*12;
Frmt((k-1):k:length(Frmt))='\';
Frmt(k:k:length(Frmt))='n';
Frmt(end-1:end)='\n';
Frmt=strrep(Frmt,'  ',' ');

szDP=size(DP);
if length(szDP)<3
    kmax=1;
else
    kmax=prod(szDP(3:end));
end
for k=1:kmax
    fprintf(fid,Frmt,DP(:,:,k));
end
