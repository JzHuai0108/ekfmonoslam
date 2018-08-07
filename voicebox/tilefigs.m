function tilefigs
%TILEFIGS tile current figures
%
% Possible future options:
%   (1) define screen region to place tiles
%   (2) include/exclude window titles (except optionally first row)
%   (3) preserve sizes and/or aspect ratios
%   (4) place disparate sizes in the most efficient way
%   (5) cascade instead of tile
%   (6) superimpose instead of tile
%   (7) cope with dual screens

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: tilefigs.m 4992 2014-08-08 08:51:12Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

starth=65; % height for start bar
tith=74; % height for window title
winb=4; % width of window border
figl=get(0,'Children');
figl=sort(figl);
nf=length(figl);
scr=get(0,'Screensize');
nfh=ceil(sqrt(nf));
nfv=ceil(nf/nfh);
hpix=floor((scr(3)-scr(1)+1)/nfh); % horizontal pixels available incl border
vpix=floor((scr(4)-scr(2)+1-starth)/nfv); % vertical pixels available incl border
for i=nf:-1:1
    figure(figl(i));
    row=1+floor((i-1)/nfh);
    col=i-(row-1)*nfh;
    %     [i hpix*(col-1)+1 scr(4)-vpix*row+1 hpix*col scr(4)-vpix*(row-1)]
    set(figl(i),'position',[hpix*(col-1)+1+winb scr(4)-vpix*row+1+winb hpix-2*winb vpix-2*winb-tith]);
end