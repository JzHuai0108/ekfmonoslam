function calib_gui(mode),

% calib_gui(mode)
%
% Runs the Camera Calibration Toolbox.
% Set mode to 1 to run the memory efficient version.
% Any other value for mode will run the normal version (see documentation)
%
% INFORMATION ABOUT THE MEMORY EFFICIENT MODE FOR THE CAMERA CALIBRATION TOOLBOX:
%
% If your calibration images are large, or if you calibrate using a lot of images, you may have experienced memory problems
% in Matlab when using the calibration toolbox (OUT OF MEMORY errors). If this is the case, you can now run the
% new memory efficient version of the toolbox that loads every image one by one without storing them all in memory.
% If you choose to run the standard version of the toolbox now, you can always switch to the other memory efficient mode
% later in case the OUT OF MEMORY error message is encountered. The two modes of operation are totally compatible.


if nargin < 1,
    
    cell_list = {};
    
    fig_number = 1;
    
    title_figure = 'Camera Calibration Toolbox - Select mode of operation:';
    
    cell_list{1,1} = {'Standard (all the images are stored in memory)','calib_gui_normal;'};
    cell_list{2,1} = {'Memory efficient (the images are loaded one by one)','calib_gui_no_read;'};
    cell_list{3,1} = {'Exit',['disp(''Bye. To run again, type calib_gui.''); close(' num2str(fig_number) ');']};
    
    
    show_window(cell_list,fig_number,title_figure);
    
    
else
    
    if ~isnumeric(mode),
        mode = str2num(mode);
    end;
    
    mode = mode(1);

    if isempty(mode),
        mode = 0;
    end;
    
    if mode == 1,
        
        calib_gui_no_read;
        
    else
        
        calib_gui_normal;
        
    end;
    
end;




%------- DO NOT EDIT ANYTHING BELOW THIS LINE -----------%

function show_window(cell_list,fig_number,title_figure,x_size,y_size,gap_x,font_name,font_size)


if ~exist('cell_list'),
    error('No description of the functions');
end;

if ~exist('fig_number'),
    fig_number = 1;
end;
if ~exist('title_figure'),
    title_figure = '';
end;
if ~exist('x_size'),
    x_size = 340;
end;
if ~exist('y_size'),
    y_size = 18.7;
end;
if ~exist('gap_x'),
    gap_x = 0;
end;
if ~exist('font_name'),
    font_name = 'clean';
end;
if ~exist('font_size'),
    font_size = 12;
end;

figure(fig_number); clf;
pos = get(fig_number,'Position');

[n_row,n_col] = size(cell_list);

fig_size_x = x_size*n_col+(n_col+1)*gap_x;
fig_size_y = y_size*n_row+(n_row+1)*gap_x;

set(fig_number,'Units','points', ...
	'BackingStore','off', ...
	'Color',[0.8 0.8 0.8], ...
	'MenuBar','none', ...
	'Resize','off', ...
	'Name',title_figure, ...
'Position',[pos(1) pos(2) fig_size_x fig_size_y], ...
'NumberTitle','off'); %,'WindowButtonMotionFcn',['figure(' num2str(fig_number) ');']);

h_mat = zeros(n_row,n_col);

posx = zeros(n_row,n_col);
posy = zeros(n_row,n_col);

for i=n_row:-1:1,
   for j = n_col:-1:1,
      posx(i,j) = gap_x+(j-1)*(x_size+gap_x);
      posy(i,j) = fig_size_y - i*(gap_x+y_size);
   end;
end;

for i=n_row:-1:1,
    for j = n_col:-1:1,
        if ~isempty(cell_list{i,j}),
            if ~isempty(cell_list{i,j}{1}) & ~isempty(cell_list{i,j}{2}),
                h_mat(i,j) = uicontrol('Parent',fig_number, ...
                    'Units','points', ...
                    'Callback',cell_list{i,j}{2}, ...
                    'ListboxTop',0, ...
                    'Position',[posx(i,j)  posy(i,j)  x_size   y_size], ...
                    'String',cell_list{i,j}{1}, ...
                    'fontsize',font_size,...
                    'fontname',font_name,...
                    'Tag','Pushbutton1');
            end;
        end;
    end;
end;

%------ END PROTECTED REGION ----------------%
