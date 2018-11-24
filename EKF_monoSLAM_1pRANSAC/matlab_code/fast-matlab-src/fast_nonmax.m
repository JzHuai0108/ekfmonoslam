% FAST_NONMAX perform non-maximal suppression on FAST features.
%
%     nonmax = FAST_NONMAX(image, threshold, FAST_CORNER_DETECT_9(image, threshold));
%     returns a list of nonmaximally suppressed corners with the X coordinate
%     in nonmax(:,1) and Y in nonmax(:,2).
%
%     If you use this in published work, please cite:
%       Fusing Points and Lines for High Performance Tracking, E. Rosten and T. Drummond, ICCV 2005
%       Machine learning for high-speed corner detection, E. Rosten and T. Drummond, ECCV 2006
%     The Bibtex entries are:
%     
%     @inproceedings{rosten_2005_tracking,
%         title       =    "Fusing points and lines for high performance tracking.",
%         author      =    "Edward Rosten and Tom Drummond",
%         year        =    "2005",
%         month       =    "October",
%         pages       =    "1508--1511",
%         volume      =    "2",
%         booktitle   =    "IEEE International Conference on Computer Vision",
%         notes       =    "Oral presentation",
%         url         =    "http://mi.eng.cam.ac.uk/~er258/work/rosten_2005_tracking.pdf"
%     }
%     
%     @inproceedings{rosten_2006_machine,
%         title       =    "Machine learning for high-speed corner detection",
%         author      =    "Edward Rosten and Tom Drummond",
%         year        =    "2006",
%         month       =    "May",
%         booktitle   =    "European Conference on Computer Vision (to appear)",
%         notes       =    "Poster presentation",
%         url         =    "http://mi.eng.cam.ac.uk/~er258/work/rosten_2006_machine.pdf"
%     }
%
%
%
%     See also FAST_CORNER_DETECT_9, FAST_CORNER_DETECT_10
%              FAST_CORNER_DETECT_11, FAST_CORNER_DETECT_12
function ret = fast_nonmax(im, barrier, c)
	sz = size(im);
	ysize = sz(1);
	
	
	% x, y
	dir=[	0,3
			1,3
			2,2
			3,1
			3,0
			3,-1
			2,-2
			1,-3
			0,-3
			-1,-3
			-2,-2
			-3,-1
			-3,0
			-3,1
			-2,2
			-1,3];

	dir = dir(:,2) + dir(:,1) * ysize;
	centres = c(:,2) + (c(:,1)-1) * ysize;

	%First pass: compute scores for all of the pixels:
	scores = zeros(sz(1) * sz(2), 1);
		
	
	for i=1:length(centres)

		p1 = im(centres(i)+dir) -( im(centres(i)) + barrier);
		pos = sum(p1 .* (p1 > 0));
		n1 = im(centres(i)) - barrier - im(centres(i)+dir);
		neg = sum(n1 .* (n1 > 0));

		scores(centres(i)) = max([pos neg]);
	end

	
	cs = zeros(length(centres),1);
	
	up = -1;
	down = 1;
	left =-ysize;
	right = ysize;

	%second pass: get local maxima
	for i=1:length(centres)
		p = centres(i);
		square=p +  [ up down left right up+left up+right down+left down+right];
		
		cs(i) = (sum(scores(p) >= scores(square)) == 8);
	end
	
	%Get the maxima positions
	maximas = centres(find(cs));

	ret = zeros(length(maximas), 2);	

	ret(:,1) = 1 + floor(maximas / ysize);
	ret(:,2) = mod(maximas, ysize);
