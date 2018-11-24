%FAST_CORNER_DETECT_11 perform an 11 point FAST corner detection.
%    corners = FAST_CORNER_DETECT_11(image, threshold) performs the detection on the image
%    and returns the X coordinates in corners(:,1) and the Y coordinares in corners(:,2).
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
%     Additional information from the generating program:
%
%     Automatically generated code
%     Parameters:
%     splat_subtree = 1
%     corner_pointers = 2
%     force_first_question = 0
%     corner_type = 11
%     barrier = 25
%     
%     Data:
%     Number of frames:    120
%     Potential features:  25786080
%     Real features:       93213
%     Questions per pixel: 2.37184
%
%
%    See also FAST_NONMAX FAST_CORNER_DETECT_9 FAST_CORNER_DETECT_10 FAST_CORNER_DETECT_11 FAST_CORNER_DETECT_12
%
function coords = fast_corner_detect_11(im, threshold)								
	sz = size(im);																
	xsize=sz(2);																
	ysize=sz(1);																
	cs = zeros(5000, 2);														
	nc = 0;																		
	for x = 4 : xsize - 3														
		for y = 4 : ysize -3													
			cb = im(y,x) + threshold;											
			c_b= im(y,x) - threshold;											
            if im(y+-3,x+0) > cb
                if im(y+3,x+1) > cb
                    if im(y+0,x+3) > cb
                        if im(y+1,x+3) > cb
                            if im(y+-2,x+2) > cb
                                if im(y+3,x+-1) > cb
                                    if im(y+-1,x+3) > cb
                                        if im(y+2,x+2) > cb
                                            if im(y+-3,x+-1) > cb
                                                if im(y+-3,x+1) > cb
                                                    if im(y+3,x+0) > cb
                                                    elseif im(y+3,x+0) < c_b
                                                        continue;
                                                    else
                                                        if im(y+-1,x+-3) > cb
                                                            if im(y+-2,x+-2) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                elseif im(y+-3,x+1) < c_b
                                                    continue;
                                                else
                                                    if im(y+0,x+-3) > cb
                                                        if im(y+2,x+-2) > cb
                                                            if im(y+3,x+0) > cb
                                                                if im(y+1,x+-3) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            elseif im(y+-3,x+-1) < c_b
                                                continue;
                                            else
                                                if im(y+2,x+-2) > cb
                                                    if im(y+-3,x+1) > cb
                                                        if im(y+3,x+0) > cb
                                                        else
                                                            continue;
                                                        end
                                                    elseif im(y+-3,x+1) < c_b
                                                        continue;
                                                    else
                                                        if im(y+0,x+-3) > cb
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        elseif im(y+2,x+2) < c_b
                                            continue;
                                        else
                                            if im(y+0,x+-3) > cb
                                                if im(y+1,x+-3) > cb
                                                    if im(y+-1,x+-3) > cb
                                                        if im(y+-3,x+1) > cb
                                                            if im(y+-3,x+-1) > cb
                                                                if im(y+-2,x+-2) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    elseif im(y+-1,x+3) < c_b
                                        continue;
                                    else
                                        if im(y+0,x+-3) > cb
                                            if im(y+-1,x+-3) > cb
                                                if im(y+1,x+-3) > cb
                                                    if im(y+2,x+-2) > cb
                                                        if im(y+-2,x+-2) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    end
                                elseif im(y+3,x+-1) < c_b
                                    if im(y+-1,x+-3) > cb
                                        if im(y+-3,x+1) > cb
                                            if im(y+-3,x+-1) > cb
                                                if im(y+-2,x+-2) > cb
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    elseif im(y+-1,x+-3) < c_b
                                        continue;
                                    else
                                        if im(y+3,x+0) > cb
                                        else
                                            continue;
                                        end
                                    end
                                else
                                    if im(y+-1,x+-3) > cb
                                        if im(y+-3,x+1) > cb
                                            if im(y+-2,x+-2) > cb
                                                if im(y+-1,x+3) > cb
                                                    if im(y+2,x+2) > cb
                                                        if im(y+-3,x+-1) > cb
                                                        else
                                                            continue;
                                                        end
                                                    elseif im(y+2,x+2) < c_b
                                                        continue;
                                                    else
                                                        if im(y+1,x+-3) > cb
                                                            if im(y+0,x+-3) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    elseif im(y+-1,x+-3) < c_b
                                        if im(y+-2,x+-2) > cb
                                            if im(y+3,x+0) > cb
                                                if im(y+2,x+-2) > cb
                                                    if im(y+2,x+2) > cb
                                                        if im(y+-3,x+1) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                elseif im(y+2,x+-2) < c_b
                                                    if im(y+-3,x+1) > cb
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    if im(y+-1,x+3) > cb
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        if im(y+3,x+0) > cb
                                            if im(y+-2,x+-2) > cb
                                                if im(y+-1,x+3) > cb
                                                    if im(y+2,x+2) > cb
                                                        if im(y+-3,x+1) > cb
                                                            if im(y+-3,x+-1) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    end
                                end
                            elseif im(y+-2,x+2) < c_b
                                if im(y+0,x+-3) > cb
                                    if im(y+3,x+-1) > cb
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                if im(y+0,x+-3) > cb
                                    if im(y+3,x+-1) > cb
                                        if im(y+2,x+-2) > cb
                                            if im(y+-1,x+-3) > cb
                                                if im(y+1,x+-3) > cb
                                                    if im(y+3,x+0) > cb
                                                        if im(y+2,x+2) > cb
                                                            if im(y+-2,x+-2) > cb
                                                            elseif im(y+-2,x+-2) < c_b
                                                                continue;
                                                            else
                                                                if im(y+-1,x+3) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        elseif im(y+2,x+2) < c_b
                                                            continue;
                                                        else
                                                            if im(y+-3,x+1) > cb
                                                                if im(y+-2,x+-2) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        elseif im(y+1,x+3) < c_b
                            if im(y+1,x+-3) > cb
                                if im(y+2,x+2) < c_b
                                    if im(y+2,x+-2) > cb
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            if im(y+0,x+-3) > cb
                                if im(y+-1,x+-3) > cb
                                    if im(y+2,x+-2) > cb
                                        if im(y+1,x+-3) > cb
                                            if im(y+-2,x+-2) > cb
                                                if im(y+-3,x+1) > cb
                                                    if im(y+-3,x+-1) > cb
                                                        if im(y+3,x+-1) > cb
                                                        elseif im(y+3,x+-1) < c_b
                                                            continue;
                                                        else
                                                            if im(y+-2,x+2) > cb
                                                                if im(y+-1,x+3) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                elseif im(y+-3,x+1) < c_b
                                                    continue;
                                                else
                                                    if im(y+2,x+2) > cb
                                                        if im(y+3,x+0) > cb
                                                            if im(y+-3,x+-1) > cb
                                                                if im(y+3,x+-1) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        end
                    elseif im(y+0,x+3) < c_b
                        if im(y+0,x+-3) > cb
                            if im(y+-1,x+-3) > cb
                                if im(y+2,x+-2) > cb
                                    if im(y+-2,x+-2) > cb
                                        if im(y+2,x+2) > cb
                                            if im(y+1,x+-3) > cb
                                                if im(y+-3,x+-1) > cb
                                                    if im(y+3,x+-1) > cb
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        elseif im(y+2,x+2) < c_b
                                            if im(y+-3,x+1) > cb
                                                if im(y+3,x+-1) > cb
                                                    if im(y+1,x+-3) > cb
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            if im(y+-3,x+1) > cb
                                                if im(y+1,x+-3) > cb
                                                    if im(y+3,x+-1) > cb
                                                        if im(y+-3,x+-1) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    else
                        if im(y+0,x+-3) > cb
                            if im(y+-1,x+-3) > cb
                                if im(y+2,x+-2) > cb
                                    if im(y+1,x+-3) > cb
                                        if im(y+-2,x+-2) > cb
                                            if im(y+3,x+-1) > cb
                                                if im(y+-3,x+1) > cb
                                                    if im(y+-3,x+-1) > cb
                                                        if im(y+3,x+0) > cb
                                                        elseif im(y+3,x+0) < c_b
                                                            continue;
                                                        else
                                                            if im(y+-1,x+3) > cb
                                                                if im(y+-2,x+2) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                elseif im(y+-3,x+1) < c_b
                                                    continue;
                                                else
                                                    if im(y+2,x+2) > cb
                                                        if im(y+-3,x+-1) > cb
                                                            if im(y+3,x+0) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    end
                elseif im(y+3,x+1) < c_b
                    if im(y+3,x+-1) > cb
                        if im(y+-1,x+3) > cb
                            if im(y+1,x+-3) > cb
                                if im(y+-2,x+-2) > cb
                                    if im(y+-1,x+-3) > cb
                                        if im(y+-2,x+2) > cb
                                            if im(y+-3,x+1) > cb
                                                if im(y+-3,x+-1) > cb
                                                    if im(y+0,x+-3) > cb
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        elseif im(y+-1,x+3) < c_b
                            if im(y+3,x+0) > cb
                                if im(y+-2,x+2) > cb
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    elseif im(y+3,x+-1) < c_b
                        if im(y+0,x+-3) > cb
                            if im(y+1,x+3) > cb
                                if im(y+1,x+-3) > cb
                                    if im(y+-1,x+3) > cb
                                        if im(y+-3,x+1) > cb
                                            if im(y+-2,x+-2) > cb
                                                if im(y+0,x+3) > cb
                                                    if im(y+-1,x+-3) > cb
                                                        if im(y+-2,x+2) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                elseif im(y+1,x+-3) < c_b
                                    continue;
                                else
                                    if im(y+2,x+2) > cb
                                        if im(y+-2,x+2) > cb
                                            if im(y+-3,x+1) > cb
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            elseif im(y+1,x+3) < c_b
                                continue;
                            else
                                if im(y+2,x+-2) > cb
                                    if im(y+0,x+3) > cb
                                        if im(y+-2,x+2) > cb
                                            if im(y+-3,x+1) > cb
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        elseif im(y+0,x+-3) < c_b
                            if im(y+0,x+3) > cb
                                if im(y+-3,x+-1) < c_b
                                    if im(y+1,x+3) < c_b
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            elseif im(y+0,x+3) < c_b
                                if im(y+-2,x+2) > cb
                                    if im(y+-1,x+-3) < c_b
                                        if im(y+-1,x+3) < c_b
                                            if im(y+2,x+2) < c_b
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                elseif im(y+-2,x+2) < c_b
                                    if im(y+1,x+3) < c_b
                                        if im(y+-1,x+3) < c_b
                                            if im(y+1,x+-3) < c_b
                                                if im(y+2,x+-2) < c_b
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    if im(y+-1,x+-3) < c_b
                                        if im(y+-1,x+3) > cb
                                            continue;
                                        elseif im(y+-1,x+3) < c_b
                                            if im(y+2,x+2) < c_b
                                                if im(y+2,x+-2) < c_b
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            if im(y+-2,x+-2) < c_b
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                end
                            else
                                if im(y+-3,x+-1) < c_b
                                else
                                    continue;
                                end
                            end
                        else
                            if im(y+-3,x+1) < c_b
                            else
                                continue;
                            end
                        end
                    else
                        if im(y+0,x+3) > cb
                            if im(y+1,x+-3) > cb
                                if im(y+2,x+-2) > cb
                                    if im(y+1,x+3) > cb
                                        if im(y+-2,x+2) > cb
                                            if im(y+-1,x+-3) > cb
                                                if im(y+3,x+0) > cb||im(y+3,x+0) < c_b
                                                    continue;
                                                else
                                                    if im(y+-3,x+1) > cb
                                                        if im(y+-2,x+-2) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    elseif im(y+1,x+3) < c_b
                                        if im(y+-3,x+-1) > cb
                                        else
                                            continue;
                                        end
                                    else
                                        if im(y+-2,x+2) > cb
                                            if im(y+-3,x+1) > cb
                                                if im(y+0,x+-3) > cb
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    end
                                elseif im(y+2,x+-2) < c_b
                                    continue;
                                else
                                    if im(y+1,x+3) > cb
                                        if im(y+-1,x+-3) > cb
                                            if im(y+-2,x+2) > cb
                                                if im(y+0,x+-3) > cb
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    end
                else
                    if im(y+0,x+-3) > cb
                        if im(y+-1,x+3) > cb
                            if im(y+1,x+3) > cb
                                if im(y+-2,x+2) > cb
                                    if im(y+-2,x+-2) > cb
                                        if im(y+0,x+3) > cb
                                            if im(y+-3,x+1) > cb
                                                if im(y+1,x+-3) > cb
                                                    if im(y+-1,x+-3) > cb
                                                        if im(y+-3,x+-1) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                elseif im(y+1,x+-3) < c_b
                                                    continue;
                                                else
                                                    if im(y+2,x+2) > cb
                                                        if im(y+-1,x+-3) > cb
                                                            if im(y+-3,x+-1) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        elseif im(y+0,x+3) < c_b
                                            continue;
                                        else
                                            if im(y+3,x+-1) > cb
                                                if im(y+-1,x+-3) > cb
                                                    if im(y+-3,x+1) > cb
                                                        if im(y+1,x+-3) > cb
                                                            if im(y+-3,x+-1) > cb
                                                                if im(y+2,x+-2) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            elseif im(y+1,x+3) < c_b
                                if im(y+3,x+-1) > cb
                                    if im(y+2,x+-2) > cb
                                        if im(y+-2,x+2) > cb
                                            if im(y+-1,x+-3) > cb
                                                if im(y+-3,x+1) > cb
                                                    if im(y+-2,x+-2) > cb
                                                        if im(y+-3,x+-1) > cb
                                                            if im(y+1,x+-3) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                elseif im(y+3,x+-1) < c_b
                                    continue;
                                else
                                    if im(y+0,x+3) > cb
                                        if im(y+2,x+-2) > cb
                                            if im(y+-2,x+2) > cb
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            else
                                if im(y+3,x+-1) > cb
                                    if im(y+2,x+-2) > cb
                                        if im(y+-1,x+-3) > cb
                                            if im(y+-2,x+2) > cb
                                                if im(y+1,x+-3) > cb
                                                    if im(y+-2,x+-2) > cb
                                                        if im(y+-3,x+1) > cb
                                                            if im(y+-3,x+-1) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                elseif im(y+3,x+-1) < c_b
                                    continue;
                                else
                                    if im(y+0,x+3) > cb
                                        if im(y+2,x+-2) > cb
                                            if im(y+-2,x+2) > cb
                                                if im(y+-1,x+-3) > cb
                                                    if im(y+-2,x+-2) > cb
                                                        if im(y+-3,x+1) > cb
                                                            if im(y+1,x+-3) > cb
                                                                if im(y+-3,x+-1) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            end
                        elseif im(y+-1,x+3) < c_b
                            if im(y+-2,x+2) > cb
                                if im(y+3,x+0) > cb
                                    if im(y+-1,x+-3) > cb
                                        if im(y+2,x+-2) > cb
                                            if im(y+-3,x+1) > cb
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            if im(y+3,x+0) > cb
                                if im(y+-2,x+2) > cb
                                    if im(y+-1,x+-3) > cb
                                        if im(y+2,x+-2) > cb
                                            if im(y+1,x+-3) > cb
                                                if im(y+-3,x+-1) > cb
                                                    if im(y+-2,x+-2) > cb
                                                        if im(y+3,x+-1) > cb
                                                            if im(y+-3,x+1) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        end
                    else
                        continue;
                    end
                end
            elseif im(y+-3,x+0) < c_b
                if im(y+3,x+0) > cb
                    if im(y+1,x+3) > cb
                        if im(y+-1,x+-3) > cb
                            if im(y+-1,x+3) > cb
                                if im(y+0,x+-3) > cb
                                    if im(y+2,x+-2) > cb
                                        if im(y+2,x+2) > cb
                                            if im(y+0,x+3) > cb
                                                if im(y+3,x+-1) > cb
                                                    if im(y+1,x+-3) > cb
                                                        if im(y+3,x+1) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            elseif im(y+-1,x+3) < c_b
                                continue;
                            else
                                if im(y+-2,x+-2) > cb
                                    if im(y+0,x+3) > cb
                                        if im(y+1,x+-3) > cb
                                            if im(y+0,x+-3) > cb
                                                if im(y+2,x+-2) > cb
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        elseif im(y+-1,x+-3) < c_b
                            if im(y+0,x+3) > cb
                                if im(y+0,x+-3) > cb
                                    if im(y+-2,x+2) > cb
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            elseif im(y+0,x+3) < c_b
                                if im(y+2,x+-2) < c_b
                                    if im(y+3,x+-1) > cb||im(y+3,x+-1) < c_b
                                        continue;
                                    else
                                        if im(y+-2,x+2) < c_b
                                        else
                                            continue;
                                        end
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            if im(y+-2,x+2) > cb
                                if im(y+0,x+-3) > cb
                                    if im(y+0,x+3) > cb
                                        if im(y+3,x+-1) > cb
                                            if im(y+-1,x+3) > cb
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        end
                    elseif im(y+1,x+3) < c_b
                        if im(y+0,x+-3) < c_b
                            if im(y+-1,x+3) < c_b
                                if im(y+1,x+-3) > cb
                                    continue;
                                elseif im(y+1,x+-3) < c_b
                                    if im(y+-1,x+-3) < c_b
                                        if im(y+-2,x+2) < c_b
                                            if im(y+0,x+3) < c_b
                                                if im(y+-3,x+1) < c_b
                                                    if im(y+-3,x+-1) < c_b
                                                        if im(y+-2,x+-2) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    if im(y+2,x+2) < c_b
                                        if im(y+-3,x+-1) < c_b
                                            if im(y+-2,x+2) < c_b
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    else
                        if im(y+2,x+-2) < c_b
                            if im(y+0,x+3) < c_b
                                if im(y+0,x+-3) < c_b
                                    if im(y+-2,x+2) < c_b
                                        if im(y+3,x+1) > cb
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    end
                elseif im(y+3,x+0) < c_b
                    if im(y+0,x+-3) > cb
                        if im(y+1,x+3) < c_b
                            if im(y+-1,x+3) < c_b
                                if im(y+-2,x+-2) > cb
                                    if im(y+2,x+-2) > cb
                                        if im(y+-3,x+-1) < c_b
                                            if im(y+3,x+-1) < c_b
                                                if im(y+-2,x+2) < c_b
                                                    if im(y+2,x+2) < c_b
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    elseif im(y+2,x+-2) < c_b
                                        if im(y+-2,x+2) < c_b
                                            if im(y+0,x+3) < c_b
                                                if im(y+2,x+2) < c_b
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        if im(y+-3,x+-1) < c_b
                                            if im(y+-2,x+2) < c_b
                                                if im(y+3,x+-1) < c_b
                                                    if im(y+2,x+2) < c_b
                                                        if im(y+0,x+3) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    end
                                elseif im(y+-2,x+-2) < c_b
                                    if im(y+0,x+3) < c_b
                                        if im(y+2,x+2) < c_b
                                            if im(y+-3,x+-1) > cb
                                                continue;
                                            elseif im(y+-3,x+-1) < c_b
                                                if im(y+-3,x+1) < c_b
                                                    if im(y+3,x+1) < c_b
                                                        if im(y+-2,x+2) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                if im(y+2,x+-2) < c_b
                                                    if im(y+3,x+-1) < c_b
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    if im(y+3,x+-1) < c_b
                                        if im(y+-3,x+-1) > cb
                                            continue;
                                        elseif im(y+-3,x+-1) < c_b
                                            if im(y+-2,x+2) < c_b
                                                if im(y+2,x+2) < c_b
                                                    if im(y+0,x+3) < c_b
                                                        if im(y+-3,x+1) < c_b
                                                            if im(y+3,x+1) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            if im(y+2,x+-2) < c_b
                                                if im(y+0,x+3) < c_b
                                                    if im(y+-2,x+2) < c_b
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    elseif im(y+0,x+-3) < c_b
                        if im(y+-1,x+-3) > cb
                            if im(y+0,x+3) < c_b
                                if im(y+3,x+-1) < c_b
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        elseif im(y+-1,x+-3) < c_b
                            if im(y+2,x+-2) > cb
                                if im(y+0,x+3) < c_b
                                    if im(y+2,x+2) < c_b
                                        if im(y+-2,x+2) < c_b
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            elseif im(y+2,x+-2) < c_b
                                if im(y+1,x+-3) > cb
                                    continue;
                                elseif im(y+1,x+-3) < c_b
                                    if im(y+-2,x+-2) > cb
                                        continue;
                                    elseif im(y+-2,x+-2) < c_b
                                        if im(y+-3,x+1) > cb
                                            if im(y+2,x+2) < c_b
                                            else
                                                continue;
                                            end
                                        elseif im(y+-3,x+1) < c_b
                                            if im(y+3,x+1) > cb
                                                continue;
                                            elseif im(y+3,x+1) < c_b
                                                if im(y+3,x+-1) > cb
                                                    continue;
                                                elseif im(y+3,x+-1) < c_b
                                                    if im(y+-3,x+-1) > cb
                                                        continue;
                                                    elseif im(y+-3,x+-1) < c_b
                                                    else
                                                        if im(y+0,x+3) < c_b
                                                            if im(y+2,x+2) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else
                                                    if im(y+0,x+3) < c_b
                                                        if im(y+-1,x+3) < c_b
                                                            if im(y+-2,x+2) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                if im(y+-2,x+2) < c_b
                                                    if im(y+3,x+-1) > cb
                                                        continue;
                                                    elseif im(y+3,x+-1) < c_b
                                                        if im(y+-3,x+-1) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        if im(y+0,x+3) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        else
                                            if im(y+2,x+2) < c_b
                                                if im(y+3,x+-1) < c_b
                                                    if im(y+3,x+1) < c_b
                                                        if im(y+-3,x+-1) > cb
                                                            continue;
                                                        elseif im(y+-3,x+-1) < c_b
                                                        else
                                                            if im(y+0,x+3) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        if im(y+0,x+3) < c_b
                                            if im(y+2,x+2) < c_b
                                                if im(y+-1,x+3) < c_b
                                                    if im(y+1,x+3) < c_b
                                                        if im(y+3,x+-1) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    end
                                else
                                    if im(y+0,x+3) < c_b
                                        if im(y+1,x+3) < c_b
                                            if im(y+-2,x+2) < c_b
                                                if im(y+2,x+2) < c_b
                                                    if im(y+-1,x+3) < c_b
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            else
                                if im(y+0,x+3) < c_b
                                    if im(y+-2,x+2) < c_b
                                        if im(y+1,x+3) < c_b
                                            if im(y+-3,x+1) < c_b
                                                if im(y+2,x+2) > cb
                                                    continue;
                                                elseif im(y+2,x+2) < c_b
                                                    if im(y+-1,x+3) < c_b
                                                        if im(y+-3,x+-1) < c_b
                                                            if im(y+-2,x+-2) > cb
                                                                continue;
                                                            elseif im(y+-2,x+-2) < c_b
                                                            else
                                                                if im(y+3,x+-1) < c_b
                                                                else
                                                                    continue;
                                                                end
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    if im(y+1,x+-3) < c_b
                                                        if im(y+-1,x+3) < c_b
                                                            if im(y+-2,x+-2) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        else
                            if im(y+0,x+3) < c_b
                                if im(y+1,x+3) < c_b
                                    if im(y+-2,x+2) < c_b
                                        if im(y+2,x+2) < c_b
                                            if im(y+-1,x+3) < c_b
                                                if im(y+3,x+-1) > cb
                                                    continue;
                                                elseif im(y+3,x+-1) < c_b
                                                    if im(y+3,x+1) < c_b
                                                        if im(y+2,x+-2) > cb
                                                            continue;
                                                        elseif im(y+2,x+-2) < c_b
                                                        else
                                                            if im(y+-3,x+-1) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    if im(y+-2,x+-2) < c_b
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        end
                    else
                        if im(y+0,x+3) < c_b
                            if im(y+1,x+3) < c_b
                                if im(y+-1,x+3) < c_b
                                    if im(y+-2,x+2) < c_b
                                        if im(y+2,x+2) < c_b
                                            if im(y+3,x+-1) > cb
                                                continue;
                                            elseif im(y+3,x+-1) < c_b
                                                if im(y+-3,x+-1) > cb
                                                    continue;
                                                elseif im(y+-3,x+-1) < c_b
                                                    if im(y+3,x+1) < c_b
                                                        if im(y+-3,x+1) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    if im(y+2,x+-2) < c_b
                                                        if im(y+3,x+1) < c_b
                                                            if im(y+-3,x+1) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                if im(y+-2,x+-2) < c_b
                                                    if im(y+-3,x+1) < c_b
                                                        if im(y+-3,x+-1) < c_b
                                                            if im(y+3,x+1) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    end
                else
                    if im(y+-1,x+-3) < c_b
                        if im(y+0,x+3) > cb
                            if im(y+-1,x+3) < c_b
                                if im(y+3,x+-1) < c_b
                                    if im(y+1,x+-3) < c_b
                                        if im(y+2,x+-2) < c_b
                                            if im(y+0,x+-3) < c_b
                                                if im(y+-2,x+2) < c_b
                                                    if im(y+-3,x+1) < c_b
                                                        if im(y+-3,x+-1) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        elseif im(y+0,x+3) < c_b
                            if im(y+1,x+3) > cb
                                if im(y+2,x+-2) < c_b
                                    if im(y+0,x+-3) < c_b
                                        if im(y+-3,x+-1) < c_b
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            elseif im(y+1,x+3) < c_b
                                if im(y+-2,x+2) < c_b
                                    if im(y+3,x+1) > cb
                                        if im(y+0,x+-3) < c_b
                                            if im(y+1,x+-3) > cb
                                                continue;
                                            elseif im(y+1,x+-3) < c_b
                                                if im(y+-1,x+3) < c_b
                                                else
                                                    continue;
                                                end
                                            else
                                                if im(y+2,x+2) < c_b
                                                else
                                                    continue;
                                                end
                                            end
                                        else
                                            continue;
                                        end
                                    elseif im(y+3,x+1) < c_b
                                        if im(y+-3,x+1) < c_b
                                            if im(y+2,x+2) > cb
                                                continue;
                                            elseif im(y+2,x+2) < c_b
                                                if im(y+-2,x+-2) < c_b
                                                    if im(y+-1,x+3) < c_b
                                                        if im(y+-3,x+-1) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                if im(y+1,x+-3) < c_b
                                                    if im(y+-2,x+-2) < c_b
                                                        if im(y+0,x+-3) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        if im(y+0,x+-3) < c_b
                                            if im(y+1,x+-3) > cb
                                                if im(y+2,x+2) < c_b
                                                else
                                                    continue;
                                                end
                                            elseif im(y+1,x+-3) < c_b
                                                if im(y+-1,x+3) < c_b
                                                    if im(y+-3,x+1) < c_b
                                                        if im(y+-2,x+-2) < c_b
                                                            if im(y+-3,x+-1) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                if im(y+2,x+2) < c_b
                                                    if im(y+-3,x+1) < c_b
                                                        if im(y+-1,x+3) < c_b
                                                            if im(y+-3,x+-1) < c_b
                                                                if im(y+-2,x+-2) < c_b
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        else
                                            continue;
                                        end
                                    end
                                else
                                    continue;
                                end
                            else
                                if im(y+2,x+-2) < c_b
                                    if im(y+0,x+-3) < c_b
                                        if im(y+-2,x+2) < c_b
                                            if im(y+1,x+-3) < c_b
                                                if im(y+-2,x+-2) < c_b
                                                    if im(y+-1,x+3) < c_b
                                                        if im(y+-3,x+1) < c_b
                                                            if im(y+-3,x+-1) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        else
                            if im(y+3,x+-1) < c_b
                                if im(y+-1,x+3) < c_b
                                    if im(y+0,x+-3) < c_b
                                        if im(y+2,x+-2) < c_b
                                            if im(y+1,x+-3) < c_b
                                                if im(y+-3,x+1) < c_b
                                                    if im(y+-2,x+2) < c_b
                                                        if im(y+-3,x+-1) < c_b
                                                            if im(y+-2,x+-2) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        end
                    else
                        continue;
                    end
                end
            else
                if im(y+1,x+-3) > cb
                    if im(y+1,x+3) > cb
                        if im(y+3,x+-1) > cb
                            if im(y+-1,x+-3) > cb
                                if im(y+-3,x+-1) > cb
                                    if im(y+2,x+-2) > cb
                                        if im(y+0,x+-3) > cb
                                            if im(y+3,x+0) > cb
                                                if im(y+2,x+2) > cb
                                                    if im(y+-2,x+-2) > cb
                                                        if im(y+3,x+1) > cb
                                                        else
                                                            continue;
                                                        end
                                                    elseif im(y+-2,x+-2) < c_b
                                                        continue;
                                                    else
                                                        if im(y+0,x+3) > cb
                                                            if im(y+-1,x+3) > cb
                                                                if im(y+3,x+1) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        elseif im(y+0,x+-3) < c_b
                                            continue;
                                        else
                                            if im(y+-3,x+1) > cb
                                                if im(y+0,x+3) > cb
                                                    if im(y+-2,x+2) > cb
                                                        if im(y+3,x+0) > cb
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                elseif im(y+-3,x+-1) < c_b
                                    if im(y+-1,x+3) > cb
                                        if im(y+3,x+0) > cb
                                            if im(y+2,x+2) > cb
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    if im(y+0,x+3) > cb
                                        if im(y+3,x+1) > cb
                                            if im(y+2,x+-2) > cb
                                                if im(y+0,x+-3) > cb
                                                    if im(y+2,x+2) > cb
                                                        if im(y+-1,x+3) > cb
                                                            if im(y+3,x+0) > cb
                                                            else
                                                                continue;
                                                            end
                                                        elseif im(y+-1,x+3) < c_b
                                                            if im(y+-2,x+-2) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            if im(y+-2,x+-2) > cb
                                                                if im(y+3,x+0) > cb
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                elseif im(y+0,x+-3) < c_b
                                                    continue;
                                                else
                                                    if im(y+-3,x+1) > cb
                                                        if im(y+-2,x+2) > cb
                                                            if im(y+2,x+2) > cb
                                                                if im(y+-1,x+3) > cb
                                                                    if im(y+3,x+0) > cb
                                                                    else
                                                                        continue;
                                                                    end
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            elseif im(y+-1,x+-3) < c_b
                                if im(y+-3,x+1) > cb
                                    if im(y+-1,x+3) > cb
                                        if im(y+0,x+3) > cb
                                            if im(y+-2,x+2) > cb
                                                if im(y+3,x+1) > cb
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                elseif im(y+-3,x+1) < c_b
                                    continue;
                                else
                                    if im(y+-2,x+2) > cb
                                        if im(y+0,x+-3) > cb
                                            if im(y+0,x+3) > cb
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                end
                            else
                                if im(y+-2,x+2) > cb
                                    if im(y+0,x+3) > cb
                                        if im(y+-3,x+1) > cb
                                            if im(y+3,x+1) > cb
                                                if im(y+-1,x+3) > cb
                                                    if im(y+2,x+-2) > cb
                                                        if im(y+2,x+2) > cb
                                                            if im(y+3,x+0) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        elseif im(y+-3,x+1) < c_b
                                            continue;
                                        else
                                            if im(y+0,x+-3) > cb
                                                if im(y+3,x+0) > cb
                                                    if im(y+-1,x+3) > cb
                                                        if im(y+2,x+-2) > cb
                                                            if im(y+2,x+2) > cb
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            end
                        else
                            continue;
                        end
                    else
                        continue;
                    end
                elseif im(y+1,x+-3) < c_b
                    if im(y+0,x+3) > cb
                        if im(y+1,x+3) < c_b
                            if im(y+-3,x+-1) < c_b
                                if im(y+-1,x+-3) < c_b
                                    if im(y+3,x+1) < c_b
                                        if im(y+-3,x+1) > cb||im(y+-3,x+1) < c_b
                                        else
                                            if im(y+-2,x+-2) < c_b
                                                if im(y+3,x+-1) < c_b
                                                    if im(y+2,x+2) < c_b
                                                        if im(y+0,x+-3) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    elseif im(y+0,x+3) < c_b
                        if im(y+3,x+-1) < c_b
                            if im(y+-1,x+-3) > cb
                                if im(y+-3,x+1) < c_b
                                    if im(y+-2,x+2) < c_b
                                        if im(y+2,x+2) < c_b
                                            if im(y+-1,x+3) < c_b
                                                if im(y+2,x+-2) < c_b
                                                    if im(y+1,x+3) < c_b
                                                        if im(y+3,x+0) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            elseif im(y+-1,x+-3) < c_b
                                if im(y+2,x+2) < c_b
                                    if im(y+2,x+-2) < c_b
                                        if im(y+0,x+-3) > cb
                                            continue;
                                        elseif im(y+0,x+-3) < c_b
                                            if im(y+3,x+0) < c_b
                                                if im(y+-1,x+3) > cb
                                                    if im(y+1,x+3) < c_b
                                                    else
                                                        continue;
                                                    end
                                                elseif im(y+-1,x+3) < c_b
                                                    if im(y+1,x+3) < c_b
                                                        if im(y+3,x+1) < c_b
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    if im(y+-2,x+-2) < c_b
                                                        if im(y+1,x+3) < c_b
                                                            if im(y+3,x+1) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            if im(y+-3,x+1) < c_b
                                                if im(y+-1,x+3) < c_b
                                                    if im(y+1,x+3) < c_b
                                                        if im(y+-2,x+2) < c_b
                                                            if im(y+3,x+0) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                if im(y+-2,x+2) < c_b
                                    if im(y+-3,x+1) > cb
                                    elseif im(y+-3,x+1) < c_b
                                        if im(y+1,x+3) < c_b
                                            if im(y+3,x+1) < c_b
                                                if im(y+2,x+-2) < c_b
                                                    if im(y+-1,x+3) < c_b
                                                        if im(y+3,x+0) < c_b
                                                            if im(y+2,x+2) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        if im(y+0,x+-3) < c_b
                                            if im(y+2,x+2) < c_b
                                                if im(y+2,x+-2) < c_b
                                                    if im(y+1,x+3) < c_b
                                                        if im(y+3,x+0) < c_b
                                                            if im(y+-1,x+3) < c_b
                                                                if im(y+3,x+1) < c_b
                                                                else
                                                                    continue;
                                                                end
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    end
                                else
                                    continue;
                                end
                            end
                        else
                            continue;
                        end
                    else
                        if im(y+-3,x+-1) < c_b
                            if im(y+1,x+3) < c_b
                                if im(y+3,x+-1) < c_b
                                    if im(y+0,x+-3) < c_b
                                        if im(y+-2,x+-2) < c_b
                                            if im(y+2,x+2) < c_b
                                                if im(y+-1,x+-3) < c_b
                                                    if im(y+2,x+-2) < c_b
                                                        if im(y+3,x+1) < c_b
                                                            if im(y+3,x+0) < c_b
                                                            else
                                                                continue;
                                                            end
                                                        else
                                                            continue;
                                                        end
                                                    else
                                                        continue;
                                                    end
                                                else
                                                    continue;
                                                end
                                            else
                                                continue;
                                            end
                                        else
                                            continue;
                                        end
                                    else
                                        continue;
                                    end
                                else
                                    continue;
                                end
                            else
                                continue;
                            end
                        else
                            continue;
                        end
                    end
                else
                    continue;
                end
            end
			nc = nc + 1;														
			if nc > length(cs)													
				cs(length(cs)*2,1) = 0;											
			end																	
			cs(nc,1) = x;														
			cs(nc,2) = y;														
		end																		
	end																			
	coords = cs([1:nc],:);														
