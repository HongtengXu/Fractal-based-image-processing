function [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS, AL, BT, GM)
    % =====================================================================
    % File name   : icbi.m
    % File Type   : m-file (script file for Matlab or Octave)
    % Begin       : 2007-10-01
    % Last Update : 2008-01-30
    % Author      : Andrea Giachetti, Nicola Asuni
    % Description : ICBI (Iteractive Curvature Based Interpolation)
    %               This function returns an enlarged image by a factor 2^N
    %               implements the enlargement methods FCBI and ICBI described in the 
    %               paper "Fast artifacts-free image interpolation", proceedings of BMVC 2008
    %               
    %
    % Copyright   : Andrea Giachetti, via 24 Maggio 17 38100 Trento, Italy
    %               Nicola Asuni, nicola.asuni@tecnick.com
    % License     : GNU GENERAL PUBLIC LICENSE v.2
    %               http://www.gnu.org/copyleft/gpl.html
    % Version     : 1.1.000
    % =====================================================================    
    %
    % DESCRIPTION
    % --------------------
    % ICBI (Iteractive Curvature Based Interpolation)
    % This function returns an enlarged image by a factor 2^N implements
    % the enlargement methods FCBI and ICBI described in the paper
    % "Curvature-Driven Image Interpolation" submitted to SIGGRAPH 2008.
    %
    % KEYWORDS
    % --------------------
    % ICBI, image, zooming, magnification, upsizing, resampling,
    % resolution enhancement, interpolation,covariance-based adaptation,
    % geometric regularity, matlab, octave.
    %
    % USAGE
    % --------------------
    % [EI] = icbi(IM)
    % [EI] = icbi(IM, ZK)
    % [EI] = icbi(IM, ZK, SZ)
    % [EI] = icbi(IM, ZK, SZ, PF)
    % [EI] = icbi(IM, ZK, SZ, PF, VR)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS, AL)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS, AL, BT)
    % [EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS, AL, BT, GM)
    %
    %
    % INPUT
    % --------------------
    % IM : Source image.
    % ZK : Power of the zoom factor (default = 1)
    %      the image enlargement on vertical and horizontal direction is
    %      2^ZK; the final image size will be (SIZE * 2^ZK) - (2^ZK - 1).
    % SZ : Number of image bits per layer (default = 8).
    % PF : Potential to be minimized (default = 1).
    % VR : Verbose mode, if true prints some information during calculation
    %      (default = false).
    % ST : Maximum number of iterations (default = 20).
    % TM : Maximum edge step (default = 100).
    % TC : Edge continuity threshold (deafult = 50).
    % SC : Stopping criterion: 1 = change under threshold, 0 = ST iterations (default = 1).
    % TS : Threshold on image change for stopping iterations (default = 100).
    % AL : Weight for Curvature Continuity energy (default = 1.0).
    % BT : Weight for Curvature enhancement energy (default = -1.0).
    % GM : Weight for Isophote smoothing energy (default = 5.0).
    %
    %
    % OUTPUT
    % --------------------
    % EI : Enlarged image.
    %
    %
    % Examples
    % --------------------
    % Please check the icbiexample.m file on how to use this function.
    %
    %
    % Notes
    % --------------------
    % This implementation is not intended to be used in a production
    % environment. The main purpose of this script is to clearly show how
    % this technique works. Better performaces could be obtained using a
    % compiled version or rewriting this technique using a low-level
    % programming language.
    %
    %
    % ---------------------------------------------------------------------
    
    
    % Some initial tests on the input arguments
    
    if (nargin < 1)
        disp('ICBI (Iteractive Curvature Based Interpolation) function.');
        disp('This function returns an enlarged image.');
        disp('Usage:');
        disp('[EI] = icbi(IM)');
        disp('[EI] = icbi(IM, ZK)');
        disp('[EI] = icbi(IM, ZK, SZ)');
        disp('[EI] = icbi(IM, ZK, SZ, PF)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS, AL)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS, AL, BT)');
        disp('[EI] = icbi(IM, ZK, SZ, PF, VR, ST, TM, TC, SC, TS, AL, BT, GM)');
        disp('Where:');
        disp('IM : Source image.');
        disp('ZK : Power of the zoom factor (default = 1); the image enlargement on vertical and horizontal direction is 2^ZK; the final image size will be (SIZE * 2^ZK) - (2^ZK - 1).');
        disp('SZ : Number of image bits per layer (default = 8).');
        disp('PF : Potential to be minimized (default = 1).');
        disp('VR : Verbose mode, if true prints some information during calculation (default = false).');
        disp('ST : Maximum number of iterations (default = 20).');
        disp('TM : Maximum edge step (default = 100).');
        disp('TC : Edge continuity threshold (deafult = 50).');
        disp('SC : Stopping criterion: 1 = change under threshold, 0 = ST iterations (default = 1).');
        disp('TS : Threshold on image change for stopping iterations (default = 100).');
        disp('AL : Weight for Curvature Continuity energy (default = 1.0).');
        disp('BT : Weight for Curvature enhancement energy (default = -1.0).');
        disp('GM : Weight for Isophote smoothing energy (default = 5.0).');
        EI = [];
        return;
    end
    
    % assign default values
    if (nargin > 13)
        error('Too many arguments');
    end
    if (nargin < 13)
        GM = 5.0;
    end
    if (nargin < 12)
        BT = -1.0;
    end
    if (nargin < 11)
        AL = 1.0;
    end
    if (nargin < 10)
        TS = 100;
    end
    if (nargin < 9)
        SC = 1;
    end
    if (nargin < 8)
        TC = 50;
    end
    if (nargin < 7)
        TM = 100;
    end
    if (nargin < 6)
        ST = 20;
    end
    if (nargin < 5)
        VR = false;
    end
    if (nargin < 4)
        PF = 1;
    end
    if (nargin < 3)
        SZ = 8;
    end
    if (nargin < 2)
        ZK = 1;
    end  
    
    % --------------------------------------------------------------------- 
    
    % check parameters
    
    if (ZK < 1)
        EI=imresize(IM,2^ZK);%error('ZK must be a positive integer');
%         EI=EI(1:end-1,1:end-1,:);
    else
        ZK = int8(ZK);
    end
    SZ = int8(SZ);
    
    IDIM = ndims(IM);
    if (IDIM == 3)
     % number of colors
        CL = size(IM,3);
    elseif (IDIM == 2)
        CL = 1;
    else
        error('Unrecognized image type, please use RGB or grayscale images');
    end
        
    % ---------------------------------------------------------------------
    
    % store time for verbose mode
    if(VR)
        t1 = cputime;
    end
    
    % get image size
    [m,n] = size(IM(:,:,1));
    % calculate final image size;
    fm = (m * double(2^ZK)) - (double(2^ZK) - 1);
    fn = (n * double(2^ZK)) - (double(2^ZK) - 1);
    
    % initialize output image
    if (SZ > 32)
        EI = zeros(fm,fn,CL,'uint64');
    elseif (SZ > 16)
        EI = zeros(fm,fn,CL,'uint32');
    elseif (SZ > 8)
        EI = zeros(fm,fn,CL,'uint16');
    else
        EI = zeros(fm,fn,CL,'uint8');
    end

    % for each image color (for each color layer)
    for CID = 1:CL
        
        if(VR)
            fprintf('\n[%8.3f sec] LAYER: %02d\n', cputime-t1, CID);
        end
        
        % convert image type to double to improve interpolation accuracy
        IMG = double(IM(:,:,CID));

        % the image is enlarged by scaling factor (2^ZK - 1) at each cycle
        for ZF = 1:ZK
            
            if(VR)
                fprintf('[%8.3f sec]    ZF: %02d\n', cputime-t1, ZF);
            end
            
            % get image size
            [m,n] = size(IMG);
            
            % size of the expanded image
            mm = 2*m - 1; % rows
            nn = 2*n - 1; % columns
            
            % initialize expanded image and support matrices with zeros
            IMGEXP = zeros(mm,nn);
            D1 = zeros(mm,nn);
            D2 = zeros(mm,nn);
            D3 = zeros(mm,nn);
            C1 = zeros(mm,nn);
            C2 = zeros(mm,nn);
            
            % copy the low resolution grid on the high resolution grid
            
            IMGEXP(1:2:end,1:2:end) = IMG;
             

            % interpolation at borders (average value of 2 neighbors)
            for i = 2:2:mm-1
                % left column
                IMGEXP(i,1) = (IMGEXP(i-1,1) + IMGEXP(i+1,1))/2;
                %right column
                IMGEXP(i,nn) = (IMGEXP(i-1,nn) + IMGEXP(i+1,nn))/2;
            end
            for i = 2:2:nn
                % top row
                IMGEXP(1,i) = (IMGEXP(1,i-1) + IMGEXP(1,i+1))/2;
                % bottom row
                IMGEXP(mm,i) = (IMGEXP(mm,i-1) + IMGEXP(mm,i+1))/2;
            end
                        
            % Calculate interpolated points in two steps.
            % s = 0 calculates on diagonal directions.
            % s = 1 calculates on vertical and horizondal directions.
            for s = 0:1
                
                if(VR)
                    fprintf('[%8.3f sec]        PHASE: %02d\n', cputime-t1, s);
                end
                
                % ---------------------------------------------------------

                % FCBI (Fast Curvature Based Interpolation)
                % We compute second order derivatives in the opposite
                % directions and interpolate the two opposite neighbors in the
                % direction where the curvature is lower.
                for i = 2:2-s:mm-s
                    for j = 2+(s*(1-mod(i,2))):2:nn-s
                        v1 = abs(IMGEXP(i-1,j-1+s)-IMGEXP(i+1,j+1-s));
                        v2 = abs(IMGEXP(i+1-s,j-1)-IMGEXP(i-1+s,j+1));
                        p1 = (IMGEXP(i-1,j-1+s)+IMGEXP(i+1,j+1-s))/2;
                        p2 = (IMGEXP(i+1-s,j-1)+IMGEXP(i-1+s,j+1))/2;
                        if( (v1 < TM) && (v2 < TM) && (i > 3-s) && (i < mm-3-s) && (j > 3-s) && (j < nn-3-s) && (abs(p1-p2) < TM) )
                            if( abs(IMGEXP(i-1-s,j-3+(2*s)) + IMGEXP(i-3+s,j-1+(2*s)) + IMGEXP(i+1+s,j+3-(2*s)) + IMGEXP(i+3-s,j+1-(2*s)) + (2 * p2) - (6 * p1)) > ...
                                    abs(IMGEXP(i-3+(2*s),j+1+s) + IMGEXP(i-1+(2*s),j+3-s) + IMGEXP(i+3-(2*s),j-1-s) + IMGEXP(i+1-(2*s),j-3+s) + (2 * p1) - (6 * p2)) )
                                IMGEXP(i,j) = p1;                    
                            else
                                IMGEXP(i,j) = p2;                    
                            end
                        else
                            if( v1 < v2)
                                IMGEXP(i,j) = p1;                   
                            else
                                IMGEXP(i,j) = p2;                   
                            end               
                        end
                    end
                end
                
                step = 4.0 / (1 + s);
                
                
                % iterative refinement
                for g = 1:ST
                    
                    diff = 0;
                    
                    if g < ST/4
                        step = 1;
                    elseif g < ST/2
                        step = 2;
                    elseif g < 3*ST/4
                        step = 2;
                    end


                    % computation of derivatives
                    for i = 4-(2*s):1:mm-3+s
                        for j = 4-(2*s)+((1-s)*mod(i,2)):2-s:nn-3+s
                            C1(i,j) = (IMGEXP(i-1+s,j-1) - IMGEXP(i+1-s,j+1))/2;
                            C2(i,j) = (IMGEXP(i+1-(2*s),j-1+s) - IMGEXP(i-1+(2*s),j+1-s))/2;
                            D1(i,j) = IMGEXP(i-1+s,j-1) + IMGEXP(i+1-s,j+1) - 2*IMGEXP(i,j);
                            D2(i,j) = IMGEXP(i+1,j-1+s) + IMGEXP(i-1,j+1-s) - 2*IMGEXP(i,j);
                            D3(i,j) = (IMGEXP(i-s,j-2+s) - IMGEXP(i-2+s,j+s) + IMGEXP(i+s,j+2-s) - IMGEXP(i+2-s,j-s))/2;
                        end
                    end
                                        
                    for i = 6-(3*s):2-s:mm-5+(3*s)
                        for j = 6+(s*(mod(i,2)-3)):2:nn-5+(3*s)
                            
                            c_1 = 1;
                            c_2 = 1;
                            c_3 = 1;
                            c_4 = 1;
                            
                            if(abs(IMGEXP(i+1-s,j+1) - IMGEXP(i,j)) > TC) 
                                c_1 = 0;
                            end
                            if(abs(IMGEXP(i-1+s,j-1) - IMGEXP(i,j)) > TC) 
                                c_2 = 0;
                            end
                            if(abs(IMGEXP(i+1,j-1+s) - IMGEXP(i,j)) > TC) 
                                c_3 = 0;
                            end
                            if(abs(IMGEXP(i-1,j+1-s) - IMGEXP(i,j)) > TC) 
                                c_4 = 0;
                            end
                            
                            EN1 = ( c_1*abs(D1(i,j)-D1(i+1-s,j+1)) + c_2*abs(D1(i,j)-D1(i-1+s,j-1)));
                            EN2 = ( c_3*abs(D1(i,j)-D1(i+1,j-1+s)) + c_4*abs(D1(i,j)-D1(i-1,j+1-s)));
                            EN3 = ( c_1*abs(D2(i,j)-D2(i+1-s,j+1)) + c_2*abs(D2(i,j)-D2(i-1+s,j-1)));
                            EN4 = ( c_3*abs(D2(i,j)-D2(i+1,j-1+s)) + c_4*abs(D2(i,j)-D2(i-1,j+1-s)));
                            EN5 = abs(IMGEXP(i-2+(2*s),j-2) + IMGEXP(i+2-(2*s),j+2) - 2*IMGEXP(i,j));
                            EN6 = abs(IMGEXP(i+2,j-2+(2*s)) + IMGEXP(i-2,j+2-(2*s)) - 2*IMGEXP(i,j));

                            EA1 = (c_1*abs(D1(i,j)-D1(i+1-s,j+1)- 3*step) + c_2*abs(D1(i,j)-D1(i-1+s,j-1)-3*step));
                            EA2 = (c_3*abs(D1(i,j)-D1(i+1,j-1+s)- 3*step) + c_4*abs(D1(i,j)-D1(i-1,j+1-s)-3*step));
                            EA3 = (c_1*abs(D2(i,j)-D2(i+1-s,j+1)- 3*step) + c_2*abs(D2(i,j)-D2(i-1+s,j-1)-3*step));
                            EA4 = (c_3*abs(D2(i,j)-D2(i+1,j-1+s)- 3*step) + c_4*abs(D2(i,j)-D2(i-1,j+1-s)-3*step));
                            EA5 = abs(IMGEXP(i-2+(2*s),j-2) + IMGEXP(i+2-(2*s),j+2) - 2*IMGEXP(i,j) -2*step);
                            EA6 = abs(IMGEXP(i+2,j-2+(2*s)) + IMGEXP(i-2,j+2-(2*s)) - 2*IMGEXP(i,j) -2*step);
                            
                            ES1 = (c_1*abs(D1(i,j)-D1(i+1-s,j+1)+3*step) + c_2*abs(D1(i,j)-D1(i-1+s,j-1)+3*step));
                            ES2 = (c_3*abs(D1(i,j)-D1(i+1,j-1+s)+3*step) + c_4*abs(D1(i,j)-D1(i-1,j+1-s)+3*step));
                            ES3 = (c_1*abs(D2(i,j)-D2(i+1-s,j+1)+3*step) + c_2*abs(D2(i,j)-D2(i-1+s,j-1)+3*step));
                            ES4 = (c_3*abs(D2(i,j)-D2(i+1,j-1+s)+3*step) + c_4*abs(D2(i,j)-D2(i-1,j+1-s)+3*step));
                            ES5 = abs(IMGEXP(i-2+(2*s),j-2) + IMGEXP(i+2-(2*s),j+2) - 2*IMGEXP(i,j) +2*step);
                            ES6 = abs(IMGEXP(i+2,j-2+(2*s)) + IMGEXP(i-2,j+2-(2*s)) - 2*IMGEXP(i,j) +2*step);
                            
                            EISO = (C1(i,j)*C1(i,j)*D2(i,j) -2*C1(i,j)*C2(i,j)*D3(i,j) + C2(i,j)*C2(i,j)*D1(i,j))/(C1(i,j)*C1(i,j)+C2(i,j)*C2(i,j));

                            if(abs(EISO) < 0.2)
                                EISO=0;
                            end

                            if(PF == 1)
                                EN = (AL * (EN1 + EN2 + EN3 + EN4)) + (BT * (EN5 + EN6));
                                EA = (AL * (EA1 + EA2 + EA3 + EA4)) + (BT * (EA5 + EA6));
                                ES = (AL * (ES1 + ES2 + ES3 + ES4)) + (BT * (ES5 + ES6));
                            elseif(PF == 2)
                                EN = (AL * (EN1 + EN2 + EN3 + EN4));
                                EA = (AL * (EA1 + EA2 + EA3 + EA4)) - (GM * sign(EISO));
                                ES = (AL * (ES1 + ES2 + ES3 + ES4)) + (GM * sign(EISO));
                            else
                                EN = (AL * (EN1 + EN2 + EN3 + EN4)) + (BT * (EN5 + EN6));
                                EA = (AL * (EA1 + EA2 + EA3 + EA4)) + (BT * (EA5 + EA6)) - (GM * sign(EISO));
                                ES = (AL * (ES1 + ES2 + ES3 + ES4)) + (BT * (ES5 + ES6)) + (GM * sign(EISO));
                            end

                            if((EN > EA) && (ES >EA))
                                IMGEXP(i,j) = IMGEXP(i,j) + step;
                                diff = diff + step;
                            elseif((EN > ES) && (EA >ES))
                                IMGEXP(i,j) = IMGEXP(i,j) - step;
                                diff = diff + step;
                            end
                        end
                    end
                    
                    if((SC==1) && (diff < TS))
                        break;
                    end
                    
                end % end of iterative refinement
 
            end % end of (for s = 0:1)
            
            % assign th expanded image to the current image
            IMG = IMGEXP;
            
        end % end of (for ZF = 1:ZK)
        
        % store this color layer to the output image
        EI(:,:,CID) = round(IMG);

    end
    
    if(VR)
        fprintf('[%8.3f sec] END\n', cputime-t1);
    end

    % === EOF =================================================================
