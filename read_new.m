function varargout = read_new(obj, varargin)
%READ Read a video file. 
%
%   VIDEO = READ(OBJ) reads in video frames from the associated file.  VIDEO
%   is an H x W x B x F matrix where H is the image frame height, W is the
%   image frame width, B is the number of bands in the image (e.g. 3 for RGB),
%   and F is the number of frames read in.  The default behavior is to read in
%   all frames unless an index is specified.  The type of data returned is 
%   always UINT8 data representing RGB24 video frames.
%
%   VIDEO = READ(...,INDEX) performs the same operation, but reads only the
%   frame(s) specified by INDEX, where the first frame number is 1.  INDEX can
%   be a single index,  or a two-element array representing an index range 
%   of the video stream.
%
%   For example:
%
%      VIDEO = READ(OBJ);           % Read in all video frames.
%      VIDEO = READ(OBJ, 1);        % Read only the first frame.
%      VIDEO = READ(OBJ, [1 10]);   % Read the first 10 frames.
%
%   If any invalid INDEX is specified, MATLAB throws an error.
%
%   Example:
%      % Construct a multimedia reader object associated with file 'xylophone.mpg' with
%      % user tag set to 'myreader1'.
%      readerobj = VideoReader('xylophone.mpg', 'tag', 'myreader1');
%
%      % Read in all video frames.
%      vidFrames = read(readerobj);
%
%      % Get the number of frames.
%      numFrames = get(readerobj, 'numberOfFrames');
%
%      % Create a MATLAB movie struct from the video frames.
%      for k = 1 : numFrames
%            mov(k).cdata = vidFrames(:,:,:,k);
%            mov(k).colormap = [];
%      end
%
%      % Create a figure
%      hf = figure; 
%      
%      % Resize figure based on the video's width and height
%      set(hf, 'position', [150 150 readerobj.Width readerobj.Height])
%
%      % Playback movie once at the video's frame rate
%      movie(hf, mov, 1, readerobj.FrameRate);
%
%   See also AUDIOVIDEO, MOVIE, VIDEOREADER, VIDEOREADER/GET, VIDEOREADER/SET, MMFILEINFO.

%    NCH DTL
%    Copyright 2005-2011 The MathWorks, Inc.
%    $Revision: 1.1.6.5 $  $Date: 2011/10/15 01:50:42 $


if length(obj) > 1
    error(message('MATLAB:audiovideo:VideoReader:nonscalar'));
end

% ensure that we pass in 1 or 2 arguments only
error(nargchk(1, 2, nargin, 'struct'));

% Verify that the index argument is of numeric type
if nargin == 2
    index = varargin{1};
    validateattributes(index, {'numeric'}, {'vector'}, 'VideoReader.read', 'index');
    index = double(index);
end

try
    if nargin == 1
        frameIndex = [1 Inf];
        % Dimensions of data returned is HxWx3xN
        videoFrames = read(getImpl(obj));
    elseif nargin == 2
        % Dimensions of data returned is HxWx3xN
        videoFrames = read(getImpl(obj), index);
        if isscalar(index)
            frameIndex = [index index];
        else
            frameIndex = index;
        end
    end
catch err
    VideoReader.handleImplException(err);
end

% Update the frameIndex only if the total number of frames in the video can
% be determined after the read operation
if( ~isempty(obj(1).NumberOfFrames) )
    frameIndex(frameIndex == Inf) = obj(1).NumberOfFrames;
end

% Check that read was complete only if the frame indices to be read have
% been accurately determined.
if ~any(frameIndex == Inf)
    checkIncompleteRead(size(videoFrames, 4), frameIndex);
end 

% Video is the output argument.
varargout{1} = videoFrames;

end

function checkIncompleteRead(actNum, index)
expNum = index(2) - index(1) + 1;
if actNum < expNum
    warning(message('MATLAB:audiovideo:VideoReader:incompleteRead', index(1), index(1)+actNum-1));
end
end
