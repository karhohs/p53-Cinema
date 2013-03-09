function [IM] = scale12to16bit(IM)
%input:
%IM = an image with bit depth 12
%
%output:
%IM = a scaled image with bit depth 16
%description:
%The class of the image is detected. Depending on the class type an image
%maybe converted into an integer format and then converted back into the
%format of the input.
%
%other notes:
%The Hamamatsu cameras for the closet scope and curtain scope create images
%with 12-bit dynamic range. However, the TIFF format that stores these
%images uses a 16-bit format. Viewing a 12-bit image in a 16-bit format on
%a computer monitor is complicated by the scaling being done at 16-bit. To
%make viewing images from the microscope easier on a computer monitor,
%without any compression or loss of data, the 12-bit data is shifted left
%4-bits to become 16-bit data.
numType = class(IM);
switch numType
    case 'double'
        IM = uint16(IM);
        IM = bitshift(IM,4);
        IM = double(IM);
    case 'uint16'
        IM = bitshift(IM,4);
end
end