function [IM] = scale12to16bit(IM)
%The Hamamatsu cameras for the closet scope and curtain scope create images
%with 12-bit dynamic range. However, the TIFF format that stores these
%images uses a 16-bit format. Viewing a 12-bit image in a 16-bit format on
%a computer monitor is complicated by the scaling being done at 16-bit. To
%make viewing images from the microscope easier on a computer monitor,
%without any compression or loss of data, the 12-bit data is shifted left
%4-bits to become 16-bit data.
IM = bitshift(IM,4);
end