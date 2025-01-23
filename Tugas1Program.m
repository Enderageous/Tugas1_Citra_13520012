classdef Tugas1Program
    methods
        function histogram = calculateGrayscaledHistogram(~, image)
            [row, col] = size(image); 
            histogram = zeros(1, 256);
            for i = 1:row
                for j = 1:col
                    intensity = image(i, j);
                    histogram(intensity + 1) = histogram(intensity + 1) + 1;
                end
            end
        end
        
        function histogram = calculateColoredHistogram(~, image)
            [row, col, channel] = size(image);
            histogram = zeros(3, 256);
            for c = 1:channel
                for i = 1:row
                    for j = 1:col
                        intensity = image(i, j, c); 
                        histogram(c, intensity + 1) = histogram(c, intensity + 1) + 1; 
                    end
                end
            end
        end

        function result = imageNegative(~, image)
            result = 255 - image; 
        end

        function result = imageBrightening(~, image, a, b)
            result = a * image + b;
        end

        function result = logTransformation(~, image, c)
            result = im2double(image);
            result = c * log(result + 1);
            result = im2uint8(result);
        end

        function result = powerTransformation(~, image, c, gamma)
            result = im2double(image);
            result = c * (result.^gamma);
            result = im2uint8(result);
        end

        function result = contrastStretching(~, image)
            rmin = min(image(:));
            rmax = max(image(:));
            result = (image - rmin) .* (255 / (rmax - rmin)); 
        end

        function result = histogramEqualization(self, image)
            imageHistogram = self.calculateGrayscaledHistogram(image);
            [row, col] = size(image);
            n = row * col;
            eqHistogram = zeros(1, 256);
            sum = 0;
            for i = 1:256
                sum = sum + imageHistogram(i);
                eqHistogram(i) = floor((double(sum) / n) * 255);
            end
            result = image;
            for i = 1:row
                for j = 1:col
                    result(i, j) = eqHistogram(image(i, j) + 1) ;
                end
            end
        end

        function result = histogramSpecification(self, image, ref)
            [row, col] = size(image);
            n = row * col;
            imageHistogram = self.calculateGrayscaledHistogram(image); 
            imageEqHist = zeros(1, 256); 
            sum = 0;
            for i = 1:256
                sum = sum + imageHistogram(i);
                imageEqHist(i) = floor((double(sum) / n) * 255);
            end
            refHistogram = self.calculateGrayscaledHistogram(ref);
            refEqHist = zeros(1, 256); 
            sum = 0;
            for i = 1:256
                sum = sum + refHistogram(i);
                refEqHist(i) = floor((double(sum) / n) * 255);
            end
            invHist = zeros(1, 256);
            for i = 1:row
                minVal = abs(imageEqHist(uint8(max(i/(row/256), 1))) - refEqHist(1));
                minJ = 0;
                for j = 1:256
                    if (abs(imageEqHist(uint8(max(i/(row/256), 1))) - refEqHist(j)) < minVal)
                        minVal = abs(imageEqHist(uint8(max(i/(row/256), 1))) - refEqHist(j));
                        minJ = j;
                    end
                end
                invHist(i) = minJ;
            end
            result = image;
            for i = 1:row
                for j = 1:col
                    result(i, j) = invHist(image(i, j) + 1) ;
                end
            end
        end
    end
end