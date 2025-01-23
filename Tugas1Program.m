classdef Tugas1Program
    methods
        function histogram = calculateGrayscaledHistogram(~, image)
            [row, col] = size(image); 
            histogram = zeros(1, 256); % Inisialisasi histogram 
            
            % Loop setiap pixel untuk mendapat intensity
            for i = 1:row
                for j = 1:col
                    intensity = image(i, j);
                    histogram(intensity + 1) = histogram(intensity + 1) + 1; % Tambah pada intensity yang sesuai
                end
            end
        end
        
        function histogram = calculateColoredHistogram(~, image)
            [row, col, channel] = size(image);
            histogram = zeros(3, 256); % Inisialisasi 3 histogram untuk tiap warna
            
            % Loop setiap pixel pada tiap channel warna untuk mendapat intensity
            for c = 1:channel
                for i = 1:row
                    for j = 1:col
                        intensity = image(i, j, c); 
                        histogram(c, intensity + 1) = histogram(c, intensity + 1) + 1; % Tambah pada intensity yang sesuai
                    end
                end
            end
        end

        function result = imageBrightening(~, image, a, b)
            result = a * image + b;
        end

        function result = imageNegative(~, image)
            result = 255 - image; % Membalikan intensitas piksel
        end

        function result = logTransformation(~, image, c)
            result = im2double(image); % Mengubah menjadi double
            result = c * log(result + 1); % Melakukan transformasi log
            result = im2uint8(result); % Mengembalikan ke uint8
        end

        function result = powerTransformation(~, image, c, gamma)
            result = im2double(image); % Mengubah menjadi double
            result = c * (result.^gamma); % Melakukan transformasi pangkat
            result = im2uint8(result); % Mengembalikan ke uint8
        end

        function result = contrastStretching(~, image)
            rmin = min(image(:)); % Mendapat intensitas minimal citra
            rmax = max(image(:)); % Mendapat intensitas maksimal citra
            result = (image - rmin) .* (255 / (rmax - rmin)); % Melakukan stretching
        end

        function result = histogramEqualization(self, image)
            imageHistogram = self.calculateGrayscaledHistogram(image); % Memperoleh histogram citra
            [row, col] = size(image);
            n = row * col;
            eqHistogram = zeros(1, 256); % Inisialisasi equalizer histogram
            
            % Mengisi equalizer histogram
            sum = 0;
            for i = 1:256
                sum = sum + imageHistogram(i);
                eqHistogram(i) = floor((double(sum) / n) * 255);
            end
            
            % Update citra sesuai dengan equalizer histogram
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
            
            imageHistogram = self.calculateGrayscaledHistogram(image); % Hitung histogram citra
            imageEqHist = zeros(1, 256); % Inisialisasi eqHist citra
            
            % Mengisi eqHist citra
            sum = 0;
            for i = 1:256
                sum = sum + imageHistogram(i);
                imageEqHist(i) = floor((double(sum) / n) * 255);
            end

            refHistogram = self.calculateGrayscaledHistogram(ref); % Hitung histogram referensi
            refEqHist = zeros(1, 256); % Inisialisasi eqHist referensi
            
            % Mengisi eqHist referensi
            sum = 0;
            for i = 1:256
                sum = sum + refHistogram(i);
                refEqHist(i) = floor((double(sum) / n) * 255);
            end
            
            % Transformasi balikan
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
            
            % Update citra
            result = image;
            for i = 1:row
                for j = 1:col
                    result(i, j) = invHist(image(i, j) + 1) ;
                end
            end
        end

        function result = conv(~, image, mask)
            [sizeMask, ~] = size(mask);
            disp(sizeMask);
            
            offset = floor(sizeMask / 2);
            padded_image = padarray(image, [offset, offset], 'replicate');
            [row, col] = size(padded_image);
            
            imageDouble = im2double(padded_image);
            ImageResult = imageDouble;
            
            for i = 1+offset:(row - offset)
                for j = 1+offset:(col - offset)
                    conv_sum = 0;
                    
                    for m = 1:sizeMask
                        for n = 1:sizeMask
                            conv_sum = conv_sum + imageDouble(i + m - 1 - offset, j + n - 1 - offset) * mask(m, n);
                        end
                    end
                    
                    ImageResult(i, j) = conv_sum;
                end
            end
            original_size_result = ImageResult(offset+1:end-offset, offset+1:end-offset, :);
            result = im2uint8(original_size_result);
        end
    end
end