#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>

class CircleSplatting {
   private:
    std::vector<uint8_t> originalImageBuffer;
    std::vector<uint8_t> drawBuffer;
    int width;
    int height;
    std::mt19937 rng;

   public:
    CircleSplatting(int w, int h) : width(w), height(h) {
        rng.seed(std::random_device()());
        const size_t bufferSize = w * h * 4;

        originalImageBuffer.resize(bufferSize, 255);
        drawBuffer.resize(bufferSize, 255);

        std::cout << "[WASM] Crated instance for " << w << "x" << h
                  << std::endl;
    }

    emscripten::val getInputBuffer() {
        return emscripten::val(emscripten::typed_memory_view(
            originalImageBuffer.size(), originalImageBuffer.data()));
    }

    void run(int numCircles, int numEpochs, int radiusMax) {
        std::memcpy(drawBuffer.data(), originalImageBuffer.data(),
                    drawBuffer.size());
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                const int idx = (y * width + x) * 4;
                if (x < numCircles) {
                    drawBuffer[idx] = 255;
                }
                if (y < numEpochs) {
                    drawBuffer[idx + 1] = 255;
                }
            }
        }
    }
    emscripten::val getDrawImageData() {
        return emscripten::val(emscripten::typed_memory_view(
            drawBuffer.size(), drawBuffer.data()));
    }
};

EMSCRIPTEN_BINDINGS(circle_splatting_module) {
    emscripten::class_<CircleSplatting>("CircleSplatting")
        .constructor<int, int>()
        .function("getInputBuffer", &CircleSplatting::getInputBuffer)
        .function("run", &CircleSplatting::run)
        .function("getDrawImageData", &CircleSplatting::getDrawImageData);
}
