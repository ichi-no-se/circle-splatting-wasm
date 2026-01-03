#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <optional>
#include <random>
#include <vector>

class CircleSplatting {
  private:
	std::vector<uint8_t> originalImageBuffer;
	std::vector<uint8_t> drawBuffer;
	int width;
	int height;

	class XorShift {
	  public:
		XorShift(uint32_t seed = 123456789) : x(seed) {}

		uint32_t next() {
			x ^= x << 13;
			x ^= x >> 17;
			x ^= x << 5;
			return x;
		}

		uint32_t next(uint32_t n) { return next() % n; }

		float nextFloat() {
			return static_cast<float>(next()) / static_cast<float>(UINT32_MAX);
		}

		float nextFloat(float min, float max) {
			return min + (max - min) * nextFloat();
		}

		void setSeed(uint32_t seed) { x = seed; }

	  private:
		uint32_t x;
	};

	XorShift rng;

	struct Vec2i {
		int x;
		int y;
		int &operator[](int i) { return *(&x + i); }
		const int &operator[](int i) const { return *(&x + i); }
	};

	struct Vec2f {
		float x;
		float y;
		float &operator[](int i) { return *(&x + i); }
		const float &operator[](int i) const { return *(&x + i); }
	};

	struct Vec3f {
		float x;
		float y;
		float z;
		float &operator[](int i) { return *(&x + i); }
		const float &operator[](int i) const { return *(&x + i); }
		Vec3f &operator+=(const Vec3f &other) {
			x += other.x;
			y += other.y;
			z += other.z;
			return *this;
		}
		Vec3f operator+(const Vec3f &other) const {
			return Vec3f{x + other.x, y + other.y, z + other.z};
		}
		Vec3f &operator-=(const Vec3f &other) {
			x -= other.x;
			y -= other.y;
			z -= other.z;
			return *this;
		}
		Vec3f operator-(const Vec3f &other) const {
			return Vec3f{x - other.x, y - other.y, z - other.z};
		}
		Vec3f operator*(float scalar) const {
			return Vec3f{x * scalar, y * scalar, z * scalar};
		}
		Vec3f operator/(float scalar) const {
			return Vec3f{x / scalar, y / scalar, z / scalar};
		}
	};

	struct Circle {
		Vec2f center;
		float radius;
		Vec3f color;
	};

	Circle shiftCircle(const Circle &circle, const float shiftX,
					   const float shiftY, const float shiftRadius) {
		Circle shiftedCircle;
		shiftedCircle.center.x = circle.center.x + shiftX;
		shiftedCircle.center.y = circle.center.y + shiftY;
		shiftedCircle.radius = circle.radius + shiftRadius;
		for (int i = 0; i < 3; i++) {
			shiftedCircle.color[i] = circle.color[i];
		}
		return shiftedCircle;
	}

	std::optional<Vec2i> calcCircleRangeY(const Circle &circle, const int x,
										  const int imgHeight) {
		const float xFloat = static_cast<float>(x);
		const float dx = xFloat - circle.center.x;
		const float radiusSq = circle.radius * circle.radius;
		const float distSq = dx * dx;
		if (distSq > radiusSq) {
			return std::nullopt;
		}
		const float dy = std::sqrt(radiusSq - distSq);
		int yStart = static_cast<int>(std::ceil(circle.center.y - dy));
		int yEnd = static_cast<int>(std::floor(circle.center.y + dy));
		yStart = std::max(0, yStart);
		yEnd = std::min(imgHeight - 1, yEnd);
		return Vec2i{yStart, yEnd};
	}

	std::optional<Vec2i> calcCircleRangeX(const Circle &circle, const int y,
										  const int imgWidth) {
		const float yFloat = static_cast<float>(y);
		const float dy = yFloat - circle.center.y;
		const float radiusSq = circle.radius * circle.radius;
		const float distSq = dy * dy;
		if (distSq > radiusSq) {
			return std::nullopt;
		}
		const float dx = std::sqrt(radiusSq - distSq);
		int xStart = static_cast<int>(std::ceil(circle.center.x - dx));
		int xEnd = static_cast<int>(std::floor(circle.center.x + dx));
		xStart = std::max(0, xStart);
		xEnd = std::min(imgWidth - 1, xEnd);
		return Vec2i{xStart, xEnd};
	}

	float calcPixelLoss(const Vec3f &accum, const int count,
						const Vec3f &imageColor) {
		if (count == 0) {
			const int EMPTY_PENALTY = 20;
			return EMPTY_PENALTY;
		}
		const Vec3f currentColor = accum / static_cast<float>(count);
		const Vec3f diff = imageColor - currentColor;
		float loss = 0.0f;
		for (int i = 0; i < 3; i++) {
			loss += diff[i] * diff[i];
		}
		return loss;
	}

	float calcPixelLossDeltaAdd(
		const int x, const int y, const Vec3f &color,
		const std::vector<std::vector<Vec3f>> &image,
		const std::vector<std::vector<Vec3f>> &canvasColorsAccum,
		const std::vector<std::vector<int>> &canvasCounts) {
		const Vec3f &imageColor = image[y][x];
		const Vec3f &currentAccum = canvasColorsAccum[y][x];
		const int currentCount = canvasCounts[y][x];
		const float currentLoss =
			calcPixelLoss(currentAccum, currentCount, imageColor);
		const float newLoss =
			calcPixelLoss(currentAccum + color, currentCount + 1, imageColor);
		return newLoss - currentLoss;
	}

	float calcPixelLossDeltaRemove(
		const int x, const int y, const Vec3f &color,
		const std::vector<std::vector<Vec3f>> &image,
		const std::vector<std::vector<Vec3f>> &canvasColorsAccum,
		const std::vector<std::vector<int>> &canvasCounts) {
		const Vec3f &imageColor = image[y][x];
		const Vec3f &currentAccum = canvasColorsAccum[y][x];
		const int currentCount = canvasCounts[y][x];
		const float currentLoss =
			calcPixelLoss(currentAccum, currentCount, imageColor);
		const float newLoss =
			calcPixelLoss(currentAccum - color, currentCount - 1, imageColor);
		return newLoss - currentLoss;
	}

	float calcLossDeltaFromRanges(
		const std::optional<Vec2i> &currentXRangeOpt,
		const std::optional<Vec2i> &newXRangeOpt, const int y,
		const Vec3f &color, const std::vector<std::vector<Vec3f>> &image,
		const std::vector<std::vector<Vec3f>> &canvasColorsAccum,
		const std::vector<std::vector<int>> &canvasCounts) {
		float lossDelta = 0.0f;
		if (currentXRangeOpt.has_value()) {
			const Vec2i currentXRange = currentXRangeOpt.value();
			if (newXRangeOpt.has_value()) {
				const Vec2i newXRange = newXRangeOpt.value();
				if (newXRange[1] < currentXRange[0] ||
					currentXRange[1] < newXRange[0]) {
					// non-overlapping
					for (int x = currentXRange[0]; x <= currentXRange[1]; ++x) {
						lossDelta += calcPixelLossDeltaRemove(
							x, y, color, image, canvasColorsAccum,
							canvasCounts);
					}
					for (int x = newXRange[0]; x <= newXRange[1]; ++x) {
						lossDelta += calcPixelLossDeltaAdd(x, y, color, image,
														   canvasColorsAccum,
														   canvasCounts);
					}
				} else if (newXRange[0] <= currentXRange[0] &&
						   currentXRange[1] <= newXRange[1]) {
					// current is inside new
					for (int x = newXRange[0]; x < currentXRange[0]; ++x) {
						lossDelta += calcPixelLossDeltaAdd(x, y, color, image,
														   canvasColorsAccum,
														   canvasCounts);
					}
					for (int x = currentXRange[1] + 1; x <= newXRange[1]; ++x) {
						lossDelta += calcPixelLossDeltaAdd(x, y, color, image,
														   canvasColorsAccum,
														   canvasCounts);
					}
				} else if (currentXRange[0] <= newXRange[0] &&
						   newXRange[1] <= currentXRange[1]) {
					// new is inside current
					for (int x = currentXRange[0]; x < newXRange[0]; ++x) {
						lossDelta += calcPixelLossDeltaRemove(
							x, y, color, image, canvasColorsAccum,
							canvasCounts);
					}
					for (int x = newXRange[1] + 1; x <= currentXRange[1]; ++x) {
						lossDelta += calcPixelLossDeltaRemove(
							x, y, color, image, canvasColorsAccum,
							canvasCounts);
					}
				} else if (currentXRange[0] <= newXRange[0]) {
					// partial overlap, current left
					for (int x = currentXRange[0]; x < newXRange[0]; ++x) {
						lossDelta += calcPixelLossDeltaRemove(
							x, y, color, image, canvasColorsAccum,
							canvasCounts);
					}
					for (int x = currentXRange[1] + 1; x <= newXRange[1]; ++x) {
						lossDelta += calcPixelLossDeltaAdd(x, y, color, image,
														   canvasColorsAccum,
														   canvasCounts);
					}
				} else {
					// partial overlap, new left
					for (int x = newXRange[0]; x < currentXRange[0]; ++x) {
						lossDelta += calcPixelLossDeltaAdd(x, y, color, image,
														   canvasColorsAccum,
														   canvasCounts);
					}
					for (int x = newXRange[1] + 1; x <= currentXRange[1]; ++x) {
						lossDelta += calcPixelLossDeltaRemove(
							x, y, color, image, canvasColorsAccum,
							canvasCounts);
					}
				}
			} else {
				for (int x = currentXRange[0]; x <= currentXRange[1]; ++x) {
					lossDelta += calcPixelLossDeltaRemove(
						x, y, color, image, canvasColorsAccum, canvasCounts);
				}
			}
		} else {
			if (newXRangeOpt.has_value()) {
				const Vec2i newXRange = newXRangeOpt.value();
				for (int x = newXRange[0]; x <= newXRange[1]; ++x) {
					lossDelta += calcPixelLossDeltaAdd(
						x, y, color, image, canvasColorsAccum, canvasCounts);
				}
			} else {
				// do nothing
			}
		}
		return lossDelta;
	}

	void renderCircles(std::vector<Circle> &circles, const float scale) {
		const int imgWidth = width;
		const int imgHeight = height;
		std::vector<std::vector<Vec3f>> canvasColorsAccum(
			imgHeight, std::vector<Vec3f>(imgWidth, Vec3f{0.0f, 0.0f, 0.0f}));
		std::vector<std::vector<int>> canvasCounts(
			imgHeight, std::vector<int>(imgWidth, 0));
		for (auto circle : circles) {
			circle.center.x *= scale;
			circle.center.y *= scale;
			circle.radius += 0.5f;
			circle.radius *= scale;
			const int yStart = std::max(
				0,
				static_cast<int>(std::ceil(circle.center.y - circle.radius)));
			const int yEnd = std::min(
				imgHeight - 1,
				static_cast<int>(std::floor(circle.center.y + circle.radius)));
			for (int y = yStart; y <= yEnd; ++y) {
				const auto xRangeOpt = calcCircleRangeX(circle, y, imgWidth);
				if (!xRangeOpt.has_value()) {
					continue;
				}
				const Vec2i xRange = xRangeOpt.value();
				canvasColorsAccum[y][xRange[0]] += circle.color;
				canvasCounts[y][xRange[0]] += 1;
				if (xRange[1] + 1 < imgWidth) {
					canvasColorsAccum[y][xRange[1] + 1] -= circle.color;
					canvasCounts[y][xRange[1] + 1] -= 1;
				}
			}
		}
		for (int y = 0; y < imgHeight; ++y) {
			for (int x = 1; x < imgWidth; ++x) {
				canvasColorsAccum[y][x] += canvasColorsAccum[y][x - 1];
				canvasCounts[y][x] += canvasCounts[y][x - 1];
			}
		}

		for (int y = 0; y < imgHeight; ++y) {
			for (int x = 0; x < imgWidth; ++x) {
				const int idx = (y * imgWidth + x) * 4;
				if (canvasCounts[y][x] > 0) {
					const Vec3f color = canvasColorsAccum[y][x] /
										static_cast<float>(canvasCounts[y][x]);
					for (int c = 0; c < 3; ++c) {
						int val =
							static_cast<int>(std::round(color[c] * 255.0f));
						val = std::clamp(val, 0, 255);
						drawBuffer[idx + c] = static_cast<uint8_t>(val);
					}
					drawBuffer[idx + 3] = 255;
				} else {
					drawBuffer[idx + 0] = 0;
					drawBuffer[idx + 1] = 0;
					drawBuffer[idx + 2] = 0;
					drawBuffer[idx + 3] = 255;
				}
			}
		}
	}

  public:
	CircleSplatting(int w, int h) : width(w), height(h) {
		rng.setSeed(static_cast<uint32_t>(std::random_device()()));
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
		const float scaleFactor =
			200.0f / static_cast<float>(std::max(width, height));
		const int newWidth = static_cast<int>(std::round(width * scaleFactor));
		const int newHeight =
			static_cast<int>(std::round(height * scaleFactor));
		std::vector<std::vector<Vec3f>> image(
			newHeight, std::vector<Vec3f>(newWidth, Vec3f{0.0f, 0.0f, 0.0f}));
		for (int y = 0; y < newHeight; ++y) {
			for (int x = 0; x < newWidth; ++x) {
				const int origX =
					std::min(width - 1, static_cast<int>(x / scaleFactor));
				const int origY =
					std::min(height - 1, static_cast<int>(y / scaleFactor));

				const int idx = (origY * width + origX) * 4;
				for (int c = 0; c < 3; ++c) {
					image[y][x][c] =
						static_cast<float>(originalImageBuffer[idx + c]) /
						255.0f;
				}
			}
		}
		std::vector<Circle> circles(numCircles);
		for (auto &circle : circles) {
			circle.center =
				Vec2f{rng.nextFloat(0.0f, static_cast<float>(newWidth - 1)),
					  rng.nextFloat(0.0f, static_cast<float>(newHeight - 1))};
			circle.radius =
				std::exp(rng.nextFloat(std::log(1.0f), std::log(radiusMax)));
			circle.color = image[static_cast<int>(circle.center.y)]
								[static_cast<int>(circle.center.x)];
		}

		for (int iter = 0; iter < numEpochs; ++iter) {
			std::vector<std::vector<Vec3f>> canvasColorsAccum(
				newHeight,
				std::vector<Vec3f>(newWidth, Vec3f{0.0f, 0.0f, 0.0f}));
			std::vector<std::vector<int>> canvasCounts(
				newHeight, std::vector<int>(newWidth, 0));
			for (auto &circle : circles) {
				const int yStart =
					std::max(0, static_cast<int>(std::ceil(circle.center.y -
														   circle.radius)));
				const int yEnd = std::min(
					newHeight - 1, static_cast<int>(std::floor(circle.center.y +
															   circle.radius)));
				for (int y = yStart; y <= yEnd; ++y) {
					const auto xRangeOpt =
						calcCircleRangeX(circle, y, newWidth);
					if (!xRangeOpt.has_value()) {
						continue;
					}
					const Vec2i xRange = xRangeOpt.value();
					canvasColorsAccum[y][xRange[0]] += circle.color;
					canvasCounts[y][xRange[0]] += 1;
					if (xRange[1] != newWidth - 1) {
						canvasColorsAccum[y][xRange[1] + 1] -= circle.color;
						canvasCounts[y][xRange[1] + 1] -= 1;
					}
				}
			}
			for (int y = 0; y < newHeight; ++y) {
				for (int x = 1; x < newWidth; ++x) {
					canvasColorsAccum[y][x] += canvasColorsAccum[y][x - 1];
					canvasCounts[y][x] += canvasCounts[y][x - 1];
				}
			}
			for (auto &circle : circles) {
				Vec3f dLdcolor = {0.0f, 0.0f, 0.0f};
				float dLdx = 0.0f;
				float dLdy = 0.0f;
				float dLdr = 0.0f;

				const auto circleShiftedX =
					shiftCircle(circle, 1.0f, 0.0f, 0.0f);
				const auto circleShiftedY =
					shiftCircle(circle, 0.0f, 1.0f, 0.0f);
				const auto circleShiftedRadius =
					shiftCircle(circle, 0.0f, 0.0f, 1.0f);

				const int yStart =
					std::max(0, static_cast<int>(std::ceil(circle.center.y -
														   circle.radius)) -
									1);
				const int yEnd = std::min(
					newHeight - 1, static_cast<int>(std::floor(circle.center.y +
															   circle.radius)) +
									   1);
				for (int y = yStart; y <= yEnd; ++y) {
					const auto xRangeOpt =
						calcCircleRangeX(circle, y, newWidth);
					const auto xRangeShiftedXOpt =
						calcCircleRangeX(circleShiftedX, y, newWidth);
					const auto xRangeShiftedYOpt =
						calcCircleRangeX(circleShiftedY, y, newWidth);
					const auto xRangeShiftedRadiusOpt =
						calcCircleRangeX(circleShiftedRadius, y, newWidth);

					dLdx += calcLossDeltaFromRanges(
						xRangeOpt, xRangeShiftedXOpt, y, circle.color, image,
						canvasColorsAccum, canvasCounts);
					dLdy += calcLossDeltaFromRanges(
						xRangeOpt, xRangeShiftedYOpt, y, circle.color, image,
						canvasColorsAccum, canvasCounts);
					dLdr += calcLossDeltaFromRanges(
						xRangeOpt, xRangeShiftedRadiusOpt, y, circle.color,
						image, canvasColorsAccum, canvasCounts);
					if (!xRangeOpt.has_value()) {
						continue;
					}
					if (iter % 10 == 0) {
						const auto xRange = xRangeOpt.value();
						for (int x = xRange[0]; x <= xRange[1]; ++x) {
							const Vec3f &imageColor = image[y][x];
							const Vec3f &currentColorAccum =
								canvasColorsAccum[y][x];
							const int currentCount = canvasCounts[y][x];
							const Vec3f currentColor =
								currentCount != 0
									? currentColorAccum /
										  static_cast<float>(currentCount)
									: Vec3f{0.0f, 0.0f, 0.0f};
							const Vec3f diff = currentColor - imageColor;
							for (int c = 0; c < 3; ++c) {
								dLdcolor[c] += diff[c] * std::abs(diff[c]) /
											   static_cast<float>(currentCount);
							}
						}
					}
				}
				// Update circle parameters
				const float learningRateColor = 1.0f;

				const float radiusSq = circle.radius * circle.radius;
				circle.color -= dLdcolor * (learningRateColor / radiusSq);
				for (int c = 0; c < 3; ++c) {
					circle.color[c] = std::clamp(circle.color[c], 0.0f, 1.0f);
				}

				const float learningRatePosition = 1.0f;
				circle.center.x -= learningRatePosition * dLdx / circle.radius;
				circle.center.y -= learningRatePosition * dLdy / circle.radius;
				circle.center.x = std::clamp(circle.center.x, 0.0f,
											 static_cast<float>(newWidth - 1));
				circle.center.y = std::clamp(circle.center.y, 0.0f,
											 static_cast<float>(newHeight - 1));
				const float learningRateRadius = 0.01f;
				circle.radius -= learningRateRadius * dLdr / circle.radius;
				const float maxRadius = std::min(newWidth, newHeight);
				circle.radius = std::clamp(circle.radius, 1.0f, maxRadius);
			}
		}
		renderCircles(circles, 1.0f / static_cast<float>(scaleFactor));
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
