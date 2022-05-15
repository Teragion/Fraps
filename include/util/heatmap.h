#ifndef HEATMAP_H
#define HEATMAP_H

#include <vector>

#include "util.h"
#include "util/color.h"

class Heatmap {
private:
    struct ColorPoint { // Internal class used to store colors at different points in the gradient.
        float r,g,b;      // Red, green and blue values of our color.
        float val;        // Position of our color along the gradient (between 0 and 1).
        ColorPoint(float red, float green, float blue, float value)
        : r(red), g(green), b(blue), val(value) {}
    };
    std::vector<ColorPoint> color;      // An array of color points in ascending value.

public:
    //-- Default constructor:
    Heatmap();

    //-- Inserts a new color point into its correct position:
    void addColorPoint(float red, float green, float blue, float value);

    //-- Inserts a new color point into its correct position:
    void clearGradient();

    //-- Places a 5 color heapmap gradient into the "color" vector:
    void createDefaultHeatMapGradient();

    //-- Inputs a (value) between 0 and 1 and outputs the (red), (green) and (blue)
    //-- values representing that position in the gradient.
    void getColorAtValue(const float value, Color& c);
};

#endif
