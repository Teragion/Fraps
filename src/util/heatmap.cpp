#include "util/heatmap.h"

//-- Default constructor:
Heatmap::Heatmap()  {  createDefaultHeatMapGradient();  }

//-- Inserts a new color point into its correct position:
void Heatmap::addColorPoint(float red, float green, float blue, float value) {
    for (unsigned int i = 0; i < color.size(); i++) {
    if (value < color[i].val) {
        color.insert(color.begin() + i, ColorPoint(red, green, blue, value));
        return;  }}
    color.push_back(ColorPoint(red, green, blue, value));
}

//-- Inserts a new color point into its correct position:
void Heatmap::clearGradient() { color.clear(); }

//-- Places a 5 color heapmap gradient into the "color" vector:
void Heatmap::createDefaultHeatMapGradient() {
    color.clear();
    color.push_back(ColorPoint(0, 0, 1,   0.0f));      // Blue.
    color.push_back(ColorPoint(0, 1, 1,   0.25f));     // Cyan.
    color.push_back(ColorPoint(0, 1, 0,   0.5f));      // Green.
    color.push_back(ColorPoint(1, 1, 0,   0.75f));     // Yellow.
    color.push_back(ColorPoint(1, 0, 0,   1.0f));      // Red.
}

//-- Inputs a (value) between 0 and 1 and outputs the (red), (green) and (blue)
//-- values representing that position in the gradient.
void Heatmap::getColorAtValue(const float value, Color& c) {
    if (color.size() == 0)
    return;

    for (unsigned int i = 0; i < color.size(); i++) {
        ColorPoint& currC = color[i];
        if (value < currC.val) {
            ColorPoint& prevC  = color[max(0u, i - 1)];
            float valueDiff    = (prevC.val - currC.val);
            float fractBetween = (valueDiff == 0) ? 0 : (value - currC.val) / valueDiff;
            c.r = (prevC.r - currC.r) * fractBetween + currC.r;
            c.g = (prevC.g - currC.g) * fractBetween + currC.g;
            c.b = (prevC.b - currC.b) * fractBetween + currC.b;
            return;
        }
    }
    c.r = color.back().r;
    c.g = color.back().g;
    c.b = color.back().b;
    return;
}
