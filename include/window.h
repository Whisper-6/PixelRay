#ifndef PIXELRAY_WINDOW_H
#define PIXELRAY_WINDOW_H

#include <SDL2/SDL.h>

// Window created by SDL
class Window {
private:
    SDL_Window* window;
public:
    SDL_Surface* surface;   // An Uint32 array, representing the screen
    Window(int width, int height, const char* title);
    ~Window();
    void flush();
};

#endif //PIXELRAY_WINDOW_H