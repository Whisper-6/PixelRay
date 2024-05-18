#include "window.h"

Window::Window(int width, int height, const char* title) {
    window = SDL_CreateWindow(
            title,
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED,
            width,
            height,
            SDL_WINDOW_SHOWN
    );
    if (!window) {
        SDL_Log("create window failed!%s\n", SDL_GetError());
        exit(-1);
    }
    surface = SDL_GetWindowSurface(window);
}

Window::~Window() {
    SDL_DestroyWindow(window);
}

void Window::flush() {
    SDL_UpdateWindowSurface(window);
}