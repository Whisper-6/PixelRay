#include <algorithm>
#include <conio.h>
#include <random>
#include <chrono>
#include <thread>
#include "render.h"

inline void normalize(Vector& vec) {
    vec = (1.0f / sqrtf(vec * vec)) * vec;
}

void addParallelogram(Point p1, Point p2, Point p3, Color col, ChunkSet& chunks) {
    Vector normal = vec3::cross(p2 - p1, p3 - p1);
    normalize(normal);
    Triangle now = {
        {p1,normal,col},
        {p2,normal,col},
        {p3,normal,col}
    };
    chunks.appendTriangle(now);
    now = {
        {p2 + p3 - p1,normal,col},
        {p2,normal,col},
        {p3,normal,col}
    };
    chunks.appendTriangle(now);
}

void gen0(ChunkSet& chunks) {
    Triangle now = {
        {make_vec(0,0,0),make_vec(0,0,1),make_col(1,0,0,1)},
        {make_vec(100,0,0),make_vec(0,0,1),make_col(0,1,0,1)},
        {make_vec(0,100,0),make_vec(0,0,1),make_col(0,0,1,1)},
    };
    chunks.appendTriangle(now);
}

void addCube(Point posM, int Len, Color col, ChunkSet& chunks) {
    addParallelogram(
        make_vec(-Len / 2, -Len / 2, Len) + posM,
        make_vec(Len / 2, -Len / 2, Len) + posM,
        make_vec(-Len / 2, Len / 2, Len) + posM,
        col,
        chunks
    );

    addParallelogram(
        make_vec(-Len / 2, -Len / 2, 0) + posM,
        make_vec(-Len / 2, -Len / 2, Len) + posM,
        make_vec(-Len / 2, Len / 2, 0) + posM,
        col,
        chunks
    );

    addParallelogram(
        make_vec(Len / 2, -Len / 2, 0) + posM,
        make_vec(Len / 2, Len / 2, 0) + posM,
        make_vec(Len / 2, -Len / 2, Len) + posM,
        col,
        chunks
    );

    addParallelogram(
        make_vec(-Len / 2, -Len / 2, 0) + posM,
        make_vec(Len / 2, -Len / 2, 0) + posM,
        make_vec(-Len / 2, -Len / 2, Len) + posM,
        col,
        chunks
    );

    addParallelogram(
        make_vec(-Len / 2, Len / 2, Len) + posM,
        make_vec(Len / 2, Len / 2, Len) + posM,
        make_vec(-Len / 2, Len / 2, 0) + posM,
        col,
        chunks
    );
}

void gen1(ChunkSet& chunks) {
    const float Len0 = 200, Len = 80;

    addParallelogram(
        make_vec(-Len0, -Len0, 0),
        make_vec(Len0, -Len0, 0),
        make_vec(-Len0, Len0, 0),
        make_col(1, 1, 1, 1),
        chunks
    );

    addCube(make_vec(0, 0, 0), Len, make_col(1, 1, 0, 0), chunks);
}

void addBall(Point posM, int Edges, float R, ChunkSet& chunks) {
    static Voxel _vox[205][205], (*vox)[205];
    vox = _vox;
    vox += Edges/4;
    const float alphaStep = 2 * pi / Edges;
    for (int i = - Edges/4; i <= Edges/4; i++)
        for (int j = 0; j < Edges; j++) {
            float beta = alphaStep * i, alpha = alphaStep * j;
            Vector e = make_vec(cos(beta) * sin(alpha), cos(beta) * cos(alpha), sin(beta));
            vox[i][j] = {
                R * e + posM,
                e,
                make_col(abs(cos(beta)),abs(cos(alpha)),1,1)
            };
        }
    for (int i = -Edges/4; i + 1 <= Edges/4; i++)
        for (int j = 0; j < Edges; j++) {
            Triangle now = {
                vox[i][j],
                vox[i + 1][j],
                vox[i + 1][(j + 1) % Edges]
            };
            chunks.appendTriangle(now);
            now = {
                vox[i][j],
                vox[i][(j + 1) % Edges],
                vox[i + 1][(j + 1) % Edges]
            };
            chunks.appendTriangle(now);
        }
    for (int j = 0; j < Edges; j++) {

    }
}


void gen2(ChunkSet& chunks) {
    const int Len0 = 320;

    addParallelogram(
        make_vec(-Len0, -Len0, 0),
        make_vec(Len0, -Len0, 0),
        make_vec(-Len0, Len0, 0),
        make_col(1, 1, 1, 1),
        chunks
    );

    addBall(make_vec(-50, -50, 50), 36, 50, chunks);

    addBall(make_vec(50, -50, 50), 36, 50, chunks);

    addBall(make_vec(-50, 50, 50), 36, 50, chunks);

    addBall(make_vec(50, 50, 50), 36, 50, chunks);
}

const int N = 288, M = 512, PixelScale = 2;

RayBoard<N, M> board;
Picture<N, M> pic;

void render1(ChunkSet& chunks, int Frames_limit, Point Mpos, Window& window) {

    /*board.setCamera(pi / 64, pi / 6, Mpos);
    board.setSun(pi / 3, pi / 4, Mpos);
    board.render(chunks, pic);
    pic.print(img, PixelScale);
    putimage(0, 0, &img);*/

    bool quit = false;
    SDL_Event e;
    int Frames = 0, sav0 = clock(), Flush_time = 0;
    while (!quit) {
        // 处理事件
        while (SDL_PollEvent(&e) != 0) {
            if (e.type == SDL_QUIT) {
                quit = true;
            }
        }

        board.setCamera(Frames*2e-2, pi / 4, Mpos);
        board.setSun(Frames*4e-2, pi / 4, Mpos);
        board.render(chunks, pic);

        Flush_time -= clock();
        pic.print(window.surface, PixelScale);
        window.flush();
        Flush_time += clock();

        Frames++;
        if (Frames == Frames_limit)
            break;
    }
    printf("cameraRender_time : %d\n", cameraRender_time);
    printf("sunRender_time : %d\n", sunRender_time);
    printf("SSAO_TIME : %d\n", SSAO_TIME);
    printf("PCSS_TIME : %d\n", PCSS_TIME);
    printf("Flush_time : %d\n", Flush_time);
    printf("Time: %d ms\n", clock() - sav0);
    printf("CPU Time: %d ms\n", cameraRender_time + sunRender_time + SSAO_TIME + PCSS_TIME + Flush_time);
    printf("Average fps: %.2f fps\n", Frames * 1e3 / (clock() - sav0));
}

const int Frames = 1024;

ChunkSet chunks;
Initializer initializer;
 
int main(int argc, char* argv[])
{
    Window window(M * PixelScale, N * PixelScale, "PixelRay");
    gen2(chunks);
    printf("Triangles : %d\n", chunks.tris.size());
    render1(chunks, Frames, make_vec(0,0,16), window);
    getch();
    return 0;
}