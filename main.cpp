#include <algorithm>
#include <conio.h>
#include <random>
#include <chrono>
#include <thread>
#include "render.h"

inline void normalize(Vector& vec) {
    vec = (1.0 / sqrtf(vec * vec)) * vec;
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

void addBall2(Point center, int Edges, float R, ChunkSet& chunks) {
    static Voxel vox[205][205];
    const float alphaStep = 2 * pi / Edges;
    for (int i = 0; i < Edges / 6; i++)
        for (int j = 0; j < Edges; j++) {
            float beta = alphaStep * i, alpha = alphaStep * j;
            Vector e = make_vec(cos(beta) * sin(alpha), cos(beta) * cos(alpha), sin(beta));
            vox[i][j] = {
                R * e,
                e,
                make_col(abs(cos(beta)),abs(cos(alpha)),1,1)
            };
        }
    for (int i = 0; i + 1 < Edges / 6; i++)
        for (int j = 0; j < Edges; j++) {
            Triangle now = {
                vox[i][j],
                vox[i + 1][j],
                vox[i + 1][(j + 1) % Edges]
            };
            chunks.appendTriangle(now);
            now.v1.normal = (-1) * now.v1.normal;
            now.v2.normal = (-1) * now.v2.normal;
            now.v3.normal = (-1) * now.v3.normal;
            now.v1.pos = now.v1.pos + now.v1.normal;
            now.v2.pos = now.v2.pos + now.v2.normal;
            now.v3.pos = now.v3.pos + now.v3.normal;
            chunks.appendTriangle(now);
            now = {
                vox[i][j],
                vox[i][(j + 1) % Edges],
                vox[i + 1][(j + 1) % Edges]
            };
            chunks.appendTriangle(now);
            now.v1.normal = (-1) * now.v1.normal;
            now.v2.normal = (-1) * now.v2.normal;
            now.v3.normal = (-1) * now.v3.normal;
            now.v1.pos = now.v1.pos + now.v1.normal;
            now.v2.pos = now.v2.pos + now.v2.normal;
            now.v3.pos = now.v3.pos + now.v3.normal;
            chunks.appendTriangle(now);
        }
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

/*
void gen3(ChunkSet& chunks) {
    const float R = 80, D = 16, H = 100;
    const int Edges = 120;
    Vector up = make_vec(0, 0, 1);
    for (float r = 0; r < R * 2; r += 1/ sqrt(2)) {
        const float alphaStep = 2 * pi / Edges;
        Vector las = make_vec(0, 1, 0);
        for (float alpha = alphaStep; alpha < 2 * pi+1e-4; alpha += alphaStep) {
            Vector e = make_vec(sin(alpha), cos(alpha), 0);
            Segment now = {
                r * las,
                r * e,
                make_col(abs(sin(alpha + r/100)),0,1,1),
                make_col(abs(sin(alpha + r/100 + alphaStep)),0,1,1),
                make_vec(0,0,1),
                make_vec(0,0,1)
            };
            chunks.addSegment(now);
            las = e;
        }
    }

    for (float r = R - 5; r < R + 0.3; r += 1 / sqrt(2)) {
        const float alphaStep = 2 * pi / Edges;
        Vector las = make_vec(0, 1, 0);
        for (float alpha = 0; alpha < 2 * pi + 1e-4; alpha += alphaStep) {
            Vector e = make_vec(sin(alpha), cos(alpha), 0);
            Segment now = {
                r * las + make_vec(0,0,R/2+0.1),
                r * e + make_vec(0,0,R/2+0.1),
                make_col(abs(sin(alpha + r / 100)),0,1,1),
                make_col(abs(sin(alpha + r / 100 + alphaStep)),0,1,1),
                make_vec(0,0,1),
                make_vec(0,0,1)
            };
            chunks.addSegment(now);
            las = e;
        }
    }

    const float betaSetp = 1 / R / sqrt(2);
    for (float beta = 0 ; beta > -pi/6; beta -= betaSetp) {
        const float alphaStep = 2 * pi / Edges;
        Vector las = make_vec(0, cos(beta), sin(beta));
        for (float alpha = 0; alpha < 2 * pi + 1e-4; alpha += alphaStep) {
            Vector e = make_vec(cos(beta) * sin(alpha), cos(beta) * cos(alpha), sin(beta));
            Vector Vr = make_vec(cos(beta) * sin(alpha), cos(beta) * cos(alpha), 0);
            Segment now = {
                R * las + make_vec(0,0,R/2),
                R * e + make_vec(0,0,R / 2),
                make_col(abs(sin(alpha + R)),abs(cos(beta)),1,1),
                make_col(abs(sin(alpha + R + alphaStep)),abs(cos(beta)),1,1),
                las,
                e
            };
            chunks.addSegment(now);
            las = e;
        }
    }

    for (float beta = 0; beta > -pi / 6; beta -= betaSetp) {
        const float alphaStep = 2 * pi / Edges;
        Vector las = make_vec(0, cos(beta), sin(beta));
        for (float alpha = 0; alpha < 2 * pi + 1e-4; alpha += alphaStep) {
            Vector e = make_vec(cos(beta) * sin(alpha), cos(beta) * cos(alpha), sin(beta));
            Vector Vr = make_vec(cos(beta) * sin(alpha), cos(beta) * cos(alpha), 0);
            Segment now = {
                (R-5) * las + make_vec(0,0,R / 2),
                (R-5) * e + make_vec(0,0,R / 2),
                make_col(abs(sin(alpha + R)),abs(cos(beta)),1,1),
                make_col(abs(sin(alpha + R + alphaStep)),abs(cos(beta)),1,1),
                (-1) * las,
                (-1) * e
            };
            chunks.addSegment(now);
            las = e;
        }
    }
}

void gen4(ChunkSet& chunks) {
    const float R = 64, D = 16, H = 100;

    for (float r = 1; r < R * 2; r += 1 / sqrt(2)) {
        const float alphaStep = pi / (int)(8 * pow(R, 0.4));
        Vector las = make_vec(0, 1, 0);
        for (float alpha = 0; alpha < 2 * pi + 1e-4; alpha += alphaStep) {
            Vector e = make_vec(sin(alpha), cos(alpha), 0);
            Segment now = {
                r * las,
                r * e,
                make_col(1,1,1,1),
                make_col(1,1,1,1),
                make_vec(0,0,1),
                make_vec(0,0,1)
            };
            chunks.addSegment(now);
            las = e;
        }
    }


   for (float r = 1; r < R * 2; r += 1 / sqrt(2)) {
        const float alphaStep = pi / (int)(8 * pow(R, 0.4));
        Vector las = make_vec(0, 0, 1);
        for (float alpha = 0; alpha < 2 * pi + 1e-4; alpha += alphaStep) {
            Vector e = make_vec(sin(alpha), 0, cos(alpha));
            Segment now = {
                r * las,
                r * e,
                make_col(1,1,1,1),
                make_col(1,1,1,1),
                make_vec(0,1,0),
                make_vec(0,1,0)
            };
            chunks.addSegment(now);
            now.normal1 = (-1) * now.normal1;
            now.normal2 = (-1) * now.normal2;
            now.pos1 = now.pos1 + now.normal1;
            now.pos2 = now.pos2 + now.normal2;
            chunks.addSegment(now);
            las = e;
        }
    }

    for (float r = 1; r < R * 2; r += 1 / sqrt(2)) {
        const float alphaStep = pi / (int)(8 * pow(R, 0.4));
        Vector las = make_vec(0, 0, 1);
        for (float alpha = 0; alpha < 2 * pi + 1e-4; alpha += alphaStep) {
            Vector e = make_vec(0, sin(alpha), cos(alpha));
            Segment now = {
                r * las,
                r * e,
                make_col(1,1,1,1),
                make_col(1,1,1,1),
                make_vec(1,0,0),
                make_vec(1,0,0)
            };
            chunks.addSegment(now);
            now.normal1 = (-1) * now.normal1;
            now.normal2 = (-1) * now.normal2;
            now.pos1 = now.pos1 + now.normal1;
            now.pos2 = now.pos2 + now.normal2;
            chunks.addSegment(now);
            las = e;
        }
    }
}*/

const int N = 320, M = 480, PixelScale = 2;

RayBoard<N, M> board;
Picture<N, M> pic;
IMAGE img(M * PixelScale, N * PixelScale);

void render1(ChunkSet& chunks, int Frames, Point Mpos) {

    board.setCamera(pi / 64, pi / 6, Mpos);
    board.setSun(pi / 3, pi / 4, Mpos);
    board.render(chunks, pic);
    pic.print(img, PixelScale);
    putimage(0, 0, &img);

    int sav0 = clock();
    int Flush_time = 0;
    for (int i = 0; i < Frames; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        board.setCamera(i*2e-2, pi / 4, Mpos);
        board.setSun(i*4e-2, pi / 4, Mpos);
        board.render(chunks, pic);
        Flush_time -= clock();
        pic.print(img, PixelScale);
        putimage(0, 0, &img);
        Flush_time += clock();
    }
    printf("cameraRender_time : %d\n", cameraRender_time);
    printf("sunRender_time : %d\n", sunRender_time);
    printf("SSAO_TIME : %d\n", SSAO_TIME);
    printf("PCSS_TIME : %d\n", PCSS_TIME);
    printf("Flush_time : %d\n", Flush_time);
    printf("时间: %d ms\n", clock() - sav0);
    printf("时间: %d ms\n", cameraRender_time + sunRender_time + SSAO_TIME + PCSS_TIME + Flush_time);
    printf("平均帧率: %.2f fps\n", Frames * 1e3 / (clock() - sav0));
}

Picture<N, M> sav;

void render2(ChunkSet& chunks, int Frames, Point Mpos) {

    int sav0 = clock();
    float total = .0;
    const int Edges = 36;
    const float betaStep = pi / Edges;
    board.setCamera(pi / 6, pi/12, Mpos);
    for (float beta = betaStep; beta < pi / 2 - 1e-4; beta += betaStep) {
        const float alphaStep = pi / max(3, (int)(Edges * cos(beta)));
        for (float alpha = 0; alpha < 2 * pi - 1e-4; alpha += alphaStep)
            total += sin(beta);
    }
    for (float beta = betaStep; beta < pi / 2 - 1e-4; beta += betaStep) {
        const float alphaStep = pi / max(3,(int)(Edges*cos(beta)));
        for (float alpha = 0; alpha < 2 * pi- 1e-4; alpha += alphaStep) {
            board.setSun(alpha, beta, Mpos);
            board.render(chunks, pic);
            add(sav, pic);
        }
    }
    for (float beta = -pi / 2 + 1e-4; beta < -1e-4; beta += betaStep) {
        const float alphaStep = pi / max(3, (int)(Edges * cos(beta)));
        for (float alpha = 0; alpha < 2 * pi - 1e-4; alpha += alphaStep) {
            board.setSun(alpha, beta, Mpos);
            board.render(chunks, pic);
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    pic.col[i][j] = 0.33f * pic.col[i][j];
            add(sav, pic);
        }
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            sav.col[i][j] = (0.5f / total) * sav.col[i][j];

        }
    board.setSun(-pi / 3 - pi / 3, pi / 6, Mpos);
    board.render(chunks, pic);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            pic.col[i][j] = 0.5 * pic.col[i][j];
    add(sav, pic);
    sav.print(img, PixelScale);
    putimage(0, 0, &img);
    printf("时间: %d ms\n", clock() - sav0);
}

const int Frames = 1024;

ChunkSet chunks;
 
int main()
{
    initgraph(M * PixelScale, N * PixelScale);
    initGamma();
    gen1(chunks);
    printf("%d\n", chunks.tris.size());
    render1(chunks, Frames, make_vec(0,0,16));
    _getch();
    closegraph();
}