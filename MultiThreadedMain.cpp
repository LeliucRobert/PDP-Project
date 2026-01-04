#include "MultithreadedMain.h"

#include <thread>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>

static double checksum_bodies(const std::vector<Body>& bodies) {
    double s = 0.0;
    for (const auto& b : bodies) {
        s += b.mass
           + b.position.x + 2.0*b.position.y
           + 3.0*b.velocity.x + 4.0*b.velocity.y;
    }
    return s;
}

static std::vector<Body> init_random_bodies(int N, unsigned seed = 42) {
    std::vector<Body> bodies(N);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> posDist(-1.0, 1.0);
    std::uniform_real_distribution<double> velDist(-0.05, 0.05);
    std::uniform_real_distribution<double> massDist(0.5, 2.0);

    for (auto& b : bodies) {
        b.mass = massDist(gen);
        b.position = Vector2D(posDist(gen), posDist(gen));
        b.velocity = Vector2D(velDist(gen), velDist(gen));
    }
    return bodies;
}

static void compute_acc_range(
    const std::vector<Body>& bodies,
    int start, int end,
    double G, double eps,
    std::vector<Vector2D>& acc
) {
    const int n = (int)bodies.size();
    const double eps2 = eps * eps;

    for (int i = start; i < end; i++) {
        Vector2D a(0.0, 0.0);
        const double xi = bodies[i].position.x;
        const double yi = bodies[i].position.y;

        for (int j = 0; j < n; j++) {
            if (j == i) continue;

            double dx = bodies[j].position.x - xi;
            double dy = bodies[j].position.y - yi;

            double dist2 = dx*dx + dy*dy + eps2;

            double invDist = 1.0 / std::sqrt(dist2);
            double invDist3 = invDist * invDist * invDist;

            double mj = bodies[j].mass;
            double s = G * mj * invDist3;

            a.x += dx * s;
            a.y += dy * s;
        }
        acc[i] = a;
    }
}

static void integrate_range(
    std::vector<Body>& bodies,
    int start, int end,
    const std::vector<Vector2D>& acc,
    double dt
) {
    for (int i = start; i < end; i++) {
        bodies[i].velocity.x += acc[i].x * dt;
        bodies[i].velocity.y += acc[i].y * dt;
        bodies[i].position.x += bodies[i].velocity.x * dt;
        bodies[i].position.y += bodies[i].velocity.y * dt;
    }
}

void MultithreadedMain::step_threaded(std::vector<Body>& bodies, double dt, double G, double eps, int threads) {
    const int n = (int)bodies.size();
    if (n == 0) return;

    int T = std::max(1, threads);
    T = std::min(T, n);

    std::vector<Vector2D> acc(n, Vector2D(0.0, 0.0));
    std::vector<std::thread> pool;
    pool.reserve(T);

    for (int t = 0; t < T; t++) {
        int start = (t * n) / T;
        int end   = ((t + 1) * n) / T;
        pool.emplace_back([&, start, end] { compute_acc_range(bodies, start, end, G, eps, acc); });
    }
    for (auto& th : pool) th.join();

    pool.clear();
    for (int t = 0; t < T; t++) {
        int start = (t * n) / T;
        int end   = ((t + 1) * n) / T;
        pool.emplace_back([&, start, end] { integrate_range(bodies, start, end, acc, dt); });
    }
    for (auto& th : pool) th.join();
}

double MultithreadedMain::run_multithreaded(std::vector<Body>& bodies, int steps, double dt, double G, double eps, int threads) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int s = 0; s < steps; s++) step_threaded(bodies, dt, G, eps, threads);
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}

void MultithreadedMain::executeMain() {
    // Simple defaults (you can later parse from argv if you want)
    int N = 1000;
    int steps = 200;
    double dt = 0.01;
    double G = 1.0;
    double eps = 1e-3;
    int threads = (int)std::thread::hardware_concurrency();
    if (threads <= 0) threads = 4;

    auto bodies = init_random_bodies(N, 42);

    double t = run_multithreaded(bodies, steps, dt, G, eps, threads);

    std::cout << "[Threads] N=" << N
              << " steps=" << steps
              << " threads=" << threads
              << " time=" << t << "s\n";
    std::cout << "Checksum: " << checksum_bodies(bodies) << "\n";
}
