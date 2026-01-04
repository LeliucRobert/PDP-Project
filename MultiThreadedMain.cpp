#include "MultithreadedMain.h"

#include <thread>
#include <vector>
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

static void computeAccelerationsRange(
    const std::vector<Body>& bodies,
    int start,
    int end,
    double G,
    double minDistance,
    std::vector<Vector2D>& accelerations
) {
    int bodyCount = (int)bodies.size();

    for (int i = start; i < end; i++)
    {
        Vector2D totalForce(0.0, 0.0);

        for (int j = 0; j < bodyCount; j++)
        {
            if (i == j) continue;

            Vector2D vectorToOtherBody = bodies[j].position - bodies[i].position;

            double distance = std::sqrt(
                vectorToOtherBody.x * vectorToOtherBody.x +
                vectorToOtherBody.y * vectorToOtherBody.y
            );

            if (distance < minDistance) distance = minDistance;

            Vector2D unitDirection = vectorToOtherBody / distance;

            // Newton: F = G * (m_i * m_j) / r^2   and   a = F / m_i
            double forceStrength =
                (G * bodies[i].mass * bodies[j].mass) / (distance * distance);

            totalForce += unitDirection * forceStrength;
        }

        accelerations[i] = totalForce / bodies[i].mass;
    }
}

static void integrateRange(
    std::vector<Body>& bodies,
    int start,
    int end,
    const std::vector<Vector2D>& accelerations,
    double dt
) {
    for (int i = start; i < end; i++)
    {
        bodies[i].velocity += accelerations[i] * dt;
        bodies[i].position += bodies[i].velocity * dt;
    }
}

void MultithreadedMain::step_threaded(
    std::vector<Body>& bodies,
    double dt,
    double G,
    double eps,
    int threads
) {
    const int n = (int)bodies.size();
    if (n == 0) return;

    int T = std::max(1, threads);
    T = std::min(T, n);

    std::vector<Vector2D> accelerations(n, Vector2D(0.0, 0.0));
    std::vector<std::thread> pool;
    pool.reserve(T);

    for (int t = 0; t < T; t++) {
        int start = (t * n) / T;
        int end   = ((t + 1) * n) / T;

        pool.emplace_back([&, start, end] {
            computeAccelerationsRange(bodies, start, end, G, eps, accelerations);
        });
    }
    for (auto& th : pool) th.join();

    pool.clear();
    for (int t = 0; t < T; t++) {
        int start = (t * n) / T;
        int end   = ((t + 1) * n) / T;

        pool.emplace_back([&, start, end] {
            integrateRange(bodies, start, end, accelerations, dt);
        });
    }
    for (auto& th : pool) th.join();
}

double MultithreadedMain::run_multithreaded(
    std::vector<Body>& bodies,
    int steps,
    double dt,
    double G,
    double eps,
    int threads
) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int s = 0; s < steps; s++) {
        step_threaded(bodies, dt, G, eps, threads);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}

void MultithreadedMain::executeMain() {
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
