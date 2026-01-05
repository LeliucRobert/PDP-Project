#include "MpiMain.h"

#include <mpi.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <array>

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

void MpiMain::step_mpi(
    std::vector<Body>& localBodies,
    std::vector<Body>& allBodies,
    double dt,
    double G,
    double eps,
    int rank,
    int size,
    int localStart
) {
    const int n = (int)allBodies.size();
    if (n == 0) return;

    // Each process computes accelerations for its local bodies
    std::vector<Vector2D> localAccelerations(localBodies.size(), Vector2D(0.0, 0.0));
    
    // Compute accelerations for local bodies by considering all bodies
    for (size_t i = 0; i < localBodies.size(); i++) {
        Vector2D totalForce(0.0, 0.0);
        
        for (size_t j = 0; j < allBodies.size(); j++) {
            // Skip self-interaction
            int globalIdx = localStart + i;
            if (globalIdx == (int)j) continue;
            
            Vector2D vectorToOtherBody = allBodies[j].position - localBodies[i].position;
            
            double distance = std::sqrt(
                vectorToOtherBody.x * vectorToOtherBody.x +
                vectorToOtherBody.y * vectorToOtherBody.y
            );
            
            if (distance < eps) distance = eps;
            
            Vector2D unitDirection = vectorToOtherBody / distance;
            
            // Newton: F = G * (m_i * m_j) / r^2   and   a = F / m_i
            double forceStrength =
                (G * localBodies[i].mass * allBodies[j].mass) / (distance * distance);
            
            totalForce += unitDirection * forceStrength;
        }
        
        localAccelerations[i] = totalForce / localBodies[i].mass;
    }
    
    // Update velocities and positions for local bodies
    for (size_t i = 0; i < localBodies.size(); i++) {
        localBodies[i].velocity += localAccelerations[i] * dt;
        localBodies[i].position += localBodies[i].velocity * dt;
    }
    
    // Gather updated local bodies back to allBodies (each process updates its portion)
    // We need to broadcast updated positions and velocities to all processes
    // Pack local body data into a buffer
    const int bodyDataSize = 5; // mass, pos.x, pos.y, vel.x, vel.y
    std::vector<double> localBuffer(localBodies.size() * bodyDataSize);
    
    for (size_t i = 0; i < localBodies.size(); i++) {
        auto data = localBodies[i].serialize();
        for (int j = 0; j < bodyDataSize; j++) {
            localBuffer[i * bodyDataSize + j] = data[j];
        }
    }
    
    // Gather all updated bodies from all processes
    std::vector<double> allBuffer(n * bodyDataSize);
    int localCount = localBodies.size() * bodyDataSize;
    std::vector<int> recvCounts(size);
    std::vector<int> displs(size);
    
    // Calculate counts and displacements
    int baseCount = n / size;
    int remainder = n % size;
    for (int i = 0; i < size; i++) {
        recvCounts[i] = (i < remainder ? baseCount + 1 : baseCount) * bodyDataSize;
        displs[i] = (i == 0 ? 0 : displs[i-1] + recvCounts[i-1]);
    }
    
    MPI_Allgatherv(
        localBuffer.data(), localCount, MPI_DOUBLE,
        allBuffer.data(), recvCounts.data(), displs.data(), MPI_DOUBLE,
        MPI_COMM_WORLD
    );
    
    // Unpack all bodies
    for (int i = 0; i < n; i++) {
        std::array<double, 5> data;
        for (int j = 0; j < bodyDataSize; j++) {
            data[j] = allBuffer[i * bodyDataSize + j];
        }
        allBodies[i] = Body::deserialize(data);
    }
}

double MpiMain::run_mpi(
    std::vector<Body>& bodies,
    int steps,
    double dt,
    double G,
    double eps,
    int rank,
    int size
) {
    const int n = (int)bodies.size();
    
    // Distribute bodies across processes
    int baseCount = n / size;
    int remainder = n % size;
    int localCount = (rank < remainder ? baseCount + 1 : baseCount);
    int localStart = rank * baseCount + std::min(rank, remainder);
    
    std::vector<Body> localBodies;
    localBodies.reserve(localCount);
    for (int i = 0; i < localCount; i++) {
        localBodies.push_back(bodies[localStart + i]);
    }
    
    // All processes need all body positions for force computation
    // So we keep a copy of all bodies
    std::vector<Body> allBodies = bodies;
    
    // Synchronize all processes
    MPI_Barrier(MPI_COMM_WORLD);
    
    auto t0 = std::chrono::high_resolution_clock::now();
    
    for (int s = 0; s < steps; s++) {
        step_mpi(localBodies, allBodies, dt, G, eps, rank, size, localStart);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    
    // Broadcast elapsed time from rank 0 (or get max time across all processes)
    double maxElapsed = elapsed;
    MPI_Allreduce(&elapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    elapsed = maxElapsed;
    
    // Gather final bodies to rank 0
    if (rank == 0) {
        const int bodyDataSize = 5;
        std::vector<double> allBuffer(n * bodyDataSize);
        std::vector<double> localBuffer(localBodies.size() * bodyDataSize);
        
        // Pack local bodies
        for (size_t i = 0; i < localBodies.size(); i++) {
            auto data = localBodies[i].serialize();
            for (int j = 0; j < bodyDataSize; j++) {
                localBuffer[i * bodyDataSize + j] = data[j];
            }
        }
        
        // Calculate receive counts and displacements
        std::vector<int> recvCounts(size);
        std::vector<int> displs(size);
        for (int i = 0; i < size; i++) {
            int procLocalCount = (i < remainder ? baseCount + 1 : baseCount);
            recvCounts[i] = procLocalCount * bodyDataSize;
            displs[i] = (i == 0 ? 0 : displs[i-1] + recvCounts[i-1]);
        }
        
        MPI_Gatherv(
            localBuffer.data(), localBodies.size() * bodyDataSize, MPI_DOUBLE,
            allBuffer.data(), recvCounts.data(), displs.data(), MPI_DOUBLE,
            0, MPI_COMM_WORLD
        );
        
        // Unpack all bodies
        for (int i = 0; i < n; i++) {
            std::array<double, 5> data;
            for (int j = 0; j < bodyDataSize; j++) {
                data[j] = allBuffer[i * bodyDataSize + j];
            }
            bodies[i] = Body::deserialize(data);
        }
    } else {
        // Other ranks send their local bodies
        const int bodyDataSize = 5;
        std::vector<double> localBuffer(localBodies.size() * bodyDataSize);
        
        for (size_t i = 0; i < localBodies.size(); i++) {
            auto data = localBodies[i].serialize();
            for (int j = 0; j < bodyDataSize; j++) {
                localBuffer[i * bodyDataSize + j] = data[j];
            }
        }
        
        MPI_Gatherv(
            localBuffer.data(), localBodies.size() * bodyDataSize, MPI_DOUBLE,
            nullptr, nullptr, nullptr, MPI_DOUBLE,
            0, MPI_COMM_WORLD
        );
    }
    
    return elapsed;
}

void MpiMain::executeMain(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int N = 1000;
    int steps = 200;
    double dt = 0.01;
    double G = 1.0;
    double eps = 1e-3;
    
    std::vector<Body> bodies;
    
    // Initialize bodies on rank 0, then broadcast
    if (rank == 0) {
        bodies = init_random_bodies(N, 42);
    } else {
        bodies.resize(N);
    }
    
    // Broadcast initial bodies to all processes
    const int bodyDataSize = 5;
    std::vector<double> buffer(N * bodyDataSize);
    
    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            auto data = bodies[i].serialize();
            for (int j = 0; j < bodyDataSize; j++) {
                buffer[i * bodyDataSize + j] = data[j];
            }
        }
    }
    
    MPI_Bcast(buffer.data(), N * bodyDataSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        for (int i = 0; i < N; i++) {
            std::array<double, 5> data;
            for (int j = 0; j < bodyDataSize; j++) {
                data[j] = buffer[i * bodyDataSize + j];
            }
            bodies[i] = Body::deserialize(data);
        }
    }
    
    double elapsed = run_mpi(bodies, steps, dt, G, eps, rank, size);
    
    if (rank == 0) {
        std::cout << "[MPI] N=" << N
                  << " steps=" << steps
                  << " processes=" << size
                  << " time=" << elapsed << "s\n";
        std::cout << "Checksum: " << checksum_bodies(bodies) << "\n";
    }
    
    MPI_Finalize();
}
