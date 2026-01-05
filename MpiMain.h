#ifndef MPIMAIN_H
#define MPIMAIN_H

#include <vector>
#include "Body.h"

class MpiMain {
public:
    void executeMain(int argc, char** argv);

private:
    void step_mpi(
        std::vector<Body>& localBodies,
        std::vector<Body>& allBodies,
        double dt,
        double G,
        double eps,
        int rank,
        int size,
        int localStart
    );
    
    double run_mpi(
        std::vector<Body>& bodies,
        int steps,
        double dt,
        double G,
        double eps,
        int rank,
        int size
    );
};

#endif
