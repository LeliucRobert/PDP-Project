#ifndef MULTITHREADEDMAIN_H
#define MULTITHREADEDMAIN_H

#include <vector>
#include "Body.h"

class MultithreadedMain {
public:
    void executeMain();

private:
    void step_threaded(std::vector<Body>& bodies, double dt, double G, double eps, int threads);
    double run_multithreaded(std::vector<Body>& bodies, int steps, double dt, double G, double eps, int threads);
};

#endif
