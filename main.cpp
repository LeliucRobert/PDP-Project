#include <iostream>
#include <cstring>
#include "MultithreadedMain.h"
// #include "MpiMain.h"

// void executeMPI(int argc, char** argv)
// {
//     MpiMain mpiVersion;
//     mpiVersion.executeMain(argc, argv);
// }

void executeMultiThreaded()
{
    MultithreadedMain multithreadedVersion;
    multithreadedVersion.executeMain();
}

int main(int argc, char** argv)
{
    if (argc > 1)
    {
        if (std::strcmp(argv[1], "mpi") == 0)
            // executeMPI(argc, argv);
            executeMultiThreaded();
        else
            executeMultiThreaded();
    }
    else
    {
        // Choose default mode:
        executeMultiThreaded(); // I recommend this as default
        // executeMPI(argc, argv);
    }

    return 0;
}
