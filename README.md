<!-- @format -->

g++ -std=c++17 -O2 \
 main.cpp \
 MultithreadedMain.cpp \
 Body.cpp \
 Vector2D.cpp \
 -o nbody_threads \
 -pthread
