cmake -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -GNinja -DCMAKE_BUILD_TYPE=Release -B out

# Building using CMake

cmake -B out

cmake --build out

./out/CmdRunner -THREADS 2 -SIM Jaccard