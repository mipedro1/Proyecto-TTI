# Proyecto TTI 
## Comando para compilar y ejecutar la aplicación principal:
g++ test/EKF_GEOS3.cpp src/*.cpp -lm -std=c++17 -o bin/main.exe
## Comando para compilar y ejecutar los test:
g++ test/tests.cpp src/*.cpp -lm -std=c++17 -o bin/tests.exe