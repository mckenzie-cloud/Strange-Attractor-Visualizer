@echo off
g++ -O2 -s -Wall -std=c++17 src\*.cpp -I Include\ -L Lib -lsfml-main -lsfml-graphics -lsfml-system -lsfml-window -o bin\run