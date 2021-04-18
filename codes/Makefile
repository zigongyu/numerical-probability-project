main: LoiExpon LoiGauss LoiGamma asa147 algo main.cpp
	g++ -o main main.cpp *.o -std=c++11
LoiExpon: LoiExpon.cpp LoiExpon.hpp
	g++ -c LoiExpon.cpp  -std=c++11

LoiGauss: LoiGauss.cpp LoiGauss.hpp
	g++ -c LoiGauss.cpp -std=c++11

LoiGamma: asa147 LoiGamma.cpp LoiGamma.hpp
	g++ -c LoiGamma.cpp  -std=c++11

asa147: asa147.cpp asa147.hpp
	g++ -c asa147.cpp -std=c++11

algo: algo.hpp algo.cpp
	g++ -c algo.cpp -std=c++11