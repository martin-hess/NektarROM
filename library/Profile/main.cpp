#include <Profile/IntWrapper.h>
#include <boost/lexical_cast.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <Profile/VariableCostObject.h>
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <iostream>
using std::cout;
using std::endl;



void VariableCostOperationAddition2(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1;
	}
	cout << "VariableCostOperationAddition2 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition2 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition3(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2;
	}
	cout << "VariableCostOperationAddition3 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition3 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition4(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3;
	}
	cout << "VariableCostOperationAddition4 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition4 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition5(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4;
	}
	cout << "VariableCostOperationAddition5 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition5 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition6(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5;
	}
	cout << "VariableCostOperationAddition6 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition6 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition7(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6;
	}
	cout << "VariableCostOperationAddition7 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition7 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition8(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7;
	}
	cout << "VariableCostOperationAddition8 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition8 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition9(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8;
	}
	cout << "VariableCostOperationAddition9 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition9 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition10(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9;
	}
	cout << "VariableCostOperationAddition10 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition10 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition11(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10;
	}
	cout << "VariableCostOperationAddition11 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition11 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition12(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	VariableCostObject m11;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11;
	}
	cout << "VariableCostOperationAddition12 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition12 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition13(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	VariableCostObject m11;
	VariableCostObject m12;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12;
	}
	cout << "VariableCostOperationAddition13 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition13 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition14(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	VariableCostObject m11;
	VariableCostObject m12;
	VariableCostObject m13;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12 + m13;
	}
	cout << "VariableCostOperationAddition14 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition14 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

int main(int argc, char** argv)
{
	if( argc != 3 ) { cout << "Usage: Profile <numTests> <problemSize>\n"; return 1; }
    cout.precision(20);
	unsigned int numTests = boost::lexical_cast<unsigned int>(argv[1]);
	unsigned int problemSize = boost::lexical_cast<unsigned int>(argv[2]);
	VariableCostObject::size = problemSize;
//VariableCostOperationAddition2(numTests, problemSize);
//VariableCostOperationAddition3(numTests, problemSize);
//VariableCostOperationAddition4(numTests, problemSize);
//VariableCostOperationAddition5(numTests, problemSize);
//VariableCostOperationAddition6(numTests, problemSize);
//VariableCostOperationAddition7(numTests, problemSize);
//VariableCostOperationAddition8(numTests, problemSize);
//VariableCostOperationAddition9(numTests, problemSize);
//VariableCostOperationAddition10(numTests, problemSize);
//VariableCostOperationAddition11(numTests, problemSize);
//VariableCostOperationAddition12(numTests, problemSize);
//VariableCostOperationAddition13(numTests, problemSize);
VariableCostOperationAddition14(numTests, problemSize);
}
