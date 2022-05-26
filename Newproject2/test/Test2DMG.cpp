#include "Multigrid/MGFactory.h"
#include "Multigrid/TestMultigrid.h"
#include <unistd.h>

using namespace std;
using iVec = Vec<int,2>;
using rVec = Vec<Real,2>;

template <class T>
using Vector = std::vector<T>;

int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    // string name = argv[1];
    string file = "../../test/Input0.json";
    int rank = -1;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int num = 3;
    if(rank == 0) {
    cout << "V-cycle test:" << endl;
    cout << "--------------------------- Test function F = \\exp(x*y) "
            "---------------------------"
         << endl;
    }
    ScalarFunction<2>* pD2Func1f = new D2Func1f;
    ScalarFunction<2>* pD2Func1F = new D2Func1F;
    ScalarFunction<2>* pD2Func1Fx = new D2Func1Fx;
    ScalarFunction<2>* pD2Func1Fy = new D2Func1Fy;
    std::array<ScalarFunction<2>*, 3> pfuncs{pD2Func1F, pD2Func1Fx, pD2Func1Fy};
    MPI_Barrier(MPI_COMM_WORLD);
    auto t = MPI_Wtime();
    TestMultigrid<2> TM(file, true);
    TM.test(pD2Func1f, pfuncs, num, false, true);
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime() - t;
    sleep(1);
    if (rank == 0) {
      cout << "Payed time " <<  t << 
      " in " << size << " processes.";
    }

    MPI_Finalize();
}
