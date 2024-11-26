#include <iostream>
#include <stdio.h>
#include "F:\soft\mpi\MicrosoftSDKs\MPI\Include\mpi.h"
/*
compile using Windows cmd:
    g++ -o program.exe program.cpp -I"path\to\mpi\include" -L"path\to\mpi\lib" -lmsmpi
run:
    mpiexec -n <number of processes needed> progarm.exe
*/
using namespace std;
/*
计算函数f(x)= 3*x*x+2*x+1的阴影面积，用户命令行输入a,b,n 信息，
输入信息由0号进程广播给其他进程，各进程计算完后由allreduce来计算总面积
PESUDO CODE Trapezoidal Rule见ppt31
get a,b,n;
h = (b-a)/n;
local_n = n/commu_sz;
local_a = a + local_n*h*my_rank;
local_b = local_a + n*h;
local_integral = Trap(local_a, local_b, local_n, h);
if(my_rank != 0){
    Send local_integral to process0;
}
else{ 
    total_integral = local_integral;
    for( proc = 1; proc<comm_sz; proc++){
        Receive local_integral from proc;
        total_integral += local_integral;
    }
}
if(my_rank == 0){
    print result;
}
*/ 

// 定义目标函数 f(x) = 3*x*x + 2*x + 1
double f(double x) {
    return 3*x*x + 2*x + 1;
}

// 使用梯形规则计算积分
double Trap(double left_end, double right_end, int trap_count, double base_len) {
    double estimate, x;
    int i;

    estimate = (f(left_end) + f(right_end)) / 2.0;
    for (i = 1; i <= trap_count - 1; i++) {
        x = left_end + i * base_len;
        estimate += f(x);
    }
    estimate = estimate * base_len;

    return estimate;
}

int main(int argc, char* argv[]) {
    int my_rank, comm_sz, n, local_n;
    double a, b, h, local_a, local_b;
    double local_int, total_int;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // 由0号进程输入并广播a, b, n的值
    if (my_rank == 0) {
        printf("Enter a, b, and n:\n");
        fflush(stdout);
        scanf("%lf %lf %d", &a, &b, &n);
    }
    MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (n <= 0) {
        if (my_rank == 0) {
            fprintf(stderr, "Error: The number of trapezoids (n) must be greater than 0.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // 计算分配给当前进程的积分区间和积分结果
    h = (b - a) / n;
    local_n = n / comm_sz;

    local_a = a + my_rank * h * local_n;
    local_b = local_a + local_n * h;
    local_int = Trap(local_a, local_b, local_n, h);

    // 使用MPI_Allreduce求总积分
    MPI_Allreduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // 由0号进程输出最终结果
    if (my_rank == 0) {
        printf("With n = %d trapezoids, our estimate\n", n);
        printf("of the integral from %f to %f = %.15e\n", a, b, total_int);
    }

    MPI_Finalize();
    return 0;
}
