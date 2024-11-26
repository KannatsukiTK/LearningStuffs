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
使用基本的MPISend和Receive完成
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
//f(x)= 3*x*x+2*x+1
double f(double x){
    return 3*x*x+2*x+1;
}

double Trap(
        double left_end,
        double right_end,
        int trap_count,
        double base_len){
    double estimate, x;
    int i;

    estimate = (f(left_end) + f(right_end)) / 2.0;
    for(i = 1; i<=trap_count-1; i++){
        x = left_end + i * base_len;
        estimate += f(x);
    }
    estimate = estimate * base_len;

    return estimate;
}

void GetInput(
        int my_rank,
        int comm_sz,
        double* a_p,
        double* b_p,
        int* n_p ){
    int dest;

    if(my_rank == 0){
        /*主进程读取stdin，并分配给其他进程*/
        printf(" Enter a, b and n:\n ");
        fflush(stdout); // 强制刷新缓冲区，立即输出提示
        scanf("%lf %lf %d", a_p, b_p, n_p);
        for(dest = 1; dest < comm_sz; dest++){
            MPI_Send(a_p, 1, MPI_DOUBLE, dest, 0 , MPI_COMM_WORLD);
            MPI_Send(b_p, 1, MPI_DOUBLE, dest, 0 , MPI_COMM_WORLD);
            MPI_Send(n_p, 1, MPI_INT, dest, 0 , MPI_COMM_WORLD);
        }
    }else{ /*my_rank != 0，接收参数*/
        MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
    }
}

int main(void){
    int my_rank, comm_sz, n, local_n;
    double a, b, h, local_a, local_b;
    double local_int, total_int;
    int source;

    // GetInput(my_rank, comm_sz, &a, &b, &n);

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    GetInput(my_rank, comm_sz, &a, &b, &n);
    /*getInput要放在rank初始化之后，它只接收abn三个参数。。。*/

    if (n <= 0) {
        if (my_rank == 0) {
            fprintf(stderr, "Error: The number of trapezoids (n) must be greater than 0.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    h = (b-a)/n;
    local_n = n / comm_sz;

    local_a = a + my_rank*h*local_n;
    local_b = local_a + local_n * h;
    local_int = Trap(local_a, local_b, local_n, h);

    if(my_rank != 0){
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else{ /*0号进程，收集其他进程的计算结果*/
        total_int = local_int;
        for(source = 1;source < comm_sz; source++){
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_int += local_int;
        }
    }
    if(my_rank == 0){
        printf("With n = %d trapezoids, our estimate\n", n);
        printf("of the integral from %f to %f = %.15e\n",
            a, b, total_int);
    }
    MPI_Finalize();
    return 0;
}
