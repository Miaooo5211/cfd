


#include<omp.h>  //omp并行，尖括号<>代表外部库文件
#include"output.h" //双引号代表该项目中的头文件

//以下是几个测试算例的头文件
#include"accuracy_test.h"
#include"boundary_layer.h"
#include"cylinder.h"
#include"riemann_problem.h"
//end

using namespace std; //默认std的命名规则

//程序运行开始，手动输入omp并行的线程数
//如果在函数名上面注释，那么函数在别处调用时，可以看到注释
void Set_omp_thread() 
{
	
	int num_thread;
	cout << "please_specify threads number for omp parallel:  ";
	cin >> num_thread;
	omp_set_num_threads(num_thread);
}

int main()
{
#pragma omp parallel
	{
#pragma omp single
		cout << "Threads used: " << omp_get_num_threads() << endl;
	}

	Set_omp_thread(); //设置omp并行线程
	make_directory_for_result(); //兼容linux的结果文件夹的临时方案

	cylinder();           
	
    return 0;
}


