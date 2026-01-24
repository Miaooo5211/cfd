
//说明指南:
//这是一个气体动理学格式的求解器，支持二维单块结构化网格模拟，目的是方便有需求的科研工作者入门之用。
//程序基本采用面向过程编写，没有什么弯弯绕，
//故建议从main函数 顺序 看起，搭配注释学习。注释采用土洋结合，以简洁方便的传递作者意图为原则。
//
//如用作学术用途，请引用任意下列文献：
//（1）X.JI, F.ZHAO, W.SHYY, & K.XU(2018).
//A family of high - order gas - kinetic schemes and its comparison with Riemann solver based high - order methods.
//Journal of Computational Physics, 356, 150 - 173.
//（2）X.JI, & K.XU(2020).Performance Enhancement for High - order Gas - kinetic Scheme Based on WENO - adaptive - order Reconstruction.
//Communication in Computational Physics, 28, 2, 539 - 590

//----原开源协议----

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


